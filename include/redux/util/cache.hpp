#ifndef REDUX_UTIL_CACHE_HPP
#define REDUX_UTIL_CACHE_HPP

#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>

#include <vector>
#include <sys/types.h> 
#include <unistd.h> 


namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */
        enum CacheOptions {
            RAM_ONLY=1,                 //!< Don't store to disk, just use as shared/singleton storage for the global scope.
            STORE_MANUALLY=2,           //!< Requires manual calls to store/restore.
            STORE_INACTIVE=4,           //!< Store when inactive for a specified time. (make sure to call check before using!!)
            STORE_ON_DELETE=8,          //!< Write to disk before removing from cache.
            DELETE_UNREFERENCED=16,     //!< Delete item if no other reference exists.
            COMPRESSED_STORAGE=32       //!< Compress buffer using zlib before writing to disk.
        };
        
        
        class Cache {

        public:
            
            typedef std::function<void(void)> void_cb;
            typedef std::function<std::string(void)> string_cb;

            struct Info {
                Info( size_t id, const std::string& name, size_t typeSize=0 );
                Info( size_t id1, size_t id2, const std::string& name, size_t typeSize=0 );
                inline bool operator<( const Info& rhs ) const { if (id1==rhs.id1) return id2 < rhs.id2; return id1 < rhs.id1; }
                std::string getInfo(void) const;
                size_t id1;
                size_t id2;
                std::string name;
                size_t size;
                mutable size_t count;
                mutable size_t totalSize;
                mutable void_cb clear;
                mutable void_cb info;
            };

            static void cleanup(void);
            static std::string getStats(void);
            static Cache& get(void);
            std::string path(void);
            void setPath(const std::string&);
            static pid_t pid(void);
            template<class T>
            static int erase(const T& entry) {
                auto s = get().getSet<T>();
                static const Info& i = get().getSetInfo<T>();
                int ret = s.second.erase(entry);
                i.count = s.second.size();
                return ret;
            }
            template<class KeyT, class T>
            static int erase(const KeyT& key) {
                auto m = get().getMap<KeyT,T>();
                static const Info& i = get().getMapInfo<KeyT,T>();
                int ret = m.second.erase(key);
                i.count = m.second.size();
                return ret;
            }
            template<class KeyT, class T>
            static T& get(const KeyT& key, const T& val=T()) {
                auto m = get().getMap<KeyT,T>();        // m.first is a unique_lock for the map in m.second
                static const Info& i = get().getMapInfo<KeyT,T>();
                auto ret = m.second.emplace(key,val);
                i.count = m.second.size();
                return ret.first->second;
            }
            template<class T>
            static const std::shared_ptr<T>& get(const std::shared_ptr<T>& entry) {
                auto s = get().getSet<std::shared_ptr<T>, PtrCompare<T>>();
                static const Info& i = get().getSetInfo<std::shared_ptr<T>, PtrCompare<T>>();
                auto ret = s.second.emplace(entry);
                i.count = s.second.size();
                return *ret.first;
            }
            template<class T>
            static const T& get(const T& entry) {
                auto s = get().getSet<T>();
                static const Info& i = get().getSetInfo<T>();
                auto ret = s.second.emplace(entry);
                i.count = s.second.size();
                return *ret.first;
            }
            template<class KeyT, class T>
            static void clear(void) {
                auto m = get().getMap<KeyT,T>();
                static const Info& i = get().getMapInfo<KeyT,T>();
                m.second.clear();
                i.count = m.second.size();
            }
            template<class T>
            static void clear(void) {
                auto s = get().getSet<T>();
                static const Info& i = get().getSetInfo<T>();
                s.second.clear();
                i.count = s.second.size();
            }
            template<class KeyT, class T>
            static size_t getID1(void) {
                auto info = get().getMapInfo<KeyT,T>();
                return info.id1;
            }
            template<class T>
            static size_t getID1(void) {
                auto info = get().getSetInfo<T>();
                return info.id1;
            }
            template<class KeyT, class T>
            static std::string getName(void) {
                auto info = get().getMapInfo<KeyT,T>();
                return info.name;
            }
            template<class T>
            static std::string getName(void) {
                auto info = get().getSetInfo<T>();
                return info.name;
            }
            template<class KeyT, class T>
            static std::string stats(void) {
                auto info = get().getMapInfo<KeyT,T>();
                return info.getInfo();
            }
            template<class T>
            static std::string stats(void) {
                auto info = get().getSetInfo<T>();
                return info.getInfo();
            }
            template<class KeyT, class T>
            static size_t size(void) {
                auto m = get().getMap<KeyT,T>();
                return m.second.size();
            }
            template<class T>
            static size_t size(void) {
                auto s = get().getSet<T>();
                return s.second.size();
            }
            /*template<class T, class KeyT, class... Args>
            static std::shared_ptr<T>& get(const KeyT& key, Args... args) {
                Cache& c = get();
                auto& m = c.getMap<KeyT,T>();
                auto ret = m.emplace(key,std::shared_ptr<T>(new T(args...)));
                return ret.first->second;
                
            }*/
            template<class T, class U=std::less<T>>
            const Info& getSetInfo(void) {
                static const Info& info = initSetInfo<T,U>();
                return info;
            }
            template<class KeyT, class T>
            const Info& getMapInfo(void) {
                static const Info& info = initMapInfo<KeyT,T>();
                return info;
            }
            template<class KeyT, class T>
            std::pair<std::unique_lock<std::mutex>,std::map<KeyT,T>&> getMap(void) {
                static std::map<KeyT,T>& m = initMap<KeyT,T>();
                static const Info& info RDX_UNUSED = getMapInfo<KeyT,T>();
                static std::mutex mtx;
                std::pair<std::unique_lock<std::mutex>,std::map<KeyT,T>&> ret(std::unique_lock<std::mutex>(mtx),m);
                return std::move(ret);
            }
            template<class T, class U=std::less<T>>
            std::pair<std::unique_lock<std::mutex>,std::set<T,U>&> getSet(void) {
                static std::set<T,U>& s = initSet<T,U>();
                static const Info& info RDX_UNUSED = getSetInfo<T,U>();
                static std::mutex mtx;
                std::pair<std::unique_lock<std::mutex>,std::set<T,U>&> ret(std::unique_lock<std::mutex>(mtx),s);
                return std::move(ret);
            }
            template<class KeyT, class T>
            void for_each( std::function<void(std::pair<const KeyT,T>&)> func ) {
                auto m = get().getMap<KeyT,T>();
                std::for_each( m.second.begin(), m.second.end(), func );
            }
            template<class T, class U=std::less<T>>
            void for_each( std::function<void(T&)> func ) {
                auto s = get().getSet<T>();
                std::for_each( s.second.begin(), s.second.end(), func );
            }

        private:
            Cache() : path_(""), pid_(getpid()) {};
            template<class KeyT, class T>
            void mapMaintenance(void) {
                auto m = getMap<KeyT,T>();
            }
            template<class T, class U>
            void setMaintenance(void) {
                auto s = getSet<T,U>();
            }
            template<class KeyT, class T>
            const Info& initMapInfo( void ) {       // called only once when a new KeyT/T pair is used.
                const std::type_info& ki = typeid(KeyT);
                const std::type_info& vi = typeid(T);
                const std::string name = "std::map<"+demangle_name(ki.name())+","+demangle_name(vi.name())+">";
                Info info( ki.hash_code(), vi.hash_code(), name, sizeof(KeyT)+sizeof(T) );
                std::unique_lock<std::mutex> lock(mtx);
                auto it = caches.emplace(info);
                if( it.second ) {   // new element inserted
                    std::function<void(void)> func = std::bind(&Cache::mapMaintenance<KeyT,T>,this);
                    std::function<void(void)> cfunc = std::bind(&Cache::clear<KeyT,T>);
                    std::function<std::string(void)> sfunc = std::bind(&Cache::stats<KeyT,T>);
                    it.first->clear = cfunc;
                    it.first->info = sfunc;
                }
                return *(it.first);
            }
            template<class T, class U>
            const Info& initSetInfo( void ) {       // called only once when a new KeyT/T pair is used.
                const std::type_info& ti = typeid(T);
                const std::type_info& ui = typeid(U);
                const std::string name = "std::set<"+demangle_name(ti.name())+">";
                Info info( ti.hash_code(), ui.hash_code(), name, sizeof(T) );
                std::unique_lock<std::mutex> lock(mtx);
                auto it = caches.emplace(info);
                if( it.second ) {   // new element inserted
                    std::function<void(void)> func = std::bind(&Cache::setMaintenance<T,U>,this);
                    std::function<void(void)> cfunc = std::bind(&Cache::clear<T>);
                    std::function<std::string(void)> sfunc = std::bind(&Cache::stats<T>);
                    it.first->clear = cfunc;
                    it.first->info = sfunc;
                }
                return it.first;
            }
            template<class KeyT, class T>
            std::map<KeyT,T>& initMap(void) {       // called only once when a new KeyT/T pair is used.
                static std::map<KeyT,T> m;
                std::function<void(void)> func = std::bind(&Cache::mapMaintenance<KeyT,T>,this);
                std::function<void(void)> cfunc = std::bind(&Cache::clear<KeyT,T>);
                std::function<std::string(void)> sfunc = std::bind(&Cache::stats<KeyT,T>);
                std::unique_lock<std::mutex> lock(mtx);
                funcs.push_back(func);
                cleanup_funcs.push_back(cfunc);
                stats_funcs.push_back(sfunc);
                return m;
            }
            template<class T, class U>
            std::set<T,U>& initSet(void) {       // called only once when a new KeyT/T pair is used.
                static std::set<T,U> s;
                std::function<void(void)> func = std::bind(&Cache::setMaintenance<T,U>,this);
                std::function<void(void)> cfunc = std::bind(&Cache::clear<T>);
                std::function<std::string(void)> sfunc = std::bind(&Cache::stats<T>);
                std::unique_lock<std::mutex> lock(mtx);
                funcs.push_back(func);
                cleanup_funcs.push_back(cfunc);
                stats_funcs.push_back(sfunc);
                return s;
            }

            std::mutex mtx;
            std::string path_;
            pid_t pid_;
            std::set<Info> caches;
            std::vector<std::function<void(void)>> funcs;
            std::vector<std::function<void(void)>> cleanup_funcs;
            std::vector<std::function<std::string(void)>> stats_funcs;
            int pollTime;
        };
        
        /*! @} */

    }


}

#endif // REDUX_UTIL_CACHE_HPP
