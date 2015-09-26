#ifndef REDUX_UTIL_CACHE_HPP
#define REDUX_UTIL_CACHE_HPP

#include "redux/util/stringutil.hpp"

#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

#include <iostream>

namespace bfs=boost::filesystem;

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
            struct ptrcomp {
                template <class T>
                bool operator() (const std::shared_ptr<T>& lhs, const std::shared_ptr<T>& rhs) const { return (*lhs)<(*rhs); }
            };
        public:
            static Cache& get(void);
            std::string path(void);
            void setPath(const std::string&);
            template<class T>
            static int erase(const T& entry) {
                auto s = get().getSet<T>();
                int ret = s.second.erase(entry);
                return ret;
            }
            template<class KeyT, class T>
            static int erase(const KeyT& key) {
                auto m = get().getMap<KeyT,T>();
                int ret = m.second.erase(key);
                return ret;
            }
            template<class KeyT, class T>
            static T& get(const KeyT& key, const T& val) {
                auto m = get().getMap<KeyT,T>();        // m.first is a unique_lock for the map in m.second
                auto ret = m.second.emplace(key,val);
                return ret.first->second;
            }
            template<class T>
            static const std::shared_ptr<T>& get(const std::shared_ptr<T>& entry) {
                auto s = get().getSet<std::shared_ptr<T>,ptrcomp>();
                auto ret = s.second.emplace(entry);
                return *ret.first;
            }
            template<class T>
            static const T& get(const T& entry) {
                auto s = get().getSet<T>();
                auto ret = s.second.emplace(entry);
                return *ret.first;
            }
            template<class KeyT, class T>
            static void clear(void) {
                auto m = get().getMap<KeyT,T>();
                m.second.clear();
            }
            template<class T>
            static void clear(void) {
                auto s = get().getSet<T>();
                s.second.clear();
            }
            /*template<class T, class KeyT, class... Args>
            static std::shared_ptr<T>& get(const KeyT& key, Args... args) {
                Cache& c = get();
                auto& m = c.getMap<KeyT,T>();
                auto ret = m.emplace(key,std::shared_ptr<T>(new T(args...)));
                std::cout << "blaha: " << hexString(ret.first->second.get()) << std::endl; 
                return ret.first->second;
                
            }*/
            template<class KeyT, class T>
            std::pair<std::unique_lock<std::mutex>,std::map<KeyT,T>&> getMap(void) {
                static std::map<KeyT,T>& m = initMap<KeyT,T>();
                static std::mutex mtx;
                std::pair<std::unique_lock<std::mutex>,std::map<KeyT,T>&> ret(std::unique_lock<std::mutex>(mtx),m);
                return std::move(ret);
            }
            template<class T, class U=std::less<T>>
            std::pair<std::unique_lock<std::mutex>,std::set<T,U>&> getSet(void) {
                static std::set<T,U>& s = initSet<T,U>();
                static std::mutex mtx;
                std::pair<std::unique_lock<std::mutex>,std::set<T,U>&> ret(std::unique_lock<std::mutex>(mtx),s);
                return std::move(ret);
            }

        private:
            Cache(){};
            template<class KeyT, class T>
            void mapMaintenance(void) {
                auto m = getMap<KeyT,T>();
                std::cout << "mapMaintenance:  &m = " << redux::util::hexString(&m.second) << "   m.sz = " << m.second.size() << std::endl;
            }
            template<class T, class U>
            void setMaintenance(void) {
                auto s = getSet<T,U>();
                std::cout << "setMaintenance:  &m = " << redux::util::hexString(&s.second) << "   s.sz = " << s.second.size() << std::endl;
            }
            template<class KeyT, class T>
            std::map<KeyT,T>& initMap(void) {       // called only once when a new KeyT/T pair is used.
                static std::map<KeyT,T> m;
                std::function<void(void)> func = std::bind(&Cache::mapMaintenance<KeyT,T>,this);
                std::unique_lock<std::mutex> lock(cacheMutex);
                funcs.push_back(func);
                return m;
            }
            template<class T, class U>
            std::set<T,U>& initSet(void) {       // called only once when a new KeyT/T pair is used.
                static std::set<T,U> s;
                std::function<void(void)> func = std::bind(&Cache::setMaintenance<T,U>,this);
                std::unique_lock<std::mutex> lock(cacheMutex);
                funcs.push_back(func);
                return s;
            }

            std::mutex cacheMutex;
            static std::string cachePath;
            std::vector<std::function<void(void)>> funcs;
            int pollTime;
        };
        
        
        class CacheItem {
        public:
            CacheItem(void);
            CacheItem(const std::string& path);
            CacheItem(const CacheItem&);
            
            /***** Methods to be overloaded *****/
            virtual size_t csize(void) const { return 0; };                //!< returns the necessary buffer size for the call to "cpack"
            virtual uint64_t cpack(char*) const { return 0; };             //!< put all relevant data into a char-array
            virtual uint64_t cunpack(const char*, bool) { return 0; };     //!< read all relevant data from the char array (inverse of the above call)
            virtual void cclear(void) {};                                  //!< free whatever memory possible, while being able to restore with "cunpack"
            /************************************/
            
            virtual void cacheTouch(void);
            virtual bool cacheLoad(bool removeAfterLoad=false);
            virtual bool cacheStore(bool clearAfterStore=false);
            void setLoaded(bool il=true);
            virtual void setPath(const std::string& path);
            std::string path(void) const { return itemPath.string(); };
            
        protected:
            bfs::path itemPath;
            int timeout;
            int lastAccessed;
            int options;
            std::mutex itemMutex;
            bool isLoaded;
            
        };
        
        
        
        /*! @} */

    }


}

#endif // REDUX_UTIL_CACHE_HPP
