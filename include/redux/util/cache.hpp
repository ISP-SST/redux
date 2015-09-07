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
                Cache& c = get();
                auto& s = c.getSet<T>();
                int ret = s.erase(entry);
                //if ( ret ) std::cout << "Cache::get<T>(): erased value " << entry << "  new sSize = " << s.size() << printArray(s,"set") << std::endl; 
                return ret;
                
            }
            template<class KeyT, class T>
            static int erase(const T& key) {
                Cache& c = get();
                auto& m = c.getMap<KeyT,T>();
                int ret = m.erase(key);
                //if ( ret ) std::cout << "Cache::get<T>(): erased value " << entry << "  new sSize = " << s.size() << printArray(s,"set") << std::endl; 
                return ret;
                
            }
            template<class KeyT, class T>
            static T& get(const KeyT& key, const T& val) {
                Cache& c = get();
                auto& m = c.getMap<KeyT,T>();
                auto ret = m.emplace(key,val);
                //std::cout << "CacheItem::get<KeyT,T>(): returning item " << hexString(ret.first->second.get()) << "   mSize = " << m.size() << std::endl; 
                return ret.first->second;
                
            }
            template<class T>
            static const std::shared_ptr<T>& get(const std::shared_ptr<T>& entry) {
                Cache& c = get();
                auto& s = c.getSet<std::shared_ptr<T>,ptrcomp>();
                auto ret = s.emplace(entry);
                //std::cout << "Cache::get<T>(): returning ptr " << hexString(ret.first->get()) << "   sSize = " << s.size() << std::endl; 
                return *ret.first;
                
            }
            template<class T>
            static const T& get(const T& entry) {
                Cache& c = get();
                auto& s = c.getSet<T>();
                auto ret = s.emplace(entry);
               // if ( ret.second ) std::cout << "Cache::get<T>(): inserted value " << *ret.first << "  new sSize = " << s.size() << printArray(s,"set") << std::endl; 
                return *ret.first;
                
            }
            template<class KeyT, class T>
            void clear(void) {
                auto& m = getMap<KeyT,T>();
                m.clear();
            }
            template<class T>
            void clear(void) {
                auto& s = getSet<T>();
                s.clear();
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
            std::map<KeyT,T>& getMap(void) {
                static std::map<KeyT,T>& m = initMap<KeyT,T>();
                return m;
            }
            template<class T, class U=std::less<T>>
            std::set<T,U>& getSet(void) {
                static std::set<T,U>& s = initSet<T,U>();
                return s;
            }

        private:
            Cache(){};
            template<class KeyT, class T>
            void mapMaintenance(std::map<KeyT,T>& m) {
                std::cout << "mapMaintenance:  &m = " << redux::util::hexString(&m) << "   m.sz = " << m.size() << std::endl;
            }
            template<class T, class U>
            void setMaintenance(std::set<T,U>& s) {
                std::cout << "setMaintenance:  &m = " << redux::util::hexString(&s) << "   s.sz = " << s.size() << std::endl;
            }
            template<class KeyT, class T>
            std::map<KeyT,T>& initMap(void) {       // called only once when a new KeyT/T pair is used.
                static std::map<KeyT,T> m;
                std::function<void(void)> func = std::bind(&Cache::mapMaintenance<KeyT,T>,this,std::ref(m));
                std::unique_lock<std::mutex> lock(cacheMutex);
                funcs.push_back(func);
                return m;
            }
            template<class T, class U>
            std::set<T,U>& initSet(void) {       // called only once when a new KeyT/T pair is used.
                static std::set<T,U> s;
                std::function<void(void)> func = std::bind(&Cache::setMaintenance<T,U>,this,std::ref(s));
                std::unique_lock<std::mutex> lock(cacheMutex);
                funcs.push_back(func);
                return s;
            }

            std::mutex cacheMutex;
            static std::string cachePath;
            std::vector<std::function<void(void)>> funcs;
            int pollTime;
        };
        
        
        
        //template<class KeyT, class T>
        class CacheItem /*: public std::enable_shared_from_this<CacheItem>*/ {
        public:
            CacheItem() : isLoaded(false) {};
            
            //template<class... Args>
            /***** Methods to be overloaded *****/
            virtual size_t csize(void) const { return 0; };                //!< returns the necessary buffer size for the call to "cpack"
            virtual uint64_t cpack(char*) const { return 0; };             //!< put all relevant data into a char-array
            virtual uint64_t cunpack(const char*, bool) { return 0; };     //!< read all relevant data from the char array (inverse of the above call)
            virtual void cclear(void) {};                                  //!< free whatever memory possible, while being able to restore with "cunpack"
            /************************************/
            
            void cacheTouch(void) { lastAccessed = 0; }
            bool cacheLoad(bool removeAfterLoad=false){
                if(!isLoaded) {
                    std::unique_lock<std::mutex> lock(itemMutex);
                    //std::cout << "CacheItem::cacheLoad() itemPath = " << itemPath << std::endl;
                    if ( bfs::exists(itemPath) ) {
                        std::ifstream in(itemPath.c_str(), std::ofstream::binary|std::ios_base::ate);
                        if( in.good() ) {
                            size_t sz = in.tellg();
                    //std::cout << "CacheItem::cacheLoad() sz = " << sz << std::endl;
                            std::unique_ptr<char[]> buf( new char[sz] );
                            in.seekg(0, std::ios_base::beg);
                            in.clear();
                            in.read(buf.get(), sz);
                            size_t psz = cunpack(buf.get(),false);
                    //std::cout << "CacheItem::cacheLoad() psz = " << psz << std::endl;
                            if(psz == sz) {
                                isLoaded = true;
                                if(removeAfterLoad) {
                                    in.close();
                                    bfs::remove(itemPath);
                                }
                                return true;
                            } else std::cout << "CacheItem::cacheLoad() " << psz << " != " << sz << " !!!" << std::endl;
                        } else std::cout << "CacheItem::cacheLoad() in.good() == false !!!" << std::endl;
                    } //else std::cout << "CacheItem::cacheLoad() " << itemPath << " doesn't exist. !!!" << std::endl;
                }
                return false;
            }
            bool cacheStore(bool clearAfterStore=false){
                if(isLoaded) {
                    std::unique_lock<std::mutex> lock(itemMutex);
                    //std::cout << "CacheItem::cacheStore() itemPath = " << itemPath << std::endl;
                    bfs::create_directories(itemPath.parent_path());
                    size_t sz = csize();
                    std::unique_ptr<char[]> buf( new char[sz] );
                    size_t psz = cpack(buf.get());
                    if(psz == sz) {
                        std::ofstream out(itemPath.c_str(), std::ofstream::binary);
                        if( out.good() ) {
                            out.write( buf.get(), sz );
                            if(clearAfterStore) {
                                cclear();
                                isLoaded = false;
                            }
                            return true;
                        }
                    }
                }
                return false;
            }
            
            void setPath(const std::string& path) {
                Cache& c = Cache::get();
                itemPath = bfs::path(c.path()) / bfs::path(path);
            }
            
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
