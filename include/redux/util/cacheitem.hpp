#ifndef REDUX_UTIL_CACHEITEM_HPP
#define REDUX_UTIL_CACHEITEM_HPP

#ifdef RDX_TRACE_CACHE
#   include "redux/util/trace.hpp"
#endif

#include <mutex>
#include <string>

#include <boost/filesystem.hpp>

namespace bfs=boost::filesystem;

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */
      
        class CacheItem
#ifdef RDX_TRACE_CACHE
            : public redux::util::TraceObject<CacheItem>
#endif
        {
        public:
            CacheItem(void);
            explicit CacheItem(const std::string& path);
            CacheItem(const CacheItem&);
            virtual ~CacheItem();
            
            /***** Methods to be overloaded *****/
            virtual size_t csize(void) const { return 0; };                //!< returns the necessary buffer size for the call to "cpack"
            virtual uint64_t cpack(char*) const { return 0; };             //!< put all relevant data into a char-array
            virtual uint64_t cunpack(const char*, bool) { return 0; };     //!< read all relevant data from the char array (inverse of the above call)
            virtual void cclear(void) {};                                  //!< free whatever memory possible, while being able to restore with "cunpack"
            /************************************/
            
            void cacheClear(void);
            void cacheRemove(void);
            virtual void cacheTouch(void);
            virtual bool cacheLoad(bool removeAfterLoad=false);
            virtual bool cacheStore(bool clearAfterStore=false);
            void setLoaded(bool il=true);
            virtual const std::string& getFullPath(void) { return fullPath.string(); }
            virtual void setPath( const std::string& path );
            std::string path(void) const { return itemPath.string(); };

            std::unique_lock<std::mutex> getLock(bool trylock=false) {
                if(trylock) return std::move( std::unique_lock<std::mutex>(itemMutex,std::try_to_lock) );
                return std::move( std::unique_lock<std::mutex>(itemMutex) );
            }

        protected:
            bfs::path itemPath;
            bfs::path fullPath;
            std::mutex itemMutex;
            bool isLoaded;
            uint64_t cachedSize;
        };
        
        
        
        /*! @} */

    }


}

#endif // REDUX_UTIL_CACHEITEM_HPP
