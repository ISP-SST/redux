#include "redux/util/cacheitem.hpp"

#ifdef DEBUG_
#   define TRACE_THREADS
#endif

#include "redux/util/cache.hpp"
#include "redux/file/fileio.hpp"
#include "redux/util/trace.hpp"

#include <fstream>
#include <memory>


using namespace redux::file;
using namespace redux::util;
using namespace std;


CacheItem::CacheItem(void) : itemPath(""), fullPath(""), isLoaded(true), cachedSize(0) {
    
}


CacheItem::CacheItem(const string& path) : itemPath(""), fullPath(""), isLoaded(true), cachedSize(0) {
    CacheItem::setPath(path);
    unique_lock<mutex> lock(itemMutex);
    if ( bfs::exists(fullPath) ) {
        bfs::remove(fullPath);      // HACK  for now we just delete pre-existing file, later we want to recover after crash...
    }
}


CacheItem::CacheItem(const CacheItem& rhs) : itemPath(rhs.itemPath), fullPath(rhs.fullPath), isLoaded(rhs.isLoaded), cachedSize(rhs.cachedSize) {
         
}


CacheItem::~CacheItem() {
    
}


void CacheItem::cacheClear(void) {

    THREAD_MARK
    unique_lock<mutex> lock(itemMutex);
    THREAD_MARK
    try {
        cclear();
        isLoaded = false;
    } catch( exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "CacheItem::cacheClear() failed: " << e.what() << endl;
        } else throw;
    }
    
}


void CacheItem::cacheRemove(void) {

    unique_lock<mutex> lock(itemMutex);
    try {
        bfs::remove(fullPath);
    } catch( exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "CacheItem::cacheRemove() failed: " << e.what() << endl;
        } else throw;
    }
    
}


void CacheItem::cacheTouch(void) {
    unique_lock<mutex> lock(itemMutex);
//     lastAccessed = 0;
}


bool CacheItem::cacheLoad(bool removeAfterLoad) {

    bool ret(false);
    try {
        THREAD_MARK
        unique_lock<mutex> lock(itemMutex);
        THREAD_MARK
        if( fullPath.empty() ) {
            return ret;
        }

        if( !isLoaded ) {
            if ( bfs::exists(fullPath) ) {
                ifstream in(fullPath.string().c_str(), ofstream::binary|ios_base::ate);
                if( in.good() ) {
                    size_t sz = in.tellg();
                    shared_ptr<char> buf = rdx_get_shared<char>(sz);
                    memset( buf.get(), 0, sz);
                    in.seekg(0, ios_base::beg);
                    in.clear();
                    in.read(buf.get(), sz);
                    size_t psz = cunpack(buf.get(),false);
                    if( psz >= sz ) {
                        isLoaded = true;
                        ret = true;
                    }
                }
            }
        }
        THREAD_MARK
        if( isLoaded && removeAfterLoad ) {
            bfs::remove(fullPath);
        }
    } catch( exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "CacheItem::cacheLoad() failed: " << e.what() << endl;
        } else throw;
    }
    
    return ret;
    
}


bool CacheItem::cacheStore(bool clearAfterStore){

    bool ret(false);
    try {
        THREAD_MARK
        unique_lock<mutex> lock(itemMutex);
        THREAD_MARK
        if(fullPath.empty()) {
            return ret;
        }

        if( isLoaded ) {
            bfs::path parent = fullPath.parent_path();
            if( !parent.empty() ) {
                bfs::create_directories( fullPath.parent_path() );
            }
            cachedSize = csize();
            shared_ptr<char> buf = rdx_get_shared<char>(cachedSize);
            size_t psz = cpack( buf.get() );
            if( psz <= cachedSize ) {
                ofstream out( fullPath.string().c_str(), ofstream::binary );
                if( out.good() ) {
                    out.write( buf.get(), psz );
                    ret = true;
                }
            }
        }
        THREAD_MARK
        if(ret && clearAfterStore) {
            cclear();
            isLoaded = false;
        }
    } catch( exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "CacheItem::cacheStore() failed: " << e.what() << endl;
        } else throw;
    }
    
    return ret;
    
}


void CacheItem::setLoaded(bool il) {
    
    unique_lock<mutex> lock(itemMutex);
    isLoaded = il;
    
}


void CacheItem::setPath(const string& path) {
    
    unique_lock<mutex> lock(itemMutex);
    itemPath = bfs::path(path);
    Cache& c = Cache::get();
    fullPath = bfs::path(c.path()) / itemPath;
    
}
        
