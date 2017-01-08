#include "redux/util/cache.hpp"


#include <iostream>
#include <memory>


using namespace redux::util;
using namespace std;


void Cache::cleanup(void) {
    for( auto& func: get().cleanup_funcs ) {
        func();
    }
}

std::string Cache::getStats(void) {
    string ret;
    for( auto& func: get().stats_funcs ) {
        ret += func() + "\n";
    }
    return ret;
}

Cache& Cache::get(void) {
    static Cache cache;
    return cache;
}


string Cache::path(void) {
    unique_lock<mutex> lock(mtx);
    return path_;
}


void Cache::setPath(const string& path) {
    unique_lock<mutex> lock(mtx);
    path_ = path;
}


pid_t Cache::pid(void) {
    return get().pid_;
}


CacheItem::CacheItem(void) : itemPath(""), fullPath(""), isLoaded(true), cachedSize(0) {
    
}


CacheItem::CacheItem(const string& path) : itemPath(""), fullPath(""), isLoaded(true), cachedSize(0) {
    setPath(path);
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

    unique_lock<mutex> lock(itemMutex);
    try {
        cclear();
        isLoaded = false;
    } catch( exception& e ) {
        cerr << "CacheItem::cacheClear() failed: " << e.what() << endl;
    } catch( ... ) {
        cerr << "CacheItem::cacheClear() failed for unknown reasons." << endl;
    }
    
}


void CacheItem::cacheRemove(void) {

    unique_lock<mutex> lock(itemMutex);
    try {
        bfs::remove(fullPath);
    } catch( exception& e ) {
        cerr << "CacheItem::cacheRemove() failed: " << e.what() << endl;
    } catch( ... ) {
        cerr << "CacheItem::cacheRemove() failed for unknown reasons." << endl;
    }
    
}


void CacheItem::cacheTouch(void) {
    unique_lock<mutex> lock(itemMutex);
//     lastAccessed = 0;
}


bool CacheItem::cacheLoad(bool removeAfterLoad) {

    bool ret(false);
    try {
        unique_lock<mutex> lock(itemMutex);
        if( fullPath.empty() ) {
            //cout << "Not loading b.c. non-existing: " << fullPath << endl;
            return ret;
        }

        if( !isLoaded ) {
            if ( bfs::exists(fullPath) ) {
                ifstream in(fullPath.string().c_str(), ofstream::binary|ios_base::ate);
                if( in.good() ) {
                    size_t sz = in.tellg();
                    unique_ptr<char[]> buf( new char[sz] );
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
        } //else cout << "Not loading b.c. loaded: " << fullPath << endl;

        if( isLoaded && removeAfterLoad ) {
            bfs::remove(fullPath);
        }
    } catch( exception& e ) {
        cerr << "CacheItem::cacheLoad() failed: " << e.what() << endl;
    } catch( ... ) {
        cerr << "CacheItem::cacheLoad() failed for unknown reasons." << endl;
    }
    
    return ret;
    
}


bool CacheItem::cacheStore(bool clearAfterStore){

    bool ret(false);
    try {
        unique_lock<mutex> lock(itemMutex);
        if(fullPath.empty()) {
            return ret;
        }

        if( isLoaded ) {
            bfs::path parent = fullPath.parent_path();
            if( !parent.empty() ) {
                bfs::create_directories( fullPath.parent_path() );
            }
            cachedSize = csize();
            unique_ptr<char[]> buf( new char[cachedSize] );
            size_t psz = cpack( buf.get() );
            if( psz <= cachedSize ) {
                ofstream out( fullPath.string().c_str(), ofstream::binary );
                if( out.good() ) {
                    out.write( buf.get(), psz );
                    ret = true;
                }
            }
        } //else cout << "Not storing b.c. not loaded: " << fullPath << endl;
        
        if(ret && clearAfterStore) {
            cclear();
            isLoaded = false;
        }
    } catch( exception& e ) {
        cerr << "CacheItem::cacheStore() failed: " << e.what() << endl;
    } catch( ... ) {
        cerr << "CacheItem::cacheStore() failed for unknown reasons." << endl;
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
        