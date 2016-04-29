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


CacheItem::CacheItem(void) : itemPath(""), isLoaded(true) {
    
}


CacheItem::CacheItem(const string& path) : itemPath(""), fullPath(""), isLoaded(true) {
    setPath(path);
    unique_lock<mutex> lock(itemMutex);
    if ( bfs::exists(fullPath) ) {
        bfs::remove(fullPath);      // HACK  for now we just delete pre-existing file, later we want to recover after crash...
    }
}


CacheItem::CacheItem(const CacheItem& rhs) : itemPath(rhs.itemPath), fullPath(rhs.fullPath), timeout(rhs.timeout), lastAccessed(rhs.lastAccessed),
     options(rhs.options), isLoaded(rhs.isLoaded) {
         
}


CacheItem::~CacheItem() {
    
}


void CacheItem::cacheRemove(void) {
    unique_lock<mutex> lock(itemMutex);
    bfs::remove(fullPath);
}


void CacheItem::cacheTouch(void) {
    unique_lock<mutex> lock(itemMutex);
    lastAccessed = 0;
}


bool CacheItem::cacheLoad(bool removeAfterLoad) {

    bool ret(false);
    if(fullPath.empty()) return ret;
    if(!isLoaded) {
        unique_lock<mutex> lock(itemMutex);
        if ( bfs::exists(fullPath) ) {
            ifstream in(fullPath.c_str(), ofstream::binary|ios_base::ate);
            if( in.good() ) {
                size_t sz = in.tellg();
                unique_ptr<char[]> buf( new char[sz] );
                in.seekg(0, ios_base::beg);
                in.clear();
                in.read(buf.get(), sz);
                size_t psz = cunpack(buf.get(),false);
                if( psz >= sz ) {
                    isLoaded = true;
                    //cout << hexString(this) << "  CacheItem::cacheLoad() " << fullPath << "   OK  remove=" << removeAfterLoad << endl;
                    ret = true;
                } else cout << hexString(this) << "  CacheItem::cacheLoad() " << psz << " < " << sz << " !!!" << endl;
            } else cout << hexString(this) << "  CacheItem::cacheLoad() in.good() == false !!!" << endl;
        } //else cout << hexString(this) << "  CacheItem::cacheLoad() " << fullPath << " doesn't exist !!!" << endl;
    } //else cout << hexString(this) << "  CacheItem::cacheLoad()  " << fullPath << "   already Loaded !!!" << endl;
    
    if(isLoaded && removeAfterLoad) {
        unique_lock<mutex> lock(itemMutex);
        bfs::remove(fullPath);
    }
    
    return ret;
}


bool CacheItem::cacheStore(bool clearAfterStore){

    bool ret(false);
    if(fullPath.empty()) return ret;
    if(isLoaded) {
        unique_lock<mutex> lock(itemMutex);
        bfs::create_directories(fullPath.parent_path());
        cachedSize = csize();
        unique_ptr<char[]> buf( new char[cachedSize] );
        size_t psz = cpack(buf.get());
        if( psz <= cachedSize ) {
            ofstream out(fullPath.c_str(), ofstream::binary);
            if( out.good() ) {
                out.write( buf.get(), cachedSize );
                //cout << hexString(this) << "  CacheItem::cacheStore() " << fullPath << "   OK  clear=" << clearAfterStore << endl;
                ret = true;
            } else cout << hexString(this) << "  CacheItem::cacheStore() out.good() == false !!!" << endl;
        } else cout << hexString(this) << "  CacheItem::cacheStore() " << psz << " > " << cachedSize << " !!!" << endl;
    } //else cout << hexString(this) << "  CacheItem::cacheStore()  " << fullPath << "   not Loaded !!!" << endl;
    
    if(clearAfterStore) {
        cclear();
        isLoaded = false;
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
        