#include "redux/util/cache.hpp"


#include <iostream>
#include <memory>


using namespace redux::util;
using namespace std;

string Cache::cachePath = "/scratch/tomas/redux/";


Cache& Cache::get(void) {
    static Cache cache;
    return cache;
}


string Cache::path(void) {
    unique_lock<mutex> lock(cacheMutex);
    return cachePath;
}


void Cache::setPath(const string& path) {
    unique_lock<mutex> lock(cacheMutex);
    cachePath = path;
}

CacheItem::CacheItem(void) : itemPath(""), isLoaded(true) {
}

CacheItem::CacheItem(const string& path) : itemPath(path), isLoaded(true) {
    unique_lock<mutex> lock(itemMutex);
    if ( bfs::exists(itemPath) ) {
        bfs::remove(itemPath);      // HACK  for now we just delete pre-existing file, later we want to recover after crash...
    }
}

CacheItem::CacheItem(const CacheItem& rhs) : itemPath(rhs.itemPath), timeout(rhs.timeout), lastAccessed(rhs.lastAccessed),
 options(rhs.options), isLoaded(rhs.isLoaded) {
}


void CacheItem::cacheTouch(void) {
    unique_lock<mutex> lock(itemMutex);
    lastAccessed = 0;
}


bool CacheItem::cacheLoad(bool removeAfterLoad) {

    bool ret(false);
    if(itemPath.empty()) return ret;
    if(!isLoaded) {
        unique_lock<mutex> lock(itemMutex);
        if ( bfs::exists(itemPath) ) {
            ifstream in(itemPath.c_str(), ofstream::binary|ios_base::ate);
            if( in.good() ) {
                size_t sz = in.tellg();
                unique_ptr<char[]> buf( new char[sz] );
                in.seekg(0, ios_base::beg);
                in.clear();
                in.read(buf.get(), sz);
                size_t psz = cunpack(buf.get(),false);
                if(psz == sz) {
                    isLoaded = true;
                    //cout << hexString(this) << "  CacheItem::cacheLoad() " << itemPath << "   OK  removed=" << removeAfterLoad << endl;
                    ret = true;
                } else cout << hexString(this) << "  CacheItem::cacheLoad() " << psz << " != " << sz << " !!!" << endl;
            } else cout << hexString(this) << "  CacheItem::cacheLoad() in.good() == false !!!" << endl;
        } //else cout << hexString(this) << "  CacheItem::cacheLoad() " << itemPath << " doesn't exist !!!" << endl;
    } //else cout << hexString(this) << "  CacheItem::cacheLoad()  " << itemPath << "   already Loaded !!!" << endl;
    
    if(isLoaded && removeAfterLoad) {
        bfs::remove(itemPath);
    }
    
    return ret;
}


bool CacheItem::cacheStore(bool clearAfterStore){

    bool ret(false);
    if(itemPath.empty()) return ret;
    if(isLoaded) {
        unique_lock<mutex> lock(itemMutex);
        bfs::create_directories(itemPath.parent_path());
        size_t sz = csize();
        unique_ptr<char[]> buf( new char[sz] );
        size_t psz = cpack(buf.get());
        if(psz == sz) {
            ofstream out(itemPath.c_str(), ofstream::binary);
            if( out.good() ) {
                out.write( buf.get(), sz );
                //cout << hexString(this) << "  CacheItem::cacheStore() " << itemPath << "   OK  cleared=" << clearAfterStore << endl;
                ret = true;
            } else cout << hexString(this) << "  CacheItem::cacheStore() out.good() == false !!!" << endl;
        } else cout << hexString(this) << "  CacheItem::cacheStore() " << psz << " != " << sz << " !!!" << endl;
    } //else cout << hexString(this) << "  CacheItem::cacheStore()  " << itemPath << "   not Loaded !!!" << endl;
    
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
    Cache& c = Cache::get();
    itemPath = bfs::path(c.path()) / bfs::path(path);
    
}
        