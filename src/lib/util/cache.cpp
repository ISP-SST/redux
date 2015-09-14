#include "redux/util/cache.hpp"


#include <iostream>
#include <memory>


using namespace redux::util;
using namespace std;

std::string Cache::cachePath = "/tmp/redux/";


Cache& Cache::get(void) {
    static Cache cache;
    return cache;
}


std::string Cache::path(void) {
    unique_lock<mutex> lock(cacheMutex);
    return cachePath;
}


void Cache::setPath(const std::string& path) {
    unique_lock<mutex> lock(cacheMutex);
    cachePath = path;
}


CacheItem::CacheItem(const CacheItem& rhs) : itemPath(rhs.itemPath), timeout(rhs.timeout), lastAccessed(rhs.lastAccessed),
 options(rhs.options), isLoaded(rhs.isLoaded) {
    
}


void CacheItem::cacheTouch(void) {
    lastAccessed = 0;
}


bool CacheItem::cacheLoad(bool removeAfterLoad) {
    if(!isLoaded) {
        std::unique_lock<std::mutex> lock(itemMutex);
        if ( bfs::exists(itemPath) ) {
            std::ifstream in(itemPath.c_str(), std::ofstream::binary|std::ios_base::ate);
            if( in.good() ) {
                size_t sz = in.tellg();
                std::unique_ptr<char[]> buf( new char[sz] );
                in.seekg(0, std::ios_base::beg);
                in.clear();
                in.read(buf.get(), sz);
                size_t psz = cunpack(buf.get(),false);
                if(psz == sz) {
                    isLoaded = true;
                    if(removeAfterLoad) {
                        in.close();
                        bfs::remove(itemPath);
                    }
                    //std::cout << "CacheItem::cacheLoad() " << itemPath << "   OK  removed=" << removeAfterLoad << std::endl;
                    return true;
                } else std::cout << "CacheItem::cacheLoad() " << psz << " != " << sz << " !!!" << std::endl;
            } else std::cout << "CacheItem::cacheLoad() in.good() == false !!!" << std::endl;
        } //else std::cout << "CacheItem::cacheLoad() " << itemPath << " doesn't exist !!!" << std::endl;
    } //else std::cout << "CacheItem::cacheLoad()  " << itemPath << "   already Loaded !!!" << std::endl;
    return false;
}


bool CacheItem::cacheStore(bool clearAfterStore){
    if(isLoaded) {
        std::unique_lock<std::mutex> lock(itemMutex);
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
                //std::cout << "CacheItem::cacheStore() " << itemPath << "   OK  cleared=" << clearAfterStore << std::endl;
                return true;
            } else std::cout << "CacheItem::cacheStore() out.good() == false !!!" << std::endl;
        } std::cout << "CacheItem::cacheStore() " << psz << " != " << sz << " !!!" << std::endl;
    } //else std::cout << "CacheItem::cacheStore()  " << itemPath << "   not Loaded !!!" << std::endl;
    return false;
}


void CacheItem::setPath(const std::string& path) {
    Cache& c = Cache::get();
    itemPath = bfs::path(c.path()) / bfs::path(path);
}
        