#include "redux/util/cache.hpp"

#include <redux/file/fileio.hpp>

#include <iostream>
#include <memory>


using namespace redux::file;
using namespace redux::util;
using namespace std;


Cache::Info::Info( size_t i, const string& n, size_t s ) : id1(i), id2(0), name(n), size(s), count(0), totalSize(0) {

}

Cache::Info::Info( size_t i1, size_t i2, const string& n, size_t s ) : id1(i1), id2(i2), name(n), size(s), count(0), totalSize(0) {

}


std::string Cache::Info::getInfo(void) const {
    return alignLeft(to_string(count),12) + alignLeft(to_string(count*size),14) + name;
    
/*              map:   return  std::string(typeid(KeyT).name()) + " -> " + std::string(typeid(T).name())
                    + "   count: " + std::to_string(m.second.size());
                set:   return std::string(typeid(T).name()) + "   count: " + std::to_string(s.second.size());
*/
}


void Cache::cleanup(void) {
    lock_guard<mutex> lock(get().mtx);
    for( auto& c: get().caches ) {
        c.clear();
    }
}

std::string Cache::getStats(void) {
    string ret;
    static const string hdr = "Items       Size          Type\n";
    lock_guard<mutex> lock(get().mtx);
    for( auto& c: get().caches ) {
        ret += c.getInfo() + "\n";
    }
    if( ret.empty() ) return ret;
    return hdr+ret;
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

