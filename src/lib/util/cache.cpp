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

