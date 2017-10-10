#include "redux/logging/logoutput.hpp"

//#include "redux/logging/logmsg.hpp"
//#include "redux/logging/log.hpp"

#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"
// #include "redux/socket.hpp"
// #include "redux/util/stringtools.hpp"
// #include "redux/threadpool.hpp"
//#include "redux/ServerSocket.hpp"

#include <algorithm>
#include <functional>
#include <thread>
//#include <stdio.h>
#include <syslog.h>
#include <string.h>

#include <boost/filesystem.hpp>

using namespace redux::logging;
using namespace redux::util;
using namespace redux;
using namespace std;

namespace bfs = boost::filesystem;


LogOutput::LogOutput( uint8_t m, unsigned int flushPeriod ) : flushPeriod(flushPeriod), async_(false), mask(m),
        itemCount(0) {

}


LogOutput::~LogOutput( )  {
    
}


void LogOutput::addItem( LogItemPtr item ) {
    
    unique_lock<mutex> lock(queueMutex);
    itemQueue.push_back(item);
    itemCount = itemQueue.size();
    if( itemCount >= flushPeriod ) {
        if(async_) {
            std::thread( std::bind(&LogOutput::flushBuffer, this) ).detach();
        } else {
            lock.unlock();
            this->flushBuffer();
        }
    }
    
}


void LogOutput::addItems( const vector<LogItemPtr>& items ) {
    
    if( items.empty() ) return;
    
    unique_lock<mutex> lock(queueMutex);
    std::move( std::begin(items), std::end(items), back_inserter(itemQueue) );
    itemCount = itemQueue.size();
    if( (flushPeriod == 0) || (itemCount > flushPeriod) ) {
        if(async_) {
            std::thread( std::bind(&LogOutput::flushBuffer, this) ).detach();
        } else {
            lock.unlock();
            this->flushBuffer();
        }
    }

}
