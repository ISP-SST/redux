#include "redux/logging/logtonetwork.hpp"

#include "redux/util/datautil.hpp"

using namespace redux::logging;
using namespace redux::util;
using namespace std;


LogToNetwork::LogToNetwork( const network::TcpConnection::Ptr& c, uint32_t id, uint8_t m, unsigned int flushPeriod) : LogOutput(m,flushPeriod), conn(c) {


}


LogToNetwork::~LogToNetwork() {
    
    if( !itemQueue.empty() ) {
        flushBuffer();
    }

}



void LogToNetwork::flushBuffer( void ) {

    unique_lock<mutex> lock( queueMutex );
    vector<LogItemPtr> sendBuffer( itemQueue.begin(), itemQueue.end() );

    itemQueue.clear();
    if( sendBuffer.empty() ) return;
    lock.unlock();

    uint64_t blockSize(0);
    for( LogItemPtr& it: sendBuffer ) {
        blockSize += it->size();
    }
    size_t totalSize = blockSize + sizeof( uint64_t );        // + blocksize
    shared_ptr<char> data( new char[totalSize], []( char* p ){ delete[] p; } );
    char* ptr = data.get();
    uint64_t count = pack( ptr, blockSize );
    for( LogItemPtr& it: sendBuffer ) {
        count += it->pack( ptr+count );
    }

    conn->lock();
    conn->asyncWrite( data, totalSize );
    conn->unlock();

}

