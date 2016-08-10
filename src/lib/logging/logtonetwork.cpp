#include "redux/logging/logtonetwork.hpp"

#include "redux/util/datautil.hpp"

using namespace redux::logging;
using namespace redux::network;
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
    
    if( itemQueue.empty() ) return;
    
    vector<LogItemPtr> sendBuffer( itemQueue.begin(), itemQueue.end() );
    itemQueue.clear();
    
    lock.unlock();
    uint64_t blockSize(0);
    for( LogItemPtr& it: sendBuffer ) {
        blockSize += it->size();
    }

    Command cmd = CMD_PUT_LOG;
    size_t totalSize = blockSize + sizeof( uint64_t ) + 1;        // + blocksize + CMD_PUT_LOG
    shared_ptr<char> data( new char[totalSize], []( char* p ){ delete[] p; } );
    char* ptr = data.get();
    uint64_t count = pack( ptr, cmd );
    count += pack( ptr+count, blockSize );
    for( LogItemPtr& it: sendBuffer ) {
        count += it->pack( ptr+count );
    }

    conn->lock();
    try {
        conn->asyncWrite( data, totalSize );
        *conn >> cmd;       // ok or err
        if( cmd == CMD_OK ) {
            conn->unlock();
            return;
        }
    } catch( ... ) {
        // TODO Narrower catch, log error message
    }
    conn->unlock();
    
    lock.lock();
    for( auto it=sendBuffer.rbegin(); it != sendBuffer.rend(); ++it ) {
        itemQueue.push_front(*it);
    }

}

