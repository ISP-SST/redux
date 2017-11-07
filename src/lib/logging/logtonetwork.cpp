#include "redux/logging/logtonetwork.hpp"

#include "redux/util/datautil.hpp"

using namespace redux::logging;
using namespace redux::network;
using namespace redux::util;
using namespace std;


LogToNetwork::LogToNetwork( boost::asio::io_service& s, const network::Host::Ptr& h, uint32_t i, uint8_t m, unsigned int flushPeriod)
    : LogOutput(m,flushPeriod), service(s), host(h), id(i) {


}


LogToNetwork::~LogToNetwork() {
    
    if( !itemQueue.empty() ) {
        flushBuffer();
    }

}


void LogToNetwork::connect(void) {

    if ( !conn ) conn = TcpConnection::newPtr( service );
    
    if ( conn ) {
        if( !conn->socket().is_open() && host ) {
            try {
                conn->connect( host->info.connectName, to_string(host->info.connectPort) );
                if( conn->socket().is_open() ) {
                    Command cmd;
                    Host::HostInfo myInfo = Host::myInfo().info;
                    myInfo.peerType &= ~Host::TP_WORKER;
                    *conn << CMD_CONNECT;
                    *conn >> cmd;
                    if( cmd == CMD_AUTH ) {
                        // implement
                    }
                    if( cmd == CMD_CFG ) {  // handshake requested
                        *conn << myInfo;
                        *conn >> host->info;
                        *conn >> cmd;       // ok or err
                    }
                    if( cmd == CMD_OK ) {
                        *conn << CMD_LOG_CONNECT << id;
                        *conn >> cmd;
                    }
                    if( cmd != CMD_OK ) throw exception();
                } else throw exception();
            } catch ( ... ) {
                conn->socket().close();
            }
        }
        try {
            auto test RDX_UNUSED = conn->socket().remote_endpoint();  // check if endpoint exists, will throw if not connected.
        } catch ( ... ) {
            conn->socket().close();
        }
    }
}


void LogToNetwork::flushBuffer( void ) {

    connect();
    if( !conn || !conn->socket().is_open() ) return;
    
    unique_lock<mutex> lock( queueMutex );
    if( itemQueue.empty() ) return;
    vector<LogItemPtr> sendBuffer( itemQueue.begin(), itemQueue.end() );
    itemQueue.clear();
    itemCount = 0;
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
            sendBuffer.clear();
        }
    } catch( ... ) {
        // TODO Narrower catch, log error message
    }
    conn->unlock();
    
    lock.lock();
    for( auto it=sendBuffer.rbegin(); it != sendBuffer.rend(); ++it ) {
        itemQueue.push_front(*it);
    }

    itemCount = itemQueue.size();

}

