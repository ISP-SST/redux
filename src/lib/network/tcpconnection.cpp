#include "redux/network/tcpconnection.hpp"

#include "redux/logger.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <thread>

namespace ba = boost::asio;

using namespace redux::util;
using namespace redux::network;
using namespace std;


#ifdef DEBUG_
//#define DBG_NET_
#endif

#define lg Logger::mlg
namespace {
    const std::string thisChannel = "net";

#ifdef DBG_NET_
    static atomic<int> connCounter(0);
#endif

}


TcpConnection::TcpConnection( boost::asio::io_service& io_service )
    : activityCallback( nullptr ), mySocket( io_service ), myService( io_service ), swapEndian(false), strand(io_service) {
#ifdef DBG_NET_
    LOG_DEBUG << "Constructing TcpConnection: (" << hexString(this) << ") new instance count = " << (connCounter.fetch_add(1)+1);
#endif
}

            
TcpConnection::~TcpConnection( void ) {
    mySocket.close();
#ifdef DBG_NET_
    LOG_DEBUG << "Destructing TcpConnection: (" << hexString(this) << ") new instance count = " << (connCounter.fetch_sub(1)-1);
#endif
}


shared_ptr<char> TcpConnection::receiveBlock( uint64_t& received ) {

    using redux::util::unpack;
    shared_ptr<char> buf;
    uint64_t blockSize(0);
    received = 0;

    char sz[sizeof(uint64_t)];
    boost::system::error_code ec;
    uint64_t count = boost::asio::read( socket(), boost::asio::buffer(sz, sizeof(uint64_t)), ec );
    if( ec ) {
        //LOG_ERR << "TcpConnection::receiveBlock(" << hexString(this) << ")  failed to receive blockSize.  error: " + ec.message();
        //cerr << "TcpConnection::receiveBlock(" << hexString(this) << ")  Failed to receive blockSize.  error: " + ec.message() << endl;
        return buf;
    }

    unpack( sz, blockSize, swapEndian );
    if( blockSize == 0 ) return buf;

    buf.reset( new char[ blockSize ], []( char * p ) { delete[] p; } );

    int64_t remain = blockSize;
    while( remain > 0 ) {
        try {
            count = boost::asio::read( socket(), boost::asio::buffer(buf.get()+received,remain), boost::asio::transfer_at_least(1) );
            received += count;
            remain -= count;
        }
        catch( const exception& ) {
            // ignore and continue.  TODO: make it safer
        }
    }
    return buf;

}


void TcpConnection::connect( string host, string service ) {
    
    if( host == "" ) host = "localhost";

    ba::ip::tcp::resolver::query query( host, service );
    ba::ip::tcp::resolver resolver( myService );
    ba::ip::tcp::resolver::iterator destination = resolver.resolve( query );
    ba::ip::tcp::resolver::iterator end ;
    ba::ip::tcp::endpoint endpoint;

    while( destination != end ) {       // TODO: try connect to each destination and return on success...
        mySocket.connect( *destination++ );
        if( mySocket.is_open() ) return;
    }
    LOG_ERR << "Connection failed";
}

void TcpConnection::idle( void ) {

    if( !activityCallback ) {
        throw invalid_argument( "TcpConnection::idle()  Attempting to idle without callback !!" );
    }

    if( mySocket.is_open() ) {
        mySocket.async_read_some( ba::null_buffers(),
                                  boost::bind( &TcpConnection::onActivity, this, shared_from_this(), ba::placeholders::error) );

    }
    
}

void TcpConnection::onActivity( Ptr connptr, const boost::system::error_code& error ) {

    if( !error ) {
        if( activityCallback ) {
            //LOG_DEBUG << "Activity on connection \"" << connptr->socket().remote_endpoint().address().to_string() << "\"";
            if( mySocket.is_open() ) {
                std::thread(activityCallback,connptr).detach();
            }
        }
    } else {
        if( ( error == ba::error::eof ) || ( error == ba::error::connection_reset ) ) {
            mySocket.close();
        } else {
            throw std::ios_base::failure( "TcpConnection::onActivity: error: " + error.message() );
        }
    }

}


TcpConnection& TcpConnection::operator<<( const Command& in ) {
    syncWrite(&in, sizeof(Command));
    return *this;
}


TcpConnection& TcpConnection::operator>>( Command& out ) {
    if( ba::read( mySocket, ba::buffer( &out, sizeof( Command ) ) ) < sizeof( Command ) ) {
        out = CMD_ERR;
        throw std::ios_base::failure( "Failed to receive command." );
    }
    return *this;
}
