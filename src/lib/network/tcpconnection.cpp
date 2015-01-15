#include "redux/network/tcpconnection.hpp"

#include "redux/logger.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

namespace ba = boost::asio;

using namespace redux::util;
using namespace redux::network;
using namespace std;


#ifdef DEBUG_
#define DBG_NET_
#endif

#define lg Logger::mlg
namespace {
    const std::string thisChannel = "net";

#ifdef DBG_NET_
    static atomic<int> connCounter(0);
#endif

}


TcpConnection::TcpConnection( boost::asio::io_service& io_service )
    : activityCallback( nullptr ), mySocket( io_service ), myService( io_service ), swapEndian(false), strand( myService ) {
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


shared_ptr<char> TcpConnection::receiveBlock( uint64_t& blockSize ) {

    using redux::util::unpack;
    shared_ptr<char> buf;
    
    try {
        char sz[sizeof( uint64_t )];
        uint64_t count = boost::asio::read( socket(), boost::asio::buffer( sz, sizeof( uint64_t ) ) );

        if( count != sizeof( uint64_t ) ) {
            throw ios_base::failure( "blockSize: Only received " + to_string( count ) + "/" + to_string( sizeof( uint64_t ) ) + " bytes." );
        }
        
        unpack( sz, blockSize, swapEndian );

        if(blockSize > 0) {
            buf.reset( new char[ blockSize ], []( char * p ) { delete[] p; } );
            count = boost::asio::read( socket(), boost::asio::buffer( buf.get(), blockSize ) );

            if( count != blockSize ) {
                throw ios_base::failure( "Only received " + to_string( count ) + "/" + to_string( blockSize ) + " bytes." );
            }
        }
    }
    catch( const exception& e ) {
        cerr << "Failed to receive datablock: " << e.what() << endl;
        //LOG_ERR << "Failed to receive datablock: " << e.what();
        blockSize = 0;
        buf.reset();
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
        mySocket.async_read_some( ba::null_buffers(), strand.wrap(boost::bind( &TcpConnection::onActivity, this, shared_from_this(),
                                                                   ba::placeholders::error ) ) );

    }
    
}

void TcpConnection::onActivity( Ptr connptr, const boost::system::error_code& error ) {

    if( !error ) {
        if( activityCallback ) {
            //LOG_DEBUG << "Activity on connection \"" << connptr->socket().remote_endpoint().address().to_string() << "\"";
            activityCallback( connptr );
        }
    } else {
        if( ( error == ba::error::eof ) || ( error == ba::error::connection_reset ) ) {
            mySocket.close();
        } else {
            throw std::ios_base::failure( "TcpConnection::onActivity: error: " + error.message() );
        }
    }

}


TcpConnection& TcpConnection::operator<<( const Command& out ) {
    writeAndCheck(out);
    return *this;
}


TcpConnection& TcpConnection::operator>>( Command& in ) {
    if( ba::read( mySocket, ba::buffer( &in, sizeof( Command ) ) ) < sizeof( Command ) ) {
        in = CMD_ERR;
        throw std::ios_base::failure( "Failed to receive command." );
    }
    return *this;
}
