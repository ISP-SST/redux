#include "redux/network/tcpconnection.hpp"

#include "redux/logger.hpp"

namespace ba = boost::asio;

using namespace redux::network;
using namespace std;

#define lg Logger::lg
namespace {
    const std::string thisChannel = "net";

}

TcpConnection::~TcpConnection( void ) {
    mySocket.close();
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
