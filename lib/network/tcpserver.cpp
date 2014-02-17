#include "redux/network/tcpserver.hpp"

#include "redux/logger.hpp"
#include "redux/util/stringutil.hpp"

using namespace redux::network;
using namespace redux::util;
using namespace std;

#define lg Logger::lg
namespace {
    const std::string thisChannel = "net";
}

TcpServer::TcpServer( boost::asio::io_service& io_service, uint16_t port )
    : acceptor( io_service, tcp::endpoint( tcp::v4(), port ) ), onConnected(nullptr) {
    acceptor.set_option(boost::asio::ip::tcp::acceptor::reuse_address(true));
    LOG_TRACE << "Binding to port: " << port;
}

void TcpServer::accept(void) {
    TcpConnection::Ptr nextConnection =
        TcpConnection::newPtr( acceptor.get_io_service() );

    acceptor.async_accept( nextConnection->socket(),
                           boost::bind( &TcpServer::onAccept, this, nextConnection,
                                        boost::asio::placeholders::error ) );
}

void TcpServer::onAccept( TcpConnection::Ptr conn,
                               const boost::system::error_code& error ) {
    accept();    // always start another accept
    if( !error ) {
        if( onConnected ) {
            LOG_DETAIL << "Accepted connection from \"" << conn->socket().remote_endpoint().address().to_string() << "\"";
            Command cmd;
            *conn >> cmd;
            if( cmd != CMD_CONNECT ) return;    // The connection is terminated when going out of scope.
            
            if( false ) {                       // TODO authentication. (key exchange ?)
                *conn << CMD_AUTH;
                // if auth fails => return, else continue and do the callback
            }

            onConnected(conn);
        }
    } else {
        LOG_ERR << "TcpServer::onAccept():  asio reports error:" << error.message();    // TODO some intelligent error handling
    }

}

