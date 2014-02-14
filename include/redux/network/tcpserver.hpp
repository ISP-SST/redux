#ifndef REDUX_NETWORK_TCPSERVER_HPP
#define REDUX_NETWORK_TCPSERVER_HPP

#include "redux/network/tcpconnection.hpp"

#include <boost/asio.hpp>

using boost::asio::ip::tcp;

namespace redux {

    namespace network {


        class TcpServer : private boost::noncopyable {

        public:

            TcpServer( boost::asio::io_service& io_service, uint16_t port );

            void accept(void);
            void setCallback( TcpConnection::callback cb = nullptr ) { onConnected = cb; };

        private:

            void onAccept( TcpConnection::ptr conn, const boost::system::error_code& error );

            tcp::acceptor acceptor;
            TcpConnection::callback onConnected;

        };

    }   // network

}   // redux

#endif // REDUX_NETWORK_TCPSERVER_HPP
