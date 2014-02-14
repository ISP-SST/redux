#ifndef REDUX_NETWORK_TCPCONNECTION_HPP
#define REDUX_NETWORK_TCPCONNECTION_HPP

#include "redux/network/protocol.hpp"

#include <boost/asio.hpp>
#include <boost/bind.hpp>

using boost::asio::ip::tcp;

namespace redux {

    namespace network {

        class TcpConnection : public std::enable_shared_from_this<TcpConnection> {

            static void writeCallback( size_t, const boost::system::error_code&, size_t );
            
        public:
            
            typedef std::shared_ptr<TcpConnection> ptr;
            typedef std::function<void( ptr )> callback;

            ~TcpConnection( void );


            static ptr newPtr( boost::asio::io_service& io_service ) {
                return ptr( new TcpConnection( io_service ) );
            }

            template <class T>
            void writeAndCheck( const T& data, size_t sz ) {
                if( mySocket.is_open() ) {
                    boost::asio::async_write( mySocket, boost::asio::buffer( data, sz ), boost::bind( &writeCallback, sz, boost::asio::placeholders::error,
                                              boost::asio::placeholders::bytes_transferred ) );
                }
            }
            tcp::socket& socket() { return mySocket; }

            void connect( std::string host, std::string service );
            void setCallback( callback cb = nullptr ) { activityCallback = cb; };
            void idle( void );
            void onActivity( ptr conn, const boost::system::error_code& error );

            TcpConnection& operator<<( const Command& );
            TcpConnection& operator>>( Command& );

        private:
            TcpConnection( boost::asio::io_service& io_service )
                : mySocket( io_service ), myService( io_service ), activityCallback( nullptr ) {

            }

            tcp::socket mySocket;
            boost::asio::io_service& myService;
            callback activityCallback;

        };

    }   // network

}   // redux

#endif // REDUX_NETWORK_TCPCONNECTION_HPP
