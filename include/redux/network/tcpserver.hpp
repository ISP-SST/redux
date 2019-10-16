#ifndef REDUX_NETWORK_TCPSERVER_HPP
#define REDUX_NETWORK_TCPSERVER_HPP

#include "redux/network/tcpconnection.hpp"
#include "redux/network/host.hpp"
#include "redux/util/datautil.hpp"

#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <boost/noncopyable.hpp>


namespace redux {

    namespace network {


        class TcpServer : private boost::noncopyable {

        public:

            TcpServer( uint16_t port, uint16_t threads );
            ~TcpServer();

            void accept(void);
            void setCallback( TcpConnection::callback cb = nullptr ) { onConnected = cb; };
            void start(void);
            void start( uint16_t port );
            void stop(void);
            unsigned short port(void);
            void cleanup(void);
            void addThread(uint16_t n=1);
            void delThread(uint16_t n=1);
            void addConnection( const Host::HostInfo&, TcpConnection::Ptr );
            void removeConnection( TcpConnection::Ptr );
            void releaseConnection( TcpConnection::Ptr );
            Host::Ptr getHost( const TcpConnection::Ptr ) const;
            TcpConnection::Ptr getConnection( Host::Ptr );
            size_t size( void ) const;
            size_t nThreads( void ) const;
            
        private:

            void onAccept( TcpConnection::Ptr conn, const boost::system::error_code& error );
            void threadLoop( void );

            std::map<TcpConnection::Ptr, Host::Ptr, redux::util::PtrCompare<TcpConnection>> connections;
            mutable std::mutex mtx;
            boost::asio::io_service ioService;
            std::shared_ptr<boost::asio::io_service::work> workLoop;
            tcp::acceptor acceptor;
            tcp::endpoint endpoint;
            TcpConnection::callback onConnected;
            boost::thread_group pool;
            uint16_t minThreads;
            std::atomic<uint16_t> nThreads_, nConnections;
            bool running,do_handshake,do_auth;

        };

    }   // network

}   // redux

#endif // REDUX_NETWORK_TCPSERVER_HPP
