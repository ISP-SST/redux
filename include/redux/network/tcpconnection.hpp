#ifndef REDUX_NETWORK_TCPCONNECTION_HPP
#define REDUX_NETWORK_TCPCONNECTION_HPP

#include "redux/network/protocol.hpp"
#include "redux/util/arrayutil.hpp"

#include <iostream>
#include <memory>
#include <mutex>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/bind/protect.hpp>

using boost::asio::ip::tcp;

namespace redux {

    namespace network {

       class TcpConnection : public std::enable_shared_from_this<TcpConnection> {

            static void writeCallback( size_t sent, const boost::system::error_code& error, size_t transferred ) {
                using namespace boost::asio;

                if( !error ) {
                    if( sent != transferred ) {
                        std::ostringstream ss;
                        std::cerr << "TcpConnection::write: only " << transferred << "/" << sent << " bytes were successfully transferred.";
                        //throw std::ios_base::failure( ss.str() );
                    }
                }
//                 else {
//                     if( ( error == error::eof ) || ( error == error::connection_reset ) ) {
//                         // TODO handle reconnects...
//                     }
//                     else {
//                         throw std::ios_base::failure( "TcpConnection::async_write: error: " + error.message() );
//                     }
//                 }
            }

            
        public:

            template <typename T>
            void asyncWrite( const std::shared_ptr<T>& data, size_t sz ) {
                if( mySocket.is_open() ) {
                    boost::asio::async_write( mySocket, boost::asio::buffer( data.get(), sz ),
                                          //    strand.wrap(
                                              boost::bind( &TcpConnection::writeCallback, sz,
                                                           boost::asio::placeholders::error,
                                                           boost::asio::placeholders::bytes_transferred )
                                          //         )
                                            );

                }
            }

            template <typename T>
            void syncWrite( const T* data, size_t sz ) {
                if( mySocket.is_open() ) {
                    boost::asio::write( mySocket, boost::asio::buffer( data, sz ) );
                }
            }

            typedef std::shared_ptr<TcpConnection> Ptr;
            typedef std::function<void(Ptr)> callback;

            ~TcpConnection( void );


            static Ptr newPtr( boost::asio::io_service& io_service ) {
                return Ptr( new TcpConnection( io_service ) );
            }

            std::shared_ptr<char> receiveBlock( uint64_t& blockSize );
            
            template <class T>
            void asyncWrite( const std::vector<T>& data ) {
                size_t sz = data.size();
                if (!sz) return;
                auto tmp = redux::util::sharedArray<T>( sz );
                memcpy(tmp.get(),data.data(),sz*sizeof(T));
                asyncWrite(tmp, sz*sizeof(T));
            }

            template <class T>
            void asyncWrite( const T& data ) {
                asyncWrite(std::make_shared<T>(data), sizeof(T) );
            }

            template <class T>
            void syncWrite( const std::vector<T>& in ) {
                size_t sz = in.size();
                if (!sz) return;
                syncWrite( in.data(), sz*sizeof(T) );
            }

            template <class T>
            void syncWrite( const T& data ) {
                syncWrite( &data, sizeof(T) );
            }

            tcp::socket& socket() { return mySocket; }
            operator bool() const { return mySocket.is_open(); };

            void connect( std::string host, std::string service );
            void setCallback( callback cb = nullptr ) { activityCallback = cb; };
            void idle( void );
            void onActivity( Ptr conn, const boost::system::error_code& error );
            void setSwapEndian(bool se) { swapEndian = se; };
            bool getSwapEndian(void) { return swapEndian; };

            TcpConnection& operator<<( const Command& );
            TcpConnection& operator>>( Command& );

       private:
            TcpConnection( boost::asio::io_service& io_service );
            TcpConnection( const TcpConnection& ) = delete;

            callback activityCallback;
            tcp::socket mySocket;
            boost::asio::io_service& myService;
            bool swapEndian;

       public:
            boost::asio::strand strand;

        };
        
    }   // network

}   // redux

#endif // REDUX_NETWORK_TCPCONNECTION_HPP
