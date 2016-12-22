#ifndef REDUX_NETWORK_TCPCONNECTION_HPP
#define REDUX_NETWORK_TCPCONNECTION_HPP

#include "redux/network/protocol.hpp"
#include "redux/util/arrayutil.hpp"

#include <iostream>
#include <memory>
#include <mutex>
#include <typeinfo>

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
                        std::cerr << "TcpConnection::write: only " << transferred << "/" << sent << " bytes were successfully transferred.";
                    }
                }
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
                std::shared_ptr<T> tmp( new T[sz], [](T* p){ delete[] p;} );
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
            callback getCallback( void ) { std::unique_lock<std::mutex> lock(mtx); return activityCallback; };
            void setCallback( callback cb = nullptr ) { std::unique_lock<std::mutex> lock(mtx); activityCallback = cb; };
            void setErrorCallback( callback cb = nullptr ) { std::unique_lock<std::mutex> lock(mtx); errorCallback = cb; };
            void idle( void );
            void onActivity( Ptr conn, const boost::system::error_code& error );
            void setSwapEndian(bool se) { swapEndian_ = se; };
            bool getSwapEndian(void) { return swapEndian_; };
            
            void lock(void) { mtx.lock(); };
            void unlock(void) { mtx.unlock(); };
            bool try_lock(void) { return mtx.try_lock(); };
            
            void sendUrgent( uint8_t c );
            void receiveUrgent( uint8_t& c );
            
            TcpConnection& operator<<( const uint8_t& );
            TcpConnection& operator>>( uint8_t& );
            TcpConnection& operator<<( const std::vector<std::string>& );
            TcpConnection& operator>>( std::vector<std::string>& );

            template <typename T>
            TcpConnection& operator<<( const T& in ) {
                syncWrite(&in, sizeof(T));
                return *this;
            }


            template <typename T>
            TcpConnection& operator>>( T& out ) {
                if( boost::asio::read( mySocket, boost::asio::buffer( &out, sizeof(T) ) ) < sizeof(T) ) {
                    out = T();
                    throw std::ios_base::failure( std::string("TcpConnection: Failed to receive ")+typeid(T).name() );
                }
                return *this;
            }

       private:
            TcpConnection( boost::asio::io_service& io_service );
            TcpConnection( const TcpConnection& ) = delete;

            callback activityCallback;
            callback errorCallback;
            tcp::socket mySocket;
            boost::asio::io_service& myService;
            bool swapEndian_;
            std::mutex mtx;

       public:
            boost::asio::strand strand;

        };
        
    }   // network

}   // redux

#endif // REDUX_NETWORK_TCPCONNECTION_HPP
