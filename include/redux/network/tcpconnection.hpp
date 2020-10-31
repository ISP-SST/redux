#ifndef REDUX_NETWORK_TCPCONNECTION_HPP
#define REDUX_NETWORK_TCPCONNECTION_HPP

#include "redux/network/protocol.hpp"
#include "redux/util/arrayutil.hpp"

#include "redux/util/trace.hpp"

#include <iostream>
#include <memory>
#include <mutex>
#include <typeinfo>

#ifndef Q_MOC_RUN
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/bind/protect.hpp>
#endif

using boost::asio::ip::tcp;

namespace redux {

    namespace network {
        class TcpServer;
        class TcpConnection : public std::enable_shared_from_this<TcpConnection>
#ifdef RDX_TRACE_NET
            ,public redux::util::TraceObject<TcpConnection>
#endif
       {

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
            void asyncWrite( const std::shared_ptr<T> data, size_t sz ) {
                if( mySocket.is_open() ) {
                    boost::asio::async_write( mySocket, boost::asio::buffer( data.get(), sz ),
                                              boost::bind( &TcpConnection::writeCallback, sz,
                                                           boost::asio::placeholders::error,
                                                           boost::asio::placeholders::bytes_transferred )
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

            uint64_t receiveN( std::shared_ptr<char> buf, uint64_t N );
            std::shared_ptr<char> receiveBlock( uint64_t& blockSize );
            
            size_t readline( std::string& line );
            void writeline( const std::string& line );
            template <class T>
            void asyncWrite( const std::vector<T>& data ) {
                size_t sz = data.size();
                if (!sz) return;
                std::shared_ptr<T> tmp = redux::util::rdx_get_shared<T>(sz);
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
            void close( void );
            callback getCallback( void ) { std::lock_guard<std::mutex> g(mtx); return activityCallback; };
            bool hasUrgentCallback( void ) const { return (urgentCallback != nullptr); };
            void setCallback( callback cb = nullptr ) { std::lock_guard<std::mutex> g(mtx); activityCallback = cb; };
            void setUrgentCallback( callback cb = nullptr ) { std::lock_guard<std::mutex> g(mtx); urgentCallback = cb; };
            void setErrorCallback( callback cb = nullptr ) { std::lock_guard<std::mutex> g(mtx); errorCallback = cb; };
            void uIdle( void );
            void urgentHandler( const boost::system::error_code& error, size_t transferred );
            void idle( void );
            void onActivity( const boost::system::error_code& error, size_t transferred );
            void setSwapEndian(bool se) { swapEndian_ = se; };
            bool getSwapEndian(void) { return swapEndian_; };
            
            void lock(void) { mtx.lock(); };
            void unlock(void) { mtx.unlock(); };
            bool try_lock(void) { return mtx.try_lock(); };
            
            void sendUrgent( uint8_t c );
            void receiveUrgent( uint8_t& c );
            uint8_t getUrgentData(void) { return urgentData; };
            
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
                if( !mySocket.is_open() || boost::asio::read( mySocket, boost::asio::buffer( &out, sizeof(T) ) ) < sizeof(T) ) {
                    out = T();
                    throw std::ios_base::failure( std::string("TcpConnection: Failed to receive ")+typeid(T).name() );
                }
                return *this;
            }

       private:
            TcpConnection( boost::asio::io_service& io_service );
            TcpConnection( const TcpConnection& ) = delete;

            callback activityCallback;
            callback urgentCallback;
            callback errorCallback;
            tcp::socket mySocket;
            boost::asio::io_service& myService;
            bool swapEndian_;
            uint8_t urgentData;
            bool urgentActive;
            uint64_t id;
            std::mutex mtx;
            
            friend class TcpServer;

       public:

        };
        
    }   // network

}   // redux

#endif // REDUX_NETWORK_TCPCONNECTION_HPP
