#include "redux/network/tcpconnection.hpp"

#include "redux/logging/logger.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <thread>

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

namespace ba = boost::asio;

using namespace redux::util;
using namespace redux::network;
using namespace std;

#ifdef DEBUG_
//#define DBG_NET_
#endif

namespace {

    set<uint64_t> connection_ids;
    mutex gmtx;

    uint64_t getID( void ) {
        lock_guard<mutex> lock( gmtx );
        uint64_t id(1);
        while( connection_ids.count( id ) ) id++;     // find first unused ID.
        connection_ids.insert( id );
        return id;
    }

    void freeID( uint64_t id ) {
        lock_guard<mutex> lock( gmtx );
        connection_ids.erase( id );
    }

    /*size_t idCount( void ) RDX_UNUSED {
        lock_guard<mutex> lock( gmtx );
        return connection_ids.size();
    }*/

    boost::asio::thread_pool& callbackPool( void ) {
        static boost::asio::thread_pool pool( std::thread::hardware_concurrency() );
        return pool;
    }

}


TcpConnection::TcpConnection( ba::io_context& ioc )
    : activityCallback( nullptr ), urgentCallback( nullptr ), errorCallback( nullptr ), mySocket( ioc ),
    ioContext( ioc ), swapEndian_(false), urgentData(0), urgentActive(false), id( getID() ) {
#ifdef DBG_NET_
    LOG_DEBUG << "Constructing TcpConnection: (" << hexString(this) << ")  ID: " << id << "/" << idCount() << ende;
#endif

}

            
TcpConnection::~TcpConnection( void ) {
    close();
#ifdef DBG_NET_
    LOG_DEBUG << "Destructing TcpConnection: (" << hexString(this) << ")  ID: " << id << "/" << idCount() << ende;
#endif
    freeID( id );
}


uint64_t TcpConnection::receiveN( std::shared_ptr<char> buf, uint64_t N ) {
    
    int64_t remain = N;
    uint64_t count(0);
    boost::system::error_code ec;
    while( remain > 0 ) {
        try {
            uint64_t this_count = boost::asio::read( mySocket, boost::asio::buffer(buf.get()+count,remain), boost::asio::transfer_at_least(1), ec );
            if( this_count == 0 ) break;
            count += this_count;
            remain -= this_count;
        }
        catch( const exception& ) {
            break;  // don't spin retrying after a genuine connection error
        }
    }
    
    if( remain ) {  // connection severed before completion, clear partial transfer.
        count = 0;
    }

    return count;
}


shared_ptr<char> TcpConnection::receiveBlock( uint64_t& received ) {

    using redux::util::unpack;
    shared_ptr<char> buf;
    uint64_t blockSize(0);
    received = 0;

    if( !mySocket.is_open() ) {
        return buf;
    }
    
    char sz[sizeof(uint64_t)];
    boost::system::error_code ec;
    uint64_t count;
    try {
        count = boost::asio::read( mySocket, boost::asio::buffer(sz, sizeof(uint64_t)), ec );
    } catch( const exception& ) {
        //LOG_ERR << "TcpConnection::receiveBlock(" << hexString(this) << ")  failed to receive blockSize.  error: " + ec.message();
        //cerr << "TcpConnection::receiveBlock(" << hexString(this) << ")  Failed to receive blockSize.  error: " + ec.message() << endl;
        return buf;
    }

    if( count < sizeof(uint64_t) ) return buf;
    unpack( sz, blockSize, swapEndian_ );
    if( blockSize == 0 ) return buf;

    buf = rdx_get_shared<char>(blockSize);
    received = receiveN( buf, blockSize );
    if( received != blockSize ) {  // connection severed before completion, clear partial transfer.
        buf.reset();
        received = 0;
    }
    
    return buf;

}


namespace {
    const string delimiter = "\r\n";
}

size_t TcpConnection::readline( string& line ) {

    line.clear();

    try {
        ba::streambuf streambuf;
        size_t bytes_transferred = ba::read_until( socket(), streambuf, delimiter );
        auto data_begin = buffers_begin( streambuf.data() );
        line = string( data_begin, data_begin + bytes_transferred - delimiter.size() );
    } catch( const exception& e ) {
        //cout << "TcpConnection::readline  e = " << e.what() << endl;
    }
    
    return !line.length();

}


void TcpConnection::writeline( const string& line ) {

    try {
        ba::write( mySocket, ba::buffer(line + delimiter) );
    } catch( const exception& e ) {
        //cout << "TcpConnection::writeline  e = " << e.what() << endl;
    }

}


void TcpConnection::connect(std::string host, std::string service) {

    if (host.empty()) host = "localhost";

    try {
        ba::ip::tcp::resolver resolver(ioContext);

        // Modern resolve: returns a range of endpoints
        auto results = resolver.resolve(host, service);

        for (const auto& entry : results) {
            try {
                mySocket.connect(entry.endpoint());
            } catch (const boost::system::system_error&) {
                mySocket.close();
            }

            if (mySocket.is_open()) return;
        }
    } catch (...) {
        // TODO: Handle error
    }

    close();

}


void TcpConnection::close( void ) {

    unique_lock<mutex> lock(mtx);
    activityCallback = nullptr;
    urgentCallback = nullptr;
    errorCallback = nullptr;
    if( mySocket.is_open() ) {
        boost::system::error_code error;
        mySocket.shutdown( ba::socket_base::shutdown_both, error );
        mySocket.close( error );
    }
}


void TcpConnection::uIdle( void ) {

    unique_lock<mutex> lock(mtx);
    if( !urgentCallback || urgentActive ) {
        return;
    }

    urgentActive = true;
    
    if( mySocket.is_open() ) {
        lock.unlock();
                
        mySocket.async_receive( ba::buffer( &urgentData, 1 ), ba::socket_base::message_out_of_band,
                                  boost::bind( &TcpConnection::urgentHandler, this,
                                               ba::placeholders::error, ba::placeholders::bytes_transferred) );
//         mySocket.async_read_some( ba::null_buffers(),
//                                   boost::bind( &TcpConnection::urgentHandler, this, ba::placeholders::error, ba::placeholders::bytes_transferred) );

    }

}


void TcpConnection::urgentHandler( const boost::system::error_code& ec, size_t transferred ) {
    
    boost::system::error_code error = ec;
    
    if( !error ) {
        auto test RDX_UNUSED = mySocket.remote_endpoint( error );
    }
    
    unique_lock<mutex> lock(mtx);
    try {
        if( !urgentActive ) return;     // prevent multiple uIdle
        urgentActive = false;
        if( mySocket.is_open() && !mySocket.at_mark() ) {
            lock.unlock();
            uIdle();
            return;
        }
    } catch(...) {
        return;
    }
    
    if( !error ) {
        if( urgentCallback ) {
            //LOG_DEBUG << "Activity on connection \"" << connptr->socket().remote_endpoint().address().to_string() << "\"";
            if( mySocket.is_open() ) {
                ba::post( callbackPool(), std::bind( urgentCallback, shared_from_this() ) );
            }
        }
    } else {
        if( errorCallback ) {
            ba::post( callbackPool(), std::bind( errorCallback, shared_from_this() ) );
        } else {
            if( ( error == ba::error::eof ) || ( error == ba::error::connection_reset ) || (  error == ba::error::operation_aborted) ) {
                lock.unlock();
                close();
            }
        }
    }

}

void TcpConnection::idle( void ) {

    unique_lock<mutex> lock(mtx);
    if( !activityCallback ) {
        return;
    }

    if( mySocket.is_open() ) {
        lock.unlock();
        mySocket.async_read_some( ba::null_buffers(),
                                  boost::bind( &TcpConnection::onActivity, this, ba::placeholders::error, ba::placeholders::bytes_transferred ) );

    }
    
}

void TcpConnection::onActivity( const boost::system::error_code& ec, size_t transferred ) {
    
    boost::system::error_code error = ec;
    
    if( !error ) {
        auto test RDX_UNUSED = mySocket.remote_endpoint( error );  // check if endpoint exists, will throw if not connected.
    }
    if( !error ) {
        uint8_t urgentData(0);
        size_t received RDX_UNUSED = socket().receive( boost::asio::buffer( &urgentData, 1 ), tcp::socket::message_peek, error );
        //if( received == 0 ) error = ba::error::eof;
    }
    
    unique_lock<mutex> lock(mtx);
    if( !error ) {
        if( activityCallback ) {
            //LOG_DEBUG << "Activity on connection \"" << connptr->socket().remote_endpoint().address().to_string() << "\"";
            if( mySocket.is_open() ) {
                ba::post( callbackPool(), std::bind( activityCallback, shared_from_this() ) );
            }
        }
    } else {
        if( errorCallback ) {
            ba::post( callbackPool(), std::bind( errorCallback, shared_from_this() ) );
        } else {
            if( ( error == ba::error::eof ) || ( error == ba::error::connection_reset ) || (  error == ba::error::operation_aborted) ) {
                lock.unlock();
                close();
            }
        }
    }

}


void TcpConnection::sendUrgent( uint8_t c ) {
    if( mySocket.is_open() ) {
        mySocket.send( ba::buffer( &c, 1 ), ba::socket_base::message_out_of_band );
    }
}

 
void TcpConnection::receiveUrgent( uint8_t& c ) {
    if( mySocket.is_open() && mySocket.at_mark() ) {
        mySocket.receive( ba::buffer( &c, 1 ), ba::socket_base::message_out_of_band );
    }
}

 
TcpConnection& TcpConnection::operator<<( const uint8_t& in ) {
    syncWrite(&in, sizeof(uint8_t));
    return *this;
}


TcpConnection& TcpConnection::operator>>( uint8_t& out ) {
    if(  !mySocket.is_open() || ba::read( mySocket, ba::buffer( &out, sizeof( uint8_t ) ) ) < sizeof( uint8_t ) ) {
        out = CMD_ERR;
        throw std::ios_base::failure( "Failed to receive command." );
    }
    return *this;
}


TcpConnection& TcpConnection::operator<<( const std::vector<std::string>& in ) {
    uint64_t inSize(0);
    if( in.size() ) {
        uint64_t messagesSize = redux::util::size( in );
        size_t totSize = messagesSize+sizeof(uint64_t);
        shared_ptr<char> tmp = rdx_get_shared<char>(totSize);
        char* ptr = tmp.get();
        ptr += redux::util::pack( ptr, messagesSize );
        redux::util::pack( ptr, in );
        syncWrite(tmp.get(), totSize);
    } else {
        syncWrite( inSize );
    }
    return *this;
}


TcpConnection& TcpConnection::operator>>( std::vector<std::string>& out ) {
    uint64_t blockSize, received;
    if( mySocket.is_open() ) {
        received = ba::read( mySocket, ba::buffer( &blockSize, sizeof(uint64_t) ) );
        if( received == sizeof(uint64_t) ) {
            if( swapEndian_ ) swapEndian( blockSize );
            if( blockSize ) {
                shared_ptr<char> buf = rdx_get_shared<char>(blockSize+1);
                char* ptr = buf.get();
                memset( ptr, 0, blockSize+1 );
                received = ba::read( mySocket, ba::buffer( ptr, blockSize ) );
                if( received ) {
                    unpack( ptr, out, swapEndian_ );
                }
            }
        }
    }
    return *this;
}

