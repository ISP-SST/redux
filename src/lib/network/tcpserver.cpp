#include "redux/network/tcpserver.hpp"

#include "redux/application.hpp"
#include "redux/logging/logger.hpp"
#include "redux/util/stringutil.hpp"

#include <boost/thread.hpp>

using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;

namespace {
    std::map<boost::thread::id,boost::thread*> thread_map;
    std::set<boost::thread::id> old_threads;
}


TcpServer::TcpServer( uint16_t port, uint16_t threads ) : ioService(), acceptor( ioService, tcp::endpoint( tcp::v4(), port ) ),
                onConnected(nullptr), minThreads(threads), nThreads(0), nConnections(0), running(true), do_handshake(true), do_auth(false) {
        
    start();
    acceptor.set_option(boost::asio::ip::tcp::acceptor::reuse_address(true));

}


void TcpServer::accept(void) {

    TcpConnection::Ptr connection = TcpConnection::newPtr( ioService );
    acceptor.async_accept( connection->socket(), boost::bind( &TcpServer::onAccept, this, connection, boost::asio::placeholders::error ) );
    
}


void TcpServer::start(void) {
    
    running = true;
    workLoop.reset( new boost::asio::io_service::work(ioService) );
    addThread( minThreads );
    accept();
    
}


void TcpServer::stop(void) {
    
    running = false;
    pool.interrupt_all();
    workLoop.reset();
    cleanup();
    
}


void TcpServer::cleanup(void) {

    lock_guard<mutex> lock(mtx);
    for( auto& t: old_threads ) {
        pool.remove_thread( thread_map[t] );
    }
    old_threads.clear();
    
    for( auto it=connections.begin(); it != connections.end(); ) {
        if( it->first && !it->first->socket().is_open() ) {
            it->second->nConnections--;
            connections.erase( it++ );          // N.B iterator is invalidated on erase, so the postfix increment is necessary.
        } else ++it;
    }
    
}


void TcpServer::addThread( uint16_t n ) {
    
    lock_guard<mutex> lock(mtx);
    while( n-- ) {
        boost::thread* t = pool.create_thread( std::bind( &TcpServer::threadLoop, this ) );
        thread_map[ t->get_id() ] = t;
    }
    
}


void TcpServer::delThread( uint16_t n ) {
    
    while( n-- ) {
        ioService.post( [](){ throw Application::ThreadExit(); } );
    }
    
}


void TcpServer::addConnection( const Host::HostInfo& remote_info, TcpConnection::Ptr conn ) {
    
    if( !conn ) {   // should never happen
        return;
    }
    
    lock_guard<mutex> clock( conn->mtx );
    lock_guard<mutex> lock( mtx );
    
    Host& hi = Host::myInfo();
    if( hi.info == remote_info ) {
        // don't add connections from own process, just append peerType
        hi.info.peerType |= remote_info.peerType;
        return;
    }
    
    auto& host = connections[conn];
    
    if( host == nullptr ) {     // new connection, possibly from a Host we already know
        conn->setSwapEndian( remote_info.littleEndian != hi.info.littleEndian );
        for( auto& p: connections ) {
            if( p.second && (p.second->info == remote_info) ) {   // Same hostname+pid, so recycle existing info.
                host = p.second;
                host->touch();
                host->nConnections++;
                return;
            }
        }

        host.reset( new Host( remote_info, conn->id ) );

    }
    
    host->touch();
    host->nConnections++;
    
}


void TcpServer::removeConnection( TcpConnection::Ptr conn ) {       //  remove from connections and close

    if( !conn ) {   // should never happen
        return;
    }
    
    releaseConnection( conn );
    
    conn->setErrorCallback(nullptr);
    conn->setCallback(nullptr);
    conn->socket().close();

}


void TcpServer::releaseConnection( TcpConnection::Ptr conn ) {  // don't close socket etc, just remove from connections

    if( !conn ) {   // should never happen
        return;
    }

    lock_guard<mutex> lock( mtx );
    
    auto connit = connections.find(conn);
    if( connit != connections.end() ) {
        auto& host = connit->second;
        host->nConnections--;
        connections.erase(connit);
    }

}


Host::Ptr TcpServer::getHost( const TcpConnection::Ptr conn ) const {
    
    lock_guard<mutex> lock( mtx );
    Host::Ptr ret;
    auto connit = connections.find( conn );
    if( connit != connections.end() ) {
        ret = connit->second;
    }
    return ret;
    
}


TcpConnection::Ptr TcpServer::getConnection( Host::Ptr h ) {
    lock_guard<mutex> lock( mtx );
    TcpConnection::Ptr ret;
    for( auto& c: connections ) {
        if( *(c.second) == *h ) {
            ret = c.first;
            break;
        }
    }
    return ret;
}


size_t TcpServer::size( void ) const {
    
    lock_guard<mutex> lock( mtx );
    return connections.size();
    
}


void TcpServer::onAccept( TcpConnection::Ptr conn,
                               const boost::system::error_code& error ) {
    accept();    // always start another accept

    if( !error ) {
        try {
            Host remote_host;
            if( do_handshake ) {
                Host& hi = Host::myInfo();
                Command cmd;
                *conn >> cmd;
                if( cmd != CMD_CONNECT ) return;    // The connection is terminated when going out of scope.
                *conn << CMD_CFG;           // request handshake
                *conn >> remote_host.info;
                *conn << hi.info;
                if( do_auth ) {                       // TODO simple authentication. (key exchange ?)
                    *conn << CMD_AUTH;
                    // if auth fails => return, else continue and do the callback
                }
                *conn << CMD_OK;           // all ok
            }
            addConnection( remote_host.info, conn );
            if( onConnected ) {
                onConnected(conn);
            }
            conn->idle();
            
        } catch( const exception& e ) {
            //LOG_TRACE << "onAccept() Failed to process new connection. Reason: " << e.what();   // only report to trace level since connect/disconnect should be quiet.
            return;
        } 
        //LOG_DETAIL << "Accepted connection from \"" << conn->socket().remote_endpoint().address().to_string() << "\"";
    } else {
        //LOG_ERR << "TcpServer::onAccept():  asio reports error:" << error.message();    // TODO some intelligent error handling
    }

}


void TcpServer::threadLoop( void ) {

    ++nThreads;
    while( running ) {
        try {
            boost::this_thread::interruption_point();
            ioService.run();
        } catch( const Application::ThreadExit& ) {
            break;
        } catch( const boost::thread_interrupted& ) {
            break;
        } catch( exception& e ) {
            cerr << "Exception in thread: " << e.what() << endl;
        } catch( ... ) {
            
        }
        usleep(100000);
    }
    --nThreads;
    lock_guard<mutex> lock(mtx);
    old_threads.insert( boost::this_thread::get_id() );

}


