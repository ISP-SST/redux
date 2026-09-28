#include "redux/network/tcpserver.hpp"

#include "redux/application.hpp"
#include "redux/logging/logger.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <boost/thread.hpp>

using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;

namespace {
    std::map<boost::thread::id,boost::thread*> thread_map;
    std::set<boost::thread::id> old_threads;

    void asyncReadTimed( TcpConnection::Ptr conn, shared_ptr<boost::asio::steady_timer> timer,
                          boost::asio::mutable_buffer buf,
                          function<void(const boost::system::error_code&, size_t)> handler ) {
        timer->async_wait( [conn]( const boost::system::error_code& ec ) {
            if( !ec ) {     // timeout
                boost::system::error_code ignore;
                conn->socket().cancel( ignore );
            }
        });
        boost::asio::async_read( conn->socket(), buf,
            [timer, handler]( const boost::system::error_code& ec, size_t transferred ) {
                timer->cancel();   // no op
                handler( ec, transferred );
            });
    }
}


TcpServer::TcpServer( uint16_t port, uint16_t threads ) : ioContext(), acceptor( ioContext ), endpoint( tcp::v4(), port ),
                onConnected(nullptr), minThreads(threads), nThreads_(0), nConnections(0), running(false), do_handshake(true), do_auth(false) {
        
    start();

}


TcpServer::~TcpServer() {

    stop();

}


void TcpServer::accept(void) {

    TcpConnection::Ptr connection = TcpConnection::newPtr( ioContext );
    acceptor.async_accept( connection->socket(), boost::bind( &TcpServer::onAccept, this, connection, boost::asio::placeholders::error ) );
    
}


void TcpServer::start(void) {
    
    try {
        if( running ) {
            stop();
        }
        running = true;
        ioContext.restart();
        workGuard.emplace(ioContext.get_executor());
        addThread( minThreads );
        acceptor.open(endpoint.protocol());
        acceptor.set_option(tcp::acceptor::reuse_address(true));
        acceptor.bind(endpoint);
        acceptor.listen();
        accept();
    } catch ( const exception& ) {
        stop();
        throw;
    }
}


void TcpServer::start( uint16_t port ) {

    if( running && (port == endpoint.port()) ) {
        return;
    }
    stop();
    endpoint.port( port );
    start();
}


void TcpServer::stop(void) {
    
    if( running ) {
        running = false;
        workGuard.reset();
        ioContext.stop();
        acceptor.close();
        pool.interrupt_all();
        pool.join_all();
    }
    cleanup();

}


unsigned short TcpServer::port(void) {
    
    lock_guard<mutex> lock(mtx);
    return endpoint.port();
    
}


void TcpServer::cleanup(void) {

    lock_guard<mutex> lock(mtx);
    for( auto& t: old_threads ) {
        boost::thread* thr = thread_map[t];
        pool.remove_thread( thr );
        thread_map.erase( t );
        if( thr->joinable() ) thr->join();
        delete thr;
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
        boost::asio::post(ioContext,  [](){ throw Application::ThreadExit(); } );
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
    conn->setUrgentCallback(nullptr);
    if( conn->socket().is_open() ) {
        conn->socket().close();
    }
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


size_t TcpServer::nThreads( void ) const {
    
    lock_guard<mutex> lock( mtx );
    return pool.size();
    
}


void TcpServer::onAccept( TcpConnection::Ptr conn,
                               const boost::system::error_code& error ) {
    accept();    // always start another accept

    if( error ) {
        //LOG_ERR << "TcpServer::onAccept():  asio reports error:" << error.message();    // TODO some intelligent error handling
        return;
    }

    if( !do_handshake ) {
        try {
            Host::HostInfo rhi;
            addConnection( rhi, conn );
            if( onConnected ) {
                onConnected(conn);
            }
            conn->idle();
        } catch( const exception& e ) {
            return;
        }
        return;
    }

    doHandshake( conn );

}


void TcpServer::doHandshake( TcpConnection::Ptr conn ) {

    auto timer = make_shared<boost::asio::steady_timer>( ioContext );
    timer->expires_after( std::chrono::seconds(5) );

    auto cmdBuf = make_shared<Command>( CMD_ERR );

    asyncReadTimed( conn, timer, boost::asio::buffer( cmdBuf.get(), sizeof(Command) ),
        [this, conn, timer, cmdBuf]( const boost::system::error_code& ec, size_t transferred ) {
            if( ec || transferred != sizeof(Command) || *cmdBuf != CMD_CONNECT ) {
                return;
            }

            try {
                *conn << CMD_CFG;           // request handshake
            } catch( const exception& e ) {
                return;
            }

            auto hdrBuf = make_shared<std::array<char, sizeof(uint64_t)+1>>();
            asyncReadTimed( conn, timer, boost::asio::buffer( hdrBuf->data(), hdrBuf->size() ),
                [this, conn, timer, hdrBuf]( const boost::system::error_code& ec2, size_t transferred2 ) {
                    if( ec2 || transferred2 != hdrBuf->size() ) return;

                    int one = 1;
                    bool swap_endian = ( *reinterpret_cast<char*>(&one) != (*hdrBuf)[0] );
                    uint64_t sz(0);
                    redux::util::unpack( hdrBuf->data()+1, sz, swap_endian );
                    if( sz == 0 || sz > (1<<20) ) return;   // HostInfo is small - reject anything absurd rather than allocate it

                    auto bodyBuf = redux::util::rdx_get_shared<char>( sz );
                    asyncReadTimed( conn, timer, boost::asio::buffer( bodyBuf.get(), sz ),
                        [this, conn, timer, bodyBuf, sz, swap_endian]( const boost::system::error_code& ec3, size_t transferred3 ) {
                            if( ec3 || transferred3 != sz ) return;

                            try {
                                Host::HostInfo rhi;
                                rhi.unpack( bodyBuf.get(), swap_endian );

                                Host& hi = Host::myInfo();
                                *conn << hi.info;
                                if( do_auth ) {                       // TODO simple authentication. (key exchange ?)
                                    *conn << CMD_AUTH;
                                    // if auth fails => return, else continue and do the callback
                                }
                                *conn << CMD_OK;           // all ok

                                addConnection( rhi, conn );
                                if( onConnected ) {
                                    onConnected(conn);
                                }
                                conn->idle();
                            } catch( const exception& e ) {
                                return;
                            }
                        });
                });
        });

}


void TcpServer::threadLoop( void ) {

    ++nThreads_;
    while( running ) {
        try {
            boost::this_thread::interruption_point();
            ioContext.run();
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
    --nThreads_;
    lock_guard<mutex> lock(mtx);
    old_threads.insert( boost::this_thread::get_id() );

}


