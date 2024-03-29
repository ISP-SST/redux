#include "redux/daemon.hpp"

#ifdef DEBUG_
#   define TRACE_THREADS
#endif

#include "redux/logging/logger.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/network/protocol.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stopwatch.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/util/trace.hpp"
#include "redux/translators.hpp"
#include "redux/image/cachedfile.hpp"
#include "redux/revision.hpp"
#include "redux/version.hpp"

#include <functional>
#include <sys/resource.h> 

#include <boost/asio/time_traits.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

using boost::algorithm::iequals;

using namespace redux::logging;
using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;

namespace {
    thread_local bool logPrint;
    std::map<boost::thread::id,boost::thread*> thread_map;
    std::set<boost::thread::id> old_threads;

    std::atomic<uint16_t> preparing_local(0);
    std::atomic<uint16_t> preparing_remote(0);

}

#ifdef DEBUG_
#define RDX_HOST_TIMEOUT 72000
#else
#define RDX_HOST_TIMEOUT 7200
#endif


Daemon::Daemon( po::variables_map& vm ) : Application( vm, LOOP ), params( vm ), jobCounter( 1 ), nQueuedJobs( 0 ),
    hostTimeout(RDX_HOST_TIMEOUT), inTransfers(-1), outTransfers(-1), myInfo(Host::myInfo()),
    max_local(5), max_remote(50), timer( ioService ), worker( *this ) {

    file::setErrorHandling( file::EH_THROW );   // we want to catch and print messages to the log.

    logger.setFlushPeriod( 100 );   // delay log-flushing until we have set local/remote logging in workerInit
    logger.setContext( myInfo.info.name + ":" +to_string(myInfo.info.pid) + " init" );

    uint16_t nThreads = params["threads"].as<uint16_t>();
    if( nThreads ) {
        myInfo.status.nThreads = myInfo.status.maxThreads = nThreads;
    }

    if( params.count("max-transfers") ) {
        uint32_t maxTransfers = params["max-transfers"].as<uint32_t>();
        inTransfers.set( maxTransfers );
        outTransfers.set( maxTransfers );
    }
    
    if( params.count("cache-dir") ) {
        auto & c = Cache::get();
        c.setPath( params["cache-dir"].as<string>() );
    }

    if( params.count("max-running") ) {
        uint32_t maxRunning = params["max-running"].as<uint32_t>();
        Job::JobPtr slask = Job::newJob( "MOMFBD" ); // create a job to trigger static Initializions, otherwise this setting might get mangled.
        momfbd::MomfbdJob::setRunningCount( Job::CountT(1, maxRunning) );
    }

}


Daemon::~Daemon( void ) {
    cleanup();
    Daemon::stop();
}


void Daemon::start_server( uint16_t& port, uint8_t tries ) {
    
    if( port < 1024 ) {
        LOG_ERR << "Daemon: Listening on a port < 1024 requires root permissions, which this program should NOT have !!!" << ende;
        return;
    }
    
    while( tries-- ) {
        try {
            if( server ) {
                server->start( port );
            } else {
                server.reset( new TcpServer( port, thread::hardware_concurrency() ) );
            }
            server->setCallback( bind( &Daemon::connected, this, std::placeholders::_1 ) );
            myInfo.info.peerType |= Host::TP_MASTER;
            myInfo.status.listenPort = server->port();
            LOG << "Started server on port " << server->port() << "." << ende;
            return;
        } catch ( const exception& e ) {
            LOG_TRACE << "Daemon: Failed to start server on port " << port << ": " << e.what() << ende;
            server.reset();
        }
        port++;
    }
}


void Daemon::stop_server( void ) {

    if( server ) {
        LOG_DEBUG << "Stopping server." << ende;
        server->stop();
        server.reset();
        myInfo.status.listenPort = 0;
    }
    
}


void Daemon::reset( void ) {

    LOG << "Resetting daemon." << ende;
    std::thread( [this](){
        std::this_thread::sleep_for(std::chrono::seconds(1));
        stop_server();
        worker.stop();
        runMode = RESET;
        logger.flushAll();
        if( myMaster.conn && myMaster.conn->socket().is_open() ) {
            *myMaster.conn << CMD_DISCONNECT;
            myMaster.conn->socket().close();
            myInfo.info.peerType &= ~Host::TP_WORKER;
        }
        ioService.stop();
        pool.interrupt_all();
    }).detach();
    
}


void Daemon::stop( void ) {

    LOG << "Stopping daemon." << ende;
    stop_server();
    worker.stop();
    runMode = EXIT;
    logger.flushAll();
    if( myMaster.conn && myMaster.conn->socket().is_open() ) {
        *myMaster.conn << CMD_DISCONNECT;
        myMaster.conn->socket().close();
        myInfo.info.peerType &= ~Host::TP_WORKER;
    }
    ioService.stop();
    pool.interrupt_all();
    
}


void Daemon::maintenance( void ) {

#ifdef DEBUG_
    LOG_TRACE << "Maintenance:  nJobs = " << jobs.size() << "  nConn = " << server->size() << "  nPeerWIP = " << peerWIP.size() << ende;
#endif
    updateLoadAvg();
    checkSwapSpace();
    cleanup();
    //checkCurrentUsage();
    

    logger.flushAll();
    updateStatus();   // TODO: use a secondary connection for auxiliary communications
    
    timer.expires_from_now( boost::posix_time::seconds( 5 ) );
    
    boost::posix_time::ptime now = boost::posix_time::second_clock::universal_time();
    boost::posix_time::time_duration elapsed = (now - myInfo.status.lastActive);
    {
        lock_guard<mutex> qlock( wip_lqueue_mtx );
        if( !wip_localqueue.empty() || (elapsed > boost::posix_time::minutes( 5 )) ) {  // kick the corpse every now and then.
            ioService.post( std::bind( &Worker::start, std::ref(worker) ));
        }
    }
    {
        lock_guard<mutex> qlock( wip_queue_mtx );
        size_t nQ = std::min<size_t>( wip_queue.size(), outTransfers.count() );
        if( nQ ) {
            ioService.post( std::bind( &Daemon::pokeSlaves, this, nQ ));
        }
    }
    
    timer.async_wait( boost::bind( &Daemon::maintenance, this ) );
    
}


void Daemon::checkSwapSpace( void ) {
    bfs::path cachePath( Cache::get().path() );
    if( !cachePath.empty() && bfs::exists(cachePath) ) {
// Disabled temporarily for test-compilation on redhat 6.6
/*        boost::system::error_code ec;
        bfs::space_info si = bfs::space(cachePath,ec);
        if( ec ) {
            LOG_ERR << "checkSwapSpace failed for path" << cachePath<< ": " << ec.message() << ende;
        } else {
            double diskFree = static_cast<double>(si.available)/si.capacity;
            double diskFreeGB = static_cast<double>(si.available)/(1<<30);
            if(false)if( diskFree < 0.05 || diskFreeGB < 100 ) {
                LOG_WARN << "Only " << (int)diskFreeGB << "Gb (" << (int)(diskFree*100)
                         << "%) free space on the cache drive (" << cachePath << ")."
                         << "\n\tYour jobs will fail if you run out of space!!" << ende;
            }
        }
*/
    }
}


void Daemon::checkCurrentUsage( void ) {
    unique_lock<mutex> lock( jobsMutex );
    for( auto& job : jobs ) {
        if( !job ) continue;
        job->updateStatus();
        continue;   // FIXME: usage check is not working.
        size_t memUsage = job->memUsage();
        size_t diskUsage = job->diskUsage();
        LOG << "Job " << job->info.id << " (" << job->info.name << "), memUsage = " << static_cast<double>(memUsage)/(1<<30)
            << " Gb, diskUsage = " << static_cast<double>(diskUsage)/(1<<30) << " Gb" << ende;
        //LLOG(job->logger) << "Job " << job->info.id << " (" << job->info.name << ") is completed, removing from queue." << ende;

    }
}


void Daemon::check_limits( void ) {

    struct rlimit rl;   // FIXME: fugly hack until the FD usage is more streamlined.
    if( getrlimit(RLIMIT_NOFILE, &rl) ) {
        LOG_WARN << "Failed to get limit on file-descriptors. errno:" << strerror(errno) << ende;
    } else if( rl.rlim_cur < rl.rlim_max ) {
        LOG_DEBUG << "Raising max open files from " << rl.rlim_cur << " to " << rl.rlim_max << ende;
        rl.rlim_cur = rl.rlim_max;
        if( setrlimit(RLIMIT_NOFILE, &rl) ) {
            LOG_ERR << "Failed to set limit on file-descriptors to max value. errno:" << strerror(errno) << ende;
        }
    }
    if( getrlimit(RLIMIT_CORE, &rl) ) {
        LOG_WARN << "Failed to get limit on coredumps. errno:" << strerror(errno) << ende;
    } else {
        if( rl.rlim_cur < rl.rlim_max ) {
            rlim_t orig = rl.rlim_cur;
            rl.rlim_cur = rl.rlim_max;
            if( setrlimit(RLIMIT_CORE, &rl) ) {
                LOG_ERR << "Failed to set limit on coredumps to max value. errno:" << strerror(errno) << ende;
            } else {
                getrlimit(RLIMIT_CORE, &rl);
                if( rl.rlim_cur != orig ) {
                    LOG_DEBUG << "Coredump limit increased from " << orig << " to " << rl.rlim_cur << ende;
                }
            }
        }
    }

}


bool Daemon::doWork( void ) {


    try {
        
        // log version info
        LOG << "Version: " << reduxCommitMessage << ende;
        
        // start the server
        uint16_t port = params["port"].as<uint16_t>();
        string master = params["master"].as<string>();
        bool dbg = params.count( "log-stdout" );    // check if -d was passed on command-line.
        
        if( master.empty() ) {      // this is a manager instance
            
            start_server( port );
            if( !server ) {
                LOG_FATAL << "The manager has to be able to listen on a port, please try with another port-number!\n"
                          << "Exiting!" << ende;
                die();
            }
            
            check_limits();
            
            std::thread( std::bind( &Daemon::prepareWork, this ) ).detach();

        } else if( dbg ) {
            uint16_t oport=port;
            start_server( port, 20 );
            if( !server ) {
                LOG_WARN << "Failed to start debug-listener on ports " << oport << "-" << (port-1) << ", giving up!" << ende;
            }
        }
        
        // start the maintenance loop
        LOG_DEBUG << "Initializing maintenance timer." << ende;
        timer.expires_from_now( boost::posix_time::seconds( 5 ) );
        timer.async_wait( boost::bind( &Daemon::maintenance, this ) );

        // Add some threads for the async work.
        addThread( thread::hardware_concurrency() );
        
        LOG_DEBUG << "Initializing worker." << ende;
        if( workerInit() ) {
            worker.start();
        }
        LOG_DEBUG << "Running the asio service." << ende;
        // 
        while( (runMode == LOOP) || preparing_local || preparing_remote ) {
            std::this_thread::sleep_for( std::chrono::milliseconds(100) );
        }
        pool.join_all();
        myInfo.info.peerType = 0;
        worker.stop();
        LOG_DETAIL << "Shutting down Daemon." << ende;
    } catch( const ResetException& e ) {
        LOG << "Hard reset requested" << ende;
        throw;
    } catch( const exception& e ) {
        LOG_FATAL << "Unhandled exception: If something got here, you forgot to catch it !!!\nI will do a hard reset now....\n   reason: " << e.what() << ende;
        throw; // Application::ResetException();
    }

    return true;

}


bool Daemon::workerInit( void ) {
    
    string master = params["master"].as<string>();
    
    if( master.empty() ) {
        logger.setContext( "master" );
        logger.setFlushPeriod(1);
    } else {
        myMaster.host.reset( new Host() );
        myMaster.host->info.connectName = master;
        myMaster.host->info.connectPort = params["port"].as<uint16_t>();
        
        connect( myMaster.host->info, myMaster.conn );
        if( !myMaster.conn->socket().is_open() ) {
            logger.addStream( cerr, Logger::getDefaultMask() );
            LOG_ERR << "Failed to connect to master at " << myMaster.host->info.connectName << ":"
                    << myMaster.host->info.connectPort << ende;
            stop();
            return false;
        }
        myMaster.conn->setUrgentCallback( bind( &Daemon::urgentHandler, this, std::placeholders::_1 ) );
        myMaster.conn->uIdle();
        
        TcpConnection::Ptr logConn;

        connect( myMaster.host->info, logConn );

        int remoteLogFlushPeriod = 5;       // TODO make this a config setting.
        logger.addNetwork( ioService, myMaster.host, 0, Logger::getDefaultMask(), remoteLogFlushPeriod );
        logger.setContext( myInfo.info.name+":"+to_string(myInfo.info.pid) );
        logger.setFlushPeriod( remoteLogFlushPeriod );

        LOG_DETAIL << "Running slave with " << myInfo.status.nThreads << " threads." << ende;
        
    }
    
    myInfo.info.peerType = Host::TP_WORKER;
    return true;
}


void Daemon::connect( network::Host::HostInfo& host, network::TcpConnection::Ptr& conn ) {
    
    if( host.connectName.empty() ) {
        LOG_ERR << "Attempting to connect without a hostname." << ende;
        return;
    }
    
    if( !conn ) {
        conn = TcpConnection::newPtr( ioService );
    }
    
    try {
        auto test RDX_UNUSED = conn->socket().remote_endpoint();  // check if endpoint exists, will throw if not connected.
        return;
    } catch ( ... ) {
        // if we get here the socket is disconnected, continue to reconnect below.
    }
    
    if( conn->socket().is_open() ) {
        conn->socket().close();
    }

#ifdef DEBUG_
    LOG_TRACE << "Attempting to connect to " << host.connectName << ":" << host.connectPort << ende;
#endif
    
    try {
        conn->connect( host.connectName, to_string(host.connectPort) );
        if( conn->socket().is_open() ) {

            Command cmd;

            myInfo.info.peerType |= Host::TP_WORKER;
            *conn << CMD_CONNECT;
            *conn >> cmd;
            if( cmd == CMD_AUTH ) {
                // implement
            }
            Host tmpHost;
            if( cmd == CMD_CFG ) {  // handshake requested
                *conn << myInfo.info;
                *conn >> tmpHost.info;
                *conn >> cmd;       // ok or err
            }
            if( cmd != CMD_OK ) {
                LOG_ERR << "Handshake with master failed  (server replied: " << cmd << ")" << ende;
                conn->socket().close();
                //myInfo.info.peerType &= ~Host::TP_WORKER;
            }
        }

    } catch ( const std::exception& e ) {
        LOG_ERR << "Failed to connect to " << host.connectName << ":" << host.connectPort
                << "  reason:" << e.what() << ende;
    } catch ( ... ) {
        LOG_WARN << "Unhandled exception when connecting: " << __LINE__ << ende;
    }
}


void Daemon::updateStatus( void ) {

    network::TcpConnection::Ptr conn = getMaster();
    //myInfo.touch();
    
    try {
        if( conn && conn->socket().is_open() ) {
#ifdef DEBUG_
            LOG_TRACE << "Sending statusupdate to server" << ende;
#endif
            size_t blockSize = myInfo.status.size();
            size_t totSize = blockSize + sizeof( size_t ) + 1;
            shared_ptr<char> buf = rdx_get_shared<char>( totSize );
            char* ptr = buf.get();
            memset( ptr, 0, totSize );

            uint64_t count = pack( ptr, CMD_STAT );
            count += pack( ptr+count, blockSize );
            count += myInfo.status.pack( ptr+count );

            conn->asyncWrite( buf, totSize );
        }
    } catch ( const std::exception& e ) {
        LOG_ERR << "Exception caught while sending statusupdate to server: " << e.what() << ende;
    } catch ( ... ) {
        LOG_ERR << "Unrecognized exception caught while sending statusupdate to server." << ende;
    }
    if( conn ) unlockMaster();
    
}


network::TcpConnection::Ptr Daemon::getMaster(void) {
    
    if ( myMaster.conn ) {
        try {
            myMaster.conn->lock();
            connect( myMaster.host->info, myMaster.conn );          // will return without doing anything if the remote endpoint exists.
            if( myMaster.conn->socket().is_open() ) {
                myMaster.conn->unlock();
                if( !myMaster.conn->hasUrgentCallback() ) {
                    myMaster.conn->setUrgentCallback( bind( &Daemon::urgentHandler, this, std::placeholders::_1 ) );
                }
                myMaster.conn->uIdle();
                myMaster.conn->lock();
                return myMaster.conn;
            }
        }
        catch( exception& e ) {
            LOG_DEBUG << "getMaster: Exception caught while getting connection to master: " << e.what() << ende;
        }
        catch( ... ) {
            LOG_DEBUG << "getMaster: Uncaught exception while getting connection to master." << ende;
        }
        myMaster.conn->unlock();        
    }

    return nullptr;
    
}


void Daemon::unlockMaster(void) {

    if ( myMaster.conn ) {
        myMaster.conn->unlock();
    }

}

        
void Daemon::connected( TcpConnection::Ptr conn ) {

    try {
        THREAD_MARK
        addConnection( conn );
        conn->setCallback( bind( &Daemon::handler, this, std::placeholders::_1 ) );
        conn->setErrorCallback( bind( &Daemon::removeConnection, this, std::placeholders::_1 ) );
        THREAD_UNMARK
        return;
    } catch( const exception& e ) {
        LOG_ERR << "connected() Failed to process new connection. Reason: " << e.what() << ende;
    } catch( ... ) {
        LOG_ERR << "Daemon::connected() Unhandled exception." << ende;
    }
    conn->socket().close();
}


void Daemon::handler( TcpConnection::Ptr conn ) {
    
    THREAD_MARK
    Command cmd = CMD_ERR;
    try {
        try {
            THREAD_MARK
            *conn >> cmd;
            THREAD_MARK
            Host::Ptr host = server->getHost( conn );
            if( !host ) throw std::runtime_error("handler(): Null host returned by server->getHost(conn)");
            host->touch();
        } catch( const exception& e ) {      // no data, disconnected or some other unrecoverable error -> close socket and return.
            server->removeConnection( conn );
            THREAD_UNMARK
            return;
        }
        THREAD_MARK
        processCommand( conn, cmd );
        THREAD_MARK
        conn->idle();
        THREAD_UNMARK
    } catch( const std::exception& e ) {      // disconnected -> close socket and return.
        LOG_ERR << "handler() Failed to process command Reason: " << e.what() << ende;
        server->removeConnection(conn);
    }
    THREAD_UNMARK

}


void Daemon::urgentHandler( TcpConnection::Ptr conn ) {
    
    uint8_t cmd = conn->getUrgentData();
    try {
        processCommand( conn, cmd, true );
        if( cmd == CMD_DIE ) {
            conn->setUrgentCallback( nullptr );
            conn->close();
            return;
        }
        conn->uIdle();
    } catch( const std::exception& e ) {
        LOG_ERR << "urgentHandler() Failed to process command Reason: " << e.what() << ende;
    } catch( ... ) {
        LOG_ERR << "urgentHandler() Unrecognized exception." << ende;
    }

}


void Daemon::processCommand( TcpConnection::Ptr conn, uint8_t cmd, bool urgent ) {

    try {

        switch( cmd ) {
            case CMD_EXIT: softExit(); break;
            case CMD_DIE: die( conn, urgent ); break;
            case CMD_WAKE: worker.start(); break;
            case CMD_RESET: reset( conn, urgent ); break;
            case CMD_FORCE_RESET: reset(); break;
            case CMD_ADD_JOB: addJobs(conn); break;
            case CMD_DEL_JOB: removeJobs(conn); break;
            case CMD_GET_WORK: sendWork(conn); break;
            case CMD_GET_JOBLIST: sendJobList(conn); break;
            case CMD_PUT_PARTS: putParts(conn); break;
            case CMD_STAT: updateHostStatus(conn); break;
            case CMD_JSTAT: sendJobStats(conn); break;
            case CMD_PSTAT: sendPeerList(conn); break;
            case CMD_LOG_CONNECT: addToLog(conn); break;
            case CMD_DISCONNECT: removeConnection(conn); break;
            case CMD_LISTEN: listen(); break;
            case CMD_DEL_SLV: ;
            case CMD_SLV_RES: resetSlaves(conn, cmd); break;
            case CMD_INTERACTIVE: interactive( conn ); break;
            default: LOG_DEBUG << "processCommand: not implemented: " << cmdToString(cmd) << " (" << (int)cmd << "," << bitString(cmd) << ")" << ende;
                removeConnection(conn);
                return;
        }
        if( !conn || !conn->socket().is_open() ) throw std::ios_base::failure("Disconnected.");

    } catch( const ios_base::failure& e ) { // severed connection, just remove it.
        removeConnection(conn);
        return;
    } catch( const exception& e ) {
        LOG_DEBUG << "processCommand(): Exception when parsing incoming command: " << cmdToString(cmd)
                  << " (" << (int)cmd << "," << bitString(cmd) << "): " << e.what() << ende;
        removeConnection(conn);
        return;
    }

    logger.flushBuffer();

}


void Daemon::addConnection( TcpConnection::Ptr conn ) {
    
    if( conn && server ) {
        Host::Ptr host = server->getHost( conn );
        if( host && (host->info.peerType & Host::TP_WORKER) ) {
            getWIP( host );
        }
    }
    
}


void Daemon::removeConnection( TcpConnection::Ptr conn ) {

    if( conn && server ) {
        Host::Ptr host = server->getHost( conn );
        if( host && (host->info.peerType & Host::TP_WORKER) ) {
            if( host->nConnections == 1 ) {
                LOG_DEBUG << "Host #" << host->id << "  (" << host->info.name << ":" << host->info.pid << ") disconnected." << ende;
            }
        }
        unique_lock<mutex> lock( peerMutex );
        peers.erase( host );
        server->removeConnection( conn );
    }

}


void Daemon::cleanup( void ) {
    
    if( server ) {
        server->cleanup();
    }
    
    vector<WorkInProgress::Ptr> timedOutWIPs;        // will clear/reset jobs when the vector goes out of scope.
    {
        unique_lock<mutex> lock( peerMutex );
        boost::posix_time::ptime now = boost::posix_time::second_clock::universal_time();
        for( auto wipit=peerWIP.begin(); wipit != peerWIP.end(); ) {
            WorkInProgress::Ptr wip = wipit->second;
            Host::Ptr host = wipit->first;
            if( wip && host ) {
                Job::JobPtr job = wip->job.lock();
                boost::posix_time::time_duration elapsed = (now - wip->workStarted);
                if( job && (elapsed > boost::posix_time::seconds( job->info.timeout )) ) {
                    LOG_DETAIL << "Work has not been completed in " << to_simple_string(elapsed) << ": " << wip->print() << ende;
                    timedOutWIPs.push_back( wip );
                    wip.reset();
                }
                
                if( host && job ) {
                    elapsed = (now - host->status.lastSeen);
                    if( elapsed > boost::posix_time::seconds( hostTimeout ) ) {
                        LOG_NOTICE << "Peer has not been active in " << to_simple_string(elapsed) << ":  " << host->info.name << ":" << host->info.pid << ende;
                        timedOutWIPs.push_back( wip );
                        wip.reset();
                    }
                } else {
                    wip.reset();
                }
            }
            if( !wip || !host ) {
                peerWIP.erase(wipit++);            // N.B iterator is invalidated on erase, so the postfix increment is necessary.
            } else ++wipit;
        }
    }
    for( auto& wip: timedOutWIPs ) failedWIP( wip );
    
    std::thread( [&](){
        vector<Job::JobPtr> deletedJobs;        // will clear/reset jobs when the vector goes out of scope.
        {
            unique_lock<mutex> lock( jobsMutex );
            for( auto& job : jobs ) {
                if( !job ) continue;
                if( job->info.step == Job::JSTEP_COMPLETED ) {
                    LOG << "Job " << job->info.id << " (" << job->info.name << ") is completed, removing from queue." << ende;
                    LLOG(job->logger) << "Job " << job->info.id << " (" << job->info.name << ") is completed, removing from queue." << ende;
                    deletedJobs.push_back( job );
                    job.reset();
                }
            }
            jobs.erase( std::remove_if(jobs.begin(), jobs.end(), [](const shared_ptr<Job>& j){ return !j; }), jobs.end() );
        }
        // here deletedJobs will be destructed, and the jobs cleaned up. This might take a while, so we do it in a detached thread.
    }).detach();
    
    cleanupThreads();
    
}


void Daemon::failedWIP( WorkInProgress::Ptr wip ) {
    
    if( wip ) {
        Job::JobPtr job = wip->job.lock();
        if( job ) {
            LOG_NOTICE << "Returning failed/unfinished part to queue: " << wip->print() << ende;
            job->failWork( wip );
            for( auto& part: wip->parts ) {
                if( part ) {
                    ioService.post( [part](){
                        part->cacheLoad();
                        part->cacheStore(true);
                    });
                }
            }
        }
        wip->parts.clear();
        wip->job.reset();
        //it->second.job.reset();
    }
}


void Daemon::die(void) {

    LOG_DEBUG << "Received exit command." << ende;
    std::thread(
        [this](){
            std::this_thread::sleep_for(std::chrono::seconds(1));
            stop();
        }).detach();

    
}


void Daemon::die( TcpConnection::Ptr& conn, bool urgent ) {
    
    if( urgent ) {
        die();
        return;
    }
    
    Host::Ptr host = myMaster.host;
    if( !host && server ) {
        host = server->getHost( conn );
    }
    
    if( !host ) return;
    const Host::HostInfo& hi = host->info;
    
    //if( hi.user == myInfo.info.user || hi.peerType == Host::TP_MASTER ) {
    if( hi.peerType == Host::TP_MASTER ) {
        LOG_DEBUG << "Received exit command from " << hi.user << "@" << hi.name << ende;
        *conn << CMD_OK;
        stop();
    } else {
        vector<string> messages;
        messages.push_back("Not permitted.");
        *conn << CMD_NOTICE << messages;
    }

}


void Daemon::softExit( void ) {
    
    if( myInfo.info.peerType == Host::TP_MASTER ) {
        LOG_ERR << "Daemon::softExit, not implemented for master yet." << ende;
    } else {
        LOG_DETAIL << "Slave will exit after the current job is completed." << ende;
        worker.exitWhenDone();
    }
}


void Daemon::reset( TcpConnection::Ptr& conn, bool urgent ) {
    
    if( urgent ) {  // received on a slave, just do it without replying
        worker.resetWhenDone();
        return;
    }
    
    Host::Ptr host = myMaster.host;
    if( !host && server ) {
        host = server->getHost( conn );
    }
    
    if( !host ) return;
    const Host::HostInfo& hi = host->info;

    uint8_t lvl(0);
    *conn >> lvl;

    //if( hi.user == myInfo.info.user || hi.peerType == Host::TP_MASTER ) {
    if( hi.peerType == Host::TP_MASTER ) {
        LOG_DEBUG << "Received reset(" << (int)lvl << ") from " << hi.user << "@" << hi.name << ende;
        switch(lvl) {
            case 0: {
                worker.stop();
                worker.start();
            }
            default: ;
        }
        //if( lvl == 5 ) throw ResetException();
        *conn << CMD_OK;
    } else {
        vector<string> messages;
        messages.push_back("Not permitted.");
        *conn << CMD_NOTICE << messages;
    }

}


void Daemon::addJobs( TcpConnection::Ptr& conn ) {

    size_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if( !blockSize ) {
        *conn << CMD_ERR;
        return;
    }

    const char* ptr = buf.get();
    uint64_t count(0);
    uint16_t nJobs(0);                  // count all jobs present in cfg-file.
    vector<uint64_t> ids( 1, 0 );       // first element in the vector will store number of successfully parsed jobs, the rest are job IDs.
    vector<string> messages;
    bool swap_endian = conn->getSwapEndian();
    try {
        while( count < blockSize ) {
            string tmpS = string( ptr+count );
            Job::JobPtr job = Job::newJob( tmpS );
            if( job ) {
                nJobs++;
                try {
                    count += job->unpack( ptr+count, swap_endian );
                    job->info.id = jobCounter;
                    if( job->info.name.empty() ) job->info.name = "job_" + to_string( job->info.id );
                    if( job->info.logFile.empty() ) {
                        job->info.logFile = job->info.name + ".log";
                    }
                    job->startLog();
                    job->printJobInfo();
                    if( !(job->info.flags&Job::CHECKED) ) {
                        Job::moveTo( job.get(), Job::JSTEP_SUBMIT );
                    }

                    job->info.submitTime = bpx::second_clock::universal_time();
                    ids.push_back( job->info.id );
                    ids[0]++;
                    jobCounter++;
                    job->stopLog();
                    lock_guard<mutex> lock( jobsMutex );
                    jobs.push_back( job );
                } catch( const job_error& e ) {
                    messages.push_back( e.what() );
                }
            } else throw invalid_argument( "Unrecognized Job tag: \"" + tmpS + "\"" );
        }
        lock_guard<mutex> lock( jobsMutex );
        std::sort( jobs.begin(), jobs.end(), [](const Job::JobPtr& a, const Job::JobPtr& b) {
            if(a->info.priority != b->info.priority) return (a->info.priority > b->info.priority);
            if(a->info.step != b->info.step) return (a->info.step > b->info.step);
            return (a->info.id < b->info.id);
        } );
    } catch( const exception& e ) {
        LOG_ERR << "addJobs: Exception caught while parsing block: " << e.what() << ende;
    } catch( ... ) {
        LOG_ERR << "addJobs: Unrecognized exception thrown when parsing cfg file." << ende;
    }

    if( count == blockSize ) {
        if( ids[0] ) LOG << "Received " << ids[0] << " jobs." << printArray(ids.data()+1,ids[0],"  IDs") << ende;
        if( nJobs > ids[0] ) {
            nJobs -= ids[0];
            LOG_WARN << "Feiled to parse " << nJobs << " jobs." << ende;
        }
        *conn << CMD_OK;           // all ok, return IDs
        conn->syncWrite(ids);
    } else {
        LOG_ERR << "addJobs: Parsing of datablock failed, count = " << count << "   blockSize = " << blockSize << "  bytes." << ende;
        *conn << CMD_ERR;
    }

    // send any messages back
    *conn << messages;

}

 
void Daemon::failJobs( const vector<size_t>& jobList ) {

    std::set<size_t> jobSet( jobList.begin(), jobList.end() );
    unique_lock<mutex> lock( jobsMutex );
    for( const auto& job: jobs ) {
        if( job && jobSet.count( job->info.id ) ) Job::moveTo( job.get(), Job::JSTATE_ERR );
    }
            
}


void Daemon::failJobs( string jobString ) {
    
    bpt::ptree tmpTree;      // just to be able to use the VectorTranslator
    tmpTree.put( "jobs", jobString );
    vector<size_t> jobList = tmpTree.get<vector<size_t>>( "jobs", vector<size_t>() );
    failJobs(jobList);
    
}


void Daemon::removeJobs( const vector<size_t>& jobList ) {

    std::thread( [&](){
        std::set<size_t> jobSet( jobList.begin(), jobList.end() );
        vector<Job::JobPtr> removedJobs;
        unique_lock<mutex> lock( jobsMutex );
        jobs.erase( std::remove_if( jobs.begin(), jobs.end(), [&](const Job::JobPtr& job) {
                    if( !job ) return true;
                    if( jobSet.count( job->info.id ) ) {
                        removedJobs.push_back( job );
                        return true;
                    }
                    return false;
                }), jobs.end() );
        // here removedJobs will be destructed, and the jobs cleaned up. This might take a while, so we do it in a detached thread.
    }).detach();
    
}


void Daemon::removeJobs( string jobString ) {
    
    bpt::ptree tmpTree;      // just to be able to use the VectorTranslator
    tmpTree.put( "jobs", jobString );
    vector<size_t> jobList = tmpTree.get<vector<size_t>>( "jobs", vector<size_t>() );
    removeJobs(jobList);
    
}


void Daemon::removeJobs( TcpConnection::Ptr& conn ) {

    size_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );
    
    Host::Ptr host = server->getHost( conn );

    if( blockSize && host ) {

        const Host::HostInfo& hi = host->info;
        
        vector<Job::JobPtr> removedJobs;
        bool done(false);
    
        //if( hi.user == myInfo.info.user || hi.peerType == Host::TP_MASTER ) {
        string jobString = string( buf.get() );
        if( iequals( jobString, "all" ) ) {
            unique_lock<mutex> lock( jobsMutex );
            if( jobs.size() ) LOG << "Clearing joblist." << ende;
            jobs.erase( std::remove_if( jobs.begin(), jobs.end(), [&](const Job::JobPtr& job) {
                    if( !job ) return true;
                    if( !job->mayBeDeleted() ) return false;
                    if( hi.user != job->info.user && hi.peerType != Host::TP_MASTER ) return false;
                    removedJobs.push_back( job );
                    return true;
                }), jobs.end() );
            done = true;
        }

        if( !done ) {
            try {   // remove by ID
                bpt::ptree tmpTree;      // just to be able to use the VectorTranslator
                tmpTree.put( "jobs", jobString );
                vector<size_t> jobList = tmpTree.get<vector<size_t>>( "jobs", vector<size_t>() );
                std::set<size_t> jobSet( jobList.begin(), jobList.end() );
                unique_lock<mutex> lock( jobsMutex );
                jobList.clear();
                jobs.erase( std::remove_if( jobs.begin(), jobs.end(), [&](const Job::JobPtr& job) {
                        if( !job ) return true;
                        if( !job->mayBeDeleted() ) return false;
                        if( hi.user != job->info.user && hi.peerType != Host::TP_MASTER ) return false;
                        if( jobSet.count( job->info.id ) ) {
                            jobList.push_back(job->info.id);
                            removedJobs.push_back( job );
                            return true;
                        }
                        return false;
                    }), jobs.end() );
                if( !jobList.empty() ) LOG << "Removed " << printArray(jobList,"jobs") << ende;
                done = true;
            }
            catch( const boost::bad_lexical_cast& e ) {
                // ignore: it just means the specified list of jobs count not be cast into unsigned integers (could also be "all" or list of names)
            }
        }
        
        if( !done ) {
            try {   // remove by name
                vector<string> jobList;
                boost::split( jobList, jobString, boost::is_any_of( "," ) );
                std::set<string> jobSet( jobList.begin(), jobList.end() );
                unique_lock<mutex> lock( jobsMutex );
                vector<size_t> deletedList;
                jobs.erase( std::remove_if( jobs.begin(), jobs.end(), [&](const Job::JobPtr& job) {
                        if( !job ) return true;
                        if( !job->mayBeDeleted() ) return false;
                        if( hi.user != job->info.user && hi.peerType != Host::TP_MASTER ) return false;
                        if( jobSet.count( job->info.name ) ) {
                            deletedList.push_back(job->info.id);
                            removedJobs.push_back( job );
                            return true;
                        }
                        return false;
                    }) , jobs.end() );
                if( !deletedList.empty() ) LOG << "Removed " << printArray( deletedList, "jobs" ) << ende;
                done = true;
            }
            catch( const std::exception& e ) {
                LOG_ERR << "Exception caught when parsing list of jobs to remove: " << e.what() << ende;
                throw;
            }
        }

        if( removedJobs.size() ) {
            std::thread( [&,removedJobs](){
                for( auto &j: removedJobs ) {
                    if( !j ) continue;
                    lock_guard<mutex> plock( peerMutex );
                    for( auto &pw: peerWIP ) {
                        if( !pw.second ) continue;
                        Job::JobPtr job = pw.second->job.lock();
                        if( !job || (j != job) ) continue;
                        Host::Ptr host = pw.first;
                        if( host ) {
                            TcpConnection::Ptr conn = server->getConnection( host );
                            if( conn ) {
                                conn->sendUrgent( CMD_RESET );
                            }
                        }
                    }
                }
                // here removedJobs will be destructed, and the jobs cleaned up. This might take a while, so we do it in a detached thread.
            }).detach();
        }
    }

}


void Daemon::pokeSlaves( size_t n ) {
    
    unique_lock<mutex> plock( peerMutex );
    vector<Host::Ptr> tmpPeers;
    std::copy( peers.begin(), peers.end(), std::back_inserter(tmpPeers) );
    plock.unlock();

    // remove slave that are already active
    tmpPeers.erase( std::remove_if( tmpPeers.begin(), tmpPeers.end(),
        []( const Host::Ptr& h ) {
            if( !h ) return true;
            if( h->status.load[1] <= 0.0 ) return true;
            return (h->status.state != Host::ST_IDLE);
    }), tmpPeers.end() );

    if( tmpPeers.empty() ) return;
    
    // sort by current load -> select the most idle computers
    std::sort( tmpPeers.begin(), tmpPeers.end(),
        [&](const Host::Ptr& a, const Host::Ptr& b ) {
        return (a->status.load[1] < b->status.load[1]);
    });
    
    string ids;
    for( auto& h: tmpPeers ) {
        if( n-- == 0 ) break;
        h->limbo();
        ids += to_string( h->id ) + " ";
    }
    
    sendToSlaves( CMD_WAKE, ids );
    
}


void Daemon::sendToSlaves( uint8_t cmd, string slvString ) {

    boost::trim( slvString );
    if( slvString.empty() ) return;
    
    string msgStr = "Killing";
    if( cmd == CMD_RESET ) msgStr = "Restarting";
    switch( cmd ) {
        case CMD_DIE:       msgStr = "die"; break;
        case CMD_EXIT:      msgStr = "exit"; break;
        case CMD_FORCE_RESET:
        case CMD_RESET:     msgStr = "reset"; break;
        case CMD_WAKE:      msgStr = "wake"; break;
        default:            msgStr = to_string((int)cmd); break;
    }
    
    if( iequals( slvString, "all" ) ) {
        try {
            unique_lock<mutex> lock( peerMutex );
            if( peers.size() ) {
                LOG_DETAIL << "Sending command \"" << msgStr << "\" to all " << peers.size() << " slaves." << ende;
                for( auto &host: peers ) {
                    if( host ) {
                        TcpConnection::Ptr conn = server->getConnection( host );
                        if( conn ) {
                            conn->sendUrgent( cmd );
                        }
                    }
                }
            }
        }
        catch( const std::exception& e ) {
            LOG_ERR << "Exception caught when parsing list of jobs to remove: " << e.what() << ende;
            throw;
        }
        return;
    }

    try {   // remove by ID
        bpt::ptree tmpTree;      // just to be able to use the VectorTranslator
        tmpTree.put( "slaves", slvString );
        vector<uint64_t> slaveList = tmpTree.get<vector<uint64_t>>( "slaves", vector<uint64_t>() );
        std::set<uint64_t> slaveSet( slaveList.begin(), slaveList.end() );
        unique_lock<mutex> lock( peerMutex );
        slaveList.clear(); 
        for( auto &host: peers ) {
            if( host && slaveSet.count( host->id ) ) {
                TcpConnection::Ptr conn = server->getConnection( host );
                if( conn ) {
                    //LOG_DEBUG << "Sending command \"" << msgStr << "\" to slave #" << host->id << " (" << host->info.name << ":" << host->info.pid << ")" << ende;
                    slaveList.push_back( host->id );
                    conn->sendUrgent( cmd );
                }
            }
        }
        if( !slaveList.empty() ) LOG_TRACE << "Sending command \"" << msgStr << "\" to " << printArray(slaveList,"slaves") << ende;
        return;
    }
    catch( const boost::bad_lexical_cast& e ) {
        // ignore: it just means the specified list of jobs count not be cast into unsigned integers (could also be "all" or list of names)
    }

    try {   // remove by name
        vector<string> slaveList;
        boost::split( slaveList, slvString, boost::is_any_of( "," ) );
        std::set<string> slaveSet( slaveList.begin(), slaveList.end() );
        unique_lock<mutex> lock( peerMutex );
        vector<uint64_t> slaveIdList;
        for( auto &host: peers ) {
            if( host && slaveSet.count( host->info.name ) ) {
                TcpConnection::Ptr conn = server->getConnection( host );
                if( conn ) {
                    slaveIdList.push_back(host->id);
                    conn->sendUrgent( cmd );
                }
            }
        }
        if( !slaveIdList.empty() ) LOG_TRACE << "Sending command \"" << msgStr << "\" to " << printArray(slaveIdList,"slaves") << ende;
        return;
    }
    catch( const std::exception& e ) {
        LOG_ERR << "Exception caught when parsing list of jobs to remove: " << e.what() << ende;
        throw;
    }

}


void Daemon::resetSlaves( TcpConnection::Ptr& conn, uint8_t cmd ) {

    uint8_t hardExit(0);
    if( cmd == CMD_DEL_SLV ) {
        *conn >> hardExit;
    }
    size_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );
    
    if( cmd == CMD_SLV_RES ) cmd = CMD_RESET;
    else if( hardExit ) cmd = CMD_DIE;
    else cmd = CMD_EXIT;

    if( blockSize ) {
        string slvString = string( buf.get() );
        sendToSlaves( cmd, slvString );
    }

}


void Daemon::interactiveCB( TcpConnection::Ptr conn ) {

    try {

        if( !conn ) {
            throw exception();
        }
        auto test RDX_UNUSED = conn->socket().remote_endpoint();  // check if endpoint exists, will throw if not connected.
        size_t blockSize;
        shared_ptr<char> buf = conn->receiveBlock( blockSize );

        if( blockSize ) {
            
            string line( buf.get() );
            if( logPrint && !line.empty() ) LOG_TRACE << "Command: \"" << line << "\"" << ende;
            string cmdStr = popword(line);
            string replyStr;
            Command replyCmd = CMD_OK;
            
            try {
                if( cmdStr == "q" || cmdStr == "quit" || cmdStr == "bye" ) {
                    replyCmd = CMD_DISCONNECT;
                } else if( cmdStr == "e" ) {
                    logPrint = !logPrint;
                    if( logPrint ) replyStr = "Printing interactive output to log.";
                } else if( cmdStr == "cf" ) {
                    auto mf = redux::util::Cache::get().getMap< redux::image::CachedFile, redux::image::Image<float> >();
                    for( auto& it: mf.second ) {
                        replyStr += it.first.filename;
                        replyStr += ": " + to_string(it.second.size()) + "\n";
                    }
                    auto mi = redux::util::Cache::get().getMap< redux::image::CachedFile, redux::image::Image<int16_t> >();
                    for( auto& it: mi.second ) {
                        replyStr += it.first.filename;
                        replyStr += "; " + to_string(it.second.size()) + "\n";
                    }
                } else if( cmdStr == "cacheclear" ) {
                    Cache::cleanup();
                } else if( cmdStr == "cacheinfo" ) {
                    replyStr = Cache::getStats();
                } else if( cmdStr == "del" ) {
                    string argStr = popword(line);
                    if( !argStr.empty() ) {
                        removeJobs(argStr);
                    }
                } else if( cmdStr == "die" ) {
                    die();
                    replyCmd = CMD_DISCONNECT;
                } else if( cmdStr == "fail" ) {
                    string argStr = popword(line);
                    if( !argStr.empty() ) {
                        failJobs(argStr);
                    }
                } else if( cmdStr == "idle" ) {
                    worker.stop();
                } else if( cmdStr == "info" ) {
                    replyStr = getLongVersionString();
                    string tmpS = Trace::getInfo();
                    if( ! tmpS.empty() ) replyStr += "\n" + tmpS;
                } else if( cmdStr == "poke" ) {
                    worker.start();
                } else if( cmdStr == "poken" ) {
                    string argStr = popword(line);
                    size_t n(0);
                    if( !argStr.empty() ) {
                        n = boost::lexical_cast<int>(argStr);
                    }
                    if( n ) pokeSlaves(n);
                } else if( cmdStr == "kill" ) {
                    if( !line.empty() ) {
                        sendToSlaves( CMD_DIE, line );
                    }
                } else if( cmdStr == "listen" ) {
                    if( !line.empty() ) {
                        sendToSlaves( CMD_LISTEN, line );
                    }
                } else if( cmdStr == "max_local" ) {
                    string argStr = popword(line);
                    if( !argStr.empty() ) {
                        max_local = boost::lexical_cast<int>(argStr);
                    }
                    replyStr = "Ok max_local "+to_string(max_local);
                } else if( cmdStr == "max_remote" ) {
                    string argStr = popword(line);
                    if( !argStr.empty() ) {
                        max_remote = boost::lexical_cast<int>(argStr);
                    }
                    replyStr = "Ok max_remote "+to_string(max_remote);
                } else if( cmdStr == "max-recv" ) {
                    string argStr = popword(line);
                    if( !argStr.empty() ) {
                        int nTransfers = boost::lexical_cast<int>(argStr);
                        inTransfers.set( nTransfers );
                    }
                    replyStr = inTransfers.getStatus();
                } else if( cmdStr == "max-send" ) {
                    string argStr = popword(line);
                    if( !argStr.empty() ) {
                        int nTransfers = boost::lexical_cast<int>(argStr);
                        outTransfers.set( nTransfers );
                    }
                    replyStr = outTransfers.getStatus();
                } else if( cmdStr == "show_threads" ) {
                    string argStr = popword(line);
                    bool showAll=false;
                    if( !argStr.empty() ) {
                        showAll = boost::lexical_cast<int>(argStr);
                    }
                    replyStr = thread_traces(showAll);
                } else if( cmdStr == "reset" ) {
                    if( !line.empty() ) {
                        sendToSlaves( CMD_RESET, line );
                    } else {
                        reset();
                        replyCmd = CMD_DISCONNECT;
                    }
                } else if( cmdStr == "freset" ) {
                    if( !line.empty() ) {
                        sendToSlaves( CMD_FORCE_RESET, line );
                    } else {
                        reset();
                        replyCmd = CMD_DISCONNECT;
                    }
                } else if(cmdStr == "segfault") {
                    ioService.post( [](){ std::this_thread::sleep_for (std::chrono::seconds(1));  // delay for reply to be sent back
                                          *(int*)(8) = 0; } );
                } else if( cmdStr == "iothreads" ) {
                    string argStr = popword(line);
                    if( !argStr.empty() ) {
                        int nThreads = boost::lexical_cast<int>(argStr) - nSysThreads();
                        if( nThreads > 0 ) {
                            addThread( nThreads );
                        } else if( nThreads < 0 ) {
                            delThread( std::abs(nThreads) );
                        }
                        std::this_thread::sleep_for( std::chrono::milliseconds( 10 ) );         // small wait so threads are started.
                        cleanupThreads();
                    }
                    replyStr = to_string( nSysThreads() );
                } else if( cmdStr == "sthreads" ) {
                    string argStr = popword(line);
                    if( !argStr.empty() ) {
                        int nThreads = boost::lexical_cast<int>(argStr) - server->nThreads();
                        if( nThreads > 0 ) {
                            server->addThread( nThreads );
                        } else if( nThreads < 0 ) {
                            server->delThread( std::abs(nThreads) );
                        }
                        std::this_thread::sleep_for( std::chrono::milliseconds( 10 ) );         // small wait so threads are started.
                        server->cleanup();
                    }
                    replyStr = to_string( server->nThreads() );
                } else if( cmdStr == "ws" ) {
                    string argStr = popword(line);
                    bool details(false);
                    if( !argStr.empty() ) {
                        details = true;
                    }
                    replyStr = workStatus( details );
                } else if( cmdStr == "threads" ) {
                    string argStr = popword(line);
                    if( !argStr.empty() ) {
                        int nThreads = boost::lexical_cast<int>(argStr);
                        myInfo.status.nThreads = nThreads;
                    }
                    replyStr = to_string( myInfo.status.nThreads );
                } else if( cmdStr == "tracestats" ) {
                    replyStr = Trace::getStats();
                } else if( cmdStr == "trace-bt" ) {
                    replyStr = Trace::getBackTraces();
                } else if( cmdStr == "trace-max-depth" ) {
                    string argStr = popword(line);
                    if( !argStr.empty() ) {
                        int tmd = boost::lexical_cast<int>(argStr);
                        Trace::setMaxDepth( tmd );
                    }
                    replyStr = to_string( Trace::maxDepth() );
                } else if( cmdStr == "version" ) {
                    replyStr = getLongVersionString();
                } else {  // unrecognized
                    if( !cmdStr.empty() ) replyStr = "Huh?";
                }
            } catch( const exception& e ) {
                replyStr = "Failed to parse line: \"" + string(buf.get()) + string("\": ") + e.what();
                LOG_ERR << replyStr << ende;
            } catch( ... ) {
                replyStr = "Failed to parse line: \"" + string(buf.get()) + "\"";
                LOG_ERR << replyStr << ende;
            }
            
            uint64_t replySize = replyStr.length() + 1;
            uint64_t totalSize = replySize + sizeof(uint64_t);
            if( blockSize <= totalSize ) {
                buf = rdx_get_shared<char>(totalSize);
            }
            char* ptr = buf.get();
            ptr += pack( ptr, replySize );
            if( !replyStr.empty() ) {
                replyCmd = CMD_INTERACTIVE;
                if( logPrint ) LOG_DEBUG << "\n" << replyStr << ende;
            }
            ptr += pack( ptr, replyCmd );
            replyStr.copy( ptr, replyStr.length() );
            conn->syncWrite( buf.get(), totalSize );
            
            if( replyCmd == CMD_DISCONNECT ) throw exception();

        } else throw exception();
    } catch( ... ) {    // just disconnect if there is any error
        removeConnection(conn);
        //logPrint = false;
        return;
    }

    conn->idle();
    
}


void Daemon::interactive( TcpConnection::Ptr& conn ) {
    
    Host::Ptr host = server->getHost( conn );
    if( host ) {
        LOG_DETAIL << "Interactive mode from: " << host->info.name << ":" << host->info.pid << ende;
        *conn << CMD_OK;
        conn->setCallback( bind( &Daemon::interactiveCB, this, std::placeholders::_1 ) );
    }
    
}


void Daemon::listen( void ) {
    
    if( myInfo.status.listenPort ) {
        stop_server();
    } else {
        myInfo.status.listenPort = params["port"].as<uint16_t>();
        start_server( myInfo.status.listenPort , 20 );
    }
}


network::Host::Ptr Daemon::getHost( network::TcpConnection::Ptr& conn ) {
    unique_lock<mutex> lock( peerMutex );
    Host::Ptr h = server->getHost( conn );
    if( !h ) {
        LOG_ERR << "Daemon::getHost: server returned a null-host, this should not happen!!" << ende;
        return h;
    }
    auto ret = peers.insert( h );
    if( !ret.second && (h->id != (*ret.first)->id) ) {     // re-connected host, fix the ID
        (*ret.first)->id = h->id;
    }
    return *ret.first;
}


WorkInProgress::Ptr Daemon::getWIP( const network::Host::Ptr& host ) {
    
    if( !host ) {
        throw runtime_error("getWIP() called for NULL host.");
    }
    
    unique_lock<mutex> lock( peerMutex );
    auto it = peerWIP.emplace( host, std::make_shared<WorkInProgress>() );
    return it.first->second;

}


void Daemon::removeWIP( const network::Host::Ptr& host ) {
    
    if( host ) {
        unique_lock<mutex> lock( peerMutex );
        peerWIP.erase( host );
    }

}


void Daemon::updateWIP( const network::Host::Ptr& host, WorkInProgress::Ptr& wip ) {
    THREAD_MARK
    if( host ) {
        THREAD_MARK
        unique_lock<mutex> lock( peerMutex );
        THREAD_MARK
        peerWIP[ host ] = wip;
    }
    
}


void Daemon::sendWork( TcpConnection::Ptr conn ) {

    shared_ptr<char> data;
    uint64_t count(0); 

    
    THREAD_MARK
    uint32_t oldJobID(0);
    *conn >> oldJobID;
    
    WorkInProgress::Ptr wip(nullptr);
    THREAD_MARK
    Host::Ptr host = getHost( conn );
    THREAD_MARK
    if( host ) {
        host->limbo();
        Semaphore::Scope ss( outTransfers, 5 ); // if we 're not allowed a transfer-slot in 5 secs, idle slave & try later.
        if( ss && getWork( wip, true ) ) {
            host->active();
            updateWIP( host, wip );
            if( wip ) {
                Job::JobPtr job = wip->job.lock();
                if( job ) {
                    wip->jobID = oldJobID;
                    uint64_t blockSize = wip->workSize() + sizeof(uint64_t);
                    host->status.statusString = alignLeft(to_string(job->info.id) + ":" + to_string(wip->parts[0]->id),8) + " ...";
                    host->active();
                    data = rdx_get_shared<char>( blockSize );
                    char* ptr = data.get()+sizeof(uint64_t);
                    count += wip->packWork( ptr+count );
                    std::thread([wip](){
                        for( auto& part: wip->parts ) {
                            part->unload();
                        }
                    }).detach();
                    wip->jobID = job->info.id;
                }
            } else {
                LOG_DETAIL << "sendWork: wip/job is NULL. This should NOT happen !!" << ende;
            }
        }
    }
    THREAD_MARK
        
    if( count ) {
        pack( data.get(), count );         // Store actual packed bytecount (something might be compressed)
        LOG_DETAIL << "Sending work to " << host->info.name << ":" << host->info.pid << "   " << wip->print()
                   << "  (size=" << count << ")" << ende;
        conn->syncWrite( data.get(), count+sizeof(uint64_t) );
    } else {
        conn->syncWrite(count);
        host->idle();
    }
    THREAD_UNMARK

}


void Daemon::putParts( TcpConnection::Ptr conn ) {
    THREAD_MARK
    Host::Ptr host = server->getHost( conn );
    THREAD_MARK
    if( host ) {
        bool endian = conn->getSwapEndian();
        string msg = "Received results from " + host->info.name + ":" + to_string(host->info.pid);
        size_t blockSize;
        shared_ptr<char> buf = conn->receiveBlock( blockSize );
        THREAD_MARK
        if( blockSize ) {
            WorkInProgress::Ptr wip = getWIP( host );
            msg += "   " + wip->print();
            //std::thread([this,wip,buf,endian,msg](){
            ioService.post([this,wip,buf,endian,msg](){
            THREAD_MARK
                THREAD_MARK
                WorkInProgress::Ptr tmpwip = getIdleWIP();
                WorkInProgress::Ptr wip_bak = getIdleWIP();
                std::shared_ptr<Job> tmpJob = wip->job.lock();
                *wip_bak = *wip;
                *tmpwip = *wip;
                wip_bak->job = tmpJob;
                tmpwip->job = tmpJob;
                wip->reset();
                THREAD_MARK
                try {
                    tmpwip->unpackWork( buf.get(), tmpJob, endian );
                    tmpwip->returnResults();   // TBD: should this step be async/by manager?
                    returnWork( tmpwip );
                    LOG_DETAIL << msg << ende;
                } catch ( exception& e ) {
                    LOG_ERR << "putParts:  exception when unpacking results: " << e.what() << ende;
                    failedWIP( wip_bak );
                }
                putIdleWIP( tmpwip );
                putIdleWIP( wip_bak );
            });
            //}).detach();
                THREAD_UNMARK
        } else {
            LOG_TRACE << "Received unexpected results. The Host/Job probably timed out." << ende;
            //throw logic_error("Received results from unexpected host. It probably timed out.");
        }
        host->idle();
    }
    THREAD_MARK
    *conn << CMD_OK;
    THREAD_UNMARK
    
}


void Daemon::sendJobList( TcpConnection::Ptr& conn ) {

    unique_lock<mutex> lock( jobsMutex );
    vector<Job::JobPtr> tmpJobs( jobs );
    lock.unlock();
    
    uint64_t blockSize(0);
    for( auto& job : tmpJobs ) {
        if(job) blockSize += job->size();
    }
    uint64_t totalSize = blockSize + sizeof( uint64_t );                    // blockSize will be sent before block
    shared_ptr<char> buf = rdx_get_shared<char>(totalSize);
    char* ptr = buf.get()+sizeof( uint64_t );
    uint64_t packedSize(0);
    for( auto& job : tmpJobs ) {
        if(job) packedSize += job->pack( ptr+packedSize );
    }
    //lock.unlock();
    if( packedSize > blockSize ) {
        string msg = "sendJobList(): Packing mismatch:  packedSize = " + to_string(packedSize) + "   blockSize = " + to_string(blockSize) + "  bytes.";
        throw length_error(msg);
        //LOG_DEBUG << msg << ende;
    }
    totalSize = packedSize + sizeof( uint64_t );
    pack( buf.get(), packedSize );                                                // store real blockSize
    conn->syncWrite( buf.get(), totalSize );

}


void Daemon::updateHostStatus( TcpConnection::Ptr& conn ) {

    size_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if( blockSize && server ) {
        try {
            Host::Ptr host = getHost( conn );
            if( host ) {
                host->status.unpack( buf.get(), conn->getSwapEndian() );
                host->touch();     // Note: lastSeen is copied over, so do a new "touch()" afterwards.
            }
        }
        catch( const exception& e ) {
            LOG_ERR << "updateStatus: Exception caught while parsing block: " << e.what() << ende;
        }
    }

}


void Daemon::sendJobStats( TcpConnection::Ptr& conn ) {

    unique_lock<mutex> lock( jobsMutex );
    vector<Job::JobPtr> tmpJobs( jobs );
    lock.unlock();
    
    uint64_t blockSize(0);
    for( auto& job : tmpJobs ) {
        if( job ) blockSize += job->info.size();
    }
    uint64_t totalSize = blockSize + sizeof( uint64_t );                    // blockSize will be sent before block
    shared_ptr<char> buf = rdx_get_shared<char>(totalSize);
    char* ptr = buf.get()+sizeof( uint64_t );
    uint64_t packedSize(0);                                                      // store real blockSize
    for( auto& job : tmpJobs ) {
        if(job) {
            // should not be necessary to lock since we use a hardcoded string-length for the statusstring.
            // auto lock = job->getLock();
            packedSize += job->info.pack( ptr+packedSize );
        }
    }

    if( packedSize > blockSize ) {
        string msg = "sendJobStats(): Packing mismatch:  packedSize = " + to_string(packedSize) + "   blockSize = " + to_string(blockSize) + "  bytes.";
        throw length_error(msg);
    }
    totalSize = packedSize + sizeof( uint64_t );
    pack( buf.get(), packedSize );                                                // store real blockSize
    conn->syncWrite( buf.get(), totalSize );

}


void Daemon::sendPeerList( TcpConnection::Ptr& conn ) {

    uint64_t blockSize = myInfo.size();

    std::set<Host::Ptr, Host::Compare> hostList;
    {
        lock_guard<mutex> lock( peerMutex );
        for( auto &host: peers ) {
            if( host && server->getConnection(host) ) { // only list hosts with an active connection
                auto ret = hostList.emplace(  host );
                if( ret.second ) {
                    blockSize +=  host->size();
                }
            }
        }
    }
    
    uint64_t totalSize = 2*blockSize + sizeof( uint64_t );                    // status updates might change packed size slightly, so add some margin
    shared_ptr<char> buf = rdx_get_shared<char>(totalSize);
    char* ptr =  buf.get()+sizeof( uint64_t );
    uint64_t packedSize = myInfo.pack( ptr );
    for( auto &host : hostList ) {
        packedSize += host->pack( ptr+packedSize );
    }
    
    if( packedSize > 2*blockSize ) {
        string msg = "sendPeerList(): Packing mismatch:  packedSize = " + to_string(packedSize) + "   blockSize = " + to_string(blockSize) + "  bytes.";
        throw length_error(msg);
    }
    totalSize = packedSize + sizeof( uint64_t );
    pack( buf.get(), packedSize );                                                // store real blockSize
    conn->syncWrite( buf.get(), totalSize );

}


void Daemon::addToLog( network::TcpConnection::Ptr& conn ) {

    try {
        uint32_t logid(0);
        *conn >> logid;

        Host::Ptr host = server->getHost( conn );
        Command ret = CMD_ERR;
        if( host ) {
            if( logid == 0 ) {
                logger.addConnection( conn, host );
                ret = CMD_OK;
            } else {
                //unique_lock<mutex> lock( jobsMutex );
                for( Job::JobPtr & job : jobs ) {
                    if( job && (job->info.id == logid) ) {
                        job->logger.addConnection( conn, host );
                        ret = CMD_OK;
                    }
                }
            }
        }
        if( ret == CMD_OK ) {
            server->releaseConnection(conn);     // FIXME: better way to separate log-connections?
        }
        *conn << ret;
    } catch ( const std::exception& e ) {
        LOG_ERR << "Daemon::addToLog(): exception: " << e.what() << ende;
        
    }
}


void Daemon::updateLoadAvg( void ) {

    static double loadAvg[3];
    static StopWatch sw;

    int ret = getloadavg( loadAvg, 3 );
    if( ret != 3 ) {
        LOG_ERR << "updateLoadAvg(): failed to get loadavg." << ende;
        myInfo.status.load[1] = 0;
        return;
    }

    myInfo.status.load[0] = sw.getLoad();       // cpu usage of this instance
    myInfo.status.load[1] = loadAvg[0];         // system-wide cpu usage 

    sw.start();         // reset stopwatch so we measure usage of the last 5s only.

}


void Daemon::addThread( uint16_t n ) {
    
    lock_guard<mutex> lock(threadMutex);
    while( n-- ) {
        boost::thread* t = pool.create_thread( boost::bind( &Daemon::threadLoop, this ) );
        thread_map[ t->get_id() ] = t;
    }
    
}


void Daemon::delThread( uint16_t n ) {
    
    while( n-- ) {
        ioService.post( [](){ throw Application::ThreadExit(); } );
    }

}


void Daemon::cleanupThreads(void) {

    lock_guard<mutex> lock(threadMutex);
    for( auto& t: old_threads ) {
        pool.remove_thread( thread_map[t] );
    }
    old_threads.clear();

}


size_t Daemon::nSysThreads( void ) {
    lock_guard<mutex> lock( threadMutex );
    return pool.size();
}


void Daemon::threadLoop( void ) {

    while( runMode == LOOP ) {
        try {
            THREAD_UNMARK
            boost::this_thread::interruption_point();
            ioService.run();
        } catch( const ThreadExit& e ) {
            break;
        } catch( job_error& e ) {
            LOG_ERR<< "Job error: " << e.what() << ende;
        } catch( const boost::thread_interrupted& ) {
            LOG_TRACE << "Daemon: Thread interrupted." << ende;
            break;
        } catch( exception& e ) {
            LOG_ERR << "Exception in thread: " << e.what() << ende;
        } catch( ... ) {
            LOG_ERR << "Unhandled exception in thread." << ende;
        }
    }

    lock_guard<mutex> lock(threadMutex);
    old_threads.insert( boost::this_thread::get_id() );


}


void Daemon::putActiveWIP( WorkInProgress::Ptr wip ) {
    
    THREAD_MARK
    lock_guard<mutex> alock( wip_active_mtx );
    THREAD_MARK
    wip_active.insert( wip );
    
}


size_t Daemon::countActiveWIP( WorkInProgress::Ptr wip ) {

    lock_guard<mutex> alock( wip_active_mtx );
    return wip_active.count( wip );

}


size_t Daemon::eraseActiveWIP( WorkInProgress::Ptr wip ) {

    lock_guard<mutex> alock( wip_active_mtx );
    return wip_active.erase( wip );

}


WorkInProgress::Ptr  Daemon::getActiveWIP( WorkInProgress::Ptr wip ) {
    WorkInProgress::Ptr ret(nullptr);
    lock_guard<mutex> alock( wip_active_mtx );
    auto it = wip_active.find( wip );
    if( it != wip_active.end() ) {
        ret = *it;
    }
    return ret;

}


void Daemon::getIdleWIP( WorkInProgress::Ptr& wip ) {

    THREAD_MARK
    lock_guard<mutex> ilock( wip_idle_mtx );
    THREAD_MARK
    if( wip_idle.size() ) {
        wip = std::move( *wip_idle.begin() );
        wip_idle.erase( wip_idle.begin() );
        if( wip ) return;
    }
    
    wip = std::make_shared<WorkInProgress>();
    
}


WorkInProgress::Ptr Daemon::getIdleWIP( void ) {
    WorkInProgress::Ptr ret(nullptr);
    THREAD_MARK
    getIdleWIP( ret );
    THREAD_MARK
    return ret;
}

void Daemon::putIdleWIP( WorkInProgress::Ptr wip ) {
    if( wip ) {
        wip->reset();
        lock_guard<mutex> ilock( wip_idle_mtx );
        wip_idle.insert( wip );
    }
}


WorkInProgress::Ptr Daemon::getLocalWIP( void ) {
    WorkInProgress::Ptr ret(nullptr);
    {
        lock_guard<mutex> qlock( wip_lqueue_mtx );
        if( wip_localqueue.size() ) {
            ret = wip_localqueue.front();
            wip_localqueue.pop_front();
        }
    }
    return ret;
}


WorkInProgress::Ptr Daemon::getRemoteWIP( void ) {
    WorkInProgress::Ptr ret(nullptr);
    {
        THREAD_MARK
        lock_guard<mutex> qlock( wip_queue_mtx );
        THREAD_MARK
        if( wip_queue.size() ) {
            std::swap(ret, wip_queue.front());
            wip_queue.pop_front();
        }
    }
    return ret;
}


bool Daemon::getWork( WorkInProgress::Ptr& wip, bool remote ) {
    THREAD_MARK
    WorkInProgress::Ptr tmp_wip;
    if( remote ) tmp_wip = getRemoteWIP();
    else tmp_wip = getLocalWIP();
    THREAD_MARK
    
    if( !tmp_wip ) {
        THREAD_MARK
        return false;
    }

    THREAD_MARK
    if( !wip ) getIdleWIP( wip );
    THREAD_MARK
    
    Job::JobPtr job = tmp_wip->job.lock();
    if( job ) {
        job->info.state.store( Job::JSTATE_ACTIVE );
    }
    tmp_wip->workStarted = boost::posix_time::second_clock::universal_time();
    for( auto& part: tmp_wip->parts ) {
        part->partStarted = tmp_wip->workStarted;
    }
    
    *wip = *tmp_wip;                //  copy part to leave the job/part unmodified
    wip->job = job;
    putActiveWIP( tmp_wip );
    THREAD_MARK

    return true;
    
}


void Daemon::returnWork( WorkInProgress::Ptr wip ) {
    
    if( wip ) {
        WorkInProgress::Ptr tmp_wip = getActiveWIP(wip);
        if( tmp_wip ) {
            eraseActiveWIP( tmp_wip );
            if( wip->hasResults ) {
                *tmp_wip = *wip;
                lock_guard<mutex> lock( wip_completed_mtx );
                if( tmp_wip->isRemote ) {
                    wip_completed.push_back( tmp_wip );
                } else {
                    wip_completed.push_front( tmp_wip );
                }
                return;
            }
            putIdleWIP( tmp_wip );
        }
        wip->hasResults = false;
    }
    
}


void Daemon::prepareLocalWork( int count ) {
    
    THREAD_MARK


    THREAD_MARK
    auto glock = Job::getGlobalLock();
    map<Job::StepID,Job::CountT> activeCounts = Job::counts;
    glock.unlock();
    THREAD_MARK

    unique_lock<mutex> jlock( jobsMutex );
    vector<Job::JobPtr> tmpJobs = jobs;        // make a local copy so we can unlock the job-list for other threads.
    jlock.unlock();
    THREAD_MARK

    tmpJobs.erase( std::remove_if( tmpJobs.begin(), tmpJobs.end(), []( const shared_ptr<Job>& j ) { return !j; }), tmpJobs.end() );
    
    THREAD_MARK
    std::sort( tmpJobs.begin(), tmpJobs.end(),[&](const Job::JobPtr& a, const Job::JobPtr& b ){
        if(a->info.step != b->info.step) return (a->info.step > b->info.step);
        if(a->info.priority != b->info.priority) return (a->info.priority > b->info.priority);
        return (a->info.id < b->info.id);
    } );

    
    WorkInProgress::Ptr wip(nullptr);
    
    try {
        while( count > 0 ) {
            THREAD_MARK
            getIdleWIP( wip );
            wip->isRemote = false;
            bool gotJob(false);
            try {
                THREAD_MARK
                for( Job::JobPtr job: tmpJobs ) {
                THREAD_MARK
                    if( job && job->getWork( wip, myInfo.status.nThreads, activeCounts ) ) {
                        wip->job = job;
                        gotJob = true;
                        break;
                    }
                }
            } catch ( ... ) { }
            THREAD_MARK

            if( !gotJob ) {
                putIdleWIP( wip );
                break;
            }
            THREAD_MARK
            
            std::thread( [this,wip](){
                THREAD_MARK
                try {
                    for( auto& part: wip->parts ) {
                        part->load();
                    }
                    THREAD_MARK
                    lock_guard<mutex> qlock( wip_queue_mtx );
                    wip_localqueue.push_back( std::move(wip) );
                    THREAD_MARK
                } catch( ... ) {
                    if( wip ) {
                        Job::JobPtr job = wip->job.lock();
                        if( job ) job->ungetWork( wip );
                    }
                    THREAD_MARK
                }
                --preparing_local;
                THREAD_UNMARK
            }).detach();
            THREAD_MARK
            --count;
        }
    } catch( ... ) { }
    THREAD_UNMARK
    preparing_local -= count;

    
}


void Daemon::prepareRemoteWork( int count ) {

    auto glock = Job::getGlobalLock();
    map<Job::StepID,Job::CountT> activeCounts = Job::counts;
    glock.unlock();

    unique_lock<mutex> jlock( jobsMutex );
    vector<Job::JobPtr> tmpJobs = jobs;        // make a local copy so we can unlock the job-list for other threads.
    jlock.unlock();

    WorkInProgress::Ptr wip(nullptr);
    
    try {
        while( count > 0 ) {
            THREAD_MARK
            getIdleWIP( wip );
            wip->isRemote = true;
            bool gotJob(false);
            try {
                THREAD_MARK
                for( Job::JobPtr job: tmpJobs ) {
                THREAD_MARK
                    if( job && job->getWork( wip, 0, activeCounts ) ) {
                        wip->job = job;
                        gotJob = true;
                        break;
                    }
                }
            } catch ( ... ) { }
            THREAD_MARK

            if( !gotJob || !wip ) {
                putIdleWIP( wip );
                break;
            }
            THREAD_MARK
            
            std::thread( [this,wip](){
                THREAD_MARK
                try {
                    for( auto& part: wip->parts ) {
                        if( part ) {
                            part->load();
                            part->prePack();
                        }
                    }
                    THREAD_MARK
                    lock_guard<mutex> qlock( wip_queue_mtx );
                    wip_queue.push_back( wip );
                } catch( ... ) {
                    if( wip ) {
                        Job::JobPtr job = wip->job.lock();
                        if( job ) job->ungetWork( wip );
                    }
                    THREAD_MARK
                }
                --preparing_remote;
                THREAD_UNMARK
            }).detach();
            THREAD_MARK
            --count;
        }
    } catch( ... ) { }
    THREAD_UNMARK
    preparing_remote -= count;
}


void Daemon::prepareWork( void ) {

    while( runMode == LOOP ) {
        try {
            THREAD_UNMARK
            std::this_thread::sleep_for( std::chrono::milliseconds(100) );      // just to avoid a busy-loop.
            THREAD_MARK
            {
                lock_guard<mutex> lock( wip_lqueue_mtx );
                int nToPrepare = 0;
                if( !preparing_local && (wip_localqueue.size()+preparing_local < max_local) ) {
                    nToPrepare = 1;
                }
                if( nToPrepare > 0 ) {
                    THREAD_MARK
                    preparing_local += nToPrepare;
                    std::thread( std::bind( &Daemon::prepareLocalWork, this, nToPrepare  ) ).detach();
                }
            }
            THREAD_MARK
            lock_guard<mutex> lock( wip_queue_mtx );
            int nToPrepare = max_remote - wip_queue.size() - preparing_remote ;
            if( nToPrepare > 0 ) {
                THREAD_MARK
                preparing_remote += nToPrepare;
                std::thread( std::bind( &Daemon::prepareRemoteWork, this, nToPrepare ) ).detach();
            }
            THREAD_MARK
        } catch(...){ }
    }
    
}



string Daemon::workStatus( bool details ) {
    
    string ret = "Work:\n";
    {
        lock_guard<mutex> alock( wip_active_mtx );
        ret += "  Active:      " + to_string(wip_active.size()) + "\n";
    }
    {
        lock_guard<mutex> ilock( wip_idle_mtx );
        ret += "  Idle:        " + to_string(wip_idle.size()) + "\n";
    }
    {
        lock_guard<mutex> qlock( wip_queue_mtx );
        ret += "  Queue:       " + to_string(wip_queue.size()) + "/" + to_string(max_remote) + "  preparing: " + to_string(preparing_remote) + "\n";
        if( details ) {
            for( const auto& w: wip_queue ) {
                if( w ) ret += "    " + w->print() + "  @(" + hexString(w.get()) + ")\n";
            }

        }
        ret += "  LocalQueue:  " + to_string(wip_localqueue.size()) + "/" + to_string(max_local) + "  preparing: " + to_string(preparing_local) + "\n";
        if( details ) {
            for( const auto& w: wip_localqueue ) {
                if( w ) ret += "    " + w->print() + "  @(" + hexString(w.get()) + ")\n";
            }
            ret += "\n";
        }
    }
    {
        lock_guard<mutex> lock( wip_completed_mtx );
        ret += "  Completed:   " + to_string(wip_completed.size());
        if( details ) {
            for( const auto& w: wip_completed ) {
                if( w ) ret += "    " + w->print() + "\n";
            }
        }
    }
    return ret;
    
}


