#include "redux/daemon.hpp"


#include "redux/logging/logger.hpp"
#include "redux/network/protocol.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stopwatch.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/translators.hpp"
#include "redux/image/cachedfile.hpp"
#include "redux/revision.hpp"

#include <functional>

#include <boost/asio/time_traits.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>


using boost::algorithm::iequals;

using namespace redux::logging;
using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;


Daemon::Daemon( po::variables_map& vm ) : Application( vm, LOOP ), params( vm ), jobCounter( 1 ), nQueuedJobs( 0 ),
    hostTimeout(3600), myInfo(Host::myInfo()), timer( ioService ), worker( *this ) {

    logger.setFlushPeriod( 100 );   // delay log-flushing until we have set local/remote logging in workerInit
    logger.setContext( myInfo.info.name + ":" +to_string(myInfo.info.pid) + " init" );

    uint16_t nThreads = params["threads"].as<uint16_t>();
    if( nThreads ) {
        myInfo.status.nThreads = myInfo.status.maxThreads = nThreads;
    }

    if( params.count("cache-dir") ) {
        auto & c = Cache::get();
        c.setPath( params["cache-dir"].as<string>() );
    }
    
    uint16_t port = params["port"].as<uint16_t>();
    string master = params["master"].as<string>();
    if( port < 1024 ) {
        LOG_FATAL << "Daemon:  using a port < 1024 requires root permissions, which this program should *not* have." << ende;
        stop();
        return;
    }
        
    if( master.empty() ) {
        try {
            server.reset( new TcpServer( ioService, port ) );
        } catch ( const exception& e ) {
            LOG_FATAL << "Daemon:  failed to start server: " << e.what() << ende;
            server.reset();
            stop();
            return;
        }
    }
    

}


Daemon::~Daemon( void ) {
    cleanup();
    stop();
}


void Daemon::serverInit( void ) {

#ifdef DEBUG_
    LOG_TRACE << "serverInit()" << ende;
#endif
    if( server ) {
        server->setCallback( bind( &Daemon::connected, this, std::placeholders::_1 ) );
        server->accept();
        myInfo.info.peerType |= Host::TP_MASTER;
        LOG_DETAIL << "Starting server on port " << params["port"].as<uint16_t>() << "." << ende;
    }
}


void Daemon::reset( void ) {
    runMode = RESET;
    ioService.stop();
    pool.interrupt_all();
}


void Daemon::stop( void ) {
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
    LOG_TRACE << "Maintenance:  nJobs = " << jobs.size() << "  nConn = " << connections.size() << "  nPeerWIP = " << peerWIP.size() << ende;
#endif
    updateLoadAvg();
    cleanup();
    logger.flushAll();
    updateStatus();   // TODO: use a secondary connection for auxiliary communications
    timer.expires_from_now( boost::posix_time::seconds( 5 ) );
    timer.async_wait( boost::bind( &Daemon::maintenance, this ) );
    
}


bool Daemon::doWork( void ) {


    try {
        // start the server
        serverInit();
        // start the maintenance loop
        LOG << "Version: " << reduxCommitMessage << ende;
        LOG_DEBUG << "Initializing maintenance timer." << ende;
        timer.expires_from_now( boost::posix_time::seconds( 5 ) );
        timer.async_wait( boost::bind( &Daemon::maintenance, this ) );
        // Add some threads for the async work.
        for( int i = 0; i < 50; ++i ) {
            pool.create_thread( [&](){
                while( runMode == LOOP ) {
                    try {
                        ioService.run();
                    } catch( job_error& e ) {
                        LOG_ERR << "Job error: " << e.what() << ende;
                    } catch( exception& e ) {
                        LOG_ERR << "Exception in thread: " << e.what() << ende;
                    } catch( ... ) {
                        LOG_ERR << "Unhandled exception in thread." << ende;
                    }
                }
            });
        }
        LOG_DEBUG << "Initializing worker." << ende;
        if( workerInit() ) {
            worker.start();
        }
        LOG_DEBUG << "Running the asio service." << ende;
        // the io_service will keep running/blocking until stop is called, then wait for the threads to make a clean exit.
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
        logger.addNetwork( logConn, 0, Logger::getDefaultMask(), remoteLogFlushPeriod );
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

    LOG_TRACE << "Attempting to connect to " << host.connectName << ":" << host.connectPort << ende;

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
            if( cmd == CMD_CFG ) {  // handshake requested
                *conn << myInfo.info;
                *conn >> host;
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
            shared_ptr<char> buf( new char[totSize], []( char* p ){ delete[] p; } );
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
            if( myMaster.conn->socket().is_open() ) return myMaster.conn;
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
        *conn << CMD_CFG;           // request handshake
        Host::HostInfo remote_info;
        *conn >> remote_info;
        *conn << myInfo.info;

//        LOG_DEBUG << "connected()  Handshake successful." << ende;
        addConnection( remote_info, conn );

        conn->setCallback( bind( &Daemon::handler, this, std::placeholders::_1 ) );
        conn->setErrorCallback( bind( &Daemon::removeConnection, this, std::placeholders::_1 ) );
        *conn << CMD_OK;           // all ok
        conn->idle();
        return;
    }
    catch( const exception& e ) {
        LOG_ERR << "connected() Failed to process new connection. Reason: " << e.what() << ende;
    } catch( ... ) {
        LOG_ERR << "Daemon::connected() Unhandled exception." << ende;
    }
    conn->socket().close();
}


void Daemon::handler( TcpConnection::Ptr conn ) {
    
    Command cmd = CMD_ERR;
    try {
        try {
            *conn >> cmd;
        } catch( ... ) {      // no data, disconnected or some other unrecoverable error -> close socket and return.
            removeConnection(conn);
            return;
        }
        connections[conn]->touch();
        processCommand( conn, cmd );
        conn->idle();
    } catch( const std::exception& e ) {      // disconnected -> close socket and return.
        LOG_ERR << "handler() Failed to process command Reason: " << e.what() << ende;
        removeConnection(conn);
        return;
    }
     
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
            case CMD_RESET: reset( conn, urgent ); break;
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
            case CMD_DEL_SLV: ;
            case CMD_SLV_RES: resetSlaves(conn, cmd); break;
            default: LOG_DEBUG << "processCommand: not implemented: " << cmdToString(cmd) << " (" << (int)cmd << "," << bitString(cmd) << ")" << ende;
                removeConnection(conn);
                return;
        }

    }
    catch( const exception& e ) {
        LOG_DEBUG << "processCommand(): Exception when parsing incoming command: " << cmdToString(cmd)
                  << " (" << (int)cmd << "," << bitString(cmd) << "): " << e.what() << ende;
        removeConnection(conn);
        return;
    }

    logger.flushBuffer();

}


void Daemon::addConnection( const Host::HostInfo& remote_info, TcpConnection::Ptr& conn ) {
    
    unique_lock<mutex> lock( peerMutex );
    if( myInfo.info == remote_info ) {
        myInfo.info.peerType |= remote_info.peerType;
        LOG_ERR << "addConnection(): Looks like this process is connecting to itself...weird. " << ende;
        //connections[conn] = myInfo;
        return;
    }
    
    auto& host = connections[conn];
    if( host == nullptr ) { // new connection
        conn->setSwapEndian( remote_info.littleEndian != myInfo.info.littleEndian );
        set<uint64_t> used_ids;
        for( auto& p: connections ) {
            if( p.second ) {
                if( p.second->info == remote_info ) {
                    host = p.second;
                    host->touch();
                    host->nConnections++;
                    return;
                }
                used_ids.insert(p.second->id);
            }
        }
        
        uint64_t id = 1; // start indexing peers with 1
        while( used_ids.find(id) != used_ids.end() ) id++;
        Host::Ptr newHost( new Host( remote_info, id ) );
        
        auto it = peerWIP.find(newHost);
        if( it != peerWIP.end() ) {
            LOG_DEBUG << "Re-connection from host with WIP: " << it->first->info.name << ":" << it->first->info.pid << "  ID=" << it->first->id << ende;
            host = it->first;
            host->touch();
            host->nConnections++;
            return;
        }
        host = newHost;

        if( host->info.peerType & Host::TP_WORKER ) {
            LOG_DEBUG << "Slave connected: " << host->info.name << ":" << host->info.pid << "  ID=" << host->id << ende;
            peerWIP.emplace( host, std::make_shared<WorkInProgress>() );
        }

    } else LOG_DEBUG << "Host reconnected: " << host->info.name << ":" << host->info.pid << "  conn=" << hexString(conn.get()) << ende;
    
    host->touch();
    host->nConnections++;

}


void Daemon::removeConnection( TcpConnection::Ptr conn ) {

    unique_lock<mutex> lock( peerMutex );
    auto connit = connections.find(conn);
    if( connit != connections.end() ) {
        Host::Ptr& host = connit->second;
        auto wipit = peerWIP.find( host );
        if( wipit != peerWIP.end() ) {
            WorkInProgress::Ptr& wip = wipit->second;
            LOG_NOTICE << "Host #" << host->id << "  (" << host->info.name << ":" << host->info.pid << ") disconnected." << ende;
            if( wip ) {
                if( wip->job && !wip->parts.empty() ) {
                    LOG_NOTICE << "Returning unfinished work to queue: " << wip->print() << ende;
                    wip->job->ungetWork( wip );
                }
                wip->parts.clear();
                wip->previousJob.reset();
                peerWIP.erase(wipit);
            }
        }
        host->nConnections--;
        host.reset();
        connections.erase(connit);
    }
    conn->setErrorCallback(nullptr);
    conn->setCallback(nullptr);
    conn->socket().close();

}


void Daemon::cleanup( void ) {

    {
        unique_lock<mutex> lock( peerMutex );
        for( auto it=connections.begin(); it != connections.end(); ) {
            if( it->first && !it->first->socket().is_open() ) {
                it->second->nConnections--;
                connections.erase( it++ );          // N.B iterator is invalidated on erase, so the postfix increment is necessary.
            } else ++it;
        }

        boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
        for( auto wipit=peerWIP.begin(); wipit != peerWIP.end(); ) {
            WorkInProgress::Ptr& wip = wipit->second;
            if( wip ) {
                Job::JobPtr& job = wip->job;
                Host::Ptr host = wipit->first;
                boost::posix_time::time_duration elapsed = (now - wip->workStarted);
                if( job && (elapsed > boost::posix_time::seconds( job->info.timeout )) ) {
                    LOG_DETAIL << "Work has not been completed in " << to_simple_string(elapsed) << ": " << wip->print() << ende;
                    failedWIP( wip );
                    //peerWIP.erase(it++);                // N.B iterator is invalidated on erase, so the postfix increment is necessary.
                    continue;
                }
                
                if( host ) {
                    elapsed = (now - host->status.lastSeen);
                    if( elapsed > boost::posix_time::seconds( hostTimeout ) ) {
                        LOG_NOTICE << "Peer has not been active in " << to_simple_string(elapsed) << ":  " << host->info.name << ":" << host->info.pid << ende;
                        failedWIP( wip );
                        for( auto it2 = begin(connections); it2 != end(connections); ++it2 ) {
                            if( (*it2->second) == (*host) ) {
                                it2->first->socket().close();
                                break;
                            }
                        }
                        peerWIP.erase(wipit++);            // N.B iterator is invalidated on erase, so the postfix increment is necessary.
                        continue;
                    }
                }
            } else {
                LOG_DEBUG << "Peer has an invalid WIP." << ende;
                failedWIP( wip );
                peerWIP.erase(wipit++);            // N.B iterator is invalidated on erase, so the postfix increment is necessary.
                continue;
            }
            ++wipit;
        }
    }
    
    {
        unique_lock<mutex> lock( jobsMutex );
        for( auto& job : jobs ) {
            if( !job ) continue;
            //job->check();
            if( job->info.step == Job::JSTEP_COMPLETED ) {
                LOG << "Job " << job->info.id << " (" << job->info.name << ") is completed, removing from queue." << ende;
                LLOG(job->logger) << "Job " << job->info.id << " (" << job->info.name << ") is completed, removing from queue." << ende;
//                 ioService.post( [job](){
//                     job->cleanup();
//                     job->logger.flushAll();
//                 });
                job.reset();
            }
        }
        jobs.erase(std::remove_if(jobs.begin(), jobs.end(), [](const shared_ptr<Job>& j){return (j == nullptr);}), jobs.end());
    }

}


void Daemon::failedWIP( WorkInProgress::Ptr wip ) {
    
    if( wip->job ) {
        LOG_NOTICE << "Returning failed/unfinished part to queue: " << wip->print() << ende;
        wip->job->failWork( wip );
        for( auto& part: wip->parts ) {
            if( part ) {
                ioService.post( [part](){
                    part->cacheLoad();
                    part->cacheStore(true);
                });
            }
        }
        wip->parts.clear();
        wip->previousJob = wip->job;
        wip->job.reset();
        //it->second.job.reset();
    }
}


void Daemon::die(void) {

    LOG_DEBUG << "Received exit command from master." << ende;
    worker.stop();
    stop();
    
}


void Daemon::die( TcpConnection::Ptr& conn, bool urgent ) {
    
    if( urgent ) {
        die();
        return;
    }
    
    unique_lock<mutex> lock( peerMutex );
    const Host::HostInfo& hi = connections[conn]->info;
    lock.unlock();
    
    //if( hi.user == myInfo.info.user || hi.peerType == Host::TP_MASTER ) {
    if( hi.peerType == Host::TP_MASTER ) {
        LOG_DEBUG << "Received exit command from " << hi.user << "@" << hi.name << ende;
        *conn << CMD_OK;
        lock.lock();
        connections.clear();
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
        reset();
        return;
    }
    
    unique_lock<mutex> lock( peerMutex );
    const Host::HostInfo& hi = connections[conn]->info;
    lock.unlock();

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
        vector<shared_ptr<Job>> newjobs;
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
                    if( !job->check() ) {
                        string msg = "Sanity check failed for \"" + tmpS + "\"-job " + job->info.name + "   ";
                        msg += "(see " + job->info.logFile + " for details)";
                        throw job_error( msg );
                    }
                    job->info.submitTime = boost::posix_time::second_clock::local_time();
                    LLOG_DETAIL(job->logger) << "Sanity check passed, adding to queue." << ende;
                    ids.push_back( job->info.id );
                    ids[0]++;
                    jobCounter++;
                    newjobs.push_back( job );
                    job->stopLog();
                } catch( const job_error& e ) {
                    messages.push_back( e.what() );
                }
            } else throw invalid_argument( "Unrecognized Job tag: \"" + tmpS + "\"" );
        }
        unique_lock<mutex> lock( jobsMutex );
        jobs.insert( jobs.end(), newjobs.begin(), newjobs.end() );
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
        *conn << CMD_OK;           // all ok, return IDs
        conn->syncWrite(ids);
    } else {
        LOG_ERR << "addJobs: Parsing of datablock failed, count = " << count << "   blockSize = " << blockSize << "  bytes." << ende;
        *conn << CMD_ERR;
    }

    // send any messages back
    *conn << messages;

}


void Daemon::removeJobs( TcpConnection::Ptr& conn ) {

    size_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if( blockSize ) {

        unique_lock<mutex> plock( peerMutex );
        const Host::HostInfo& hi = connections[conn]->info;
        plock.unlock();
        
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
                throw e;
            }
        }

        if( removedJobs.size() ) {
            for( auto &j: removedJobs ) {
                if( !j ) continue;
                plock.lock();
                for( auto &pw: peerWIP ) {
                    if( !pw.second || !pw.second->job ) continue;
                    if( j != pw.second->job ) continue;
                    Host::Ptr host = pw.first;
                    for( auto& hconn: connections ) {
                        if( (*hconn.second) == (*host) ) {
                            hconn.first->sendUrgent( CMD_RESET );
                        }
                    }
                }
                plock.unlock();
            }
        }

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

    string msgStr = "Killing";
    if( cmd == CMD_RESET ) msgStr = "Restarting";
    
    if( blockSize ) {

        string slvString = string( buf.get() );
       if( iequals( slvString, "all" ) ) {
            unique_lock<mutex> lock( peerMutex );
            if( peerWIP.size() ) {
                LOG << msgStr << " " << peerWIP.size() << " slaves." << ende;
                for( auto &pw: peerWIP ) {
                    Host::Ptr host = pw.first;
                    for( auto& hconn: connections ) {
                        if( (*hconn.second) == (*host) ) {
                            hconn.first->sendUrgent( cmd );
                        }
                    }
                }
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
            for( auto &pw: peerWIP ) {
                Host::Ptr host = pw.first;
                if( slaveSet.count( host->id ) ) {
                    LOG << msgStr << " slave #" << host->id << " (" << host->info.name << ":" << host->info.pid << ")" << ende;
                    slaveList.push_back(host->id);
                    for( auto& hconn: connections ) {
                        if( (*hconn.second) == (*host) ) {
                            hconn.first->sendUrgent( cmd );
                        }
                    }
                }
            }
            if( !slaveList.empty() ) LOG << msgStr << " " << printArray(slaveList,"slaves") << ende;
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
            for( auto &pw: peerWIP ) {
                Host::Ptr host = pw.first;
                if( slaveSet.count( host->info.name ) ) {
                    LOG << msgStr << " slave #" << host->id << " (" << host->info.name << ":" << host->info.pid << ")" << ende;
                    slaveIdList.push_back(host->id);
                    for( auto& hconn: connections ) {
                        if( (*hconn.second) == (*host) ) {
                            hconn.first->sendUrgent( cmd );
                        }
                    }
                }
            }
            if( !slaveIdList.empty() ) LOG << msgStr << " " << printArray(slaveIdList,"slaves") << ende;
            return;
        }
        catch( const std::exception& e ) {
            LOG_ERR << "Exception caught when parsing list of jobs to remove: " << e.what() << ende;
            throw e;
        }


    }

}


bool Daemon::getWork( WorkInProgress::Ptr& wip, uint8_t nThreads ) {

    bool gotJob(false);
    bool newJob(false);

    unique_lock<mutex> lock( jobsMutex );
    vector<Job::JobPtr> tmpJobs = jobs;        // make a local copy so we can unlock the job-list for other threads.
    lock.unlock();
    map<size_t,map<uint16_t,uint16_t>> activeCounts;
    for( Job::JobPtr& job: tmpJobs ) {
        if( job ) {
            size_t tid = job->getTypeID();
            activeCounts[tid][job->info.step]++;
        }
    }

    for( Job::JobPtr& job: tmpJobs ) {
        if( job && job->getWork( wip, nThreads, activeCounts[job->getTypeID()] ) ) {
            newJob = (job.get() != wip->job.get());
            wip->job = job;
            gotJob = true;
            break;
        }
    }
    
    if( newJob ) {
        wip->job->info.state.store( Job::JSTATE_ACTIVE );
        wip->workStarted = boost::posix_time::second_clock::local_time();
        for( auto& part: wip->parts ) {
            part->cacheLoad(false);
            part->partStarted = wip->workStarted;
        }
    }
    
    return gotJob;
}


void Daemon::sendWork( TcpConnection::Ptr& conn ) {

     shared_ptr<char> data;
     uint64_t count(0); 
     
    {
        unique_lock<mutex> lock( peerMutex );
        
        auto wipit = peerWIP.find(connections[conn]);
        if( wipit == peerWIP.end() ) {
            auto rit = peerWIP.emplace( connections[conn], std::make_shared<WorkInProgress>() );
            if( !rit.second ) {
                LOG_ERR << "Failed to insert new connection into peerWIP" << ende;
            }
            wipit = rit.first;
        }
        
        WorkInProgress::Ptr wip = wipit->second;
        Host::Ptr host = wipit->first;
        lock.unlock();

        if( wip->job && wip->parts.size() ) {   // parts should have been cleared when results returned.
            wip->job->ungetWork( wip );
            //wip->job->failWork( wip );
            wip->resetParts();
        }
        
        uint64_t blockSize = 0;
        wip->isRemote = true;
        wip->job.reset();
        wip->parts.clear();
        wip->nParts = 0;

        if( getWork( wip, host->status.nThreads ) ) {
            auto jlock = wip->job->getLock();
            for( auto& part: wip->parts ) {
                part->cacheLoad();
            }
            blockSize += wip->workSize();
            //LLOG_DETAIL(wip->job->getLogger()) << "Sending work to " << host->info.name << ":" << host->info.pid << "   " << wip->print() << ende;
            LOG_DETAIL << "Sending work to " << host->info.name << ":" << host->info.pid << "   " << wip->print() << ende;
            host->status.statusString = alignLeft(to_string(wip->job->info.id) + ":" + to_string(wip->parts[0]->id),8) + " ...";
            host->active();
        } else {
            wip->previousJob.reset();
            host->idle();
        }
        data.reset( new char[blockSize+sizeof(uint64_t)], []( char* p ){ delete[] p; } );
        char* ptr = data.get()+sizeof(uint64_t);
        if( wip->job ) {
            auto jlock = wip->job->getLock();
            count += wip->packWork( ptr+count );
            for( auto& part: wip->parts ) {
                part->cacheClear();
            }
            wip->previousJob = wip->job;
        }
    }
    
    if( count ) {
        pack( data.get(), count );         // Store actual packed bytecount (something might be compressed)
        conn->syncWrite( data.get(), count+sizeof(uint64_t) );
    } else {
        conn->syncWrite(count);
    }
    
}


void Daemon::putParts( TcpConnection::Ptr& conn ) {

    size_t blockSize;
    Command reply = CMD_ERR;            // return err unless everything seems fine.
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if( blockSize == 0 ) {
        LOG_ERR << "Received empty result." << ende;
    } else {
        unique_lock<mutex> lock( peerMutex );
        Host::Ptr host = connections[conn];
        if( !host ) {
            LOG_ERR << "This connection is not listed in connections: " << hexString(conn.get()) << "  ignoring." << ende;
            return;
        }
        auto wipit = peerWIP.find(host);
        if( wipit == peerWIP.end() ) {
            LOG_ERR << "This slave is not listed in peerWIP: " << host->info.name << ":" << host->info.pid << ende;
        } else {
            
            if( wipit->second->job ) {
                shared_ptr<WorkInProgress> tmpwip( std::make_shared<WorkInProgress>(*wipit->second) );
                bool endian = conn->getSwapEndian();
                wipit->first->status.state = Host::ST_IDLE;
                wipit->second->resetParts();
                string msg = "Received results from " + wipit->first->info.name + ":" + to_string(wipit->first->info.pid);
                std::thread([this,buf,tmpwip,endian,msg](){
                    auto oldParts = tmpwip->parts;     // so we can restore the patch-data in case of a failure below.
                    auto oldJob = tmpwip->job;
                    try {
                        tmpwip->unpackWork( buf.get(), endian );
                        LOG_DETAIL << msg << "   " + tmpwip->print() << ende;
//                        LOG_DETAIL << "  Job = " << hexString(tmpwip->job.get()) << ende;
                        //LLOG_DETAIL(wip->job->getLogger()) << msg << "   " + wip->print() << ende;
                        tmpwip->returnResults();
                    } catch ( exception& e ) {
                        LOG_ERR << "putParts:  exception when unpacking results: " << e.what() << ende;
                        tmpwip->parts = std::move(oldParts);
                        tmpwip->job = oldJob;
                        failedWIP(tmpwip);
                    }
                }).detach();
                reply = CMD_OK;         // all ok
            } else {
                LOG_TRACE << "Received results from unexpected host. It probably timed out." << ende;
            }
            
        }
    }
    *conn << reply;

}


void Daemon::sendJobList( TcpConnection::Ptr& conn ) {

    uint64_t blockSize(0);
    unique_lock<mutex>( jobsMutex );
    for( auto & job : jobs ) {
        blockSize += job->size();
    }
    uint64_t totalSize = blockSize + sizeof( uint64_t );                    // blockSize will be sent before block
    shared_ptr<char> buf( new char[totalSize], []( char* p ){ delete[] p; } );
    char* ptr = buf.get()+sizeof( uint64_t );
    uint64_t packedSize(0);
    for( auto & job : jobs ) {
        packedSize += job->pack( ptr+packedSize );
    }
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

    if( blockSize ) {
        try {
            unique_lock<mutex> lock( peerMutex );
            connections[conn]->status.unpack( buf.get(), conn->getSwapEndian() );
            connections[conn]->touch();     // Note: lastSeen is copied over, so do a new "touch()" afterwards.
        }
        catch( const exception& e ) {
            LOG_ERR << "updateStatus: Exception caught while parsing block: " << e.what() << ende;
        }
    }

}


void Daemon::sendJobStats( TcpConnection::Ptr& conn ) {

    uint64_t blockSize(0);
    unique_lock<mutex>( jobsMutex );
    for( auto & job : jobs ) {
        if(job) {
            blockSize += job->info.size();
        }
    }
    uint64_t totalSize = blockSize + sizeof( uint64_t );                    // blockSize will be sent before block
    shared_ptr<char> buf( new char[totalSize], []( char* p ){ delete[] p; } );
    char* ptr = buf.get()+sizeof( uint64_t );
    uint64_t packedSize(0);                                                      // store real blockSize
    for( auto & job : jobs ) {
        if(job) {
            // should not be necessary to lock since we use a hardcoded string-length for the statusstring.
            // auto lock = job->getLock();
            packedSize += job->info.pack( ptr+packedSize );
        }
    }
    if( packedSize > blockSize ) {
        string msg = "sendJobStats(): Packing mismatch:  packedSize = " + to_string(packedSize) + "   blockSize = " + to_string(blockSize) + "  bytes.";
        throw length_error(msg);
        //LOG_DEBUG << msg << ende;
    }
    totalSize = packedSize + sizeof( uint64_t );
    pack( buf.get(), packedSize );                                                // store real blockSize
    conn->syncWrite( buf.get(), totalSize );

}


void Daemon::sendPeerList( TcpConnection::Ptr& conn ) {

    uint64_t blockSize = myInfo.size();

    unique_lock<mutex> lock( peerMutex );
    std::set<Host::Ptr, Host::Compare> hostList;
    for( auto& conn : connections ) {
        if( peerWIP.find(conn.second) != peerWIP.end() ) {
            auto ret = hostList.emplace( conn.second );
            if( ret.second ) {
                blockSize += conn.second->size();
            }
        }
    }
    lock.unlock();
    
    uint64_t totalSize = 2*blockSize + sizeof( uint64_t );                    // status updates might change packed size slightly, so add a margin
    shared_ptr<char> buf( new char[totalSize], []( char* p ){ delete[] p; } );
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

        Command ret = CMD_ERR;
        if( logid == 0 ) {
            logger.addConnection( conn, connections[conn] );
            ret = CMD_OK;
        } else {
            unique_lock<mutex>( jobsMutex );
            for( Job::JobPtr & job : jobs ) {
                if( job && (job->info.id == logid) ) {
                    job->logger.addConnection( conn, connections[conn] );
                    ret = CMD_OK;
                }
            }
        }
        if( ret == CMD_OK ) {
            unique_lock<mutex> lock( peerMutex );
            auto it = connections.find(conn);
            if( it != connections.end() ) {     // Let the Logger class deal with this connection from now on.
                it->second->nConnections--;
                it->second.reset();
                connections.erase(it);
            }
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
