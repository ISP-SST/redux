#include "redux/daemon.hpp"


#include "redux/logging/logger.hpp"
#include "redux/network/protocol.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/translators.hpp"

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


namespace {

    typedef boost::asio::time_traits<boost::posix_time::ptime> time_traits_t;

}


Daemon::Daemon( po::variables_map& vm ) : Application( vm, LOOP ), params( vm ), jobCounter( 0 ), nQueuedJobs( 0 ),
    hostTimeout(600), myInfo(Host::myInfo()), timer( ioService ), worker( *this ) {

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
}


void Daemon::stop( void ) {
    runMode = EXIT;
    if( myMaster.conn && myMaster.conn->socket().is_open() ) {
        *myMaster.conn << CMD_DISCONNECT;
        myMaster.conn->socket().close();
        myInfo.info.peerType &= ~Host::TP_WORKER;
    }
    ioService.stop();
    worker.stop();
}


void Daemon::maintenance( void ) {

#ifdef DEBUG_
      LOG_TRACE << "Maintenance:  nJobs = " << jobs.size() << "  nConn = " << connections.size() << "  nPeerWIP = " << peerWIP.size() << ende;
#endif

    updateLoadAvg();
    cleanup();
    logger.flushAll();
    updateStatus();   // TODO: use a secondary connection for auxiliary communications
    timer.expires_at( time_traits_t::now() + boost::posix_time::seconds( 5 ) );
    timer.async_wait( boost::bind( &Daemon::maintenance, this ) );
    
}


bool Daemon::doWork( void ) {


    try {
        // start the server
        serverInit();
        // start the maintenance loop
        LOG_DEBUG << "Initializing maintenance timer." << ende;
        timer.expires_at( time_traits_t::now() + boost::posix_time::seconds( 5 ) );
        timer.async_wait( boost::bind( &Daemon::maintenance, this ) );
        // Add some threads for the async work.
        for( std::size_t i = 0; i < 50; ++i ) {
            shared_ptr<thread> t( new thread( boost::bind( &boost::asio::io_service::run, &ioService ) ) );
            threads.push_back( t );
        }
        LOG_DEBUG << "Initializing worker." << ende;
        workerInit();
        worker.init();
        LOG_DEBUG << "Running the asio service." << ende;
        ioService.run();
        // the io_service will keep running/blocking until stop is called, then wait for the threads to make a clean exit.
        LOG_DETAIL << "Waiting for all threads to terminate." << ende;
        for( auto & t : threads ) {
            t->join();
        }
        myInfo.info.peerType = 0;
        threads.clear();
        worker.stop();
        LOG_DETAIL << "Shutting down Daemon." << ende;
    }
    catch( const exception& e ) {
        LOG_FATAL << "Unhandled exception: If something got here, you forgot to catch it !!!\nI will do a hard reset now....\n   reason: " << e.what() << ende;
        throw; // Application::ResetException();
    }

    return true;

}


void Daemon::workerInit( void ) {
    
    string master = params["master"].as<string>();
    
    if( master.empty() ) {
        logger.setContext( "master" );
        logger.setFlushPeriod(1);
    } else {
        myMaster.host.reset( new Host() );
        myMaster.host->info.connectName = master;
        myMaster.host->info.connectPort = params["port"].as<uint16_t>();
       
        connect( myMaster.host->info, myMaster.conn );
        TcpConnection::Ptr logConn;
        connect( myMaster.host->info, logConn );
        int remoteLogFlushPeriod = 5;       // TODO make this a config setting.
        logger.addNetwork( logConn, 0, Logger::getDefaultMask(), remoteLogFlushPeriod );
        logger.setContext( myInfo.info.name );
        logger.setFlushPeriod( remoteLogFlushPeriod );
        LOG_DETAIL << "Running slave with " << myInfo.status.nThreads << " threads." << ende;
        
    }
    
    myInfo.info.peerType = Host::TP_WORKER;

}


void Daemon::connect( network::Host::HostInfo& host, network::TcpConnection::Ptr& conn ) {
    
    if( host.connectName.empty() ) {
        return;
    }
    
    if( !conn ) {
        conn = TcpConnection::newPtr( ioService );
    }
    
    try {
        auto test RDX_UNUSED = conn->socket().remote_endpoint();  // check if endpoint exists, if not we need to re-establish connection.
        return;
    } catch ( ... ) { }

    if( conn->socket().is_open() ) {
        conn->socket().close();
    }

    LOG_TRACE << "Attempting to connect to " << host.connectName << ":" << host.connectPort << ende;

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

}


void Daemon::updateStatus( void ) {

    network::TcpConnection::Ptr conn = getMaster();
    //myInfo.touch();
    
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

        conn->setCallback( bind( &Daemon::activity, this, std::placeholders::_1 ) );
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


void Daemon::activity( TcpConnection::Ptr conn ) {

    
    Command cmd = CMD_ERR;
    try {
        *conn >> cmd;
    }
    catch( const boost::exception& ) {      // disconnected -> close socket and return.
        removeConnection(conn);
        return;
    }

    //LOG_TRACE << "activity():  received cmd = " << ( int )cmd << "  (" << bitString( cmd ) << ")" << ende;

    try {
        
        connections[conn]->touch();

        switch( cmd ) {
            case CMD_DIE: die(conn); break;
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
            default: LOG_DEBUG << "Daemon: Unrecognized command: " << ( int )cmd << "  " << bitString( cmd ) << ende;
                removeConnection(conn);
                return;
        }

    }
    catch( const exception& e ) {
        LOG_DEBUG << "activity(): Exception when parsing incoming command: " << e.what() << ende;
        removeConnection(conn);
        return;
    }

    conn->idle();
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

        if( host->info.peerType & (Host::TP_WORKER|Host::TP_MASTER) ) {
            LOG_DEBUG << "Slave connected: " << host->info.name << ":" << host->info.pid << "  ID=" << host->id << ende;
            peerWIP[host];
        }

    } else LOG_DEBUG << "Host reconnected: " << host->info.name << ":" << host->info.pid << "  conn=" << hexString(conn.get()) << ende;
    
    host->touch();
    host->nConnections++;

}


void Daemon::removeConnection( TcpConnection::Ptr conn ) {

    unique_lock<mutex> lock( peerMutex );
    auto it = connections.find(conn);
    if( it != connections.end() ) {
        auto wip = peerWIP.find( it->second );
        if( wip != peerWIP.end() ) {
            LOG_NOTICE << "Host #" << it->second->id << "  (" << it->second->info.name << ":" << it->second->info.pid << ") disconnected." << ende;
            LOG_NOTICE << "    nConnections: " << it->second->nConnections << "    ptrCnt: " << it->second.use_count() << ende;
            if( it->second->nConnections==0 ) {
                if( wip->second.job ) {
                    LOG_NOTICE << "Returning unfinished work to queue: " << wip->second.print() << ende;
                    wip->second.job->failWork( wip->second );
                }
                for( auto& part: wip->second.parts ) {
                    if( part ) {
                        part->cacheLoad();
                        part->cacheStore(true);
                    }
                }
                wip->second.parts.clear();
                wip->second.previousJob.reset();
                peerWIP.erase(wip);
            }
        }
        it->second->nConnections--;
        it->second.reset();
        connections.erase(it);
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
        for( auto it=peerWIP.begin(); it != peerWIP.end(); ) {
            Job::JobPtr& job = it->second.job;
            Host::Ptr host = it->first;
            boost::posix_time::time_duration elapsed = (now - it->second.workStarted);
            if( job && (elapsed > boost::posix_time::seconds( job->info.timeout )) ) {
                LOG_DETAIL << "Work has not been completed in " << to_simple_string(elapsed) << ": " << it->second.print() << ende;
                failedWIP( it->second );
                peerWIP.erase(it++);                // N.B iterator is invalidated on erase, so the postfix increment is necessary.
                continue;
            }
            
            if( host ) {
                elapsed = (now - host->status.lastSeen);
                if( elapsed > boost::posix_time::seconds( hostTimeout ) ) {
                    LOG_NOTICE << "Peer has not been active in " << to_simple_string(elapsed) << ":  " << host->info.name << ":" << host->info.pid << ende;
                    failedWIP( it->second );
                    for( auto it2 = begin(connections); it2 != end(connections); ++it2 ) {
                        if( (*it2->second) == (*host) ) {
                            it2->first->socket().close();
                            break;
                        }
                    }
                    peerWIP.erase(it++);            // N.B iterator is invalidated on erase, so the postfix increment is necessary.
                    continue;
                }
            }
            
            ++it;
        }
    }
    
    {
        unique_lock<mutex> lock( jobsMutex );
        for( auto& job : jobs ) {
            if( job && (job->info.step == Job::JSTEP_COMPLETED) ) {
                LOG << "Job " << job->info.id << " (" << job->info.name << ") is completed, removing from queue." << ende;
                LLOG(job->logger) << "Job " << job->info.id << " (" << job->info.name << ") is completed, removing from queue." << ende;
                ioService.post( [job](){
                    job->cleanup();
                    job->logger.flushAll();
                });
                job.reset();
            }
        }
        jobs.erase(std::remove_if(jobs.begin(), jobs.end(), [](const shared_ptr<Job>& j){return (j == nullptr);}), jobs.end());
    }

}


void Daemon::failedWIP( WorkInProgress& wip ) {
    
    if( wip.job ) {
        LOG_NOTICE << "Returning failed/unfinished part to queue: " << wip.print() << ende;
        wip.job->failWork( wip );
        for( auto& part: wip.parts ) {
            if( part ) {
                ioService.post( [part](){
                    part->cacheLoad();
                    part->cacheStore(true);
                });
            }
        }
        wip.parts.clear();
        wip.previousJob = std::move(wip.job);
        //it->second.job.reset();
    }
}


void Daemon::die( TcpConnection::Ptr& conn ) {
    LOG_DEBUG << "Received exit command." << ende;
    unique_lock<mutex> lock( peerMutex );
    connections.clear();
    conn->socket().close();
    stop();
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
                    job->info.id = ++jobCounter;
                    if( job->info.name.empty() ) job->info.name = "job_" + to_string( job->info.id );
                    if( job->info.logFile.empty() ) {
                        job->info.logFile = job->info.name + ".log";
                    }
                    job->startLog();
                    if( !job->check() ) throw job_check_failed( "Sanity check failed for \"" + tmpS + "\"-job " + job->info.name );

                    job->info.submitTime = boost::posix_time::second_clock::local_time();
                    LLOG_DETAIL(job->logger) << "Sanity check passed, adding to queue." << ende;
                    ids.push_back( job->info.id );
                    ids[0]++;
                    newjobs.push_back( job );
                } catch( const job_check_failed& e ) {
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
    uint64_t messagesSize(0);
    if( messages.size() ) {
        messagesSize = redux::util::size( messages );
        size_t totSize = messagesSize+sizeof(uint64_t);
        shared_ptr<char> tmp( new char[totSize], []( char* p ){ delete[] p; } );
        char* ptr = tmp.get();
        ptr += redux::util::pack( ptr, messagesSize );
        redux::util::pack( ptr, messages );
        conn->syncWrite(tmp.get(), totSize);
    } else {
        conn->syncWrite(messagesSize);
    }



}


void Daemon::removeJobs( TcpConnection::Ptr& conn ) {

    size_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if( blockSize ) {

        string jobString = string( buf.get() );
        if( iequals( jobString, "all" ) ) {
            unique_lock<mutex> lock( jobsMutex );
            if( jobs.size() ) LOG << "Clearing joblist." << ende;
            jobs.clear();
            return;
        }

        try {   // remove by ID
            bpt::ptree tmpTree;      // just to be able to use the VectorTranslator
            tmpTree.put( "jobs", jobString );
            vector<size_t> jobList = tmpTree.get<vector<size_t>>( "jobs", vector<size_t>() );
            unique_lock<mutex> lock( jobsMutex );
            vector<size_t> deleted_ids; 
            for( auto & jobId : jobList ) {
                jobs.erase( std::remove_if( jobs.begin(), jobs.end(), [jobId, &deleted_ids](const Job::JobPtr& j) {
                        if( !j ) return true;
                        if( j->info.id == jobId ) {
                            deleted_ids.push_back(jobId);
                            return true;
                        }
                        return false;
                    }), jobs.end() );
            }
            if( !deleted_ids.empty() ) LOG << "Removed " << printArray(deleted_ids,"jobs") << ende;
            return;
        }
        catch( const boost::bad_lexical_cast& e ) {
            // ignore: it just means the specified list of jobs count not be cast into unsigned integers (could also be "all" or list of names)
        }

        try {   // remove by name
            vector<string> jobList;
            boost::split( jobList, jobString, boost::is_any_of( "," ) );
            unique_lock<mutex> lock( jobsMutex );
            vector<string> deletedList;
            for( auto & jobId : jobList ) {
                jobs.erase( std::remove_if( jobs.begin(), jobs.end(), [jobId, &deletedList](const Job::JobPtr& j) {
                    if( !j ) return true;
                    if( j->info.name == jobId ) {
                        deletedList.push_back(jobId);
                        return true;
                    }
                    return false;
                }), jobs.end() );
            }
            if( !deletedList.empty() ) LOG << "Removed " << printArray(deletedList,"jobs") << ende;
            return;
        }
        catch( const std::exception& e ) {
            LOG_ERR << "Exception caught when parsing list of jobs to remove: " << e.what() << ende;
            throw e;
        }


    }

}

/*Job::JobPtr Daemon::selectJob( bool localRequest ) {
    for( auto & job : jobs )
        if( it->info.step.load() < Job::JSTEP_RUNNING )
	    return job;
}*/

bool Daemon::getWork( WorkInProgress& wip, uint8_t nThreads ) {

    wip.previousJob = wip.job;
    wip.job.reset();

    unique_lock<mutex> lock( jobsMutex );
    
    if( wip.isRemote ) {              // remote worker
        for( Job::JobPtr& job : jobs ) {
            if( job && (job->info.step & Job::StepUserMask) &&
                job->check() && job->getWork( wip, nThreads, false ) ) {
                wip.job = job;
                wip.workStarted = boost::posix_time::second_clock::local_time();
                for( auto& part: wip.parts ) {
                    part->partStarted = boost::posix_time::second_clock::local_time();
                }
                break;
            }
        }
    } else {                            // local worker
        uint16_t nActive(0);
        for( Job::JobPtr& job : jobs ) {
            if( job && job->active() ) nActive++; // count all user defined steps 
        }
        static int maxActiveJobs(3);            // FIXME: configurable
        bool mayStartNewJob = (nActive < maxActiveJobs);
        for( Job::JobPtr& job : jobs ) {
            if( job && job->check() && job->getWork( wip, nThreads, mayStartNewJob ) ) {
                wip.job = job;
                wip.workStarted = boost::posix_time::second_clock::local_time();
                for( auto& part: wip.parts ) {
                    part->partStarted = boost::posix_time::second_clock::local_time();
                }
                break;
            }
        }
    }
    lock.unlock();
    
    if( wip.job ) {

        wip.job->info.state.store( Job::JSTATE_ACTIVE );
        for( auto& part: wip.parts ) {
            part->cacheLoad(false);
        }
        
    }
    
    return bool(wip.job);
}


void Daemon::sendWork( TcpConnection::Ptr& conn ) {

     shared_ptr<char> data;
     uint64_t count(0); 
     
    {
        unique_lock<mutex> lock( peerMutex );
        
        auto it = peerWIP.find(connections[conn]);
        if( it == peerWIP.end() ) {
            auto rit = peerWIP.emplace( connections[conn], WorkInProgress() );
            if( !rit.second ) {
                LOG_ERR << "Failed to insert new connection into peerWIP" << ende;
            }
            it = rit.first;
        }
        lock.unlock();      // this should be safe as long as the "it" entry in peerWip is not erased while the below is executing.
        
        WorkInProgress& wip = it->second;
        auto host = it->first;

        uint64_t blockSize = 0;
        wip.isRemote = true;
        if( getWork( wip, host->status.nThreads ) ) {
            blockSize += wip.workSize();
            //LLOG_DETAIL(wip.job->getLogger()) << "Sending work to " << host->info.name << ":" << host->info.pid << "   " << wip.print() << ende;
            LOG_DETAIL << "Sending work to " << host->info.name << ":" << host->info.pid << "   " << wip.print() << ende;
            host->status.state = Host::ST_ACTIVE;
        }

        data.reset( new char[blockSize+sizeof(uint64_t)], []( char* p ){ delete[] p; } );
        char* ptr = data.get()+sizeof(uint64_t);
        if(wip.job) {

            count += wip.packWork( ptr+count );
            for(auto& part: wip.parts) {
                part->cacheStore(true);
            }
            wip.previousJob = wip.job;
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
        auto it = peerWIP.find(host);
        if( it == peerWIP.end() ) {
            LOG_ERR << "This slave is not listed in peerWIP: " << host->info.name << ":" << host->info.pid << ende;
        } else {

            shared_ptr<WorkInProgress> wip( new WorkInProgress() );
            std::swap( it->second, *wip );

            if( wip->job ) {
                it->second.job = wip->job;
                bool endian = conn->getSwapEndian();
                it->first->status.state = Host::ST_IDLE;
                string msg = "Received results from " + it->first->info.name + ":" + to_string(it->first->info.pid);
                std::thread([this,buf,wip,endian,msg](){
                    try {
                        wip->unpackWork( buf.get(), endian );
                        LOG_DETAIL << msg << "   " + wip->print() << ende;
                        //LLOG_DETAIL(wip->job->getLogger()) << msg << "   " + wip->print() << ende;
                        wip->returnResults();
                    } catch ( exception& e ) {
                        LOG_ERR << "putParts:  exception when unpacking results: " << e.what() << ende;
                        failedWIP(*wip);
                    }
                }).detach();
                reply = CMD_OK;         // all ok
            } else {
                LOG_DETAIL << "Received results from unexpected host. It probably timed out." << ende;
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
            //job->setProgressString();
            blockSize += job->info.size();
        }
    }
    uint64_t totalSize = blockSize + sizeof( uint64_t );                    // blockSize will be sent before block
    shared_ptr<char> buf( new char[totalSize], []( char* p ){ delete[] p; } );
    char* ptr = buf.get()+sizeof( uint64_t );
    uint64_t packedSize(0);                                                      // store real blockSize
    for( auto & job : jobs ) {
        if(job) {
            auto lock = job->getLock();
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
        
        
    }
}


void Daemon::updateLoadAvg( void ) {
    //LOG_TRACE << "updateLoadAvg()" << ende;

    static double loadAvg[3];

    int ret = getloadavg( loadAvg, 3 );
    if( ret != 3 ) {
        LOG_ERR << "updateLoadAvg(): failed to get loadavg." << ende;
        myInfo.status.loadAvg = 0;
        return;
    }

    myInfo.status.loadAvg = loadAvg[0]; // / myInfo.info.nCores * 100.0;
    //LOG_TRACE << "updateLoadAvg()  = " << myInfo.status.loadAvg << ende;
}
