#include "redux/daemon.hpp"


#include "redux/network/protocol.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/logger.hpp"
#include "redux/translators.hpp"

#include <functional>

#include <boost/asio/time_traits.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>


using boost::algorithm::iequals;

using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;


#define lg Logger::mlg
#define logChannel "deamon"

namespace {

    typedef boost::asio::time_traits<boost::posix_time::ptime> time_traits_t;

}


Daemon::Daemon( po::variables_map& vm ) : Application( vm, LOOP ), params( vm ), jobCounter( 0 ), nQueuedJobs( 0 ),
    hostTimeout(3600), timer( ioService ), worker( *this ) {

    myInfo.reset( new Host() );

    uint16_t nThreads = params["threads"].as<uint16_t>();
    if( nThreads ) {
        myInfo->status.nThreads = myInfo->status.maxThreads = nThreads;
    }

    if( params.count("cache-dir") ) {
        auto & c = Cache::get();
        c.setPath( params["cache-dir"].as<string>() );
    }
    
    uint16_t port = params["port"].as<uint16_t>();
    string master = params["master"].as<string>();
    if( port < 1024 ) {
        LOG_CRITICAL << "Daemon:  using a port < 1024 requires root permissions, which this program should *not* have.";
        stop();
        return;
    }
        
    if( master.empty() ) {
        try {
            server.reset( new TcpServer( ioService, port ) );
        } catch ( const exception& e ) {
            LOG_CRITICAL << "Daemon:  failed to start server: " << e.what();
            server.reset();
            stop();
            return;
        }
    }
    
    LOG << "Daemon:  using " << myInfo->status.nThreads << " threads.";
    

}


Daemon::~Daemon( void ) {
    cleanup();
    stop();
}


void Daemon::serverInit( void ) {

    LOG_TRACE << "serverInit()";
    if( server ) {
        server->setCallback( bind( &Daemon::connected, this, std::placeholders::_1 ) );
        server->accept();
        myInfo->info.peerType |= Host::TP_MASTER;
        LOG_DETAIL << "Starting server on port " << params["port"].as<uint16_t>() << ".";
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
        myInfo->info.peerType &= ~Host::TP_WORKER;
    }
    ioService.stop();
    worker.stop();
}


void Daemon::maintenance( void ) {

//     LOG_DEBUG << "Maintenance:   nJobs = " << jobs.size() << "  nConn = " << connections.size()
//     << "  nPeers = " << peers.size() << " nPeerWIP = " << peerWIP.size() << " nThreads = " << threads.size();

    updateLoadAvg();
    cleanup();
    static int cnt(0);
    if( (++cnt % 6) == 0 ) {
        updateStatus();   // TODO: use a secondary connection for auxiliary communications
        cnt = 0;
    }
    timer.expires_at( time_traits_t::now() + boost::posix_time::seconds( 5 ) );
    timer.async_wait( boost::bind( &Daemon::maintenance, this ) );
}


bool Daemon::doWork( void ) {


    try {
        // start the server
        serverInit();
        // start the maintenance loop
        LOG_TRACE << "Initializing maintenance timer.";
        timer.expires_at( time_traits_t::now() + boost::posix_time::seconds( 5 ) );
        timer.async_wait( boost::bind( &Daemon::maintenance, this ) );
        // Add some threads for the async work.
        for( std::size_t i = 0; i < 50; ++i ) {
            shared_ptr<thread> t( new thread( boost::bind( &boost::asio::io_service::run, &ioService ) ) );
            threads.push_back( t );
        }
        LOG_TRACE << "Initializing worker.";
        workerInit();
        worker.init();
        LOG_TRACE << "Running the asio service.";
        ioService.run();
        // the io_service will keep running/blocking until stop is called, then wait for the threads to make a clean exit.
        LOG_TRACE << "Waiting for all threads to terminate.";
        for( auto & t : threads ) {
            t->join();
        }
        myInfo->info.peerType = 0;
        threads.clear();
        worker.stop();
    }
    catch( const exception& e ) {
        LOG_CRITICAL << "Unhandled exception: If something got here, you forgot to catch it !!!\nI will do a hard reset now....\n   reason: " << e.what();
        throw; // Application::ResetException();
    }

    return true;

}


void Daemon::workerInit( void ) {
    
    string master = params["master"].as<string>();
    
    if( master.empty() ) {
        LOG_DEBUG << "Running local worker only.";
    } else {
        LOG_DEBUG << "Running slave.";
        myMaster.host.reset( new Host() );
        myMaster.conn = TcpConnection::newPtr( ioService );
        connect();
    }
    
}


void Daemon::connect(void) {

    uint16_t port = params["port"].as<uint16_t>();
    string master = params["master"].as<string>();
    
    if( master.empty() ) {
        return;
    }
    
    try {
        auto test UNUSED = myMaster.conn->socket().remote_endpoint();  // check if endpoint exists, if not we need to re-establish connection.
        return;
    } catch ( ... ) { }

    if( myMaster.conn->socket().is_open() ) {
        myMaster.conn->socket().close();
    }

    LOG_TRACE << "Attempting to connect to master at " << master << ":" << port;

    myMaster.conn->connect( master, to_string(port) );
    if( myMaster.conn->socket().is_open() ) {

        Command cmd;

        myInfo->info.peerType |= Host::TP_WORKER;
        *myMaster.conn << CMD_CONNECT;
        *myMaster.conn >> cmd;
        if( cmd == CMD_AUTH ) {
            // implement
        }
        if( cmd == CMD_CFG ) {  // handshake requested
            *myMaster.conn << myInfo->info;
            *myMaster.conn >> myMaster.host->info;
            *myMaster.conn >> cmd;       // ok or err
        }
        if( cmd != CMD_OK ) {
            LOG_ERR << "Handshake with master failed  (server replied: " << cmd << ")";
            myMaster.conn->socket().close();
            myInfo->info.peerType &= ~Host::TP_WORKER;
        } else LOG_TRACE << "Connected.";

    }

}


void Daemon::updateStatus( void ) {

    network::TcpConnection::Ptr conn = getMaster();
    //myInfo->touch();
    
    if( conn && conn->socket().is_open() ) {
#ifdef DEBUG_
        LOG_TRACE << "Sending statusupdate to server";
#endif
        size_t blockSize = myInfo->status.size();
        size_t totSize = blockSize + sizeof( size_t ) + 1;
        shared_ptr<char> buf( new char[totSize], []( char* p ){ delete[] p; } );
        char* ptr = buf.get();
        memset( ptr, 0, totSize );

        uint64_t count = pack( ptr, CMD_STAT );
        count += pack( ptr+count, blockSize );
        count += myInfo->status.pack( ptr+count );

        conn->asyncWrite( buf, totSize );
    }
    
    if( conn ) unlockMaster();
    
}


network::TcpConnection::Ptr Daemon::getMaster(void) {
    
    if ( myMaster.conn ) {
        try {
            myMaster.conn->lock();
            connect();          // will return without doing anything if the remote endpoint exists.
            if( myMaster.conn->socket().is_open() ) return myMaster.conn;
        }
        catch( exception& e ) {
            LOG_DEBUG << "getMaster: Exception caught while getting connection to master: " << e.what();
        }
        catch( ... ) {
            LOG_DEBUG << "getMaster: Uncaught exception while getting connection to master.";
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
        *conn << myInfo->info;

        //LOG_DEBUG << "connected()  Handshake successful.";
        addConnection( remote_info, conn );

        conn->setCallback( bind( &Daemon::activity, this, std::placeholders::_1 ) );
        conn->setErrorCallback( bind( &Daemon::removeConnection, this, std::placeholders::_1 ) );
        *conn << CMD_OK;           // all ok
        conn->idle();
        return;
    }
    catch( const exception& e ) {
        LOG_ERR << "connected() Failed to process new connection. Reason: " << e.what();
    } catch( ... ) {
        LOG_ERR << "Daemon::connected() Unhandled exception.";
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

    //LOG_TRACE << "activity():  received cmd = " << ( int )cmd << "  (" << bitString( cmd ) << ")";

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
            case CMD_DISCONNECT: removeConnection(conn); break;
            default: LOG_DEBUG << "Daemon: Unrecognized command: " << ( int )cmd << "  " << bitString( cmd );
                removeConnection(conn);
                return;
        }

    }
    catch( const exception& e ) {
        LOG_DEBUG << "activity(): Exception when parsing incoming command: " << e.what();
        removeConnection(conn);
        return;
    }

    conn->idle();

}


void Daemon::addConnection( const Host::HostInfo& remote_info, TcpConnection::Ptr& conn ) {
    
    unique_lock<mutex> lock( peerMutex );
    if( myInfo->info == remote_info ) {
        myInfo->info.peerType |= remote_info.peerType;
        connections[conn] = myInfo;
        return;
    }
    
    auto& host = connections[conn];
    if( host == nullptr ) { // new connection
        set<uint64_t> used_ids;
        for( auto& p: connections ) if(p.second) used_ids.insert(p.second->id);
        uint64_t id = 1; // start indexing peers with 1
        while( used_ids.find(id) != used_ids.end() ) id++;
        Host* newHost = new Host( remote_info, id );
        for( auto& h: peerWIP ) {
            if( h.first && (*newHost) == (*h.first) ) {
                delete newHost;
                host = h.first;
                host->id = id;
                id = 0;
                break;
            }
        }
        
        if( id ) {
            host.reset( newHost );
            if( host->info.peerType & (Host::TP_WORKER|Host::TP_MASTER) ) {  // filter out the "ui" connections which will connect/disconnect silently
                LOG_DEBUG << "Slave connected: " << host->info.name << ":" << host->info.pid << "  ID=" << host->id;
            }
            //peers.insert( make_pair( id, host ) );
        } else LOG_DEBUG << "New connection from existing host: " << host->info.name << ":" << host->info.pid << "  ID=" << host->id;
        conn->setSwapEndian(remote_info.littleEndian != myInfo->info.littleEndian);
    }
    
    host->touch();
    host->nConnections++;

}


void Daemon::removeConnection( TcpConnection::Ptr conn ) {

    unique_lock<mutex> lock( peerMutex );
    auto it = connections.find(conn);
    if( it != connections.end()) {
        if( it->second && (it->second->info.peerType & (Host::TP_WORKER|Host::TP_MASTER))) {  // filter out the "ui" connections which will connect/disconnect silently
            it->second->nConnections--;
            LOG_DEBUG << "Host #" << it->second->id << "  (" << it->second->info.name << ":" << it->second->info.pid << ") disconnected.";
        }
        connections.erase(it);
    }
    conn->setErrorCallback(nullptr);
    conn->setCallback(nullptr);
    conn->socket().close();
}


// TODO a timeout to allow a peer to reconnect before removing it (which would clear stats/id when/if it reconnects).
void Daemon::cleanup( void ) {

    {
        unique_lock<mutex> lock( peerMutex );
        for( auto it=connections.begin(); it != connections.end(); ) {
            if( it->first && !it->first->socket().is_open() ) {
                it->second->nConnections--;
                connections.erase( it++ );          // N.B iterator is invalidated on erase, so the postfix increment is necessary.
            } else ++it;
        }

    /*    for( auto it=peers.begin(); it != peers.end(); ) {
            if( it->second->nConnections < 1 ) {
                /-*auto wip = peerWIP.find( it->second );
                if( wip != peerWIP.end() ) {
                    if( wip->second.job ) {
                        LOG_DETAIL << "Peer #" << it->first << " disconnected, returning unfinished work to queue: " << wip->second.print();
                        wip->second.job->ungetWork( wip->second );
                        for( auto& part: wip->second.parts ) {
                            part->cacheLoad();
                            part->cacheStore(true);
                        }
                    }
                    wip->second.parts.clear();
                    wip->second.previousJob.reset();
                    peerWIP.erase(wip);
                }*-/
                peers.erase( it++ );
            } else ++it;
        }
     */   

        boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
        for( auto it=peerWIP.begin(); it != peerWIP.end(); ) {
            Job::JobPtr job = it->second.job;
            Host::Ptr host = it->first;
            boost::posix_time::time_duration elapsed = (now - it->second.workStarted);
            if( job && (elapsed > boost::posix_time::seconds( job->info.timeout )) ) {
                LOG_DETAIL << "Work has not been completed in " << to_simple_string(elapsed) << ": " << it->second.print();
                LOG_WARN << "Returning unfinished work to queue: " << it->second.print();
                job->failWork( it->second );
                for( auto& part: it->second.parts ) {
                    if( part ) {
                        part->cacheLoad();
                        part->cacheStore(true);
                        //LOG_TRACE << "Returned part: #" << part->id << "   step=" << (int)part->step;
                    }
                }
                it->second.parts.clear();
                it->second.previousJob.reset();
                peerWIP.erase(it++);
                continue;
            }
            
            if( host && host->nConnections > 0 ) {
                elapsed = (now - host->status.lastSeen);
                if( elapsed > boost::posix_time::seconds( hostTimeout ) ) {
                    LOG_DETAIL << "Peer has not been active in " << to_simple_string(elapsed) << ":  " << host->info.name << ":" << host->info.pid << "   nConn = " << host->nConnections;
                    for( auto it2 = begin(connections); it2 != end(connections); ++it2 ) {
                        if( (*it2->second) == (*host) ) {
                            it2->first->socket().close();
                            break;
                        }
                    }
                }
            }

            ++it;
        }
    }
    
    {
        unique_lock<mutex> lock( jobsMutex );
        for( auto& job : jobs ) {
            if( job && (job->info.step == Job::JSTEP_COMPLETED) ) {
                LOG_DEBUG << "Job " << job->info.id << " (" << job->info.name << ") is completed, removing from queue.";
                job.reset();
                //job->cleanup();
            }
        }
        jobs.erase(std::remove_if(jobs.begin(), jobs.end(), [](const shared_ptr<Job>& j){return (j == nullptr);}), jobs.end());
    }

}


void Daemon::die( TcpConnection::Ptr& conn ) {
    LOG_DEBUG << "Received exit command.";
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
                        job->setLogChannel(logChannel);
                    } else {
                        job->startLog();
                        job->setLogChannel( job->getLogChannel() );
                    }
                    if( !job->check() ) throw job_check_failed( "Sanity check failed for \"" + tmpS + "\"-job " + job->info.name );
                    job->info.submitTime = boost::posix_time::second_clock::local_time();
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
        LOG_ERR << "addJobs: Exception caught while parsing block: " << e.what();
    } catch( ... ) {
        LOG_ERR << "addJobs: Unrecognized exception thrown when parsing cfg file.";
    }

    if( count == blockSize ) {
        if( ids[0] ) LOG << "Received " << ids[0] << " jobs." << printArray(ids.data()+1,ids[0],"  IDs");
        *conn << CMD_OK;           // all ok, return IDs
        conn->syncWrite(ids);
    } else {
        LOG_ERR << "addJobs: Parsing of datablock failed, count = " << count << "   blockSize = " << blockSize << "  bytes.";
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
            if( jobs.size() ) LOG << "Clearing joblist.";
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
            if( !deleted_ids.empty() ) LOG << "Removed " << printArray(deleted_ids,"jobs");
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
            if( !deletedList.empty() ) LOG << "Removed " << printArray(deletedList,"jobs");
            return;
        }
        catch( const std::exception& e ) {
            LOG_ERR << "Exception caught when parsing list of jobs to remove: " << e.what();
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
                LOG_ERR << "Failed to insert new connection into peerWIP";
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
            LOGC_DETAIL(wip.job->getLogChannel()) << "Sending work to " << host->info.name << ":" << host->info.pid << "   " << wip.print();
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
        LOG_ERR << "Received empty result.";
    } else {
        unique_lock<mutex> lock( peerMutex );
        auto it = peerWIP.find(connections[conn]);
        if( it == peerWIP.end() ) {
            LOG_ERR << "Received results from unexpected host. (probably timed out)";
        } else {

            shared_ptr<WorkInProgress> wip( new WorkInProgress() );
            std::swap( it->second, *wip );

            it->second.job = wip->job;
            bool endian = conn->getSwapEndian();
            it->first->status.state = Host::ST_IDLE;
            string msg = "Received results from " + it->first->info.name + ":" + to_string(it->first->info.pid);
            std::thread([buf,wip,endian,msg](){
                try {
                    wip->unpackWork( buf.get(), endian );
                    LOGC_DETAIL(wip->job->getLogChannel()) << msg << "   " + wip->print();
                    wip->returnResults();
                } catch ( exception& e ) {
                    LOG_ERR << "putParts:  exception when unpacking results: " << e.what();
                }
            }).detach();
            reply = CMD_OK;         // all ok
            //lock.release();
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
        //LOG_DEBUG << msg;
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
            LOG_ERR << "updateStatus: Exception caught while parsing block: " << e.what();
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
        //LOG_DEBUG << msg;
    }
    totalSize = packedSize + sizeof( uint64_t );
    pack( buf.get(), packedSize );                                                // store real blockSize
    conn->syncWrite( buf.get(), totalSize );

}


void Daemon::sendPeerList( TcpConnection::Ptr& conn ) {

    uint64_t blockSize = myInfo->size();
    unique_lock<mutex> lock( peerMutex );

    for( auto & peer : connections ) {
        if(peer.second && peer.second->info.peerType&Host::TP_WORKER) blockSize += peer.second->size();
    }
    uint64_t totalSize = blockSize + sizeof( uint64_t );                    // blockSize will be sent before block
    shared_ptr<char> buf( new char[totalSize], []( char* p ){ delete[] p; } );
    char* ptr =  buf.get()+sizeof( uint64_t );
    uint64_t packedSize = myInfo->pack( ptr );
    for( auto & peer : connections ) {
        if(peer.second && peer.second->info.peerType&Host::TP_WORKER) packedSize += peer.second->pack( ptr+packedSize );
    }
    if( packedSize > blockSize ) {
        string msg = "sendPeerList(): Packing mismatch:  packedSize = " + to_string(packedSize) + "   blockSize = " + to_string(blockSize) + "  bytes.";
        throw length_error(msg);
        //LOG_DEBUG << msg;
    }
    totalSize = packedSize + sizeof( uint64_t );
    pack( buf.get(), packedSize );                                                // store real blockSize
    conn->syncWrite( buf.get(), totalSize );

}


void Daemon::updateLoadAvg( void ) {
    //LOG_TRACE << "updateLoadAvg()";

    static double loadAvg[3];

    int ret = getloadavg( loadAvg, 3 );
    if( ret != 3 ) {
        LOG_ERR << "updateLoadAvg(): failed to get loadavg.";
        myInfo->status.loadAvg = 0;
        return;
    }

    myInfo->status.loadAvg = loadAvg[0]; // / myInfo->info.nCores * 100.0;
    //LOG_TRACE << "updateLoadAvg()  = " << myInfo->status.loadAvg;
}
