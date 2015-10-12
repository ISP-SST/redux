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
namespace {
    const string thisChannel = "deamon";

    mutex jobMutex;

    typedef boost::asio::time_traits<boost::posix_time::ptime> time_traits_t;

}


Daemon::Daemon( po::variables_map& vm ) : Application( vm, LOOP ), master(""), port(0), params( vm ), jobCounter( 0 ), nQueuedJobs( 0 ),
    hostTimeout(1000), timer( ioService ), worker( *this ) {

    myInfo.reset( new Host() );

    if( params["master"].as<string>() == "" ) {
        server.reset( new TcpServer( ioService, params["port"].as<uint16_t>() ) );
    } else {    // we are connecting to a remote host
        master = params["master"].as<string>();
        port = params["port"].as<uint16_t>();
    }


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
    ioService.stop();
}


void Daemon::maintenance( void ) {

//     LOG_DEBUG << "Maintenance:   nJobs = " << jobs.size() << "  nConn = " << connections.size()
//     << "  nPeers = " << peers.size() << " nPeerWIP = " << peerWIP.size() << " nThreads = " << threads.size();

    updateLoadAvg();
    cleanup();
    worker.updateStatus();

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
        for( std::size_t i = 0; i < 5; ++i ) {
            shared_ptr<thread> t( new thread( boost::bind( &boost::asio::io_service::run, &ioService ) ) );
            threads.push_back( t );
        }
        LOG_TRACE << "Initializing worker.";
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


void Daemon::connected( TcpConnection::Ptr conn ) {

    try {
        *conn << CMD_CFG;           // request handshake
        Host::HostInfo remote_info;
        *conn >> remote_info;
        *conn << myInfo->info;

        LOG_DEBUG << "connected()  Handshake successful.";
        addConnection( remote_info, conn );

        conn->setCallback( bind( &Daemon::activity, this, std::placeholders::_1 ) );
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
            case CMD_ADD_JOB: conn->strand.post( boost::bind( &Daemon::addJobs, this, conn ) ); break;
            case CMD_DEL_JOB: conn->strand.post( boost::bind( &Daemon::removeJobs, this, conn ) ); break;
            case CMD_GET_WORK: conn->strand.post( boost::bind( &Daemon::sendWork, this, conn ) ); break;
            case CMD_GET_JOBLIST: conn->strand.post( boost::bind( &Daemon::sendJobList, this, conn ) ); break;
            case CMD_PUT_PARTS: conn->strand.post( boost::bind( &Daemon::putParts, this, conn ) ); break;
            case CMD_STAT: conn->strand.post( boost::bind( &Daemon::updateHostStatus, this, conn ) ); break;
            case CMD_JSTAT: conn->strand.post( boost::bind( &Daemon::sendJobStats, this, conn ) ); break;
            case CMD_PSTAT: conn->strand.post( boost::bind( &Daemon::sendPeerList, this, conn ) ); break;
            case CMD_DISCONNECT: conn->socket().close(); break;
            default: LOG_DETAIL << "Daemon: Unrecognized command: " << ( int )cmd << "  " << bitString( cmd );
        }

    }
    catch( const exception& e ) {
        LOG_ERR << "activity(): Exception when parsing incoming command: " << e.what();
    }

    conn->strand.post( boost::bind( &TcpConnection::idle, conn ) );

}


void Daemon::addConnection( const Host::HostInfo& remote_info, TcpConnection::Ptr& conn ) {
    
    unique_lock<mutex> lock( peerMutex );
    if( myInfo->info == remote_info ) {
        myInfo->info.peerType |= remote_info.peerType;
        connections[conn] = myInfo;
        return;
    }
    
    auto& host = connections[conn];
    if( host == nullptr ) { // not found
        uint64_t id = 1; // start indexing peers with 1
        while( peers.find( id ) != peers.end() ) id++;
        host.reset( new Host( remote_info, id ) );
        host->nConnections++;
        peers.insert( make_pair( id, host ) );
        conn->setSwapEndian(remote_info.littleEndian != myInfo->info.littleEndian);
    } else {
        host->touch();
        host->nConnections++;
    }


}


void Daemon::removeConnection( TcpConnection::Ptr& conn ) {
    
    unique_lock<mutex> lock( peerMutex );
    conn->socket().close();
    auto& host = connections[conn];
    if( host ) {
        LOG << "removeConnection(): nConn = " << host->nConnections;
        if(--host->nConnections == 0) {
            connections.erase(conn);
        }
    } else connections.erase(conn);
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

        for( auto it=peers.begin(); it != peers.end(); ) {
            if( it->second->nConnections < 1 ) {
                auto wip = peerWIP.find( it->second );
                if( wip != peerWIP.end() ) {
                    if( wip->second.job ) {
                        LOG_DETAIL << "Peer #" << it->first << " disconnected, returning unfinished work to queue: " << wip->second.print();
                        wip->second.job->ungetWork( wip->second );
                    }
                    peerWIP.erase(wip);
                }
                peers.erase( it++ );
            } else ++it;
        }
        
        boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
        for( auto it=peerWIP.begin(); it != peerWIP.end(); ) {
            Job::JobPtr job = it->second.job;
            if( job ) {
                bool timedOut(false);
                boost::posix_time::time_duration elapsed = (now - it->first->status.lastSeen);
                if( elapsed > boost::posix_time::seconds( hostTimeout ) ) {
                    LOG_DETAIL << "Peer has not been active in " << to_simple_string(elapsed) << ", returning unfinished work to queue: " << it->second.print();
                    it->second.job->ungetWork( it->second );
                    timedOut = true;
                }
                elapsed = (now - it->second.workStarted);
                if( elapsed > boost::posix_time::seconds( it->second.job->info.timeout ) ) {
                    LOG_DETAIL << "Work has not been completed in " << to_simple_string(elapsed) << ", returning unfinished work to queue: " << it->second.print();
                    it->second.job->ungetWork( it->second );
                    timedOut = true;
                }
                
                if( timedOut ) {
                    peerWIP.erase(it++);
                    continue;
                }
            }
            ++it;
        }
    }
    
/*    {
        unique_lock<mutex> lock( jobMutex );

        for( auto& job : jobs ) {
            if( job->info.step >= 32 ) {
                LOG_DEBUG << "Job " << job->info.id << " (" << job->info.name << ") is completed, removing from queue.  uc=" << job.use_count() << endl;
                job.reset();
            }
        }
        jobs.erase(std::remove_if(jobs.begin(), jobs.end(), [](const shared_ptr<Job>& j){return (j == nullptr);}), jobs.end());
    }
*/
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
                    if( !job->check() ) throw job_check_failed( "Sanity check failed for \"" + tmpS + "\"-job #" + to_string(nJobs) );
                    unique_lock<mutex> lock( jobMutex );
                    job->info.id = ++jobCounter;
                    job->info.name = "job_" + to_string( job->info.id );
                    job->info.submitTime = boost::posix_time::second_clock::local_time();
                    ids.push_back( jobCounter );
                    ids[0]++;
                    jobs.push_back( job );
                } catch( const job_check_failed& e ) {
                    messages.push_back( e.what() );
                }
            } else throw invalid_argument( "Unrecognized Job tag: \"" + tmpS + "\"" );
        }
    } catch( const exception& e ) {
        LOG_ERR << "addJobs: Exception caught while parsing block: " << e.what();
    } catch( ... ) {
        LOG_ERR << "addJobs: Unrecognized exception thrown when parsing cfg file.";
    }

    if( count == blockSize ) {
        if( ids[0] ) LOG << "Received " << ids[0] << " jobs." << printArray(ids.data()+1,ids[0],"  IDs");
        *conn << CMD_OK;           // all ok, return IDs
        conn->writeAndCheck( ids );
    } else {
        LOG_ERR << "addJobs: Parsing of datablock failed, count = " << count << "   blockSize = " << blockSize << "  bytes.";
        *conn << CMD_ERR;
    }

    // send any messages back
    uint64_t messagesSize(0);
    if( messages.size() ) {
        messagesSize = redux::util::size( messages );
        shared_ptr<char> tmp = sharedArray<char>(messagesSize+sizeof(uint64_t));
        char* ptr = tmp.get();
        ptr += redux::util::pack( ptr, messagesSize );
        redux::util::pack( ptr, messages );
        conn->writeAndCheck( tmp, messagesSize+sizeof(uint64_t) );
    } else {
        conn->writeAndCheck( messagesSize );
    }



}


void Daemon::removeJobs( TcpConnection::Ptr& conn ) {

    size_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if( blockSize ) {

        string jobString = string( buf.get() );
        if( iequals( jobString, "all" ) ) {
            unique_lock<mutex> lock( jobMutex );
            if( jobs.size() ) LOG << "Clearing joblist.";
            jobs.clear();
            stop();
            return;
        }

        try {
            bpt::ptree tmpTree;      // just to be able to use the VectorTranslator
            tmpTree.put( "jobs", jobString );
            vector<size_t> jobList = tmpTree.get<vector<size_t>>( "jobs", vector<size_t>() );
            int cnt = 0;
            unique_lock<mutex> lock( jobMutex );
            vector<size_t> ids; 
            for( auto & jobId : jobList ) {
                for( auto it2 = jobs.begin(); it2 != jobs.end(); ) {
                    if( ( *it2 )->info.id == jobId ) {
                        ids.push_back(jobId);
                        jobs.erase( it2++ );
                        cnt++;
                    } else {
                        ++it2;
		    }
                }
            }
            if( cnt ) LOG << "Removed " << printArray(ids,"jobs");
            return;
        }
        catch( const boost::bad_lexical_cast& e ) {
            // ignore
	}

        try {
            vector<string> jobList;
            boost::split( jobList, jobString, boost::is_any_of( "," ) );
            int cnt = 0;
            unique_lock<mutex> lock( jobMutex );
            vector<string> deletedList;
            for( auto & jobId : jobList ) {
                for( auto it2 = jobs.begin(); it2 != jobs.end(); ) {
                    if( ( *it2 )->info.name == jobId ) {
                        deletedList.push_back(jobId);
                        jobs.erase( it2++ );
                        cnt++;
                    } else {
                        ++it2;
		    }
                }
            }
            if( cnt ) LOG << "Removed " << printArray(deletedList,"jobs");
            return;
        }
        catch( const boost::bad_lexical_cast& e ) {
            // ignore
	}


    }

}

/*Job::JobPtr Daemon::selectJob( bool localRequest ) {
    for( auto & job : jobs )
        if( it->info.step.load() < Job::JSTEP_RUNNING )
	    return job;
}*/

bool Daemon::getWork( WorkInProgress& wip, uint8_t nThreads ) {

    unique_lock<mutex> lock( jobMutex );
    wip.previousJob = wip.job;
    wip.job.reset();
    // TODO: sort by priority
    for( auto& job : jobs ) {
        if( job->check() && job->getWork( wip, nThreads ) ) {
            wip.job = job;
            for( auto& part: wip.parts ) {
                part->cacheLoad(false);
            }
            wip.workStarted = boost::posix_time::second_clock::local_time();
            break;
        }
    }
    /*if( wip.connection ) {                // remote worker
        for( auto & it : jobs ) {
            uint8_t step = it->info.step.load();
            if( step == Job::JSTEP_QUEUED ) {
                if( it->getWork( wip, nThreads ) ) {
                    ret = true;
                }
            }
            else if( step == Job::JSTEP_RUNNING ) {
                if( it->getWork( wip, nThreads ) ) {
                    ret = true;
                }
            }
            if ( ret ) {
                job = it;
                break;
            }
        }
    }
    else {                          // local worker
        for( auto & it : jobs ) {
            uint8_t step = it->info.step.load();
            if( (nQueuedJobs < 2) && (step < Job::JSTEP_QUEUED) ) {
                it->getWork( wip, nThreads );
                ret = true;
            } else if( step & ( Job::JSTEP_QUEUED | Job::JSTEP_RUNNING | Job::JSTEP_POSTPROCESS ) ) {
                if( it->getWork( wip, nThreads ) || ( step != Job::JSTEP_RUNNING ) ) {
                    ret = true;
                }
            }
            if ( ret ) {
                job = it;
                break;
            }
        }
    }*/
    
    if( wip.job ) {
        //cout << "getWork("<<(int)nThreads<<"): Handing out " << wip.parts.size() << " part(s) from job #" << wip.job->info.id << endl;
        wip.job->info.state.store( Job::JSTATE_ACTIVE );
    }
    
    return bool(wip.job);
}


void Daemon::returnResults( WorkInProgress& wip ) {
    
    wip.job->returnResults( wip );
    for( auto& part : wip.parts ) {
        part->cacheStore(true);       // store and free some resources. (if implemented for this Part-derivative)
    }
    wip.parts.clear();

}


void Daemon::sendWork( TcpConnection::Ptr& conn ) {

    auto h = connections[conn];
    
    map<Host::Ptr, WorkInProgress>::iterator it = peerWIP.find(h);
    if( it == peerWIP.end() ) {
        auto rit = peerWIP.emplace( connections[conn], WorkInProgress( conn ) );
        if( !rit.second ) {
            LOG_ERR << "Failed to insert new connection into peerWIP";
        }
        it = rit.first;
    }
    WorkInProgress& wip = it->second;
    auto host = it->first;

    uint64_t blockSize = 0;
    bool includeJob = false;
    shared_ptr<char> data;
    if( getWork( wip, host->status.nThreads ) ) {
        blockSize += wip.workSize();
       // if(!previousJob || *wip.job != *previousJob) {
       //     includeJob = true;
       // }
        LOG_DEBUG << "Sending work" << (includeJob?"(new job) ":"") << " to host:pid = " << host->info.name << ":" << host->info.pid;
       // blockSize = wip.size( includeJob );
    }

    data = sharedArray<char>(blockSize+sizeof(uint64_t));
    char* ptr = data.get()+sizeof(uint64_t);
    uint64_t count(0); 
    if(wip.job) {
        count += wip.packWork( ptr+count );
        for(auto& part: wip.parts) {
            part->cacheStore(true);
        }
        wip.previousJob = wip.job;
    }
    pack( data.get(), count );         // Store actual packed bytecount (something might be compressed)
    conn->writeAndCheck( data, count+sizeof(uint64_t) );

    
/*    auto buf = sharedArray<char>( blockSize + sizeof( size_t ) );
    char* ptr = buf.get();
    uint64_t count = sizeof( size_t );
    if( blockSize > 0 ) {
        count += wip.pack( ptr+count, previousJob );
    }
*/

}


void Daemon::putParts( TcpConnection::Ptr& conn ) {

    size_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    map<Host::Ptr, WorkInProgress>::iterator it = peerWIP.find(connections[conn]);
    if( it == peerWIP.end() ) {
        LOG_ERR << "Received results from unknown host.";
    }

    
    WorkInProgress& wip = it->second;
    if( blockSize ) {
        wip.unpackWork( buf.get(), conn->getSwapEndian() );
        returnResults( wip );
    } else {
        LOG_DEBUG << "putParts():  EMPTY   blockSize = " << blockSize;
    }
    
    *conn << CMD_OK;         // all ok

}


void Daemon::sendJobList( TcpConnection::Ptr& conn ) {

    uint64_t blockSize = 0;
    unique_lock<mutex>( jobMutex );
    for( auto & job : jobs ) {
        blockSize += job->size();
    }
    auto buf = sharedArray<char>( blockSize + sizeof( uint64_t ) );         // blockSize will be sent before block
    char* ptr = buf.get();
    uint64_t count = pack( ptr, blockSize );                                // store blockSize
    blockSize += sizeof( uint64_t );                                        // ...and add sizeof(blockSize) to make verification & write correct below.
    for( auto & job : jobs ) {
        count += job->pack( ptr+count );
    }
    if( count != blockSize ) {
        LOG_ERR << "sendJobList(): Mismatch when packing joblist:  count = " << count << "   blockSize = " << blockSize << "  bytes.";
        return;
    }
    conn->writeAndCheck( buf, blockSize );

}


void Daemon::updateHostStatus( TcpConnection::Ptr& conn ) {

    size_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if( blockSize ) {
        try {
            connections[conn]->status.unpack( buf.get(), conn->getSwapEndian() );
        }
        catch( const exception& e ) {
            LOG_ERR << "updateStatus: Exception caught while parsing block: " << e.what();
        }
    }

}


void Daemon::sendJobStats( TcpConnection::Ptr& conn ) {

    size_t blockSize = 0;
    unique_lock<mutex>( jobMutex );
    for( auto & job : jobs ) {
        blockSize += job->info.size();
    }
    size_t totalSize = blockSize + sizeof( size_t );
    auto buf = sharedArray<char>( totalSize );
    char* ptr = buf.get();
    uint64_t count = pack( ptr, blockSize );
    for( auto & job : jobs ) {
        count += job->info.pack( ptr+count );
    }
    if( count != totalSize ) {
        LOG_ERR << "sendJobStats(): Packing of job infos failed:  count = " << count << "  totalSize = " << totalSize << " bytes.";
        return;
    }

    conn->writeAndCheck( buf, totalSize );

}


void Daemon::sendPeerList( TcpConnection::Ptr& conn ) {

    size_t blockSize = myInfo->size();
    for( auto & peer : peers ) {
        blockSize += peer.second->size();
    }
    size_t totalSize = blockSize + sizeof( size_t );
    auto buf = sharedArray<char>( totalSize );
    char* ptr =  buf.get();
    uint64_t count = pack(ptr, blockSize );
    count += myInfo->pack( ptr+count );
    for( auto & peer : peers ) {
        count += peer.second->pack( ptr+count );
    }
    if( count != totalSize ) {
        LOG_ERR << "Packing of peers infos failed:  count = " << count << "   totalSize = " << totalSize << "  bytes.";
        return;
    }

    conn->writeAndCheck( buf, totalSize );

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

    myInfo->status.loadAvg = loadAvg[0] / myInfo->info.nCores * 100.0;
    //LOG_TRACE << "updateLoadAvg()  = " << myInfo->status.loadAvg;
}
