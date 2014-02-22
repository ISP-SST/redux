#include "redux/daemon.hpp"

#include "redux/network/peer.hpp"
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


Daemon::Daemon( po::variables_map& vm ) : Application( vm, LOOP ), master( "" ), params( vm ), jobCounter( 0 ), nQueuedJobs( 0 ), timer( ioService ), worker( *this ) {

    myInfo.reset( new Peer() );

    if( params["master"].as<string>() == "" ) {
        server.reset( new TcpServer( ioService, params["port"].as<uint16_t>() ) );
    }


}


Daemon::~Daemon( void ) {
    stop();
}


void Daemon::serverInit( void ) {

    LOG << "serverInit()";
    if( server ) {
        server->setCallback( bind( &Daemon::connected, this, std::placeholders::_1 ) );
        server->accept();
        myInfo->host.peerType |= Peer::PEER_MASTER;
        LOG_DEBUG << "Starting server on port " << params["port"].as<uint16_t>() << ".";
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

    LOG_TRACE << "Maintenance";

    updateLoadAvg();
    cleanupPeers();

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
        for( auto & it : threads ) {
            it->join();
        }
        myInfo->host.peerType = 0;
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
        Peer::HostInfo peerhi;
        *conn >> peerhi;
        *conn << myInfo->host;

        LOG_DEBUG << "connected()  Handshake successful.";
        Peer::Ptr& peer = addOrGetPeer( peerhi, conn );
        peer->lastSeen = boost::posix_time::second_clock::local_time();

        conn->setCallback( bind( &Daemon::activity, this, std::placeholders::_1 ) );
        *conn << CMD_OK;           // all ok
        conn->idle();
    }
    catch( const exception& e ) {
        LOG_ERR << "connected() Failed to process new connection.";
    }

}


void Daemon::activity( TcpConnection::Ptr conn ) {

    Command cmd = CMD_ERR;
    try {
        *conn >> cmd;
    }
    catch( const boost::exception& ) {      // disconnected -> close socket and return.
        conn->socket().close();
        return;
    }

    LOG_TRACE << "activity():  received cmd = " << ( int )cmd << "  (" << bitString( cmd ) << ")";
    
    try {
        Peer::Ptr peer = getPeer( conn );
        peer->lastSeen = boost::posix_time::second_clock::local_time();
        
        switch( cmd ) {
            case CMD_ADD_JOB: conn->strand.post( boost::bind( &Daemon::addJobs, this, peer ) ); break;
            case CMD_DEL_JOB: conn->strand.post( boost::bind( &Daemon::removeJobs, this, peer ) ); break;
            case CMD_GET_WORK: conn->strand.post( boost::bind( &Daemon::sendWork, this, peer ) ); break;
            case CMD_GET_JOBLIST: conn->strand.post( boost::bind( &Daemon::sendJobList, this, conn ) ); break;
            case CMD_PUT_PARTS: conn->strand.post( boost::bind( &Daemon::putParts, this, peer ) ); break;
            case CMD_STAT: conn->strand.post( boost::bind( &Daemon::updatePeerStatus, this, peer ) ); break;
            case CMD_JSTAT: conn->strand.post( boost::bind( &Daemon::sendJobStats, this, conn ) ); break;
            case CMD_PSTAT: conn->strand.post( boost::bind( &Daemon::sendPeerList, this, conn ) ); break;
            case CMD_DISCONNECT: conn->socket().close(); break;
            default: LOG_DETAIL << "Daemon: Unrecognized command: " << (int)cmd << "  " << bitString(cmd);
        }

    }
    catch( const exception& e ) {
        LOG_ERR << "activity(): Exception when parsing incoming command: " << e.what();
    }

    conn->strand.post( boost::bind( &TcpConnection::idle, conn ) );

}


Peer::Ptr& Daemon::addOrGetPeer( const Peer::HostInfo& phi, TcpConnection::Ptr& conn ) {
    unique_lock<mutex> lock( peerMutex );
    if( myInfo->host == phi ) {
        myInfo->host.peerType |= phi.peerType;
        myInfo->conn = conn;
        return myInfo;
    }
    for( auto & it : peers ) {
        if( it.second->host == phi ) {
            myInfo->host.peerType = phi.peerType;
            it.second->conn = conn;      // re-connect, replace connection
            return it.second;
        }
    }
    // not found
    size_t id = 1; // start indexing peers with 1
    while( peers.find( id ) != peers.end() ) id++;
    Peer::Ptr p( new Peer( phi, conn, id ) );
    auto ret = peers.insert( make_pair( id, p ) );
    if( !ret.second ) {
        throw invalid_argument( "Failed to insert peer." );
    }
    return ret.first->second;
}


Peer::Ptr& Daemon::getPeer( const TcpConnection::Ptr& conn ) {
    unique_lock<mutex> lock( peerMutex );
    if( myInfo->conn == conn ) {
        return myInfo;
    }
    for( auto & it : peers ) {
        if( it.second->conn == conn ) {
            return it.second;
        }
    }
    // not found
    throw invalid_argument( "Failed to get peer." );

}

// TODO a timeout to allow a peer to reconnect before "hard-dropping" it.
void Daemon::cleanupPeers( void ) {

    unique_lock<mutex> lock( peerMutex );

    for( auto & it : peers ) {
        if( it.second && it.second->conn && !it.second->conn->socket().is_open() ) {
            auto wip = peerWIP.find( it.second );
            if( wip != peerWIP.end() ) {
                if( wip->second.job ) {
                    LOG_DETAIL << "Peer disconnected, returning work to queue: " << wip->second.print();
                    wip->second.job->ungetParts( wip->second );
                }
            }
            peers.erase( it.first );
        }
    }

}


void Daemon::addJobs( Peer::Ptr& peer ) {

    size_t blockSize;
    bool swap_endian;
    shared_ptr<char> buf = peer->receiveBlock( blockSize, swap_endian );
    if( !blockSize ) {
        *peer->conn << CMD_ERR;
        return;
    }

    const char* cptr = buf.get();
    char* end = buf.get() + blockSize;
    vector<size_t> ids( 1, 0 );
    try {
        while( cptr < end ) {
            string tmpS = string( cptr );
            Job::JobPtr job = Job::newJob( tmpS );
            if( job ) {
                cptr = job->unpack( cptr, swap_endian );
                unique_lock<mutex> lock( jobMutex );
                job->info.id = ++jobCounter;
                job->info.step.store( Job::JSTEP_RECEIVED );
                job->info.name = "job_" + to_string( job->info.id );
                job->info.submitTime = boost::posix_time::second_clock::local_time();
                ids.push_back( jobCounter );
                ids[0]++;
                jobs.push_back( job );
            }
            else throw invalid_argument( "Unrecognized Job tag: \"" + tmpS + "\"" );
        }
    }
    catch( const exception& e ) {
        LOG_ERR << "addJobs: Exception caught while parsing block: " << e.what();
    }

    if( cptr == end ) {
        LOG << "Received " << ids[0] << " jobs.";
        *peer->conn << CMD_OK;           // all ok, return IDs
        if( swap_endian ) swapEndian( ids.data(), ids.size() );
        peer->conn->writeAndCheck( ids );
    }
    else {
        LOG_ERR << "addJobs: Parsing of datablock failed, there was a missmatch of " << ( cptr - end ) << " bytes.";
        *peer->conn << CMD_ERR;
    }

}


void Daemon::removeJobs( Peer::Ptr& peer ) {

    size_t blockSize;
    bool swap_endian;
    shared_ptr<char> buf = peer->receiveBlock( blockSize, swap_endian );

    if( blockSize ) {

        string jobString = string( buf.get() );
        if( iequals( jobString, "all" ) ) {
            unique_lock<mutex> lock( jobMutex );
            if( jobs.size() ) LOG << "Clearing joblist.";
            jobs.clear();
            return;
        }

        try {
            bpt::ptree tmpTree;      // just to be able to use the VectorTranslator
            tmpTree.put( "jobs", jobString );
            vector<size_t> jobList = tmpTree.get<vector<size_t>>( "jobs", vector<size_t>() );
            int cnt = 0;
            unique_lock<mutex> lock( jobMutex );
            for( auto & it : jobList ) {
                for( auto it2 = jobs.begin(); it2 < jobs.end(); ++it2 ) {
                    if( ( *it2 )->info.id == it ) {
                        jobs.erase( it2 );
                        cnt++;
                    }
                }
            }
            if( cnt ) LOG << "Removed jobs: " << cnt << " jobs by id.";
            return;
        }
        catch( ... ) {}      // catch and ignore bad_lexical_cast

        try {
            vector<string> jobList;
            boost::split( jobList, jobString, boost::is_any_of( "," ) );
            int cnt = 0;
            unique_lock<mutex> lock( jobMutex );
            for( auto & it : jobList ) {
                for( auto it2 = jobs.begin(); it2 < jobs.end(); ++it2 ) {
                    if( ( *it2 )->info.name == it ) {
                        jobs.erase( it2 );
                        cnt++;
                    }
                }
            }
            if( cnt ) LOG << "Removed " << cnt << " jobs by name";
            return;
        }
        catch( ... ) {}      // catch and ignore bad_lexical_cast


    }

}

Job::JobPtr Daemon::selectJob( bool localRequest ) {
    Job::JobPtr job;
    for( auto & it : jobs ) {
        if( it->info.step.load() < Job::JSTEP_RUNNING ) {
            job = it;
            break;
        }
    }
    return job;
}

bool Daemon::getWork( WorkInProgress& wip ) {

    if( wip.job ) {        // job already checked out. lost ??      TODO: deal with it...
        wip.job.reset();
    }

    unique_lock<mutex> lock( jobMutex );
    bool ret = false;
    if( wip.peer ) {                // remote worker
        for( auto & it : jobs ) {
            uint8_t step = it->info.step.load();
            if( step == Job::JSTEP_QUEUED ) {
                wip.job = it;
                it->getParts( wip );
                ret = true;
                break;
            }
            else if( step == Job::JSTEP_RUNNING ) {
                wip.job = it;
                if( it->getParts( wip ) ) {
                    ret = true;
                    break;
                }
            }
        }
    }
    else {                          // local worker
        if( nQueuedJobs < 2 ) {     // see if another job can be prepared
            for( auto & it : jobs ) {
                if( it->info.step.load() < Job::JSTEP_QUEUED ) {
                    wip.job = it;
                    it->getParts( wip );
                    ret = true;
                    break;
                }
            }
        }
        for( auto & it : jobs ) {
            uint8_t step = it->info.step.load();
            if( step & ( Job::JSTEP_QUEUED | Job::JSTEP_RUNNING | Job::JSTEP_POSTPROCESS ) ) {
                wip.job = it;
                if( it->getParts( wip ) || ( step != Job::JSTEP_RUNNING ) ) {
                    ret = true;
                    break;
                }
            }
        }
        if(ret) wip.peer = myInfo;
    }
    if( wip.job ) {
        wip.job->info.state.store( Job::JSTATE_ACTIVE );
    }
    return ret;
}


void Daemon::sendWork( Peer::Ptr& peer ) {

    auto it = peerWIP.insert( std::pair<Peer::Ptr, WorkInProgress>( peer, WorkInProgress( peer ) ) );
    WorkInProgress& pw = it.first->second;

    size_t blockSize = 0;
    bool includeJob = true;
    if( getWork( pw ) ) {
        blockSize = pw.size( includeJob );
    }

    size_t totalSize = blockSize + sizeof( size_t );
    auto buf = sharedArray<char>( totalSize );
    char* ptr = pack( buf.get(), blockSize );

    if( blockSize > 0 ) {
        ptr = pw.pack( ptr, includeJob );
        if( pw.job ) pw.job->info.step.store( Job::JSTEP_RUNNING );
    }

    peer->conn->writeAndCheck( buf, totalSize );

}


void Daemon::putParts( Peer::Ptr& peer ) {

    size_t blockSize;
    bool swap_endian;
    shared_ptr<char> buf = peer->receiveBlock( blockSize, swap_endian );

    auto it = peerWIP.insert( std::pair<Peer::Ptr, WorkInProgress>( peer, WorkInProgress( peer ) ) );
    WorkInProgress& wip = it.first->second;
    if( blockSize ) {
        LOG_DEBUG << "putParts():  blockSize = " << blockSize;
        LOG_DEBUG << "putParts():  wip = " << wip.print();
        wip.unpack( buf.get(), false );
        wip.job->returnParts( wip );
        LOG_DEBUG << "putParts(): wip2 = " << wip.print();
    }
    else LOG_DEBUG << "putParts():  EMPTY   blockSize = " << blockSize;

        *(peer->conn) << CMD_OK;           // all ok

}


void Daemon::sendJobList( TcpConnection::Ptr& conn ) {

    size_t blockSize = 0;
    unique_lock<mutex>( jobMutex );
    for( auto & it : jobs ) {
        blockSize += it->size();
    }
    size_t totalSize = blockSize + sizeof( size_t );
    auto buf = sharedArray<char>( totalSize );
    char* ptr = pack( buf.get(), blockSize );
    for( auto & it : jobs ) {
        ptr = it->pack( ptr );
    }
    if( buf.get() + totalSize != ptr ) {
        LOG_ERR << "Packing of job infos failed:  there is a mismatch of " << ( ptr - buf.get() - totalSize ) << " bytes.";
        return;
    }

    conn->writeAndCheck( buf, totalSize );

}


void Daemon::updatePeerStatus( Peer::Ptr& peer ) {

    size_t blockSize;
    bool swap_endian;
    shared_ptr<char> buf = peer->receiveBlock( blockSize, swap_endian );

    if( !blockSize ) return;

    try {
        peer->stat.unpack( buf.get(), swap_endian );
    }
    catch( const exception& e ) {
        LOG_ERR << "updateStatus: Exception caught while parsing block: " << e.what();
    }

}


void Daemon::sendJobStats( TcpConnection::Ptr& conn ) {

    size_t blockSize = 0;
    unique_lock<mutex>( jobMutex );
    for( auto & it : jobs ) {
        blockSize += it->info.size();
    }
    size_t totalSize = blockSize + sizeof( size_t );
    auto buf = sharedArray<char>( totalSize );
    char* ptr = pack( buf.get(), blockSize );
    for( auto & it : jobs ) {
        ptr = it->info.pack( ptr );
    }
    if( buf.get() + totalSize != ptr ) {
        LOG_ERR << "Packing of job infos failed:  there is a mismatch of " << ( ptr - buf.get() - totalSize ) << " bytes.";
        return;
    }

    conn->writeAndCheck( buf, totalSize );

}


void Daemon::sendPeerList( TcpConnection::Ptr& conn ) {

    size_t blockSize = myInfo->size();
    for( auto & it : peers ) {
        blockSize += it.second->size();
    }
    size_t totalSize = blockSize + sizeof( size_t );
    auto buf = sharedArray<char>( totalSize );
    char* ptr = pack( buf.get(), blockSize );
    ptr = myInfo->pack( ptr );
    for( auto & it : peers ) {
        ptr = it.second->pack( ptr );
    }
    if( buf.get() + totalSize != ptr ) {
        LOG_ERR << "Packing of peers infos failed:  there is a mismatch of " << ( ptr - buf.get() - totalSize ) << " bytes.";
        return;
    }

    conn->writeAndCheck( buf, totalSize );

}


void Daemon::updateLoadAvg( void ) {
    LOG_TRACE << "updateLoadAvg()";

    static double loadAvg[3];

    int ret = getloadavg( loadAvg, 3 );
    if( ret != 3 ) {
        LOG_ERR << "updateLoadAvg(): failed to get loadavg.";
        myInfo->stat.loadAvg = 0;
        return;
    }

    myInfo->stat.loadAvg = loadAvg[0] / myInfo->host.nCores * 100.0;
}
