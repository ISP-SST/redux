#include "redux/daemon.hpp"

#include "redux/network/peer.hpp"
#include "redux/network/protocol.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stringutil.hpp"

#include <functional>

#include <boost/date_time/posix_time/posix_time.hpp>

using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;


#define lg Logger::lg
namespace {
    const string thisChannel = "deamon";
}

Daemon::Daemon( po::variables_map& vm ) : Application( vm, LOOP ), master( "" ), params(vm), jobCounter(0), worker(*this) {

    if( params["master"].as<string>() == "" ) {
        server.reset( new TcpServer( ioService, params["port"].as<uint16_t>() ) );
    }
    

}

Daemon::~Daemon( void ) {
    stop();
}

void Daemon::serverInit( void ) {
    if(server) {
        server->setCallback( bind( &Daemon::connected, this, std::placeholders::_1 ) );
        server->accept();
        LOG_DEBUG << "Starting server on port " << params["port"].as<uint16_t>() << ", listening with 5 threads.";
        for( std::size_t i = 0; i < 5; ++i ) {
            shared_ptr<thread> t( new thread( boost::bind( &boost::asio::io_service::run, &ioService ) ) );
            threads.push_back( t );
        }
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
    while( runMode == LOOP ) {
        sleep( 1 );
        cleanupPeers();
    }
    LOG_DETAIL << "exiting maintenance";
}

bool Daemon::doWork( void ) {

    try {
        // start the server
        serverInit();
        // start the maintenance loop
        shared_ptr<thread> t( new thread( boost::bind( &Daemon::maintenance, this ) ) );
        threads.push_back( t );
        //sleep(1);
        worker.init();
        // the io_service will keep running until stopped, then release the threads, so we just let the main thread wait on the join.
        for( auto & it : threads ) {
            it->join();
        }
        threads.clear();
    }
    catch( const exception& e ) {
        LOG_ERR << "doWork()  unhandled exception: " << e.what();
    }

    return true;

}

void Daemon::connected( TcpConnection::ptr conn ) {

    LOG_TRACE << "connected()";

    try {
        *conn << CMD_CFG;           // request handshake
        Peer::HostInfo peerhi;
        *conn >> peerhi;
        *conn << myInfo.host;

        Peer::ptr& peer = addOrGetPeer( peerhi, conn );
        peer->lastSeen = boost::posix_time::second_clock::local_time();

        conn->setCallback( bind( &Daemon::activity, this, std::placeholders::_1 ) );
        *conn << CMD_OK;           // all ok
        conn->idle();
    }
    catch( const exception& e ) {
        LOG_ERR << "connected() Failed to process new connection.";
    }

}

void Daemon::activity( TcpConnection::ptr conn ) {

    Command cmd = CMD_ERR;
    try {
        *conn >> cmd;
    }
    catch( const boost::exception& ) {      // disconnected -> close socket and return.
        conn->socket().close();
        return;
    }

    LOG_TRACE << "activity():  received cmd = " << ( int )cmd << "  (" << bitString( cmd ) << ")";

    switch( cmd ) {
        case CMD_NEW_JOB: addJobs( getPeer( conn ) ); break;
        case CMD_JSTAT: sendJobList( conn ); break;
        case CMD_PSTAT: sendPeerList( conn ); break;
        case CMD_DISCONNECT: *conn << CMD_OK; conn->socket().close(); break;
        default: LOG_DETAIL << "Daemon: Unrecognized command.";
    }

    // conn->idle();

//     try {
//         *conn << CMD_CFG;           // request handshake
//         Peer::HostInfo peerhi;
//         *conn >> peerhi;
//         *conn << myInfo.host;
//
//         Peer& peer = addOrGetPeer( peerhi, conn );
//
//         conn->setCallback( bind( &Daemon::activity, this, std::placeholders::_1 ) );
//
//         cout << "me: " << endl;
//         myInfo.host.print();
//         cout << "peer: " << endl;
//         peer.host.print();
//         cout << "peers.size(): " << peers.size() << endl;
//     }
//     catch( const exception& e ) {
//         LOG_ERR << "connected() Failed to process new connection.";
//     }

}

Peer::ptr& Daemon::addOrGetPeer( const Peer::HostInfo& phi, TcpConnection::ptr& conn ) {
    unique_lock<mutex> lock( peerMutex );
    for( auto & it : peers ) {
        if( it.second->host == phi ) {
            it.second->conn = conn;      // re-connect, replace connection
            return it.second;
        }
    }
    // not found
    size_t id = 1; // start indexing peers with 1
    while( peers.find( id ) != peers.end() ) id++;
    Peer::ptr p( new Peer( phi, conn, id ) );
    auto ret = peers.insert( make_pair( id, p ) );
    if( !ret.second ) {
        throw invalid_argument( "Failed to insert peer." );
    }
    return ret.first->second;
}

Peer::ptr& Daemon::getPeer( const TcpConnection::ptr& conn ) {
    unique_lock<mutex> lock( peerMutex );
    for( auto & it : peers ) {
        if( it.second->conn == conn ) {
            return it.second;
        }
    }
    // not found
    throw invalid_argument( "Failed to get peer." );

}

void Daemon::cleanupPeers( void ) {

    unique_lock<mutex> lock( peerMutex );
    map<size_t, Peer::ptr>::iterator it = peers.begin();
    while( it != peers.end() ) {
        if( !it->second->conn->socket().is_open() ) {
            peers.erase( it++ );
        }
        else {
            ++it;
        }
    }

}

void Daemon::addJobs( Peer::ptr& peer ) {

    TcpConnection::ptr& conn = peer->conn;

    size_t sz = sizeof( size_t );
    unique_ptr<char[]> buf( new char[ sz ] );

    size_t count = boost::asio::read( conn->socket(), boost::asio::buffer( buf.get(), sz ) );
    if( count != sz ) {
        LOG_ERR << "addJobs: Failed to receive blocksize.";
        return;
    }

    bool swapNeeded = ( peer->host.littleEndian != myInfo.host.littleEndian );
    char* ptr = buf.get();
    sz = *reinterpret_cast<size_t*>( ptr );
    if( swapNeeded ) {
        swapEndian( sz );
    }

    buf.reset( new char[ sz ] );
    ptr = buf.get();
    char* end = ptr + sz;
    const char* cptr = nullptr;
    vector<size_t> ids( 1, 0 );
    try {
        count = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, sz ) );
        if( count != sz ) {
            LOG_ERR << "addJobs: Failed to receive data block.";
        }
        else {
            cptr = ptr;
            while( cptr < end ) {
                string tmpS = string( cptr );
                Job::JobPtr job = Job::newJob( tmpS );
                if( job ) {
                    ids[0]++;
                    cptr = job->unpack( cptr, swapNeeded );
                    job->setID( ++jobCounter );
                    job->info.submitTime = boost::posix_time::second_clock::local_time();
                    ids.push_back( jobCounter );
                    unique_lock<mutex> lock( jobMutex );
                    jobs.push_back( job );
                }
                else throw invalid_argument( "Unrecognized Job tag: \"" + tmpS + "\"" );
            }
        }
    }
    catch( const exception& e ) {
        LOG_ERR << "addJobs: Exception caught while parsing block: " << e.what();
    }
    if( cptr == end ) {
        LOG << "Received " << ids[0] << " jobs.";
        *conn << CMD_OK;           // all ok, return IDs
        if( swapNeeded ) swapEndian( ids.data(), ids.size() );
        conn->writeAndCheck( reinterpret_cast<char*>( ids.data() ), ids.size()*sizeof( size_t ) );
    }
    else {
        LOG_ERR << "addJobs: Parsing of datablock failed, there was a missmatch of " << ( cptr - end ) << " bytes.";
        *conn << CMD_ERR;
    }
   
}

void Daemon::sendJobList( TcpConnection::ptr& conn ) {

    *conn << CMD_OK;

    size_t sz = sizeof( size_t );
    for( auto & it : jobs ) {
        sz += it->info.size();
    }
    unique_ptr<char[]> buf( new char[ sz ] );
    char* ptr = buf.get();
    *reinterpret_cast<size_t*>( ptr ) = sz - sizeof( size_t );;
    ptr += sizeof( size_t );
    for( auto & it : jobs ) {
        ptr = it->info.pack( ptr );
    }
    if( buf.get() + sz != ptr ) {
        LOG_ERR << "Packing of job infos failed:  there is a mismatch of " << ( ptr - buf.get() - sz ) << " bytes.";
        return;
    }

    conn->writeAndCheck( buf.get(), sz );

}

void Daemon::sendPeerList( TcpConnection::ptr& conn ) {

    *conn << CMD_OK;

    size_t sz = sizeof( size_t );
    for( auto & it : peers ) {
        sz += it.second->size();
    }
    unique_ptr<char[]> buf( new char[ sz ] );
    char* ptr = buf.get();
    *reinterpret_cast<size_t*>( ptr ) = sz - sizeof( size_t );;
    ptr += sizeof( size_t );
    for( auto & it : peers ) {
        ptr = it.second->pack( ptr );
    }
    if( buf.get() + sz != ptr ) {
        LOG_ERR << "Packing of peers infos failed:  there is a mismatch of " << ( ptr - buf.get() - sz ) << " bytes.";
        return;
    }

    conn->writeAndCheck( buf.get(), sz );

}

