#include "redux/worker.hpp"

#include "redux/daemon.hpp"
#include "redux/logger.hpp"
#include "redux/network/protocol.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"


using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {
    const std::string thisChannel = "worker";

    typedef boost::asio::time_traits<boost::posix_time::ptime> time_traits_t;
}


Worker::Worker( Daemon& d ) : strand(d.ioService), runTimer( d.ioService ), daemon( d ) {

}


Worker::~Worker( void ) {

}


void Worker::init( void ) {

    master.reset( new Peer() );
    master->conn = TcpConnection::newPtr( daemon.ioService );
    connect();
    runTimer.expires_at( time_traits_t::now() + boost::posix_time::seconds( 1 ) );
    runTimer.async_wait(strand.wrap(boost::bind(&Worker::run, this)));

}

void Worker::connect( void ) {

    master->conn->connect( daemon.params["master"].as<string>(), to_string( daemon.params["port"].as<uint16_t>() ) );
    if( master->conn->socket().is_open() ) {
        Command cmd;

        daemon.myInfo->host.peerType |= Peer::PEER_WORKER;
        *master->conn << CMD_CONNECT;
        *master->conn >> cmd;
        if( cmd == CMD_AUTH ) {
            // implement
        }
        if( cmd == CMD_CFG ) {  // handshake requested
            *master->conn << daemon.myInfo->host;
            *master->conn >> master->host;
            *master->conn >> cmd;       // ok or err
        }
        if( cmd != CMD_OK ) {
            LOG_ERR << "Handshake with master failed  (server replied: " << cmd << ")";
            master->conn->socket().close();
            daemon.myInfo->host.peerType &= ~Peer::PEER_WORKER;
        }

    }

}

void Worker::stop( void ) {
    if( master->conn->socket().is_open() ) {
        *master->conn << CMD_DISCONNECT;
        master->conn->socket().close();
    }
    ioService.stop();
    daemon.myInfo->host.peerType &= ~Peer::PEER_WORKER;

}


void Worker::updateStatus( void ) {

    if( !master->conn->socket().is_open() ) {
        connect();
    }

    if( master->conn->socket().is_open() ) {
        size_t blockSize = daemon.myInfo->stat.size();
        size_t totSize = blockSize + sizeof( size_t ) + 1;
        auto buf = sharedArray<char>(totSize);
        char* ptr = buf.get();
        memset( ptr, 0, totSize );

        uint64_t count = pack( ptr, CMD_STAT );
        count += pack( ptr+count, blockSize );
        count += daemon.myInfo->stat.pack( ptr+count );

        master->conn->writeAndCheck( buf, totSize );
    }
    else {
        LOG_WARN << "No connection to master.";
    }

}

bool Worker::fetchWork( void ) {
    LOG_DEBUG << "fetchWork() ";
    try {

        if( !master->conn->socket().is_open() ) {
            connect();
        }

        if( !master->conn->socket().is_open() ) {
            LOG_DEBUG << "fetchWork()  no socket, returning.";
            return false;
        }

        *master->conn << CMD_GET_WORK;

        size_t blockSize;
        bool swap_endian;
        auto buf = master->receiveBlock( blockSize, swap_endian );               // reply

        if( !blockSize ) return false;

        const char* ptr = buf.get();

        uint64_t count = wip.unpack( ptr, swap_endian );

        if( count != blockSize ) {
            throw invalid_argument( "Parsing of datablock failed, there was a missmatch of " + to_string( count - blockSize ) + " bytes." );

        }

        wip.peer = master;

        return true;
    }
    catch( const exception& e ) {
        LOG_ERR << "fetchWork: Exception caught while fetching job " << e.what() << endl;
    }
    catch( ... ) {
        LOG_ERR << "fetchWork: Unrecognized exception caught while fetching job." << endl;
    }

    return false;
}



bool Worker::getWork( void ) {

    Job::JobPtr lastJob = wip.job;

    if( wip.peer && (wip.peer != daemon.myInfo)) {            // remote work: return parts.
        returnWork();
    }

    wip.peer.reset();

    if( daemon.getWork( wip ) || fetchWork() ) {    // first check for local work, then remote
        if( !lastJob || (wip.job && *(wip.job) != *lastJob) ) {
            //LOG_DETAIL << "Different job !!!   " + wip.print();
            wip.job->init();
            if( lastJob ) {
                lastJob->cleanup();
            }
        } //else LOG_DETAIL << "Same job !!!   " + wip.print();
        LOG_DETAIL << "Got work: " + wip.print();
        return true;
    }
    else LOG_TRACE << "No work";

    return false;

}



void Worker::returnWork( void ) {

    if( wip.parts.size() ) {

        try {

            if( !wip.peer->conn->socket().is_open() ) {
                return;     // TODO handle reconnects
            }

            LOG_DETAIL << "Returning result: " + wip.print();

            bool includeJob = false;
            size_t blockSize = wip.size(includeJob);

            size_t totalSize = blockSize + sizeof( size_t ) + 1;
            auto buf = sharedArray<char>( totalSize );        // + blocksize + cmd

            char* ptr = buf.get();
            uint64_t count = pack( ptr, CMD_PUT_PARTS );
            count += pack( ptr+count, blockSize );
            if(blockSize) {
                count += wip.pack( ptr+count, includeJob );

            }

            master->conn->writeAndCheck( buf, totalSize );

            Command cmd = CMD_ERR;

            *(master->conn) >> cmd;

        }
        catch( const exception& e ) {
            LOG_ERR << "getJob: Exception caught while returning work: " << e.what() << endl;
        }
        catch( ... ) {
            LOG_ERR << "getJob: Unrecognized exception caught while returning work." << endl;
        }

    }

}


void Worker::run( void ) {

    static int sleepS(1);

    while( getWork() ) {
        sleepS = 1;
        while( wip.job && wip.job->run( wip, ioService, threadPool ) ) ;
    }
    
    runTimer.expires_at(time_traits_t::now() + boost::posix_time::seconds(sleepS));
    runTimer.async_wait( strand.wrap(boost::bind( &Worker::run, this )) );
    
    if( sleepS < 4 ) {
        sleepS <<= 1;
    }

}
