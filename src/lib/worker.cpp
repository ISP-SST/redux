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

    peer.reset( new Host() );
    connection = TcpConnection::newPtr( daemon.ioService );
    if( daemon.port < 1024 ) {
        LOG_DEBUG << "Running local worker only.";
    } else {
        connect();
    }
    runTimer.expires_at( time_traits_t::now() + boost::posix_time::seconds( 1 ) );
    runTimer.async_wait(strand.wrap(boost::bind(&Worker::run, this)));

}

void Worker::connect( void ) {
    
    if( daemon.port < 1024 ) {
        return;
    }
    
    try {
        auto test UNUSED = connection->socket().remote_endpoint();  // check if endpoint exists, if not we need to re-establish connection.
        return;
    } catch ( ... ) { }

    if( connection->socket().is_open() ) {
        connection->socket().close();
    }

    LOG_DEBUG << "Attempting to connect to master at " << daemon.master << ":" << daemon.port;

    connection->connect( daemon.params["master"].as<string>(), to_string( daemon.params["port"].as<uint16_t>() ) );
    if( connection->socket().is_open() ) {

        Command cmd;

        daemon.myInfo->info.peerType |= Host::TP_WORKER;
        *connection << CMD_CONNECT;
        *connection >> cmd;
        if( cmd == CMD_AUTH ) {
            // implement
        }
        if( cmd == CMD_CFG ) {  // handshake requested
            *connection << daemon.myInfo->info;
            *connection >> peer->info;
            *connection >> cmd;       // ok or err
        }
        if( cmd != CMD_OK ) {
            LOG_ERR << "Handshake with master failed  (server replied: " << cmd << ")";
            connection->socket().close();
            daemon.myInfo->info.peerType &= ~Host::TP_WORKER;
        } else LOG_DEBUG << "Connected.";

    }

}

void Worker::stop( void ) {
    if( connection->socket().is_open() ) {
        *connection << CMD_DISCONNECT;
        connection->socket().close();
    }
    ioService.stop();
    daemon.myInfo->info.peerType &= ~Host::TP_WORKER;

}


void Worker::updateStatus( void ) {

    if( !connection->socket().is_open() ) {
        connect();
    }

    if( connection->socket().is_open() ) {
#ifdef DEBUG_
        LOG_TRACE << "Sending statusupdate to server";
#endif
        size_t blockSize = daemon.myInfo->status.size();
        size_t totSize = blockSize + sizeof( size_t ) + 1;
        auto buf = sharedArray<char>(totSize);
        char* ptr = buf.get();
        memset( ptr, 0, totSize );

        uint64_t count = pack( ptr, CMD_STAT );
        count += pack( ptr+count, blockSize );
        count += daemon.myInfo->status.pack( ptr+count );

        connection->writeAndCheck( buf, totSize );
    }
    else {
        if( daemon.port > 1024 ) {
            LOG_WARN << "No connection to master:  " << daemon.master << ":" << daemon.port;
        }
    }

}

bool Worker::fetchWork( void ) {

    try {

        if( !connection->socket().is_open() ) {
            connect();
        }

        if( !connection->socket().is_open() ) {
            return false;
        }

        LOG_DEBUG << "Requesting remote work.";
        *connection << CMD_GET_WORK;

        size_t blockSize;
        shared_ptr<char> buf = connection->receiveBlock( blockSize );               // reply
        LOG_TRACE << "Received " << blockSize << " bytes.";

        if( !blockSize ) return false;

        const char* ptr = buf.get();
        uint64_t count = wip.unpackWork( ptr, connection->getSwapEndian() );

        if( count != blockSize ) {
            throw invalid_argument( "Failed to unpack data, blockSize=" + to_string( blockSize ) + "  unpacked=" + to_string( count ) );
        }
        wip.connection = connection;

        return true;
    }
    catch( const exception& e ) {
        LOG_ERR << "fetchWork: Exception caught while fetching job: " << e.what() << endl;
    }
    catch( ... ) {
        LOG_ERR << "fetchWork: Unrecognized exception caught while fetching job." << endl;
    }

    return false;
}



bool Worker::getWork( void ) {

    if( wip.connection ) {            // remote work: return parts.
        returnWork();
    } else if ( wip.parts.size() ) {
        daemon.returnResults( wip );
    }

    if( daemon.getWork( wip, daemon.myInfo->status.nThreads ) || fetchWork() ) {    // first check for local work, then remote
        if( wip.job && (!wip.previousJob || *(wip.job) != *(wip.previousJob)) ) {
            LOG_DEBUG << "Initializing new job: " + wip.print();
            wip.job->init();
            if( wip.previousJob ) {
                wip.previousJob->cleanup();
            }
            wip.previousJob = wip.job;
        } else LOG_DEBUG << "Starting new part of the same job: " + wip.print();
        for( auto& part: wip.parts ) {
            part->cacheLoad(false);               // load data for local jobs, but don't delete the storage
        }
        return true;
    }
    
#ifdef DEBUG_
    LOG_TRACE << "No work available.";
#endif

    wip.job.reset();
    wip.previousJob.reset();
    wip.connection.reset();
    wip.parts.clear();
    
    return false;

}



void Worker::returnWork( void ) {

    if( wip.parts.size() ) {

        try {

            if( !wip.connection->socket().is_open() ) {
                return;     // TODO handle reconnects
            }

            LOG_DETAIL << "Returning result: " + wip.print();
            uint64_t blockSize = wip.workSize();
            size_t totalSize = blockSize + sizeof( uint64_t ) + 1;        // + blocksize + cmd
            
            shared_ptr<char> data = sharedArray<char>( totalSize );
            char* ptr = data.get();
            
            uint64_t count = pack( ptr, CMD_PUT_PARTS );
            count += pack( ptr+count, blockSize );
           if(blockSize) {
                count += wip.packWork(ptr+count);

            }

            connection->writeAndCheck( data, totalSize );

            Command cmd = CMD_ERR;

            *(connection) >> cmd;

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

 //   LOG_TRACE << "run:   nWipParts = " << wip.parts.size() << "  conn = " << hexString(wip.connection.get()) << "  job = " << hexString(wip.job.get());
    while( getWork() ) {
        //cout << "Worker::run()  got work, calling run..." << endl;
        sleepS = 1;
        while( wip.job && wip.job->run( wip, ioService, daemon.myInfo->status.nThreads ) ) ;
    }

    runTimer.expires_at(time_traits_t::now() + boost::posix_time::seconds(sleepS));
    runTimer.async_wait( strand.wrap(boost::bind( &Worker::run, this )) );
    
    if( sleepS < 4 ) {
        sleepS <<= 1;
    }

}
