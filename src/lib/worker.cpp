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
#define logChannel "worker"

namespace {
    typedef boost::asio::time_traits<boost::posix_time::ptime> time_traits_t;
}


Worker::Worker( Daemon& d ) : strand(d.ioService), runTimer( d.ioService ), daemon( d ) {

}


Worker::~Worker( void ) {

}


void Worker::init( void ) {

    runTimer.expires_at( time_traits_t::now() + boost::posix_time::seconds( 1 ) );
    runTimer.async_wait(strand.wrap(boost::bind(&Worker::run, this)));

}


void Worker::stop( void ) {

    ioService.stop();

}


bool Worker::fetchWork( void ) {

    bool ret = false;
    network::TcpConnection::Ptr conn;
    
    try {

        if( (conn = daemon.getMaster()) ) {

            *conn << CMD_GET_WORK;

            size_t blockSize;
            shared_ptr<char> buf = conn->receiveBlock( blockSize );               // reply

            if( blockSize ) {

                const char* ptr = buf.get();
                uint64_t count = wip.unpackWork( ptr, conn->getSwapEndian() );

                if( count != blockSize ) {
                    throw invalid_argument( "Failed to unpack data, blockSize=" + to_string( blockSize ) + "  unpacked=" + to_string( count ) );
                }
                wip.isRemote = true;

                LOG_TRACE << "Received work: " << wip.print();
                
                ret = true;
                
            }

        }
        
    }
    catch( const exception& e ) {
        LOG_ERR << "fetchWork: Exception caught while fetching job: " << e.what();
        if( conn ) conn->socket().close();
        ret = false;
    }
    catch( ... ) {
        LOG_ERR << "fetchWork: Unrecognized exception caught while fetching job.";
        if( conn ) conn->socket().close();
        ret = false;
    }
    
    if(conn) daemon.unlockMaster();
    
    return ret;
}



bool Worker::getWork( void ) {

    if( wip.isRemote ) {            // remote work: return parts.
        returnWork();
        int count(0);
        while( wip.hasResults && count++ < 5 ) {
            LOG_DEBUG << "Failed to return data, trying again in 5 seconds.";
            runTimer.expires_at(time_traits_t::now() + boost::posix_time::seconds(5));
            runTimer.wait();
            returnWork();
        }
        daemon.myInfo->active();
    } else if ( wip.hasResults ) {
        wip.returnResults();
        daemon.myInfo->active();
    }

    if( wip.hasResults ) {
        LOG_WARN << "Failed to return data, this part will be discarded.";
    }
    
    wip.isRemote = wip.hasResults = false;
    
    if( daemon.getWork( wip, daemon.myInfo->status.nThreads ) || fetchWork() ) {    // first check for local work, then remote
        if( wip.job && (!wip.previousJob || *(wip.job) != *(wip.previousJob)) ) {
            LOG_DEBUG << "Initializing new job: " + wip.print();
            wip.job->init();
            wip.previousJob = wip.job;
        } else LOG_DEBUG << "Starting new part of the same job: " + wip.print();
        for( auto& part: wip.parts ) {
            part->cacheLoad(false);               // load data for local jobs, but don't delete the storage
        }
        daemon.myInfo->active();
        return true;
    }
    
#ifdef DEBUG_
    LOG_TRACE << "No work available.";
#endif

    wip.job.reset();
    wip.previousJob.reset();
    wip.parts.clear();
    
    daemon.myInfo->status.state = Host::ST_IDLE;
    return false;

}



void Worker::returnWork( void ) {

    if( wip.hasResults ) {
        
        try {
            
            network::TcpConnection::Ptr conn = daemon.getMaster();

            if( conn && conn->socket().is_open() ) {

                LOG_DEBUG << "Returning result: " + wip.print();
                uint64_t blockSize = wip.workSize();
                size_t totalSize = blockSize + sizeof( uint64_t ) + 1;        // + blocksize + cmd
                
                shared_ptr<char> data = sharedArray<char>( totalSize );
                char* ptr = data.get();
                
                uint64_t count = pack( ptr, CMD_PUT_PARTS );
                count += pack( ptr+count, blockSize );
                if(blockSize) {
                    count += wip.packWork(ptr+count);

                }

                conn->asyncWrite( data, totalSize );

                Command cmd = CMD_ERR;
                *(conn) >> cmd;
                
                if( cmd == CMD_OK ) {
                    wip.hasResults = false;
                } 
                
            }
            
            if( conn ) daemon.unlockMaster();
            
        }
        catch( const exception& e ) {
            LOG_ERR << "getJob: Exception caught while returning work: " << e.what();
        }
        catch( ... ) {
            LOG_ERR << "getJob: Unrecognized exception caught while returning work.";
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
