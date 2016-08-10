#include "redux/worker.hpp"

#include "redux/daemon.hpp"
#include "redux/network/protocol.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

using namespace redux::logging;
using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;


namespace {
    typedef boost::asio::time_traits<boost::posix_time::ptime> time_traits_t;
}


Worker::Worker( Daemon& d ) : strand(d.ioService), runTimer( d.ioService ), daemon( d ), myInfo(Host::myInfo()) {

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
                
                LLOG_TRACE(daemon.logger) << "Received work: " << wip.print() << ende;
                
                ret = true;
                
            }

        }
        
    }
    catch( const exception& e ) {
        LLOG_ERR(daemon.logger) << "fetchWork: Exception caught while fetching job: " << e.what() << ende;
        if( conn ) conn->socket().close();
        ret = false;
    }
    catch( ... ) {
        LLOG_ERR(daemon.logger) << "fetchWork: Unrecognized exception caught while fetching job." << ende;
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
            LLOG_DEBUG(daemon.logger) << "Failed to return data, trying again in 5 seconds." << ende;
            runTimer.expires_at(time_traits_t::now() + boost::posix_time::seconds(5));
            runTimer.wait();
            returnWork();
        }
        myInfo.active();
    } else if ( wip.hasResults ) {
        wip.returnResults();
        myInfo.active();
    }

    if( wip.hasResults ) {
        LLOG_WARN(daemon.logger) << "Failed to return data, this part will be discarded." << ende;
    }
    
    wip.isRemote = wip.hasResults = false;
    
    if( daemon.getWork( wip, myInfo.status.nThreads ) || fetchWork() ) {    // first check for local work, then remote
        if( wip.job && (!wip.previousJob || *(wip.job) != *(wip.previousJob)) ) {
            if( wip.previousJob ) {
                wip.previousJob->logger.flushAll();
            }
            wip.job->logger.setLevel( wip.job->info.verbosity );
            if( wip.isRemote ) {
                if( daemon.params.count( "log-stdout" ) ) { // -d flag was passed on cmd-line
                    wip.job->logger.addLogger( daemon.logger );
                }
                TcpConnection::Ptr logConn;
                daemon.connect( daemon.myMaster.host->info, logConn );
                wip.job->logger.addNetwork( logConn, wip.job->info.id, 0, 5 );   // TODO make flushPeriod a config setting.
            }
            wip.job->init();
            wip.previousJob = wip.job;
        }
        for( auto& part: wip.parts ) {
            part->cacheLoad(false);               // load data for local jobs, but don't delete the storage
        }
        myInfo.active();
        myInfo.status.statusString = "...";
        return true;
    }
    
#ifdef DEBUG_
    LLOG_TRACE(daemon.logger) << "No work available." << ende;
#endif

    wip.job.reset();
    wip.previousJob.reset();
    wip.parts.clear();
    
    myInfo.status.state = Host::ST_IDLE;
    myInfo.status.statusString = "idle";
    return false;

}



void Worker::returnWork( void ) {

    if( wip.hasResults ) {
        
        try {
            
            network::TcpConnection::Ptr conn = daemon.getMaster();

            if( conn && conn->socket().is_open() ) {

                LLOG_DEBUG(daemon.logger) << "Returning result: " + wip.print() << ende;
                uint64_t blockSize = wip.workSize();
                size_t totalSize = blockSize + sizeof( uint64_t ) + 1;        // + blocksize + cmd
                
                shared_ptr<char> data( new char[totalSize], []( char* p ){ delete[] p; } );
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
            LLOG_ERR(daemon.logger) << "getJob: Exception caught while returning work: " << e.what() << ende;
        }
        catch( ... ) {
            LLOG_ERR(daemon.logger) << "getJob: Unrecognized exception caught while returning work." << ende;
        }

    }

}


void Worker::run( void ) {

    static int sleepS(1);

 //   LOG_TRACE << "run:   nWipParts = " << wip.parts.size() << "  conn = " << hexString(wip.connection.get()) << "  job = " << hexString(wip.job.get());
    while( getWork() ) {
        //cout << "Worker::run()  got work, calling run..." << endl;
        sleepS = 1;
        try {
            while( wip.job && wip.job->run( wip, ioService, myInfo.status.nThreads ) ) ;
            if( wip.job ) wip.job->logger.flushAll(); 
        }
        catch( const exception& e ) {
            LLOG_ERR(daemon.logger) << "Worker: Exception caught while processing job: " << e.what() << ende;
        }
        catch( ... ) {
            LLOG_ERR(daemon.logger) << "Worker: Unrecognized exception caught while processing job." << ende;
        }
    }

    runTimer.expires_at(time_traits_t::now() + boost::posix_time::seconds(sleepS));
    runTimer.async_wait( strand.wrap(boost::bind( &Worker::run, this )) );
    
    if( sleepS < 4 ) {
        sleepS <<= 1;
    }

}
