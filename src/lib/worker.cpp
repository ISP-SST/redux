#include "redux/worker.hpp"

#ifdef DEBUG_
#   define TRACE_THREADS
#endif

#include "redux/daemon.hpp"
#include "redux/network/protocol.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/util/trace.hpp"

using namespace redux::logging;
using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;


Worker::Worker( Daemon& d ) : running_(false), stopped_(true), exitWhenDone_(false), resetWhenDone_(false),
    wip(nullptr), daemon( d ), myInfo(Host::myInfo()) {

}


Worker::~Worker( void ) {

}


void Worker::start( void ) {

    if( !stopped_ ) return;
    stopped_ = false;

    if( !wip ) {
        wip.reset( new WorkInProgress() );
    } else {
        wip->reset();
    }
    
    daemon.ioService.post( std::bind(&Worker::run, this) );

    myInfo.touch();
    myInfo.active();
    
}


void Worker::stop( void ) {

    running_ = false;

}


void Worker::exitWhenDone( void ) {
    
    if( running_ ) {
        exitWhenDone_ = true;
        stop();
    }
    
}


void Worker::resetWhenDone( void ) {
    
    if( running_ ) {
        resetWhenDone_ = true;
        stop();
    }
    
}


void Worker::done( void ) {

    running_ = false;
    wip->reset();
    myInfo.touch();
    myInfo.idle();

    if( exitWhenDone_ ) {
        daemon.stop();
    }

    if( resetWhenDone_ ) {
        daemon.reset();
    }

    stopped_ = true;
}


bool Worker::fetchWork( void ) {

    bool ret = false;
    network::TcpConnection::Ptr conn;
    string msg;

    try {

        if( (conn = daemon.getMaster()) ) {

            *conn << CMD_GET_WORK;
            *conn << wip->jobID;

            size_t blockSize;
            shared_ptr<char> buf = conn->receiveBlock( blockSize );               // reply

            if( blockSize ) {

                const char* ptr = buf.get();
                uint64_t count = wip->unpackWork( ptr, currentJob, conn->getSwapEndian() );

                if( count != blockSize ) {
                    throw invalid_argument( "Failed to unpack data, blockSize=" + to_string( blockSize ) + "  unpacked=" + to_string( count ) );
                }
                wip->isRemote = true;
                
                LLOG_TRACE(daemon.logger) << "Received work: " << wip->print() << ende;
                
                ret = true;

            } else {
                // This just means there is no work available at the moment.
            }

        }
        
    }
    catch( const exception& e ) {
        msg = "fetchWork: Exception caught while fetching job: ";
        msg += e.what();
    }
    catch( ... ) {
        msg = "fetchWork: Unrecognized exception caught while fetching job.";
    }
    
    if( !msg.empty() ) {
        try {   // only log if the connection was not severed (otherwise the manager will get spammed by messages on restart)
            auto test RDX_UNUSED = conn->socket().remote_endpoint();  // check if endpoint exists, will throw if not connected.
            LLOG_ERR(daemon.logger) << "fetchWork: Unrecognized exception caught while fetching job." << ende;
        } catch ( ... ) {}
        ret = false;
        if( conn ) conn->socket().close();
    }
    
    if(conn) daemon.unlockMaster();
    
    return ret;
}



bool Worker::getWork( void ) {

    myInfo.active();
    Job::JobPtr thisJob = wip->job.lock();
    if( wip->isRemote ) {            // remote work: return parts.
        returnWork();
        int count(0);
        while( wip->hasResults && count++ < 5 ) {
            LLOG_DEBUG(daemon.logger) << "Failed to return data, trying again in 5 seconds." << ende;
            std::this_thread::sleep_for( std::chrono::seconds(5) );
            returnWork();
        }
    } else if( thisJob ) {
        daemon.returnWork( wip );
    }

    if( wip->hasResults ) {
        LLOG_WARN(daemon.logger) << "Failed to return data, this part will be discarded." << ende;
    }
    
    myInfo.limbo();
    wip->isRemote = wip->hasResults = false;
    wip->resetParts();
    if( thisJob ) {
        thisJob->logger.flushAll();
        wip->jobID = thisJob->info.id;
    }
    
    try {
        
        boost::this_thread::interruption_point();
        if( running_ && !exitWhenDone_ && !resetWhenDone_ ) {
            if( daemon.getWork( wip, false ) || fetchWork() ) {    // first check for local work, then remote
                myInfo.active();
                myInfo.status.statusString = "...";
                thisJob = wip->job.lock();
                if(thisJob && (wip->jobID != thisJob->info.id)) {                                  // initialize if it is a new job.
                    thisJob->logger.setLevel( thisJob->info.verbosity );
                    if( wip->isRemote ) {
                        TcpConnection::Ptr logConn;
                        daemon.connect( daemon.myMaster.host->info, logConn );
                        thisJob->logger.addNetwork( daemon.ioService, daemon.myMaster.host, thisJob->info.id, 0, 5 );   // TODO make flushPeriod a config setting.
                    }
                    thisJob->init();
                    wip->jobID = thisJob->info.id;
                    currentJob = thisJob;
                }
                THREAD_MARK
                if( !wip->isRemote ) {
                    for( auto& part: wip->parts ) {
                        part->cacheLoad(false);         // load data for local jobs, but don't delete the storage
                    }
                }
                THREAD_UNMARK
                return true;
            }
        }
    } catch (const std::exception& e) {
        LLOG_TRACE(daemon.logger) << "Failed to get new work: reson: " << e.what() << ende;
    }
    
    boost::this_thread::interruption_point();
    currentJob.reset();
    THREAD_UNMARK;
    return false;

}



void Worker::returnWork( void ) {

    if( wip->hasResults ) {

        network::TcpConnection::Ptr conn = daemon.getMaster();
    
        try {


            if( conn && conn->socket().is_open() ) {

                LLOG_TRACE(daemon.logger) << "Returning result: " + wip->print() << ende;

                uint64_t blockSize = wip->workSize();       // might be slightly bigger than needed due to compression.
                size_t totalSize = blockSize + sizeof( uint64_t ) + 1;        // + blocksize + cmd
                
                shared_ptr<char> data = rdx_get_shared<char>(totalSize);
                char* ptr = data.get() + sizeof( uint64_t ) + 1;
                uint64_t count(0);
                if( blockSize ) {
                    count = wip->packWork( ptr );
                }
                count += pack( data.get()+1, count );
                count += pack( data.get(), CMD_PUT_PARTS );

                conn->asyncWrite( data, count );

                Command cmd = CMD_ERR;
                *(conn) >> cmd;
                
                if( cmd == CMD_OK ) {
                    wip->hasResults = false;
                } 
                
            }
            
        }
        catch( const exception& e ) {
            LLOG_ERR(daemon.logger) << "getJob: Exception caught while returning work: " << e.what() << ende;
        }
        catch( ... ) {
            LLOG_ERR(daemon.logger) << "getJob: Unrecognized exception caught while returning work." << ende;
        }

        if( conn ) daemon.unlockMaster();
            
    }

}


void Worker::run( void ) {

    running_ = true;

    if( wip ) {
        // LOG_TRACE << "run:   nWipParts = " << wip->parts.size() << "  conn = " << hexString(wip->connection.get()) << "  job = " << hexString(wip->job.get());
        while( getWork() ) {
            try {
                Job::JobPtr thisJob = wip->job.lock();
                while( thisJob && thisJob->run( wip, myInfo.status.nThreads ) ) ;
                if( thisJob ) thisJob->logger.flushAll(); 
            }
            catch( const boost::thread_interrupted& ) {
                LLOG_DEBUG(daemon.logger) << "Worker: Job interrupted."  << ende;
                throw;
            }
            catch( const exception& e ) {
                LLOG_ERR(daemon.logger) << "Worker: Exception caught while processing job: " << e.what() << ende;
            }
            catch( ... ) {
                LLOG_ERR(daemon.logger) << "Worker: Unrecognized exception caught while processing job." << ende;
            }
        }
    }

    done();
    
}
