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
    std::map<boost::thread::id,boost::thread*> tmap;
}

Worker::Worker( Daemon& d ) : runTimer( d.ioService ), running_(false), exitWhenDone_(false),
    wip(nullptr), daemon( d ), myInfo(Host::myInfo()) {

}


Worker::~Worker( void ) {

}


void Worker::threadLoop( void ) {

    while( true ) {
        try {
            boost::this_thread::interruption_point();
            ioService.run();
        } catch( const Application::ThreadExit& e ) {
            break;
        } catch( job_error& e ) {
            LLOG_ERR(daemon.logger) << "Job error: " << e.what() << ende;
        } catch( const boost::thread_interrupted& ) {
            LLOG_TRACE(daemon.logger) << "Worker: Thread interrupted." << ende;
            break;
        } catch( exception& e ) {
            LLOG_ERR(daemon.logger) << "Exception in thread: " << e.what() << ende;
        } catch( ... ) {
            LLOG_ERR(daemon.logger) << "Unhandled exception in thread." << ende;
        }
    }
return; // FIXME it gets stuck on reset/kill
    pool.remove_thread( tmap[boost::this_thread::get_id()] );

}


void Worker::start( void ) {

    running_ = true;
    wip.reset( new WorkInProgress() );
    ioService.reset();
    workLoop.reset( new boost::asio::io_service::work(ioService) );

    uint16_t nThreads = myInfo.status.nThreads;
    myInfo.status.nThreads = 0;
    //boost::thread* t = pool.create_thread( std::bind( &Worker::threadLoop, this ) ); // add a timer thread which will not "count" as a thread
    //tmap[ t->get_id() ] = t;
    setThreads(nThreads);
    
    runTimer.expires_from_now( boost::posix_time::seconds( 5 ) );
    runTimer.async_wait( std::bind( &Worker::run, this, placeholders::_1 ) );
    myInfo.touch();

}


void Worker::stop( void ) {

    running_ = false;
    pool.interrupt_all();
    workLoop.reset();
    ioService.stop();
    runTimer.cancel();
    pool.join_all();
    wip.reset();

}


void Worker::setThreads( int nThreads ) {

    int diff = nThreads - myInfo.status.nThreads;
    myInfo.status.nThreads += diff;
    if( diff > 0 ) {
        for( int i=0; i < diff; ++i ) {
            boost::thread* t = pool.create_thread( std::bind( &Worker::threadLoop, this ) );
            tmap[ t->get_id() ] = t;
        }
    } else while( diff++ < 0 ) {
        ioService.post( [](){ throw Application::ThreadExit(); } );     // Make diff threads throw exception to exit.
    }
    
}


bool Worker::fetchWork( void ) {

    bool ret = false;
    network::TcpConnection::Ptr conn;
    string msg;
    
    try {

        if( (conn = daemon.getMaster()) ) {

            *conn << CMD_GET_WORK;

            size_t blockSize;
            shared_ptr<char> buf = conn->receiveBlock( blockSize );               // reply

            if( blockSize ) {

                const char* ptr = buf.get();
                uint64_t count = wip->unpackWork( ptr, conn->getSwapEndian() );

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

    if( wip->isRemote ) {            // remote work: return parts.
        returnWork();
        int count(0);
        while( wip->hasResults && count++ < 5 ) {
            LLOG_DEBUG(daemon.logger) << "Failed to return data, trying again in 5 seconds." << ende;
            runTimer.expires_from_now(boost::posix_time::seconds(5));
            runTimer.wait();
            returnWork();
        }
        myInfo.active();
    } else if ( wip->hasResults ) {
        wip->returnResults();
        myInfo.active();
    }

    if( wip->hasResults ) {
        LLOG_WARN(daemon.logger) << "Failed to return data, this part will be discarded." << ende;
    }
    
    try {
        
        wip->isRemote = wip->hasResults = false;
        if( wip->job ) {
            wip->job->logger.flushAll();
        }
    
        boost::this_thread::interruption_point();
        if( running_ ) {
            if( daemon.getWork( wip, myInfo.status.nThreads ) || fetchWork() ) {    // first check for local work, then remote
                if( wip->previousJob.expired() ) {                                  // initialize if it is a new job.
                    wip->job->logger.setLevel( wip->job->info.verbosity );
                    if( wip->isRemote ) {
                        if( daemon.params.count( "log-stdout" ) ) { // -d flag was passed on cmd-line
                            wip->job->logger.addLogger( daemon.logger );
                        }
                        TcpConnection::Ptr logConn;
                        daemon.connect( daemon.myMaster.host->info, logConn );
                        wip->job->logger.addNetwork( daemon.ioService, daemon.myMaster.host, wip->job->info.id, 0, 5 );   // TODO make flushPeriod a config setting.
                    }
                    wip->job->init();
                    wip->previousJob = wip->job;
                }
                for( auto& part: wip->parts ) {
                    part->cacheLoad(false);               // load data for local jobs, but don't delete the storage
                }
                myInfo.active();
                myInfo.status.statusString = "...";
                return true;
            }
        }
    } catch (const std::exception& e) {
        LLOG_TRACE(daemon.logger) << "Failed to get new work: reson: " << e.what() << ende;
    }
    
#ifdef DEBUG_
    LLOG_TRACE(daemon.logger) << "No work available." << ende;
#endif

    wip->job.reset();
    wip->parts.clear();
    
    myInfo.idle();
    boost::this_thread::interruption_point();
    return false;

}



void Worker::returnWork( void ) {

    if( wip->hasResults ) {

        network::TcpConnection::Ptr conn = daemon.getMaster();
    
        try {


            if( conn && conn->socket().is_open() ) {

                LLOG_DEBUG(daemon.logger) << "Returning result: " + wip->print() << ende;

                uint64_t blockSize = wip->workSize();
                size_t totalSize = blockSize + sizeof( uint64_t ) + 1;        // + blocksize + cmd
                
                shared_ptr<char> data( new char[totalSize], []( char* p ){ delete[] p; } );
                char* ptr = data.get();
                
                uint64_t count = pack( ptr, CMD_PUT_PARTS );
                count += pack( ptr+count, blockSize );
                if(blockSize) {
                    count += wip->packWork(ptr+count);

                }

                conn->asyncWrite( data, totalSize );

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


void Worker::run( const boost::system::error_code& error ) {

    static int sleepS(1);
    if( error == boost::asio::error::operation_aborted ) {
        sleepS=1;
        return;
    }
 //   LOG_TRACE << "run:   nWipParts = " << wip->parts.size() << "  conn = " << hexString(wip->connection.get()) << "  job = " << hexString(wip->job.get());
    while( getWork() ) {
        sleepS = 1;
        try {
            while( wip->job && wip->job->run( wip, ioService, myInfo.status.nThreads ) ) ;
            if( wip->job ) wip->job->logger.flushAll(); 
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

    if( running_ ) {
        boost::this_thread::interruption_point();
        runTimer.expires_from_now( boost::posix_time::seconds(sleepS) );
        runTimer.async_wait( std::bind( &Worker::run, this, placeholders::_1 ) );
        if( sleepS < 4 ) {
            sleepS <<= 1;
        }
    } else if( exitWhenDone_ ) {
        workLoop.reset();
        ioService.stop();
        daemon.stop();
    }

}
