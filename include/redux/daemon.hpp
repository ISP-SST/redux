#ifndef REDUX_REDUXAPP_HPP
#define REDUX_REDUXAPP_HPP

#include "redux/application.hpp"
#include "redux/job.hpp"
#include "redux/work.hpp"
#include "redux/worker.hpp"
#include "redux/network/host.hpp"
#include "redux/network/tcpserver.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/semaphore.hpp"

#include <mutex>
#include <thread>
#include <vector>

#include <boost/asio.hpp>

namespace redux {

    class Daemon : public redux::Application {

    public:

        explicit Daemon ( po::variables_map& vm );
        virtual ~Daemon ( void );
        
        void reset(void);
        void stop(void);
        
    private:
        
        void start_server( uint16_t&, uint8_t tries=1 );
        void stop_server( void );
        void maintenance( void );
        void checkSwapSpace( void );
        void checkCurrentUsage( void );
        void check_limits( void );
        bool doWork(void);
        
        bool workerInit( void );
        void connect( network::Host::HostInfo& host, network::TcpConnection::Ptr& conn );
        void updateStatus(void);
        network::TcpConnection::Ptr getMaster(void);
        void unlockMaster(void);
        
        void connected( network::TcpConnection::Ptr );
        void handler( network::TcpConnection::Ptr );
        void urgentHandler( network::TcpConnection::Ptr );
        void processCommand( network::TcpConnection::Ptr, uint8_t, bool urgent=false);
        
        void addConnection( network::TcpConnection::Ptr );
        void removeConnection( network::TcpConnection::Ptr );
        void cleanup(void);
        void failedWIP( WorkInProgress::Ptr wip );
        void die(void);
        void die( network::TcpConnection::Ptr&, bool urgent=false );
        void softExit(void);
        void reset( network::TcpConnection::Ptr&, bool urgent=false );
        void addJobs( network::TcpConnection::Ptr& );
        void failJobs( const std::vector<size_t>& );
        void failJobs( std::string );
        void removeJobs( const std::vector<size_t>& );
        void removeJobs( std::string );
        void removeJobs( network::TcpConnection::Ptr& );
        void interactiveCB( network::TcpConnection::Ptr );
        void interactive( network::TcpConnection::Ptr& );
        void listen( void );
        void pokeSlaves( size_t n );
        void sendToSlaves( uint8_t cmd, std::string slvString );
        void resetSlaves( network::TcpConnection::Ptr&, uint8_t );
        //Job::JobPtr selectJob(bool);
        network::Host::Ptr getHost( network::TcpConnection::Ptr& );
        WorkInProgress::Ptr getWIP( const network::Host::Ptr& );
        void removeWIP( const network::Host::Ptr& );
        void updateWIP( const network::Host::Ptr&, WorkInProgress::Ptr& );
        void sendWork( network::TcpConnection::Ptr );
        void putParts( network::TcpConnection::Ptr );
        void sendJobList( network::TcpConnection::Ptr& );
        void updateHostStatus( network::TcpConnection::Ptr& );
        void sendJobStats( network::TcpConnection::Ptr& );
        void sendPeerList( network::TcpConnection::Ptr& );
        void addToLog( network::TcpConnection::Ptr& );
        void updateLoadAvg( void );
        
        void putActiveWIP( WorkInProgress::Ptr );
        size_t countActiveWIP( WorkInProgress::Ptr );
        size_t eraseActiveWIP( WorkInProgress::Ptr );
        WorkInProgress::Ptr  getActiveWIP( WorkInProgress::Ptr );
        void getIdleWIP( WorkInProgress::Ptr& );
        WorkInProgress::Ptr getIdleWIP( void );
        void putIdleWIP( WorkInProgress::Ptr );
        WorkInProgress::Ptr getLocalWIP( void );
        WorkInProgress::Ptr getRemoteWIP( void );
        bool getWork( WorkInProgress::Ptr&, bool remote=false );
        void returnWork( WorkInProgress::Ptr );
        void prepareLocalWork( int n=1 );
        void prepareRemoteWork( int n=1 );
        void prepareWork( void );
        std::string workStatus( bool details=false );
        
        const po::variables_map& params;
        
        std::mutex jobsMutex;
        std::vector<Job::JobPtr> jobs;
        size_t jobCounter;
        uint16_t nQueuedJobs;
        uint32_t hostTimeout;
        util::Semaphore inTransfers, outTransfers;

        network::Host& myInfo;
        
        std::mutex peerMutex;
        std::set<network::Host::Ptr, network::Host::Compare> peers;     // list of slaves (i.e. will not include connections from e.g. rdx_stat)
        
        std::mutex threadMutex;
        void addThread( uint16_t n );
        void delThread( uint16_t n );
        void cleanupThreads( void );
        size_t nSysThreads( void );
        void threadLoop( void );
        
        std::mutex wip_active_mtx, wip_idle_mtx, wip_queue_mtx, wip_lqueue_mtx, wip_completed_mtx;
        std::atomic<uint16_t>  max_local, max_remote;
        std::set<WorkInProgress::Ptr, redux::util::ObjCompare<WorkInProgress>> wip_active;
        std::map<network::Host::Ptr, WorkInProgress::Ptr, network::Host::Compare> peerWIP;
        std::set<WorkInProgress::Ptr> wip_idle;
        std::deque<WorkInProgress::Ptr> wip_queue;
        std::deque<WorkInProgress::Ptr> wip_localqueue;
        std::deque<WorkInProgress::Ptr> wip_completed;
        
        struct {
            network::TcpConnection::Ptr conn;
            network::Host::Ptr host;
        } myMaster;
        
        boost::asio::io_service ioService;
        boost::thread_group pool;
        boost::asio::deadline_timer timer;
        std::unique_ptr<network::TcpServer> server;
        
        Worker worker;
        
        friend class network::TcpServer;
        friend class Worker;
        
    };


}

#endif // REDUX_REDUXAPP_HPP
