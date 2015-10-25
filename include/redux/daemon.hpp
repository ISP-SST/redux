#ifndef REDUX_REDUXAPP_HPP
#define REDUX_REDUXAPP_HPP

#include "redux/application.hpp"
#include "redux/job.hpp"
#include "redux/worker.hpp"
#include "redux/network/host.hpp"
#include "redux/network/tcpserver.hpp"

#include <mutex>
#include <thread>
#include <vector>

#include <boost/asio.hpp>

namespace redux {

    struct WorkInProgress;

    class Daemon : public redux::Application {

    public:

        Daemon ( po::variables_map& vm );
        virtual ~Daemon ( void );
        
        void reset(void);
        void stop(void);
        
    private:
        
        void serverInit( void );
        void maintenance( void );
        bool doWork(void);
        
        void connected(network::TcpConnection::Ptr);
        void activity( network::TcpConnection::Ptr );
        void addConnection(const network::Host::HostInfo&, network::TcpConnection::Ptr&);
        void removeConnection(network::TcpConnection::Ptr&);
        void cleanup(void);
        void die( network::TcpConnection::Ptr& );
        void addJobs( network::TcpConnection::Ptr& );
        void removeJobs( network::TcpConnection::Ptr& );
        //Job::JobPtr selectJob(bool);
        bool getWork( WorkInProgress&, uint8_t nThreads = 1);
        void returnResults( WorkInProgress& );
        void sendWork( network::TcpConnection::Ptr& );
        void putParts( network::TcpConnection::Ptr& );
        void sendJobList( network::TcpConnection::Ptr& );
        void updateHostStatus( network::TcpConnection::Ptr& );
        void sendJobStats( network::TcpConnection::Ptr& );
        void sendPeerList( network::TcpConnection::Ptr& );
        void updateLoadAvg( void );
        
        std::string master;
        uint16_t port;
        
        const po::variables_map& params;
        
        std::mutex jobMutex;
        size_t jobCounter;
        std::vector<Job::JobPtr> jobs;
        uint16_t nQueuedJobs;
        uint32_t hostTimeout;
        
        std::mutex peerMutex;
        network::Host::Ptr myInfo;
        std::map<network::TcpConnection::Ptr, network::Host::Ptr> connections;
        std::map<size_t, network::Host::Ptr> peers;
        std::map<network::Host::Ptr, WorkInProgress, network::Host::Compare> peerWIP;
        
        boost::asio::io_service ioService;
        boost::asio::deadline_timer timer;
        std::unique_ptr<network::TcpServer> server;
        
        Worker worker;
        std::vector<std::shared_ptr<std::thread> > threads;
        
        friend class network::TcpServer;
        friend class Worker;
        
    };


}

#endif // REDUX_REDUXAPP_HPP
