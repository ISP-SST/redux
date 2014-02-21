#ifndef REDUX_REDUXAPP_HPP
#define REDUX_REDUXAPP_HPP

#include "redux/application.hpp"
#include "redux/job.hpp"
#include "redux/worker.hpp"
#include "redux/network/peer.hpp"
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
        network::Peer::Ptr& addOrGetPeer(const network::Peer::HostInfo&, network::TcpConnection::Ptr&);
        network::Peer::Ptr& getPeer(const network::TcpConnection::Ptr&);
        void cleanupPeers(void);
        void addJobs( network::Peer::Ptr& );
        void removeJobs( network::Peer::Ptr& );
        Job::JobPtr selectJob(bool);
        bool getWork( WorkInProgress& );
        void sendWork( network::Peer::Ptr& );
        void putParts( network::Peer::Ptr& );
        void sendJobList( network::TcpConnection::Ptr& );
        void updatePeerStatus( network::Peer::Ptr& );
        void sendJobStats( network::TcpConnection::Ptr& );
        void sendPeerList( network::TcpConnection::Ptr& );
        void updateLoadAvg( void );
        
        std::string master;
        uint16_t port;
        
        const po::variables_map& params;
        
        std::mutex jobMutex;
        size_t jobCounter;
        std::vector<Job::JobPtr> jobs;
        uint8_t nQueuedJobs;
        
        std::mutex peerMutex;
        network::Peer::Ptr myInfo;
        std::map<size_t, network::Peer::Ptr> peers;
        std::map<network::Peer::Ptr, WorkInProgress> peerWIP;
        
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
