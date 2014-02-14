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
        
        void connected(network::TcpConnection::ptr);
        void activity( network::TcpConnection::ptr );
        network::Peer::ptr& addOrGetPeer(const network::Peer::HostInfo&, network::TcpConnection::ptr&);
        network::Peer::ptr& getPeer(const network::TcpConnection::ptr&);
        void cleanupPeers(void);
        void addJobs( network::Peer::ptr& );
        void sendJobList( network::TcpConnection::ptr& );
        void sendPeerList( network::TcpConnection::ptr& );
        
        std::string master;
        uint16_t port;
        
        const po::variables_map& params;
        
        std::mutex jobMutex;
        size_t jobCounter;
        std::vector<Job::JobPtr> jobs;
        
        std::mutex peerMutex;
        network::Peer myInfo;
        std::map<size_t, network::Peer::ptr> peers;
        
        boost::asio::io_service ioService;
        std::unique_ptr<network::TcpServer> server;
        
        Worker worker;
        std::vector<std::shared_ptr<std::thread> > threads;
        
        friend class network::TcpServer;
        friend class Worker;
        
    };


}

#endif // REDUX_REDUXAPP_HPP
