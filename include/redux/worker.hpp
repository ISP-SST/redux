#ifndef REDUX_WORKER_HPP
#define REDUX_WORKER_HPP

#include "redux/network/peer.hpp"
#include "redux/job.hpp"
#include "redux/work.hpp"

#include <thread>

#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/program_options.hpp>


namespace po = boost::program_options;

using boost::asio::ip::tcp;

namespace redux {

    class Daemon;
    class Worker : private boost::noncopyable {

    public:

        Worker( Daemon& );
        ~Worker( void );

        void init( void );
        void connect( void );
        void stop( void );

        void updateStatus(void);
    private:

        bool fetchWork(void);
        bool getWork(void);
        
        bool getJob(void);
        bool getParts(void);
        void returnWork(void);
        void returnJob(void);
        void returnParts(void);
        
        void run(void);
        void maintenance(void);
        
        boost::asio::io_service ioService;
        boost::asio::strand strand;
        boost::asio::deadline_timer maintenanceTimer, runTimer;

        std::vector<std::shared_ptr<std::thread> > threads;
        
        WorkInProgress wip;

        Daemon& daemon;
        network::Peer::Ptr master;
    };


}   // redux

#endif // REDUX_WORKER_HPP
