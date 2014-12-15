#ifndef REDUX_WORKER_HPP
#define REDUX_WORKER_HPP

#include "redux/network/host.hpp"
#include "redux/job.hpp"
#include "redux/work.hpp"

#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/thread/thread.hpp>

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
        
        boost::asio::io_service ioService;

        boost::asio::strand strand;
        boost::asio::deadline_timer runTimer;
        
        WorkInProgress wip;

        Daemon& daemon;
        network::Host::Ptr peer;
        network::TcpConnection::Ptr connection;
        
    };


}   // redux

#endif // REDUX_WORKER_HPP
