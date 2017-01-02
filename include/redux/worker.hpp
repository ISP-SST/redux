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

        explicit Worker( Daemon& );
        ~Worker( void );

        void connect( void );
        void start( void );
        void stop( void );
        void exitWhenDone( void ) { running_ = false; exitWhenDone_ = true; };

    private:

        bool fetchWork(void);
        
        bool getWork(void);
        void returnWork(void);
        void returnJob(void);
        void returnResults(void);
        
        void run( const boost::system::error_code& );
        
        boost::asio::io_service ioService;

        std::shared_ptr<boost::asio::io_service::work> workLoop;
        boost::thread_group pool;
        boost::asio::strand strand;
        boost::asio::deadline_timer runTimer;
        bool running_;
        bool exitWhenDone_;
        
        WorkInProgress::Ptr wip;

        Daemon& daemon;
        network::Host& myInfo;

        
    };


}   // redux

#endif // REDUX_WORKER_HPP
