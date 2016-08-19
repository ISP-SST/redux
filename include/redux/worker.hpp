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

        void init( void );
        void connect( void );
        void stop( void );


    private:

        bool fetchWork(void);
        
        bool getWork(void);
        void returnWork(void);
        void returnJob(void);
        void returnResults(void);
        
        void run(void);
        
        boost::asio::io_service ioService;

        std::shared_ptr<boost::asio::io_service::work> workLoop;
        boost::thread_group pool;
        boost::asio::strand strand;
        boost::asio::deadline_timer runTimer;
        
        WorkInProgress::Ptr wip;

        Daemon& daemon;
        network::Host& myInfo;

        
    };


}   // redux

#endif // REDUX_WORKER_HPP
