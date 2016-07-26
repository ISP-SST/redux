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

        boost::asio::strand strand;
        boost::asio::deadline_timer runTimer;
        
        WorkInProgress wip;

        Daemon& daemon;

        
    };


}   // redux

#endif // REDUX_WORKER_HPP
