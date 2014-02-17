#ifndef REDUX_WORKER_HPP
#define REDUX_WORKER_HPP

#include "redux/network/peer.hpp"

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
        void stop( void );

        void updateStatus(void);
    private:

//         boost::asio::io_service ioService;
//         const po::variables_map& params;
        network::TcpConnection::Ptr conn;
        Daemon& daemon;
        network::Peer::HostInfo master;
    };


}   // redux

#endif // REDUX_WORKER_HPP
