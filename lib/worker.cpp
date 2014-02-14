#include "redux/worker.hpp"

#include "redux/daemon.hpp"
#include "redux/logger.hpp"

using namespace redux::network;
using namespace redux;
using namespace std;

#define lg Logger::lg
namespace {
    const std::string thisChannel = "worker";

}


Worker::Worker( Daemon& rd ) : daemon( rd ) {

}


Worker::~Worker( void ) {

}


void Worker::init( void ) {
    string portstring = to_string( daemon.params["port"].as<uint16_t>() );
    LOG_DEBUG << "init port = " << portstring << "  master = \"" << daemon.params["master"].as<string>() << "\"";
    conn = TcpConnection::newPtr( daemon.ioService );
    conn->connect( daemon.params["master"].as<string>(), portstring );
    if( conn->socket().is_open() ) {
        Command cmd;

        *conn << CMD_CONNECT;
        *conn >> cmd;
        if( cmd == CMD_AUTH ) {
            // implement
        }
        if( cmd == CMD_CFG ) {  // handshake requested
             *conn << daemon.myInfo.host;
             *conn >> master;
             *conn >> cmd;       // ok or err
        }
        if( cmd != CMD_OK ) {
            LOG_ERR << "Handshake with master failed  (server replied: " << cmd << ")";
            conn->socket().close();
        }
    }

}

