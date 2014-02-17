#include "redux/worker.hpp"

#include "redux/daemon.hpp"
#include "redux/logger.hpp"
#include "redux/network/protocol.hpp"
#include "redux/util/datautil.hpp"

using namespace redux::network;
using namespace redux::util;
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

        daemon.myInfo->host.peerType |= Peer::PEER_WORKER;
        *conn << CMD_CONNECT;
        *conn >> cmd;
        if( cmd == CMD_AUTH ) {
            // implement
        }
        if( cmd == CMD_CFG ) {  // handshake requested
            *conn << daemon.myInfo->host;
            *conn >> master;
            *conn >> cmd;       // ok or err
        }
        if( cmd != CMD_OK ) {
            LOG_ERR << "Handshake with master failed  (server replied: " << cmd << ")";
            conn->socket().close();
            daemon.myInfo->host.peerType &= ~Peer::PEER_WORKER;
        } else {
            conn->setCallback( bind( &Daemon::activity, &daemon, std::placeholders::_1 ) );
        }
    }

}

void Worker::stop( void ) {
    if( conn->socket().is_open() ) {
        *conn << CMD_DISCONNECT;
        conn->socket().close();
    }
    daemon.myInfo->host.peerType &= ~Peer::PEER_WORKER;
}


void Worker::updateStatus( void ) {

    if( !conn->socket().is_open() ) {
        init();
    }

    if( conn->socket().is_open() ) {
        LOG_TRACE << "Updating workerStatus...";
        daemon.myInfo->stat.nThreads++;
        size_t blockSize = daemon.myInfo->stat.size();
        size_t totSize = blockSize + sizeof( size_t ) + 1;
        std::unique_ptr<char[]> buf( new char[totSize] );
        char* ptr = buf.get();
        memset( ptr, 0, totSize );

        ptr = pack( ptr, CMD_STAT );
        ptr = pack( ptr, blockSize );
        ptr = daemon.myInfo->stat.pack( ptr );

        conn->writeAndCheck( buf.get(), totSize );
    } else {
        LOG_WARN << "No connection to master.";
    }

}

