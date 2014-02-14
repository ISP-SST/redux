#include "redux/application.hpp"
#include "redux/job.hpp"
#include "redux/logger.hpp"
#include "redux/network/tcpconnection.hpp"
#include "redux/network/peer.hpp"
#include "redux/util/endian.hpp"

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

using namespace redux::util;
using namespace redux::network;
using namespace redux;

using namespace std;

#define lg Logger::lg
namespace {

    const string thisChannel = "jstat";

    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "Program Options" );
        options.add_options()
        ( "master,m", bpo::value<string>()->default_value( "localhost" ), "Hostname/IP of the master" )
        ( "port,p", bpo::value<string>()->default_value( "30000" ), "Port to use, either when connecting to the master." )
        ( "jobs,j", "Get joblist" )
        ( "slaves,s", "Get slavelist" )
        ( "time,t", bpo::value<int>()->implicit_value( 1 ), "Loop and display list every (n) seconds" )
        ;

        return options;
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static map<string, string> vmap;
        if( vmap.empty() ) {
            vmap["REDUX_VERBOSITY"] = "verbosity";  // For debugging this might be convenient.
            vmap["REDUX_PORT"] = "port";            // This means the environment variable REDUX_PORT will override the
            // default value of 30000 specified above.
        }
        map<string, string>::const_iterator ci = vmap.find( envName );
        if( ci == vmap.end() ) {
            return "";
        }
        else {
            return ci->second;
        }
    }
}

void printJobList( TcpConnection::ptr conn, bool swapNeeded ) {

    Command cmd = CMD_JSTAT;
    *conn << cmd;
    *conn >> cmd;
    if( cmd != CMD_OK ) {
        LOG_ERR << "Failure while requesting joblist  (server replied = " << cmd << ")";
        return;
    }
    
    size_t sz = sizeof( size_t );
    unique_ptr<char[]> buf( new char[ sz ] );

    LOG_DEBUG << "printJobList(): about to receive " << sz << " bytes.";
    size_t count = boost::asio::read( conn->socket(), boost::asio::buffer( buf.get(), sz ) );
    if( count != sz ) {
        LOG_ERR << "printJobList: Failed to receive blocksize.";
        return;
    }

    sz = *reinterpret_cast<size_t*>( buf.get() );
    if( swapNeeded ) {
        swapEndian( sz );
    }

    buf.reset( new char[ sz ] );
    char* ptr = buf.get();
    char* end = ptr+sz;
    const char* cptr=nullptr;
    try {
        LOG_DEBUG << "printJobList(): about to receive " << sz << " bytes.";
        count = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, sz ) );
        if( count != sz ) {
            LOG_ERR << "printJobList: Failed to receive data block.";
        } else {
            cptr = ptr;
            Job::Info info;
            cout << info.printHeader() << endl;
            while( cptr < end ) {
                cptr = info.unpack(cptr,swapNeeded);
                cout << info.print() << endl;
            }
        }
    } catch ( const exception& e) {
        LOG_ERR << "printJobList: Exception caught while parsing block: " << e.what();
    }
    if( cptr != end ) {
        LOG_ERR << "printJobList: Parsing of datablock failed, there was a missmatch of " << (cptr-end) << " bytes.";
    }

}

void printPeerList( TcpConnection::ptr conn, bool swapNeeded ) {

    Command cmd = CMD_PSTAT;
    *conn << cmd;
    *conn >> cmd;
    if( cmd != CMD_OK ) {
        LOG_ERR << "Failure while requesting peerlist  (server replied = " << cmd << ")";
        return;
    }
    
    size_t sz = sizeof( size_t );
    unique_ptr<char[]> buf( new char[ sz ] );

    LOG_DEBUG << "printPeerList(): about to receive " << sz << " bytes.";
    size_t count = boost::asio::read( conn->socket(), boost::asio::buffer( buf.get(), sz ) );
    if( count != sz ) {
        LOG_ERR << "printPeerList: Failed to receive blocksize.";
        return;
    }

    sz = *reinterpret_cast<size_t*>( buf.get() );
    if( swapNeeded ) {
        swapEndian( sz );
    }

    buf.reset( new char[ sz ] );
    char* ptr = buf.get();
    char* end = ptr+sz;
    const char* cptr=nullptr;
    try {
        LOG_DEBUG << "printPeerList(): about to receive " << sz << " bytes.";
        count = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, sz ) );
        if( count != sz ) {
            LOG_ERR << "printPeerList: Failed to receive data block.";
        } else {
            cptr = ptr;
            Peer peer;
            cout << peer.printHeader() << endl;
            while( cptr < end ) {
                cptr = peer.unpack(cptr,swapNeeded);
                cout << peer.print() << endl;
            }
        }
    } catch ( const exception& e) {
        LOG_ERR << "printPeerList: Exception caught while parsing block: " << e.what();
    }
    if( cptr != end ) {
        LOG_ERR << "printPeerList: Parsing of datablock failed, there was a missmatch of " << (cptr-end) << " bytes.";
    }
}

int main( int argc, char *argv[] ) {

    bpo::variables_map vm;
    bpo::options_description programOptions = getOptions();

    bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions );

    // load matched environment variables according to the getOptionName() above.
    bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );
    vm.notify();

    if( !vm.count( "jobs" ) && !vm.count( "slaves" ) ) {
        cout << allOptions << endl;
        return EXIT_SUCCESS;
    }

    bool loop = vm.count( "time" );
    try {
        Logger logger( vm );
        boost::asio::io_service ioservice;
        TcpConnection::ptr conn = TcpConnection::newPtr( ioservice );
        conn->connect( vm["master"].as<string>(), vm["port"].as<string>() );

        if( conn->socket().is_open() ) {
            Command cmd;
            Peer::HostInfo me, master;
            *conn << CMD_CONNECT;
            *conn >> cmd;
            if( cmd == CMD_AUTH ) {
                // implement
            }
            if( cmd == CMD_CFG ) {  // handshake requested
                *conn << me;
                *conn >> master;
                *conn >> cmd;       // ok or err
            }
            if( cmd != CMD_OK ) {
                LOG_ERR << "Handshake with server failed.";
                return EXIT_FAILURE;
            }
            bool swapNeeded = ( me.littleEndian != master.littleEndian );
            if( vm.count( "jobs" ) ) printJobList( conn, swapNeeded );
            if( vm.count( "slaves" ) ) printPeerList( conn, swapNeeded );
            while( loop ) {
                sleep(vm["time"].as<int>());
                if( vm.count( "jobs" ) ) printJobList( conn, swapNeeded );
                if( vm.count( "slaves" ) ) printPeerList( conn, swapNeeded );
            }
        }
    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }
    return EXIT_SUCCESS;

}

