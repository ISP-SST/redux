#include "redux/application.hpp"
#include "redux/job.hpp"
#include "redux/logger.hpp"
#include "redux/network/tcpconnection.hpp"
#include "redux/network/host.hpp"
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

void printJobList( TcpConnection::Ptr conn ) {

    uint8_t cmd = CMD_JSTAT;
    boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));

    size_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if ( !blockSize ) return;

    const char* ptr = buf.get();
    uint64_t count(0);
    try {
        Job::Info info;
        cout << info.printHeader() << endl;
        while( count < blockSize ) {
            count += info.unpack(ptr+count,conn->getSwapEndian());
            cout << info.print() << endl;
        }
    } catch ( const exception& e) {
        LOG_ERR << "printJobList: Exception caught while parsing block: " << e.what();
    }
    if( count != blockSize ) {
        LOG_ERR << "printJobList: Parsing of datablock failed,  count = " << count << "   blockSize = " << blockSize << "  bytes.";
    }

}

void printPeerList( TcpConnection::Ptr conn ) {

    uint8_t cmd = CMD_PSTAT;
    boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));

    size_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if ( !blockSize ) return;

    const char* ptr = buf.get();
    uint64_t count(0);
    try {
        Host peer;
        cout << peer.printHeader() << endl;
        while( count < blockSize ) {
            count += peer.unpack(ptr+count,conn->getSwapEndian());
            cout << peer.print() << endl;
        }
    } catch ( const exception& e) {
        LOG_ERR << "printPeerList: Exception caught while parsing block: " << e.what();
    }
    if( count != blockSize ) {
        LOG_ERR << "printPeerList: Parsing of datablock failed, count = " << count << "  blockSize = " << blockSize << "  bytes.";
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
        TcpConnection::Ptr conn = TcpConnection::newPtr( ioservice );
        conn->connect( vm["master"].as<string>(), vm["port"].as<string>() );

        if( conn->socket().is_open() ) {
            Host::HostInfo me;
            Host master;
            uint8_t cmd = CMD_CONNECT;
            boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));
            boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));
            if( cmd == CMD_AUTH ) {
                // implement
            }
            if( cmd == CMD_CFG ) {  // handshake requested
                *conn << me;
                *conn >> master.info;
                boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));       // ok or err
            }
            if( cmd != CMD_OK ) {
                LOG_ERR << "Handshake with server failed.";
                return EXIT_FAILURE;
            }

            if( vm.count( "jobs" ) ) printJobList( conn );
            if( vm.count( "slaves" ) ) printPeerList( conn );
            while( loop ) {
                sleep(vm["time"].as<int>());
                if( vm.count( "jobs" ) ) printJobList( conn );
                if( vm.count( "slaves" ) ) printPeerList( conn );
            }
        }
    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }
    return EXIT_SUCCESS;

}

