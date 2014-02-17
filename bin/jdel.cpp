#include "redux/application.hpp"
#include "redux/logger.hpp"

#include "redux/debugjob.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/translators.hpp"
#include "redux/network/tcpconnection.hpp"
#include "redux/network/peer.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <iostream>

#include <boost/asio.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/serialization/vector.hpp>       // for serializing a vector

using namespace redux::util;
using namespace redux::network;
using namespace redux;

using namespace std;

namespace bpo = boost::program_options;
namespace bpt = boost::property_tree;

#define lg Logger::lg
namespace {

    const string thisChannel = "jdel";

    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "Program Options" );
        options.add_options()
        ( "master,m", bpo::value<string>()->default_value( "localhost" ), "Hostname/IP of the master" )
        ( "port,p", bpo::value<string>()->default_value( "30000" ), "Port to use, either when connecting to the master." )
        ( "jobs,j", po::value< vector<string> >()->composing(), "Job(s) to delete, as job-ID(s) or job-name(s)."
        "All additional arguments on the command-line will be interpreted as jobs to delete."
        "Comma (or space) separated list of IDs/names is accepted. IDs may also be specified as ranges (e.g 4-17)." )
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

void printJobList( TcpConnection::Ptr conn, bool swap_endian ) {

    Command cmd = CMD_JSTAT;
    *conn << cmd;
    *conn >> cmd;
    if( cmd != CMD_OK ) {
        LOG_ERR << "Failure while requesting joblist  (server replied = " << cmd << ")";
        return;
    }
    
    size_t sz = sizeof( size_t );
    unique_ptr<char[]> buf( new char[ sz ] );

    LOG_DEBUG << "sendJobList(): about to receive " << sz << " bytes.";
    size_t count = boost::asio::read( conn->socket(), boost::asio::buffer( buf.get(), sz ) );
    if( count != sz ) {
        LOG_ERR << "printJobList: Failed to receive blocksize.";
        return;
    }

    sz = *reinterpret_cast<size_t*>( buf.get() );
    if( swap_endian ) {
        swapEndian( sz );
    }

    buf.reset( new char[ sz ] );
    char* ptr = buf.get();
    char* end = ptr+sz;
    const char* cptr=nullptr;
    try {
        LOG_DEBUG << "sendJobList(): about to receive " << sz << " bytes.";
        count = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, sz ) );
        if( count != sz ) {
            LOG_ERR << "printJobList: Failed to receive data block.";
        } else {
            cptr = ptr;
            Job::Info info;
            cout << info.printHeader() << endl;
            while( cptr < end ) {
                cptr = info.unpack(cptr,swap_endian);
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

void printSlaveList( TcpConnection::Ptr conn, bool swap_endian ) {
}

int main( int argc, char *argv[] ) {

    bpo::variables_map vm;
    bpo::options_description programOptions = getOptions();

    bpo::positional_options_description pod;
    pod.add("jobs", -1);
    
    bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions, &pod );

    // load matched environment variables according to the getOptionName() above.
    bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );
    vm.notify();

    if( !vm.count( "jobs" ) && !vm.count( "all" ) ) {
        cout << allOptions << endl;
        return EXIT_SUCCESS;
    }

    try {
        Logger logger( vm );
        boost::asio::io_service ioservice;
        TcpConnection::Ptr conn = TcpConnection::newPtr( ioservice );
        string jobString;        
        
        if( vm.count( "jobs" ) ) {
            jobString = boost::algorithm::join(vm["jobs"].as<vector<string>>(), ",");
        }
        
        if( jobString.empty() ) {
            cout << allOptions << endl;
            return EXIT_SUCCESS;
        }
        
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
            size_t stringSize = jobString.length()+1;
            size_t totalSize = stringSize + sizeof( size_t ) + 1;
            unique_ptr<char[]> buf( new char[ totalSize ] );
            char* ptr = buf.get();
            ptr = pack(ptr,CMD_DEL_JOB);
            ptr = pack(ptr,stringSize);
            ptr = pack(ptr,jobString);
            conn->writeAndCheck( buf.get(), totalSize );
        }
    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }
    return EXIT_SUCCESS;

}

