#include "redux/application.hpp"
#include "redux/logger.hpp"

#include "redux/debugjob.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/translators.hpp"
#include "redux/network/tcpconnection.hpp"
#include "redux/network/host.hpp"
#include "redux/util/arrayutil.hpp"
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


int main( int argc, char *argv[] ) {

    bpo::variables_map vm;
    bpo::options_description programOptions = getOptions();

    bpo::positional_options_description pod;
    pod.add("jobs", -1);
    
    bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions, &pod );

    // load matched environment variables according to the getOptionName() above.
    bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );
    vm.notify();

//     if( !vm.count( "jobs" ) && !vm.count( "all" ) ) {
//         return EXIT_SUCCESS;
//     }
// 
    try {
        string jobString;        
        if( vm.count( "jobs" ) ) {
            jobString = boost::algorithm::join(vm["jobs"].as<vector<string>>(), ",");
        }
        
        if( jobString.empty() ) {
            return EXIT_SUCCESS;
        }
        
        Logger logger( vm );
        boost::asio::io_service ioservice;
        auto conn = TcpConnection::newPtr( ioservice );
        
        conn->connect( vm["master"].as<string>(), vm["port"].as<string>() );

        if( conn->socket().is_open() ) {
            Host::HostInfo me, master;
            
            uint8_t cmd = CMD_CONNECT;
            boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));
            boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));
            if( cmd == CMD_AUTH ) {
                // implement
            }
            if( cmd == CMD_CFG ) {  // handshake requested
                *conn << me;
                *conn >> master;
                boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));       // ok or err
            }
            if( cmd != CMD_OK ) {
                LOG_ERR << "Handshake with server failed.";
                return EXIT_FAILURE;
            }

            size_t stringSize = jobString.length()+1;
            size_t totalSize = stringSize + sizeof( size_t ) + 1;
            auto buf = sharedArray<char>( totalSize );
            char* ptr = buf.get();
            uint64_t count = pack(ptr,CMD_DEL_JOB);
            count += pack(ptr+count,stringSize);
            count += pack(ptr+count,jobString);
            boost::asio::write(conn->socket(),boost::asio::buffer(buf.get(),totalSize));
        }
    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }
    return EXIT_SUCCESS;

}

