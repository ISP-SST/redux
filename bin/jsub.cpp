#include "redux/application.hpp"
#include "redux/logger.hpp"
#include "redux/debugjob.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/network/tcpconnection.hpp"
#include "redux/network/peer.hpp"
#include "redux/util/stringutil.hpp"

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/program_options.hpp>
namespace bpo = boost::program_options;
namespace bpt = boost::property_tree;

using namespace redux::util;
using namespace redux::network;
using namespace redux;

using namespace std;


#define lg Logger::lg
namespace {

    const string thisChannel = "jsub";

    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "Program Options" );
        options.add_options()
        ( "master,m", bpo::value<string>()->default_value( "localhost" ),
          "Hostname/IP of a master to connect to."
          " The environment variable REDUX_MASTER can be used to override the default value." )
        ( "port,p", bpo::value<string>()->default_value( "30000" ),
          "Port to use when connecting to a master."
          " The environment variable REDUX_PORT can be used to override the default value." )
        ( "force,f", "Overwrite output file if exists" )
        ( "swap,s", "swap mode: write compressed data to swap file instead of keeping it in memory (useful for large problems)" )
        ( "config,c", bpo::value<string>()->default_value( "momfbd.cfg" ), "Configuration file to process." )
        ( "simx", "x coordinate[s] of subimages to restore" )
        ( "simy", "y coordinate[s] of subimages to restore" )
        ( "nimages,n", bpo::value<int>(), "Image numbers" )
        ( "sequence", "sequence number to insert in filename template." )
        ( "print,P", "(debug) print the parsed configuration to console and exit without uploading." )
        ( "output-files,o", "Comma separated list of output file base names."
          "File names are applied to the objects in the order they are found in the config file."
          "If insufficient names are provided, a default name will be created for the remaining objects."
          "Excess names will generate a warning but are ignored otherwise. The names are base names only,"
          "an appropriate suffix will be attached (.fits/.f0)." )
        ;

        return options;
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static map<string, string> vmap;
        if( vmap.empty() ) {
            vmap["REDUX_VERBOSITY"] = "verbosity";  // For debugging this might be convenient.
            vmap["REDUX_MASTER"] = "master";        // If it exists, it will override the default value (localhost) above
            vmap["REDUX_PORT"] = "port";            // If it exists, it will override the default value (30000) above
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

void uploadJobs(TcpConnection::Ptr conn, vector<Job::JobPtr>& jobs) {
    
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
        return;
    }
    // all ok, upload jobs.
    size_t sz = 0;
    for( auto &it: jobs ) {
        sz += it->size();
    }
    unique_ptr<char[]> buf( new char[ sz+sizeof(size_t)+1 ] );
    char* ptr = buf.get();
    *ptr++ = CMD_ADD_JOB;
    *reinterpret_cast<size_t*>( ptr ) = sz;
    ptr += sizeof(size_t);
    sz += sizeof(size_t)+1;
    for( auto &it: jobs ) {
        ptr = it->pack(ptr);
    }
    if( buf.get() + sz != ptr ) {
        LOG_ERR << "Packing of jobs failed:  there is a mismatch of " << (ptr-buf.get()-sz) << " bytes.";
        return;
    }

    conn->writeAndCheck( buf.get(), sz );
    *conn >> cmd;
    if( cmd != CMD_OK ) {
        LOG_ERR << "Failure while sending jobs  (server reply = " << cmd << ")";
        return;
    }
    ptr = buf.get();
    size_t* idPtr = reinterpret_cast<size_t*>( ptr );
    size_t count = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, sizeof(size_t) ) );
    if( count == sizeof(size_t) ) {
        count = boost::asio::read( conn->socket(), boost::asio::buffer( idPtr+1, (*idPtr)*sizeof(size_t) ) );
    } else *idPtr = 0;
    
    LOG << "Upload of " << *idPtr << " jobs completed successfully";
    if( *idPtr ) {
        LOG << printArray(idPtr+1,*idPtr,"Job IDs ");
    }

}

int main( int argc, char *argv[] ) {

    bpo::variables_map vm;
    bpo::options_description programOptions = getOptions();

    bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions );

    // load matched environment variables according to the getOptionName() above.
    bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );
    vm.notify();

    try {
        Logger logger( vm );
        bpt::ptree momfbd;
        bpt::read_info( vm["config"].as<string>(), momfbd );
        vector<Job::JobPtr> jobs = Job::parseTree( vm, momfbd );

        if( vm.count( "print" ) ) {     // dump configuration to console and exit
            bpt::ptree dump;
            for( auto & it : jobs ) {
                it->getPropertyTree( &dump );
            }
            bpt::write_info( cout<<endl, dump );
            return EXIT_SUCCESS;
        }

        boost::asio::io_service ioservice;
        TcpConnection::Ptr conn = TcpConnection::newPtr(ioservice);
        conn->connect( vm["master"].as<string>(), vm["port"].as<string>() );

        if( conn->socket().is_open() ) {
            uploadJobs(conn, jobs);
        }

    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }

    return EXIT_SUCCESS;

}

