#include "redux/application.hpp"

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


namespace {

    const string logChannel = "rdx_del";

    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "Program Options" );
        options.add_options()
        ( "master,m", bpo::value<string>()->default_value( "localhost" ),
          "Hostname/IP of a master to connect to."
          " The environment variable RDX_HOST can be used to override the default value." )
        ( "port,p", bpo::value<string>()->default_value( "30000" ), "Port to use, either when connecting to the master." )
        ( "user", bpo::value<string>(), "Username (default is to use current user)." )
        ( "force,f", "Bypass safeguards" )
        ( "ids,i", bpo::value< vector<string> >()->composing(), "ID(s) to delete/reset, as ID(s) or name(s)."
        "All additional arguments on the command-line will be interpreted as IDs/names to delete."
        "Comma (or space) separated list of IDs/names is accepted. IDs may also be specified as ranges (e.g 4-17)." )
        ( "kill,k", "Hard-kill slaves (default is to exit after patch is completed)." )
        ( "restart,r", "Restart slaves instead of exiting." )
        ( "slaves,s", "Interpret IDs as slaves, instead of jobs." )
        ;

        return options;
        
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static map<string, string> vmap;
        if( vmap.empty() ) {
            vmap["RDX_VERBOSITY"] = "verbosity";  // For debugging this might be convenient.
            vmap["RDX_HOST"] = "master";        // If it exists, it will override the default value (localhost) above
            vmap["RDX_PORT"] = "port";            // This means the environment variable RDX_PORT will override the
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
    pod.add("ids", -1);
    
    bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions, &pod );

    // load matched environment variables according to the getOptionName() above.
    bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );

#if BOOST_VERSION > 104800  // TODO check which version notify appears in
    vm.notify();
#endif

    try {
        string idString;        
        if( vm.count( "ids" ) ) {
            idString = boost::algorithm::join(vm["ids"].as<vector<string>>(), ",");
        }
        
        if( idString.empty() ) {
            return EXIT_SUCCESS;
        }
        
        boost::asio::io_service ioservice;
        auto conn = TcpConnection::newPtr( ioservice );
        conn->connect( vm["master"].as<string>(), vm["port"].as<string>() );

        if( conn->socket().is_open() ) {
            Host& me = Host::myInfo();
            Host::HostInfo master;
            
            if( vm.count ("user") ) {
                me.info.user = vm["user"].as<string>();
            }
            
            uint8_t cmd = CMD_CONNECT;
            boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));
            boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));
            if( cmd == CMD_AUTH ) {
                // implement
            }
            if( cmd == CMD_CFG ) {  // handshake requested
                if( vm.count("force") ) me.info.peerType = Host::TP_MASTER;
                *conn << me.info;
                *conn >> master;
                boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));       // ok or err
            }
            if( cmd != CMD_OK ) {
                cerr << "Handshake with server failed." << endl;
                return EXIT_FAILURE;
            }

            size_t stringSize = idString.length()+1;
            size_t totalSize = stringSize + sizeof( size_t ) + 2;
            auto buf = sharedArray<char>( totalSize );
            char* ptr = buf.get();
            
            if( vm.count( "slaves" ) ) {
                if( vm.count( "restart" ) ) {
                    cmd = CMD_SLV_RES;
                } else cmd = CMD_DEL_SLV;
            } else cmd = CMD_DEL_JOB;
            uint64_t count = pack(ptr,cmd);
            if( cmd == CMD_DEL_SLV  ) {
                uint8_t hardExit(0);
                if( vm.count( "kill" ) ) hardExit = 1;
                count += pack( ptr+count, hardExit );
            }
            count += pack(ptr+count,stringSize);
            count += pack(ptr+count,idString);
            boost::asio::write( conn->socket(), boost::asio::buffer(buf.get(), count));
        } else {
            cout << "Connection failed: " << vm["master"].as<string>() << ":" << vm["port"].as<string>() << endl;
        }
    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }
    return EXIT_SUCCESS;

}

