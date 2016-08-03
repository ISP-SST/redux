#include "redux/daemon.hpp"

#include "redux/logging/logger.hpp"
#include "redux/debugjob.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/network/protocol.hpp"

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

using namespace redux::logging;
using namespace redux;
using namespace std;


namespace {

    const string logChannel = "reduxd";
    
    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "REDUXd Options" );
        options.add_options()
        ( "cache-dir,C", po::value<string>()->default_value( "" ), "path where to store auxiliary data."
          " The environment variable RDX_CACHEDIR will be used instead of the default if this option"
          " is not specified." )
        ( "master,m", po::value<string>()->default_value( "" ), "name or ip of master."
          " If left blank, this instance will start as a master." )
        ( "port,p", bpo::value<uint16_t>()->default_value( 30000 ), "Port to listen on, or connect to."
          " The environment variable RDX_PORT will be used instead of the default if this option"
          " is not specified." )
        ( "threads,t", po::value<uint16_t>()->default_value( 0 ), "max number of threads to use.")
        ;

        return options;
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static map<string, string> vmap;
        if( vmap.empty() ) {
            vmap["RDX_CACHEDIR"] = "cache-dir";
            vmap["RDX_VERBOSITY"] = "verbosity";  // For debugging this might be convenient.
            vmap["RDX_PORT"] = "port";            // This means the environment variable RDX_PORT will override the
                                                  // default value of 30000 specified above.
        }
        map<string, string>::const_iterator ci = vmap.find( envName );
        if( ci == vmap.end() ) {
            return "";
        } else {
            return ci->second;
        }
    }
}


int main( int argc, char *argv[] ) {

    bpo::variables_map vm;
    bpo::options_description programOptions = getOptions();

    bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions );

    // load matched environment variables according to the getOptionName() above.
    bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );
    
#if BOOST_VERSION > 104800  // TODO check which version notify appears in
    vm.notify();
#endif

    try {
        while( true ) {
            try {
                Daemon daemon( vm );
                int res = daemon.run();
                return res;
            }
            catch( Application::KillException ) {
                cerr << "Application terminated by kill request." << endl;
                break;
            }
            catch( Application::ResetException ) {
                cout << "Application restarting." << endl;
                continue;
            }
        }
    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }

    return EXIT_SUCCESS;

}

