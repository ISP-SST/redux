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
          " The environment variable RDX_CACHEDIR will be used as default if it is defined." )
        ( "master,m", po::value<string>()->default_value( "" ), "name or ip of master."
          " If left blank, this instance will start as a master." )
        ( "port,p", bpo::value<uint16_t>()->default_value( 30000 ), "Port to listen on, or connect to."
          " The environment variable RDX_PORT will be used as default if it is defined." )
        ( "threads,t", po::value<uint16_t>()->default_value( 0 ), "max number of threads to use.")
        ( "max-transfers,T", po::value<uint32_t>()->implicit_value( 20 ), "max simultaneous data transfers.")
        ( "foreground,F", "Do not detach/background process.")
        ;

        return options;
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static map<string, string> vmap;
        if( vmap.empty() ) {
            vmap["RDX_CACHEDIR"] = "cache-dir";
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
                if( (vm.count( "log-stdout" ) == 0) && (vm.count( "foreground" ) == 0) ) {
                    if( daemon( 1, 0 ) ) {
                        throw runtime_error( string("Failed to background process: ") + strerror( errno ) );
                    }

                }
                Daemon daemon( vm );
                return daemon.run();
            }
            catch( Application::KillException ) {
                break;
            }
            catch( Application::ResetException ) {
                continue;
            }
        }
    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }

    return EXIT_SUCCESS;

}

