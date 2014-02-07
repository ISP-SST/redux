#include "redux/application.hpp"
#include "redux/logger.hpp"

#include "redux/debugjob.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/serialization.hpp"              // for serializing std::shared_ptr

#include <iostream>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/serialization/vector.hpp>       // for serializing a vector

using namespace redux::util;
using namespace redux;

using namespace std;

namespace bpo = boost::program_options;
namespace bpt = boost::property_tree;

#define lg Logger::lg
namespace {

    const string thisChannel = "jsub";

    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "Program Options" );
        options.add_options()
        ( "master,m", bpo::value<string>()->implicit_value( "" ),
          "Hostname/IP of a master to connect to (i.e. this instance is a \"slave\")." )
        ( "port,p", bpo::value<uint16_t>()->default_value( 30000 ),
          "Port to use, either when connecting to a master, or to bind to if this instance is a master."
          " The environment variable REDUX_PORT will be used instead of the default if this option"
          " is not specified." )
        ( "force,f", "Overwrite output file if exists" )
        ( "swap,s", "swap mode: write compressed data to swap file instead of keeping it in memory (useful for large problems)" )
        ( "config,c", bpo::value<string>()->default_value( "momfbd.cfg" ), "Configuration file to process." )
        ( "simx", "x coordinate[s] of subimages to restore" )
        ( "simy", "y coordinate[s] of subimages to restore" )
        ( "nimages,n", bpo::value<int>(), "Image numbers" )
        ( "sequence", "sequence number to insert in filename template." )
        ( "print,P", "(debug) print configuration to console." )
        ( "serialize", "(debug) print serialized data to console." )
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

    bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions );

    // load matched environment variables according to the getOptionName() above.
    bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );
    vm.notify();

    try {
        Logger logger( vm );
        bpt::ptree momfbd;
        bpt::read_info( vm["config"].as<string>(), momfbd );
        vector<Job::JobPtr> jobs = Job::parseTree( vm, momfbd );

        if( vm.count( "print" ) ) {
            bpt::ptree dump;     // dump configuration to console
            for( auto & it : jobs ) {
                it->getProperties( &dump );
            }
            bpt::write_info( cout << endl, dump );
        }

        if( vm.count( "serialize" ) ) {
            // Serialize the data and print it (text_archive just for testing, it will be a binary archive)
            std::ostringstream archive_stream;
            boost::archive::text_oarchive archive( archive_stream );
            archive << jobs;
            cout << "[" << archive_stream.str() << "]\n";
        }


    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }

    return EXIT_FAILURE;

}

