#include "redux/application.hpp"

#include "redux/logger.hpp"
#include "redux/revision.hpp"
#include "redux/version.hpp"

#include <vector>

#include <boost/filesystem.hpp>
#include <boost/property_tree/info_parser.hpp>

namespace fs = boost::filesystem;
using namespace redux;
using namespace std;

string Application::executableName;

#define lg Logger::mlg
namespace {
    const string thisChannel = "app";
}

void Application::getOptions( po::options_description& options, const string& name ) {

    po::options_description general( "General Options" );
    general.add_options()
    ( "version,V", "Print version information and quit." )
    ( "copyright", "Print copyright information and quit." )
    ( "help,h", "Show command line options and quit." )
    ( "tutorial", "Print full tutorial and quit." )
    ( "sample", "Produce a sample config script on standard output, that uses" " various application features." )
    ;

    po::options_description config( "Configuration" );
    config.add_options()
    ( "settings", "Configuration file to load settings from." )
    ( "aname", po::value<string>()->default_value( name ),
      "Name used to identify this application instance." )
    ( "unique,u", "Prevent any other application-based process with this flag from"
      " using the same name." )
    ;

    options.add( general ).add( config ).add( Logger::getOptions( name ) );
}

pair<string, string> Application::appCmdParser( const string& s ) {
    // Nothing implemented for the Application class, so just return the customParser from the Logger.
    return Logger::customParser( s );
}

po::options_description& Application::parseCmdLine( int argc, const char* const argv[],  po::variables_map& vm,
                                                    po::options_description* programOptions,
                                                    po::positional_options_description *positionalOptions,
                                                    parserFunction custom_parser ) {
    static po::options_description all;
    if( programOptions ) {
        all.add( *programOptions );
    }

    Application::executableName = fs::path( argv[0] ).filename().string();
    getOptions( all, Application::executableName );
    
    po::command_line_parser parser( argc, argv );
    parser.options( all );
    parser.extra_parser( appCmdParser );

    if( positionalOptions ){
        parser.allow_unregistered().positional(*positionalOptions);
    }
    po::store( parser.run(), vm );

    if( custom_parser ) {
        custom_parser( all, vm );
    }

    // If e.g. --help was specified, just dump output and exit.
    checkGeneralOptions( all, vm );

    vm.notify();

    return all;

}


void Application::checkGeneralOptions( po::options_description& desc, po::variables_map& vm ) {

    if( vm.count( "help" ) ) {
        cout << desc << endl;
        exit( 0 );
    }
    if( vm.count( "version" ) ) {
        if( vm.count( "verbose" ) ) {
            cout << "Version:  " << getLongVersionString() << endl;
            cout << "Commited: " << reduxCommitTime << endl;
            cout << "Compiled: " << reduxBuildTime << endl;
        } else {
            cout << getVersionString() << endl;
        }
        exit( 0 );
    }
    if( vm.count( "copyright" ) ) {
        cout << "Not implemented\n";
        exit( 0 );
    }
    if( vm.count( "tutorial" ) ) {
        cout << "Not implemented\n";
        exit( 0 );
    }
    if( vm.count( "sample" ) ) {
        cout << "Not implemented\n";
        exit( 0 );
    }

}


Application::Application( po::variables_map& vm, RunMode rm ) : runMode(rm), returnValue(0), logger(vm) {

    if( vm.count("settings") ) {
        settingsFile = vm["settings"].as<string>();
        if( fs::is_regular(settingsFile) ) {
            LOG_DETAIL << "Loading file \"" <<  settingsFile << "\"";
            pt::read_info( settingsFile, propTree );
        } else {
            LOG_ERR << "Failed to load file \"" << settingsFile << "\", starting with default settings.";
            //throw KillException();
        }
    }
    applicationName = vm["aname"].as<string>();

}


Application::~Application( void ) {

}


string Application::getName( void ) const {
    return applicationName;
}


int Application::run( void ) {

    while( doWork() && !runMode ) {       // keep calling doWork() until it returns false, or shouldStop is raised

    }
    
    switch(runMode) {
        case RESET: throw ResetException();           // this will cause a complete reset (i.e. creation of a new Application instance)
        //case EXIT:  throw KillException();            // exits the loop in main()
        default: ;
    }

    // No work left.
    return returnValue;
}


