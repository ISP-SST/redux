#include "redux/application.hpp"

#include "redux/logger.hpp"
#include "redux/version.hpp"

#include <vector>

#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;
using namespace redux;
using namespace std;

string Application::executableName;

#define lg Logger::lg
namespace {
    const string thisChannel = "main";
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
    ( "name", po::value<string>()->default_value( name ),
      "Name used to identify this application instance." )
    ( "unique,u", "Prevent any other application-based process with this flag from"
      " using the same name." )
    ;

    options.add( general ).add( config ).add( Logger::getOptions( name ) );
}

pair<string, string> Application::customParser( const string& s ) {
    // Nothing implemented for the Application class, so just return the customParser from the Logger.
    return Logger::customParser( s );
}

po::options_description& Application::parseCmdLine( int argc, const char* const argv[],  po::variables_map& vm,
                                                    po::options_description* programOptions, parserFunction parser ) {
    static po::options_description all;
    if( programOptions ) {
        all.add( *programOptions );
    }

    Application::executableName = fs::path( argv[0] ).filename().string();
    getOptions( all, Application::executableName );

    po::store( po::command_line_parser( argc, argv )
               .options( all )
               .extra_parser( customParser )
               .run(), vm );

    if( parser ) {
        parser( all, vm );
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
        cout << "Version: " << getVersionString() << endl;
        cout << getLongVersionString() << endl;
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


Application::Application( po::variables_map& vm ) : m_ShouldStop( false ), m_ShouldRestart( true ) {

    m_Name = vm["name"].as<string>();

}


Application::~Application( void ) {

}


void Application::reset( void ) {
    m_ShouldStop = true;
}


void Application::kill( void ) {
    m_ShouldRestart = false;
    m_ShouldStop = true;
}


string Application::getName( void ) const {
    return m_Name;
}


int Application::run( void ) {

    while( !m_ShouldStop && dispatch() ) ;

    if( m_ShouldStop ) {
        if( m_ShouldRestart ) {
            throw ResetException();
        }
        else {
            throw KillException();
        }
    }

    // end of execution.
    return 0;
}


bool Application::dispatch( void ) {
    static size_t count = 0;
    try {
        //usleep( 100000 );
        BOOST_LOG_CHANNEL_SEV( lg, "loop", static_cast<severity_level>( count % 8 ) ) << "tic";
        if( ( count + 1 ) % 30 == 0 ) {
            m_ShouldRestart = false;
        }
        if( ++count % 10 == 0 ) {
            m_ShouldStop = true;
            return false;
        }
        return true; //eventQueue.dispatch();
    }
    catch( ... ) {
        LOG_DEBUG << "Event dispatch interrupted";
    }

    return true;
}

