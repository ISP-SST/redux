#include "redux/logger.hpp"

#include <fstream>
#include <iostream>

#include <boost/log/utility/setup/common_attributes.hpp>

using namespace redux;
using namespace std;

redux::util::BoundValue<int> Logger::defaultSeverity( sev_normal, sev_trace, sev_critical );
Logger::logger_t Logger::lg;
Logger::logger_mt Logger::mlg;

#define lg Logger::lg
namespace {
    const string thisChannel = "log";
}

pair<string, string> Logger::customParser( const string& s ) { // custom parser to handle multiple -q/-v flags (e.g. -vvvv)

    if( s.find( "-v" ) == 0 || s.find( "--verbose" ) == 0 ) { //
        int count = std::count( s.begin(), s.end(), 'v' );
        while( count-- ) --defaultSeverity;
    }
    else if( s.find( "-q" ) == 0 || s.find( "--quiet" ) == 0 ) { //
        int count = std::count( s.begin(), s.end(), 'q' );
        while( count-- ) ++defaultSeverity;
    }

    return make_pair( string(), string() );                 // no need to return anything, we handle the verbosity directly.
}


bpo::options_description Logger::getOptions( const string& application_name ) {

    bpo::options_description logging( "Logging Options" );
    logging.add_options()
    ( "verbosity", bpo::value< int >(), "Specify verbosity level" )
    ( "verbose,v", bpo::value<vector<string>>()->implicit_value( vector<string>( 1, "1" ), "" )
      ->composing(), "More output. (ignored if --verbosity is specified)" )
    ( "quiet,q", bpo::value<vector<string>>()->implicit_value( vector<string>( 1, "-1" ), "" )
      ->composing(), "Less output. (ignored if --verbosity is specified)" )

    ( "log-file,L", bpo::value< vector<string> >()->default_value( vector<string>( 1, "" ), "" )
      ->composing(),
      "Print logevents to file. If no log-file name is"
      " specified the log will be written to ./<name>.log,"
      " where <name> is the name given to this application instance." )
    ( "log-stdout,d", "Debug mode. Will write output from all channels to stdout. --log-file can not be used together with this option." )
    ;

    return logging;
}



Logger::Logger( bpo::variables_map& vm ) {

    boost::log::add_common_attributes();

    if( vm.count( "verbosity" ) > 0 ) {         // if --verbosity N is specified, use it.
        defaultSeverity = vm["verbosity"].as<int>();
    }
        
    if( vm.count( "log-stdout" ) ) {
        boost::shared_ptr<LogSink> sink( new StreamSink( cout, defaultSeverity ) );
        logSinks.push_back( sink );
    }
    else if( vm.count( "log-file" ) ) {
        vector<string> logfiles = vm["log-file"].as<vector<string>>();
        bool hasDefault = false;
        for( auto & filename : logfiles ) {
            if( filename == "" ) {
                if( hasDefault ) {
                    continue;
                }
                filename = vm["name"].as<string>() + ".log";
                hasDefault = true;
            }
            addFileLog( filename );
        }
    }

    LOG_DEBUG << "Logging started.";

}

Logger::~Logger( void ) {
    // LOG_DEBUG( lg, thisChannel ) << "Shutting down the logging facility.";
}

void Logger::addFileLog( const string& file ) {

    boost::shared_ptr<LogSink> sink( new FileSink( file, defaultSeverity ) );
    logSinks.push_back( sink );

}

void Logger::addNullLog( void ) {

    boost::shared_ptr<LogSink> sink( new LogSink(sev_none) );
    logSinks.push_back( sink );

}
