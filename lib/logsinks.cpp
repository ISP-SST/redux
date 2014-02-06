#include "redux/logsinks.hpp"

#include "redux/logger.hpp"
#include "redux/util/stringutil.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/utility/empty_deleter.hpp>

namespace blog = boost::log;
namespace expr = boost::log::expressions;

using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

    const string thisChannel = "log";

    const map<string, severity_level> sevMap = {
        {"trace", sev_trace},
        {"debug", sev_debug},
        {"detail", sev_detail},
        {"normal", sev_normal},
        {"warning", sev_warning},
        {"error", sev_error},
        {"critical", sev_critical},
        {"1", sev_trace},
        {"2", sev_debug},
        {"3", sev_detail},
        {"4", sev_normal},
        {"5", sev_warning},
        {"6", sev_error},
        {"7", sev_critical}
    };


    static const string severity_tags[] = {
        "",
        "[T]",
        "[D]",
        "[d]",
        "",         // don't output severity-tag for normal messages.
        "[W]",
        "[E]",
        "[C]"
    };


    BOOST_LOG_ATTRIBUTE_KEYWORD( severity, "Severity", severity_level )
    BOOST_LOG_ATTRIBUTE_KEYWORD( channel, "Channel", std::string )

    struct severity_tag;

    // The operator is used when putting the severity level to log
    blog::formatting_ostream& operator<<( blog::formatting_ostream& strm,
                                          blog::to_log_manip< severity_level, severity_tag > const& manip ) {

        size_t index = static_cast< size_t >( manip.get() );
        if( index <= sev_critical && index > 0 ) {
            strm << alignLeft( severity_tags[index], 4 );
        }
        else {
            strm << string( 4, ' ' );
        }

        return strm;
    }


    severity_level parseSeverity( string s ) {
        severity_level ret = sev_none;
        auto it = sevMap.find( s );
        if( it != sevMap.end() ) {
            ret = it->second;
        }
        return ret;
    }


    blog::formatter log_format = expr::stream
                                 << expr::format_date_time<boost::posix_time::ptime>( "TimeStamp", "%Y-%m-%d %H:%M:%S.%f" ) << " "
                                 << expr::attr< severity_level, severity_tag >( "Severity" ).or_default( sev_normal )
                                 << expr::if_( expr::has_attr( channel ) ) [ expr::stream << "(" << channel << ") " ]
                                 << expr::message;

}


LogSink::LogSink( const severity_level& sev ) : minSeverity( sev ), truncate( false ) {

}


void LogSink::parseFilter( string f ) {

    channels.clear();

    auto parts = split( f, ":" );
    for( auto & it : parts ) {
        severity_level tmp = parseSeverity( lowercase( it ) );
        if( tmp ) {
            minSeverity = tmp;
        }
        else {
            channels.push_back( it );
        }
    }

}


FileSink::FileSink( const string& s, const severity_level& sev ) : LogSink( sev ) {
    parse( s );
    init();
    LOG_DEBUG << "Started logging to file \"" << fileName << "\".";
}


FileSink::~FileSink( void ) {
    LOG_DEBUG << "Stopping file log: \"" << fileName << "\".";
    blog::core::get()->remove_sink( sink );
    if( sink.get() ) {
        sink->stop();
        sink->flush();
        sink.reset();
    }
}


void FileSink::parse( const string& s ) {

    size_t pos = s.rfind( ":" );
    if( pos == string::npos ) {
        fileName = s;
    }
    else {
        fileName = s.substr( pos + 1 );
        parseFilter( s.substr( 0, pos ) );
    }

}


void FileSink::init( void ) {

    // Create a backend and initialize it with a stream
    boost::shared_ptr< bsinks::text_ostream_backend > backend =
        boost::make_shared< bsinks::text_ostream_backend >();

    ios_base::openmode mode = ios_base::out;
    if( truncate ) {
        mode |= ios_base::trunc;
    }
    else {
        mode |= ios_base::app;
    }
    backend->add_stream( boost::shared_ptr< std::ofstream >( new ofstream( fileName, mode ) ) );

    // Wrap it into the frontend and register in the core
    sink.reset( new sink_type( backend ) );

    if( channels.size() ) {         // this logfile has specific "channels" specified.
        typedef expr::channel_severity_filter_actor< std::string, severity_level > min_severity_filter;
        min_severity_filter min_severity = expr::channel_severity_filter( channel, severity );
        for( auto & it : channels ) {
            min_severity[it] = minSeverity;
        }
        sink->set_filter( min_severity );
    }
    else {
        sink->set_filter( expr::attr< severity_level >( "Severity" ).or_default( sev_normal ) >= minSeverity );
    }

    sink->set_formatter( log_format );

    blog::core::get()->add_sink( sink );

}


StreamSink::StreamSink( std::ostream& s, int sev ) : strm( s ) {
    minSeverity = static_cast<severity_level>( sev );
    init();
    LOG_DEBUG << "Started logging to stream.";
}


StreamSink::~StreamSink( void ) {
    LOG_DEBUG << "Stopping stream log.";
    blog::core::get()->remove_sink( sink );
    if( sink.get() ) {
        sink->stop();
        sink->flush();
        sink.reset();
    }
}


void StreamSink::init( void ) {

    // Create a backend and initialize it with a stream
    boost::shared_ptr< bsinks::text_ostream_backend > backend =
        boost::make_shared< bsinks::text_ostream_backend >();

    backend->add_stream( boost::shared_ptr< std::ostream >( &strm, blog::empty_deleter() ) );

    // Wrap it into the frontend and register in the core
    sink.reset( new sink_type( backend ) );
    sink->set_filter( expr::attr< severity_level >( "Severity" ).or_default( sev_normal ) >= minSeverity );

    sink->set_formatter( log_format );

    blog::core::get()->add_sink( sink );

}


