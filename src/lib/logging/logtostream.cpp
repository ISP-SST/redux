#include "redux/logging/logtostream.hpp"

#include "redux/util/stringutil.hpp"

#include <boost/date_time/posix_time/time_formatters.hpp>
#include <boost/date_time/c_local_time_adjustor.hpp>

using namespace redux::logging;
using namespace redux::util;
using namespace std;
namespace bfs = boost::filesystem;
namespace bpx = boost::posix_time;

namespace {

    std::map<uint8_t,string> level_tags = {
        { 0, "" },
        { LOG_MASK_FATAL,   "[F]" },
        { LOG_MASK_ERROR,   "[E]" },
        { LOG_MASK_WARNING, "[W]" },
        { LOG_MASK_NOTICE,  "[n]" },
        { LOG_MASK_NORMAL,  "   " },   // empty tag for normal messages.
        { LOG_MASK_DETAIL,  "[d]" },
        { LOG_MASK_DEBUG,   "[D]" },
        { LOG_MASK_TRACE,   "[T]" }
    };
    
    std::map<uint8_t,string> color_level_tags = {
        { 0, "" },
        { LOG_MASK_FATAL,   "\033[31m[F]\033[0m" },
        { LOG_MASK_ERROR,   "\033[91m[E]\033[0m" },
        { LOG_MASK_WARNING, "\033[33m[W]\033[0m" },
        { LOG_MASK_NOTICE,  "\033[93m[n]\033[0m" },
        { LOG_MASK_NORMAL,  "   " },   // empty tag for normal messages.
        { LOG_MASK_DETAIL,  "\033[32m[d]\033[0m" },
        { LOG_MASK_DEBUG,   "\033[36m[D]\033[0m" },
        { LOG_MASK_TRACE,   "\033[95m[T]\033[0m" }
    };

}


LogToStream::LogToStream( ostream &os, uint8_t m, unsigned int flushPeriod) : LogOutput(m,flushPeriod), out(os),
    color(true), localtime(true) {

}


LogToStream::~LogToStream() {
    
    if( !itemQueue.empty() ) {
        flushBuffer();
    }

}


void LogToStream::writeFormatted( const LogItem &i ) {
    
    typedef boost::date_time::c_local_adjustor<bpx::ptime> local_adj;

    if( out.good() ) {
        boost::io::ios_all_saver settings(out);
        if( localtime ) {
            out << to_iso_extended_string( local_adj::utc_to_local(i.entry.getTime()) ) << " ";
        } else {
            out << to_iso_extended_string( i.entry.getTime() ) << " ";
        }
        if( color ) {
            out << color_level_tags[ i.entry.getMask() ];
        } else {
            out << level_tags[ i.entry.getMask() ];
        }

        if( ! i.context.empty() ) {
            out << " (" << i.context << ")";
        }
        out << " " << i.entry.getMessage() << endl;
    }
}



void LogToStream::flushBuffer( void ) {
    
    unique_lock<mutex> lock( queueMutex );
    vector<LogItem> tmpQueue;
    tmpQueue.reserve( itemQueue.size() );
    for( const auto& i : itemQueue ) {
        if( i ) tmpQueue.push_back(*i); 
    }
    itemQueue.clear();
    itemCount = 0;

    for( LogItem& it: tmpQueue ) {
        writeFormatted( it );
    }

}

