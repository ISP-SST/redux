#include "redux/logging/logtostream.hpp"

#include "redux/util/stringutil.hpp"

#include <boost/date_time/posix_time/time_formatters.hpp>
#include <boost/date_time/c_local_time_adjustor.hpp>

using namespace redux::logging;
using namespace redux::util;
using namespace std;

namespace bpx = boost::posix_time;

namespace {

    const char* level_tags[] = { "", "[F]", "[E]", "[W]", "[n]", "   ", "[d]", "[D]", "[T]" };
    const char* color_level_tags[] = {
        "",
        "\033[31m[F]\033[0m",
        "\033[91m[E]\033[0m",
        "\033[33m[W]\033[0m",
        "\033[93m[n]\033[0m",
        "   ",
        "\033[32m[d]\033[0m",
        "\033[36m[D]\033[0m",
        "\033[95m[T]\033[0m"
    };

}


LogToStream::LogToStream( ostream &os, uint8_t m, unsigned int flushPeriod) : LogOutput(m,flushPeriod), out(os),
    color(true), localtime(true) {

}


LogToStream::~LogToStream() {
    
    if( !itemQueue.empty() ) {
        LogToStream::flushBuffer();
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
        uint8_t m = i.entry.getMask();
        if( m ) {
            int pos = ffs(m);
            if( color ) {
                out << color_level_tags[ pos ];
            } else {
                out << level_tags[ pos ];
            }
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

