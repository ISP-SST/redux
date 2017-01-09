#include "redux/logging/logentry.hpp"

#include "redux/util/convert.hpp"
#include "redux/util/datautil.hpp"

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>
using namespace std;
using namespace redux::logging;
using namespace redux::util;
using namespace boost::posix_time;

uint64_t LogEntry::size( void ) const {
    uint64_t sz = sizeof(time_t) + sizeof(int32_t);   // ptime is converted and transferred as time_t + int32_t (
    sz += message.length() + 2;
    return sz;
}


uint64_t LogEntry::pack( char* ptr ) const {
    using redux::util::pack;
    *ptr = mask;
    uint64_t count = 1;     // mask
    time_t timestamp = redux::util::to_time_t( entryTime );
    count += pack( ptr+count, timestamp );
    int32_t micros = entryTime.time_of_day().fractional_seconds();
    count += pack( ptr+count, micros );
    count += pack( ptr+count, message );
    return count;
}


uint64_t LogEntry::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    mask = *ptr;
    uint64_t count = 1;     // endianflag
    time_t timestamp;
    count += unpack( ptr+count, timestamp, swap_endian );
    int32_t micros(0);
    count += unpack( ptr+count, micros, swap_endian );
    entryTime = boost::posix_time::from_time_t( timestamp );
    entryTime += time_duration(0, 0, 0, micros);
    count += unpack( ptr+count, message );
    return count;
}


std::ostream &operator<<( std::ostream &os, const redux::logging::LogEntry &le ) {
    
    os << le.getMask() << '|' << to_iso_extended_string(le.getTime()) << '|' << le.getMessage();
    
    return os;
}
