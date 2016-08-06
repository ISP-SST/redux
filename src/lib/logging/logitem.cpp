#include "redux/logging/logitem.hpp"

#include "redux/logging/logger.hpp"
#include "redux/util/datautil.hpp"

using namespace redux::logging;

void LogItem::endEntry(void) {
    entry.finalize();
    if( logger) {
        logger->append( *this );
    }
}


uint64_t LogItem::size( void ) const {
    uint64_t sz = entry.size() + context.length() + 1;
    return sz;
}


uint64_t LogItem::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = entry.pack(ptr);
    count += pack(ptr+count, context );
    return count;
}


uint64_t LogItem::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = entry.unpack(ptr, swap_endian);
    count += unpack(ptr+count, context, swap_endian);
    return count;
}
