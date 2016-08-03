#include "redux/logging/logitem.hpp"

#include "redux/logging/logger.hpp"

using namespace redux::logging;

void LogItem::endEntry(void) {
    entry.finalize();
    if( logger) {
        logger->append( *this );
    }
}


uint64_t LogItem::size( void ) const {
    uint64_t sz = entry.size() + context.size();
    return sz;
}


uint64_t LogItem::pack( char* ptr ) const {
    uint64_t count = entry.pack(ptr);
    count += context.pack(ptr+count);
    return count;
}


uint64_t LogItem::unpack( const char* ptr, bool swap_endian ) {
    uint64_t count = entry.unpack(ptr, swap_endian);
    count += context.unpack(ptr+count, swap_endian);
    return count;
}
