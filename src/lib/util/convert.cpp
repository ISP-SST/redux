#include "redux/util/convert.hpp"

time_t redux::util::to_time_t( const boost::posix_time::ptime& pt ) {
    
    if( pt.is_special() ) {
        return time_t(-1);
    }

    using namespace boost::posix_time;
    static ptime epoch( boost::gregorian::date( 1970, 1, 1 ) );
    time_duration::sec_type x = (pt - epoch).total_seconds();

    return time_t(x);
    
}
