#ifndef REDUX_UTIL_CONVERT_HPP
#define REDUX_UTIL_CONVERT_HPP

#include <ctime>
#include <limits>

#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/numeric/conversion/cast.hpp>

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */


        /*!  @file      convert.hpp
         *   @brief     Some basic tools for converting types and values.
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         */


        /*!
         *  @brief      Cast with min/max boundaries.
         *  @param      src Input
         *  @param      min min
         *  @param      max max
         */
        template<typename Target, typename Source>
        Target bound_cast( Source src,
                           Target min = std::numeric_limits<Target>::lowest(),
                           Target max = std::numeric_limits<Target>::max() ) {

            if( min > max ) {
                std::swap( min, max );
            }

            try {
                Target tmp = boost::numeric_cast<Target>( src );
                return std::min( std::max( tmp, min ), max );
            }
            catch( const boost::numeric::negative_overflow & ) {
                return min;
            }
            catch( const boost::numeric::positive_overflow & ) {
                return max;
            }
        }

        /*! Convert a boost ptime into std time_t
         */
        time_t to_time_t(const boost::posix_time::ptime& pt);
        
        
        /*! @} */

    }

}

#endif // REDUX_UTIL_CONVERT_HPP
