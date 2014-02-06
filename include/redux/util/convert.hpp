#ifndef REDUX_UTIL_CONVERT_HPP
#define REDUX_UTIL_CONVERT_HPP

#include <limits>

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
                           Target min = std::numeric_limits<Target>::min(),
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

        /*! @} */

    }

}

#endif // REDUX_UTIL_CONVERT_HPP
