#ifndef REDUX_MATH_HELPERS_HPP
#define REDUX_MATH_HELPERS_HPP

#include <limits>
#include <utility>      // std::swap

#include <boost/numeric/conversion/cast.hpp>

/**
* @file
* Some general purpose functions.
*
* @author Tomas Hillberg
*/


namespace redux {

    namespace math {


        /**
         * Almost eq. function that can be used to compare floats.
         * Ref: http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
         * maxUlps = max no. floting point possitions (tic's) away,
         *   must be in range (maxUlps > 0 && maxUlps < 4 * 1024 * 1024)
         */
        bool almostEqual( float A, float B, int maxUlps = 1 );

        /*!
         *  @brief      Cast with min/max boundaries.
         *  @param      src Input
         *  @param      min min
         *  @param      max max
         */
        template<typename Target, typename Source>
        Target bound_cast(Source src,
                          Target min=std::numeric_limits<Target>::min(),
                          Target max=std::numeric_limits<Target>::max()) {
            
            if ( min > max ) {
                std::swap( min, max );
            }

            try {
                Target tmp = boost::numeric_cast<Target>(src);
                return std::min(std::max(tmp,min),max);
            }
            catch (const boost::numeric::negative_overflow &) {
                return min;
            }
            catch (const boost::numeric::positive_overflow &) {
                return max;
            }
        }
        
        
    }
}

#endif // REDUX_MATH_HELPERS_HPP
