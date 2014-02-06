#ifndef REDUX_UTIL_ENDIAN_HPP
#define REDUX_UTIL_ENDIAN_HPP

//#include "redux/util/datautil.hpp"
#include <cstdlib>
#include <utility>

#define REDUX_LITTLE_ENDIAN 1234
#define REDUX_BIG_ENDIAN    4321

// REDUX_BYTE_ORDER defined from cmake as one of the above values, if not -> panic
#ifndef REDUX_BYTE_ORDER
//#error REDUX_BYTE_ORDER is undefined, it should be set from cmake.
#define REDUX_BYTE_ORDER REDUX_LITTLE_ENDIAN
#endif

namespace redux {

    namespace util {

        /*!  @ingroup util
         *  @{
         */

        /*!  @file      endian.hpp
         *   @brief     byte-swapping
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         *   @todo      re-implement using ifdefs & compiler-specific swapping functions.
         */

        
        /*! @fn void swapEndian(T& x)
         *  @brief Reverse the byte-order of x
         *  @param x input
         */
        template <class T> void swapEndian(T& x) {

            size_t size = sizeof(T);

            if(size > 1) {
                size_t mid = size >> 1;
                size--;
                char* p = reinterpret_cast<char*>(&x);

                for(size_t j = 0; j < mid; ++j) {
                    std::swap(p[j], p[size - j]);
                }
            }

        }

        /*! @fn void swapEndian(T* x, size_t n = 1)
         *  @brief Reverse the byte-order of x[n]
         *  @param x input
         *  @param n number of elements to iterate over.
         */
        template <class T> void swapEndian(T* x, size_t n = 1) {

            size_t size = sizeof(T);
            if(size > 1) {
                while(n--) {
                    swapEndian(*x++);
                }
            }

        }


        /*! @} */

        
    } // namespace util

} // namespace redux

#endif // REDUX_UTIL_ENDIAN_HPP
