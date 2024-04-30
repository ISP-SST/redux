#ifndef REDUX_UTIL_ENDIAN_HPP
#define REDUX_UTIL_ENDIAN_HPP

#include <cstdlib>
#include <utility>

#define RDX_LITTLE_ENDIAN 1234
#define RDX_BIG_ENDIAN    4321

// RDX_BYTE_ORDER defined from cmake as one of the above values, if not -> panic
#ifndef RDX_BYTE_ORDER
//#error RDX_BYTE_ORDER is undefined, it should be set from cmake.
#define RDX_BYTE_ORDER RDX_LITTLE_ENDIAN
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

            size_t nBytes = sizeof(T);
            if(nBytes > 1) {
                size_t mid = nBytes >> 1;
                nBytes--;
                char* p = reinterpret_cast<char*>(&x);
                for(size_t j = 0; j < mid; ++j) {
                    std::swap(p[j], p[nBytes - j]);
                }
            }

        }

        /*! @fn void swapEndian(T* x, size_t n = 1)
         *  @brief Reverse the byte-order of x[n]
         *  @param x input
         *  @param n number of elements to iterate over.
         */
        template <class T> void swapEndian(T* x, size_t n = 1) {
            size_t nBytes = sizeof(T);
            if(nBytes > 1) {
                while(n--) {
                    swapEndian(*x++);
                }
            }
        }


        /*! @} */

        
    } // namespace util

} // namespace redux

#endif // REDUX_UTIL_ENDIAN_HPP
