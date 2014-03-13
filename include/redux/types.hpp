#ifndef REDUX_TYPES_HPP
#define REDUX_TYPES_HPP

#include "redux/util/datautil.hpp"

#include <cstdint>

namespace redux {


        /*!  @ingroup redux
         *  @{
         */
        
        /*!  @file      types.hpp
         *   @brief     Generic type definitions, commonly used all over the library
         *   @name      Types
         * 
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2014
         */

        /*!
        * A general 2D point
        */
        template<class T> struct PointType {
            PointType ( T yy=0, T xx=0 ) : x ( xx ), y ( yy ) {}
            size_t size(void) const { return 2*sizeof(T); };
            char* pack(char* ptr) const {
                ptr = redux::util::pack(ptr,x);
                return redux::util::pack(ptr,y);
            }
            const char* unpack(const char* ptr, bool swap_endian=false) {
                ptr = redux::util::unpack(ptr,x,swap_endian);
                return redux::util::unpack(ptr,y,swap_endian);
            }
            T x,y;
        };
        typedef PointType<double> PointD;
        typedef PointType<float> PointF;
        typedef PointType<int> PointI;
        typedef PointType<uint32_t> Point;

        /*!
        * A general 2D rectangle
        */
        template<class T> struct RegionType {
            RegionType ( T y1=0, T x1=0, T y2=0, T x2=0 ) : first(x1,y1), last(x2,y2) {}
            size_t size(void) const { return first.size()+last.size(); };
            char* pack(char* ptr) const {
                ptr = first.pack(ptr);
                return ptr = last.pack(ptr);
            }
            const char* unpack(const char* ptr, bool swap_endian=false) {
                ptr = first.unpack(ptr,swap_endian);
                return ptr = last.unpack(ptr,swap_endian);
            }
            PointType<T> first,last;
        };
        typedef RegionType<double> RegionD;
        typedef RegionType<float> RegionF;
        typedef RegionType<int> RegionI;
        typedef RegionType<uint32_t> Region;

        
        /* @} */

}


#endif // REDUX_TYPES_HPP


