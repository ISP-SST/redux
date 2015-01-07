#ifndef REDUX_TYPES_HPP
#define REDUX_TYPES_HPP

#include "redux/util/datautil.hpp"

#include <cstdint>
#include <complex>

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
            static size_t size(void) { return 2*sizeof(T); };
            uint64_t pack(char* ptr) const {
                uint64_t count = redux::util::pack(ptr,x);
                count += redux::util::pack(ptr+count,y);
                return count;
            }
            uint64_t unpack(const char* ptr, bool swap_endian=false) {
                uint64_t count = redux::util::unpack(ptr,x,swap_endian);
                count += redux::util::unpack(ptr+count,y,swap_endian);
                return count;
            }
            template <typename U> bool operator==(const PointType<U>& rhs) const { return (x == rhs.x && y == rhs.y); }
            template <typename U> bool operator!=(const PointType<U>& rhs) const { return !(*this == rhs); }
            template <typename U> bool operator<(const PointType<U>& rhs) const { if (x==rhs.x) return (y < rhs.y); return (x < rhs.x); }
            T x,y;
        };
        typedef PointType<double> PointD;
        typedef PointType<float> PointF;
        typedef PointType<int> PointI;
        typedef PointType<uint16_t> Point16;
        typedef PointType<uint32_t> Point;

        /*!
        * A general 2D rectangle
        */
        template<class T> struct RegionType {
            RegionType ( T y1=0, T x1=0, T y2=0, T x2=0 ) : first(x1,y1), last(x2,y2) {}
            static size_t size(void) { return 4*sizeof(T); };
            uint64_t pack(char* ptr) const {
                uint64_t count = first.pack(ptr);
                count += last.pack(ptr+count);
                return count;
            }
            uint64_t unpack(const char* ptr, bool swap_endian=false) {
                uint64_t count = first.unpack(ptr,swap_endian);
                count += last.unpack(ptr+count,swap_endian);
                return count;
            }
            PointType<T> first,last;
        };
        typedef RegionType<double> RegionD;
        typedef RegionType<float> RegionF;
        typedef RegionType<int> RegionI;
        typedef RegionType<uint32_t> Region;

        
        typedef std::complex<double> complex_t;
        
        /* @} */

}


#endif // REDUX_TYPES_HPP


