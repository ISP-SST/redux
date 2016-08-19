#ifndef REDUX_UTIL_REGION_HPP
#define REDUX_UTIL_REGION_HPP

#include <iostream>

#include "redux/util/point.hpp"

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */

        /*!  @file      region.hpp
         *   @brief     A general 2D rectangle
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2016
         */

        template<class T> struct RegionType {
            RegionType ( T firsty, T firstx, T lasty, T lastx ) : first(firsty,firstx), last(lasty,lastx) {}
            RegionType ( T lasty, T lastx ) : first(0,0), last(lasty,lastx) {}
            RegionType ( void ) : first(0,0), last(0,0) {}
            RegionType ( RegionType<T>&& rhs ) : first(std::move(rhs.first)), last(std::move(rhs.last)) {}
            RegionType ( const RegionType<T>& rhs ) : first(rhs.first), last(rhs.last) {}
            template <typename U> RegionType ( const RegionType<U>& rhs ) : first(rhs.first), last(rhs.last) {}
            static inline uint64_t size(void) { return 4*sizeof(T); };
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
            void shiftX(T x) { first.x += x; last.x += x; }
            void shiftY(T y) { first.y += y; last.y += y; }
            template <typename U> PointI outside(const RegionType<U>& rhs) const {
                PointI res(0,0);
                if(rhs.first.x<first.x) res.x = rhs.first.x-first.x;
                else if(rhs.last.x>last.x) res.x = rhs.last.x-last.x;
                if(rhs.first.y<first.y) res.y = rhs.first.y-first.y;
                else if(rhs.last.y>last.y) res.y = rhs.last.y-last.y;
                return std::move(res);
            }
            template <typename U> void grow(const RegionType<U>& rhs) { first=first.min(rhs.first); last=last.max(rhs.last); }
            template <typename U> void grow(const PointType<U>& rhs) { first-=rhs; last+=rhs; }
            void grow(T rhs) { PointType<T> tmp(rhs,rhs); first-=tmp; last+=tmp; }
            template <typename U> RegionType<T> grown(const PointType<U>& rhs) const { RegionType<T> res(*this); res.grow(rhs); return std::move(res); }
            RegionType<T> grown(T rhs) const { RegionType<T> res(*this); res.grow(PointType<T>(rhs,rhs)); return std::move(res); }
            template <typename U> void restrict(const RegionType<U>& rhs) { first=first.max(rhs.first); last=last.min(rhs.last); }
            operator std::string() const { return (std::string)first + "->" + (std::string)last; }
            template <typename U> RegionType<T>& operator+=(const PointType<U>& rhs) { first += rhs; last += rhs; return *this; }
            template <typename U> RegionType<T>& operator-=(const PointType<U>& rhs) { first -= rhs; last -= rhs; return *this; }
            template <typename U> RegionType<T> operator+(const PointType<U>& rhs) const { RegionType<T> res(*this); return std::move(res+=rhs); }
            template <typename U> RegionType<T> operator-(const PointType<U>& rhs) const { RegionType<T> res(*this); return std::move(res-=rhs); }
            RegionType<T>& operator=(const RegionType<T>& rhs) { first = rhs.first; last = rhs.last; return *this; }
            template <typename U> RegionType<T>& operator=(const RegionType<U>& rhs) { first = rhs.first; last = rhs.last; return *this; }
            RegionType<T>& operator=(T rhs) { first = rhs; last = rhs; return *this; }
            template <typename U> bool operator==(const RegionType<U>& rhs) const { return (first == rhs.first && last == rhs.last); }
            bool operator==(T rhs) const { return (first == rhs && last == rhs); }
            template <typename U> bool operator!=(const RegionType<U>& rhs) const { return !(*this == rhs); }
            bool operator!=(T rhs) const { return !(*this == rhs); }
            template <typename U> bool operator<(const RegionType<U>& rhs) const { if (first==rhs.first) return (last < rhs.last); return (first < rhs.first); }
            PointType<T> first,last;
        };
        typedef RegionType<double> RegionD;
        typedef RegionType<float> RegionF;
        typedef RegionType<int> RegionI;
        typedef RegionType<uint16_t> Region16;
        typedef RegionType<uint32_t> Region;
        template<typename T>
        std::ostream& operator<<(std::ostream& os, const RegionType<T>& rt) {
            os << (std::string)rt;
            return os;
        }


        /*! @} */

    }

}

#endif // REDUX_UTIL_REGION_HPP
