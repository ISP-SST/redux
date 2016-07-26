#ifndef REDUX_UTIL_POINT_HPP
#define REDUX_UTIL_POINT_HPP

#include <cstdio>
#include <iostream>

#include "redux/util/datautil.hpp"

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */

        /*!  @file      point.hpp
         *   @brief     A general 2D point
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2016
         */

        template<typename T> struct PointType {
            explicit PointType ( T yy=0, T xx=0 ) : x(xx), y(yy) {}
            PointType ( PointType<T>&& rhs ) : x(std::move(rhs.x)), y(std::move(rhs.y)) {}
            PointType ( const PointType<T>& rhs ) : x(rhs.x), y(rhs.y) {}
            template <typename U> PointType ( const PointType<U>& rhs ) : x(rhs.x), y(rhs.y) {}
            static uint64_t size(void) { return 2*sizeof(T); };
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
            template <typename U> PointType<T> max(const PointType<U>& rhs) {
                PointType<T> tmp(*this);
                tmp.x=std::max<T>(tmp.x,rhs.x);
                tmp.y=std::max<T>(tmp.y,rhs.y);
                return std::move(tmp); }
            template <typename U> PointType<T> min(const PointType<U>& rhs) {
                PointType<T> tmp(*this);
                tmp.x=std::min<T>(tmp.x,rhs.x);
                tmp.y=std::min<T>(tmp.y,rhs.y);
                return std::move(tmp); }
            //operator std::string() const { std::ostringstream out; out << "(" << y << "," << x << ")"; return out.str(); }
            operator std::string() const { return "(" + std::to_string(y) + "," + std::to_string(x) + ")"; }
            template <typename U> PointType<T>& operator+=(const PointType<U>& rhs) { x += rhs.x; y += rhs.y; return *this; }
            template <typename U> PointType<T>& operator-=(const PointType<U>& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
            template <typename U> PointType<T> operator+(const PointType<U>& rhs) const { PointType<T> res(*this); return std::move(res+=rhs); }
            template <typename U> PointType<T> operator-(const PointType<U>& rhs) const { PointType<T> res(*this); return std::move(res-=rhs); }
            PointType<T>& operator=(const PointType<T>& rhs) { x = rhs.x; y = rhs.y; return *this; }
            template <typename U> PointType<T>& operator=(const PointType<U>& rhs) { x = rhs.x; y = rhs.y; return *this; }
            PointType<T>& operator=(const T& rhs) { x = rhs; y = rhs; return *this; }
            PointType<T>& operator+=(const T& rhs) { x += rhs; y += rhs; return *this; }
            PointType<T>& operator-=(const T& rhs) { x -= rhs; y -= rhs; return *this; }
            PointType<T>& operator*=(const T& rhs) { x *= rhs; y *= rhs; return *this; }
            PointType<T> operator+(const T& rhs) const { PointType<T> tmp(*this); tmp += rhs; return std::move(tmp); }
            PointType<T> operator-(const T& rhs) const { PointType<T> tmp(*this); tmp -= rhs; return std::move(tmp); }
            PointType<T> operator*(const T& rhs) const { PointType<T> tmp(*this); tmp *= rhs; return std::move(tmp); }
            template <typename U> bool operator==(const PointType<U>& rhs) const { return (x == rhs.x && y == rhs.y); }
            bool operator==(T rhs) const { return (x == rhs && y == rhs); }
            template <typename U> bool operator!=(const PointType<U>& rhs) const { return !(*this == rhs); }
            bool operator!=(T rhs) const { return !(*this == rhs); }
            template <typename U> bool operator<(const PointType<U>& rhs) const { if (y==rhs.y) return (x < rhs.x); return (y < rhs.y); }
            T x,y;
        };
        typedef PointType<double> PointD;
        typedef PointType<float> PointF;
        typedef PointType<int> PointI;
        typedef PointType<uint16_t> Point16;
        typedef PointType<uint32_t> Point;
        template<typename T>
        std::ostream& operator<<(std::ostream& os, const PointType<T>& pt) {
            os << (std::string)pt;
            return os;
        }
        

        /*! @} */

    }

}

#endif // REDUX_UTIL_POINT_HPP
