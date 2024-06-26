#ifndef REDUX_UTIL_POINT_HPP
#define REDUX_UTIL_POINT_HPP

#include <cstdio>
#include <sstream>
#include <math.h>

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
            static inline uint64_t size(void) { return 2*sizeof(T); };
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
                return tmp; }
            template <typename U> PointType<T> min(const PointType<U>& rhs) {
                PointType<T> tmp(*this);
                tmp.x=std::min<T>(tmp.x,rhs.x);
                tmp.y=std::min<T>(tmp.y,rhs.y);
                return tmp; }
            double avg( void ) const { return static_cast<double>(x+y)*0.5; }
            T max( void ) const { return std::max<T>(x,y); }
            T min( void ) const { return std::min<T>(x,y); }
            PointType<T> Abs( void ) const { return PointType<T>(redux::util::Abs<T>(y),redux::util::Abs<T>(x)); }
            T maxAbs( void ) const { return std::max<T>(redux::util::Abs<T>(x),redux::util::Abs<T>(y)); }
            PointType<long int> round(void) const { return PointType<long int>( lround(y), lround(x) ); }
            PointType<double> remainder(void) const { PointType<double> rem(*this); rem -= round(); return rem; }
            operator std::string() const { std::ostringstream out; out << "(" << y << "," << x << ")"; return out.str(); }
            template <typename U> PointType<T>& operator+=(const PointType<U>& rhs) { x += rhs.x; y += rhs.y; return *this; }
            template <typename U> PointType<T>& operator-=(const PointType<U>& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
            template <typename U> PointType<T> operator+(const PointType<U>& rhs) const { PointType<T> res(*this); return res+=rhs; }
            template <typename U> PointType<T> operator-(const PointType<U>& rhs) const { PointType<T> res(*this); return res-=rhs; }
            template <typename U> PointType<T>& operator*=(const PointType<U>& rhs) { x *= rhs.x; y *= rhs.y; return *this; }
            template <typename U> PointType<T>& operator/=(const PointType<U>& rhs) { x /= rhs.x; y /= rhs.y; return *this; }
            template <typename U> PointType<T> operator*(const PointType<U>& rhs) const { PointType<T> res(*this); return res*=rhs; }
            template <typename U> PointType<T> operator/(const PointType<U>& rhs) const { PointType<T> res(*this); return res/=rhs; }
            PointType<T>& operator=(const PointType<T>& rhs) { x = rhs.x; y = rhs.y; return *this; }
            template <typename U> PointType<T>& operator=(const PointType<U>& rhs) { x = rhs.x; y = rhs.y; return *this; }
            PointType<T>& operator=(const T& rhs) { x = rhs; y = rhs; return *this; }
            template <typename U> PointType<T>& operator+=(const U& rhs) { x += rhs; y += rhs; return *this; }
            template <typename U> PointType<T>& operator-=(const U& rhs) { x -= rhs; y -= rhs; return *this; }
            template <typename U> PointType<T>& operator*=(const U& rhs) { x *= rhs; y *= rhs; return *this; }
            template <typename U> PointType<T>& operator/=(const U& rhs) { x /= rhs; y /= rhs; return *this; }
            template <typename U> PointType<T> operator+(const U& rhs) const { PointType<T> tmp(*this); tmp += rhs; return tmp; }
            template <typename U> PointType<T> operator-(const U& rhs) const { PointType<T> tmp(*this); tmp -= rhs; return tmp; }
            template <typename U> PointType<T> operator*(const U& rhs) const { PointType<T> tmp(*this); tmp *= rhs; return tmp; }
            template <typename U> PointType<T> operator/(const U& rhs) const { PointType<T> tmp(*this); tmp /= rhs; return tmp; }
            PointType<T> operator-(void) const { PointType<T> tmp(-y,-x); return tmp; }
            PointType<T> operator*(const T& rhs) const { PointType<T> tmp(*this); tmp *= rhs; return tmp; }
            PointType<T> operator/(const T& rhs) const { PointType<T> tmp(*this); tmp /= rhs; return tmp; }
            template <typename U> bool operator==(const PointType<U>& rhs) const { return (x == rhs.x && y == rhs.y); }
            bool operator==(T rhs) const { return (x == rhs && y == rhs); }
            template <typename U> bool operator!=(const PointType<U>& rhs) const { return !(*this == rhs); }
            bool operator!=(T rhs) const { return !(*this == rhs); }
            template <typename U> bool operator<(const PointType<U>& rhs) const { if (y==rhs.y) return (x < rhs.x); return (y < rhs.y); }
            template <typename U> bool operator>(const PointType<U>& rhs) const { if (y==rhs.y) return (x > rhs.x); return (y > rhs.y); }
            template <typename U> bool operator<(const U& rhs) const { if (y<rhs) return (x < rhs); return false; }
            template <typename U> bool operator>(const U& rhs) const { if (y>rhs) return (x > rhs); return false; }
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
        
        template<typename T>
        PointType<T> operator+( T a, const PointType<T>& pt ) {
            return PointType<T>( a+pt.y, a+pt.x );
        }
        template<typename T>
        PointType<T> operator-( T a, const PointType<T>& pt ) {
            return PointType<T>( a-pt.y, a-pt.x );
        }
        template<typename T>
        PointType<T> operator*( T a, const PointType<T>& pt ) {
            return PointType<T>( a*pt.y, a*pt.x );
        }
        template<typename T>
        PointType<T> operator/( T a, const PointType<T>& pt ) {
            return PointType<T>( a/pt.y, a/pt.x );
        }
        
        template <typename T>
        PointD pointWarp( std::vector<T> P, std::vector<T> Q, const PointD& point ) {
            if( P.size() != 16 || Q.size() != 16 ) return point;
            PointD ret(0,0);
            for( size_t j(0); j < 4; ++j ) {
                for( size_t k(0); k < 4; ++k ) {
                    ret.x += P.at(4*j+k) * pow(point.x,j) * pow(point.y,k);
                    ret.y += Q.at(4*j+k) * pow(point.x,j) * pow(point.y,k);
                }
            }
            return ret;
        }

        /*! @} */

    }

}

#endif // REDUX_UTIL_POINT_HPP
