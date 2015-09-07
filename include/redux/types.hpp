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
        template<typename T> struct PointType {
            PointType ( T yy=0, T xx=0 ) : x(xx), y(yy) {}
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
            template <typename U> PointType<T> max(const PointType<U>& rhs) { PointType<T> tmp(*this); tmp.x=std::max<T>(tmp.x,rhs.x); tmp.y=std::max<T>(tmp.y,rhs.y); return tmp; }
            template <typename U> PointType<T> min(const PointType<U>& rhs) { PointType<T> tmp(*this); tmp.x=std::min<T>(tmp.x,rhs.x); tmp.y=std::min<T>(tmp.y,rhs.y); return tmp; }
            operator std::string() const { return "("+std::to_string(y)+","+std::to_string(x)+")"; }
            template <typename U> PointType<T>& operator+=(const PointType<U>& rhs) { x += rhs.x; y += rhs.y; return *this; }
            template <typename U> PointType<T>& operator-=(const PointType<U>& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
            template <typename U> PointType<T> operator+(const PointType<U>& rhs) const { PointType<T> res(*this); return std::move(res+=rhs); }
            template <typename U> PointType<T> operator-(const PointType<U>& rhs) const { PointType<T> res(*this); return std::move(res-=rhs); }
            PointType<T>& operator=(const PointType<T>& rhs) { x = rhs.x; y = rhs.y; return *this; }
            template <typename U> PointType<T>& operator=(const PointType<U>& rhs) { x = rhs.x; y = rhs.x; return *this; }
            template <typename U> PointType<T>& operator=(const U& rhs) { x = rhs; y = rhs; return *this; }
            template <typename U> PointType<T>& operator+=(const U& rhs) { x += rhs; y += rhs; return *this; }
            template <typename U> PointType<T>& operator-=(const U& rhs) { x -= rhs; y -= rhs; return *this; }
            template <typename U> PointType<T> operator+(const U& rhs) const { PointType<T> tmp(*this); tmp += rhs; return tmp; }
            template <typename U> PointType<T> operator-(const U& rhs) const { PointType<T> tmp(*this); tmp -= rhs; return tmp; }
            template <typename U> bool operator==(const PointType<U>& rhs) const { return (x == rhs.x && y == rhs.y); }
            template <typename U> bool operator!=(const PointType<U>& rhs) const { return !(*this == rhs); }
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
        
        /*!
        * A general 2D rectangle
        */
        template<class T> struct RegionType {
            RegionType ( T firsty, T firstx, T lasty, T lastx ) : first(firsty,firstx), last(lasty,lastx) {}
            RegionType ( T lasty, T lastx ) : first(0,0), last(lasty,lastx) {}
            RegionType ( void ) : first(0,0), last(0,0) {}
            RegionType ( RegionType<T>&& rhs ) : first(std::move(rhs.first)), last(std::move(rhs.last)) {}
            RegionType ( const RegionType<T>& rhs ) : first(rhs.first), last(rhs.last) {}
            template <typename U> RegionType ( const RegionType<U>& rhs ) : first(rhs.first), last(rhs.last) {}
            static uint64_t size(void) { return 4*sizeof(T); };
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
            template <typename U> void shiftX(const U& x) { first.x += x; last.x += x; }
            template <typename U> void shiftY(const U& y) { first.y += y; last.y += y; }
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
            template <typename U> void grow(const U& rhs) { PointType<U> tmp(rhs,rhs); first-=tmp; last+=tmp; }
            template <typename U> RegionType<T> grown(const PointType<U>& rhs) const { RegionType<T> res(*this); res.grow(rhs); return std::move(res); }
            template <typename U> RegionType<T> grown(const U& rhs) const { RegionType<T> res(*this); res.grow(PointType<U>(rhs,rhs)); return std::move(res); }
            template <typename U> void restrict(const RegionType<U>& rhs) { first=first.max(rhs.first); last=last.min(rhs.last); }
            operator std::string() const { return (std::string)first + " -> " + (std::string)last; }
            template <typename U> RegionType<T>& operator+=(const PointType<U>& rhs) { first += rhs; last += rhs; return *this; }
            template <typename U> RegionType<T>& operator-=(const PointType<U>& rhs) { first -= rhs; last -= rhs; return *this; }
            template <typename U> RegionType<T> operator+(const PointType<U>& rhs) const { RegionType<T> res(*this); return std::move(res+=rhs); }
            template <typename U> RegionType<T> operator-(const PointType<U>& rhs) const { RegionType<T> res(*this); return std::move(res-=rhs); }
            template <typename U> RegionType<T>& operator=(const RegionType<U>& rhs) { first = rhs.first; last = rhs.last; return *this; }
            template <typename U> RegionType<T>& operator=(const U& rhs) { first = rhs; last = rhs; return *this; }
            template <typename U> bool operator==(const RegionType<U>& rhs) const { return (first == rhs.first && last == rhs.last); }
            template <typename U> bool operator!=(const RegionType<U>& rhs) const { return !(*this == rhs); }
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

        
        typedef std::complex<double> complex_t;
        
        /* @} */

}


#endif // REDUX_TYPES_HPP


