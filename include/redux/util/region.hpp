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
            explicit RegionType ( PointType<T> f, PointType<T> l ) : first(f), last(l) {}
            RegionType ( T lasty, T lastx ) : first(0,0), last(lasty,lastx) {}
            explicit RegionType ( PointType<T> l ) : first(0,0), last(l) {}
            RegionType ( void ) : first(0,0), last(0,0) {}
            RegionType ( RegionType<T>&& rhs ) : first(std::move(rhs.first)), last(std::move(rhs.last)) {}
            RegionType ( const RegionType<T>& rhs ) : first(rhs.first), last(rhs.last) {}
            template <typename U> explicit RegionType ( const RegionType<U>& rhs ) : first(rhs.first), last(rhs.last) {}
            static inline uint64_t size(void) { return 2*PointType<T>::size(); };
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
            void normalize( void ) {
                if( first.y > last.y ) std::swap( first.y, last.y );
                if( first.x > last.x ) std::swap( first.x, last.x );
            }
            template <typename U> void shift( const PointType<U>& s ) { first += s; last += s; }
            void shiftX(T x) { first.x += x; last.x += x; }
            void shiftY(T y) { first.y += y; last.y += y; }
            template <typename U> void shift(U y, U x) { shift( PointType<U>(y,x) ); }
            template <typename U> void shift(U s) { first += s; last += s; }
            template <typename U> bool isInside( const PointType<U> rhs ) const {
                if( (rhs.x > (U)last.x) || (rhs.x < (U)first.x) ||
                    (rhs.y > (U)last.y) || (rhs.y < (U)first.y) ) return false;
                return true;
            }
            inline bool isNormal( void ) const {
                return !( (first.x > last.x) || (first.y > last.y) );
            }
            template <typename U> PointI outside(const RegionType<U>& rhs) const {
                PointI res(0,0);
                if(rhs.first.x<first.x) res.x = rhs.first.x-first.x;
                else if(rhs.last.x>last.x) res.x = rhs.last.x-last.x;
                if(rhs.first.y<first.y) res.y = rhs.first.y-first.y;
                else if(rhs.last.y>last.y) res.y = rhs.last.y-last.y;
                return res;
            }
            template <typename U> RegionType<T> intersection( RegionType<U> rhs ) const {
                RegionType<T> tmp(*this);
                tmp.normalize();
                rhs.normalize();
                if( tmp.isInside(rhs.first) ) {
                    if( tmp.isInside(rhs.last) ) return rhs;
                    tmp.first = rhs.first;
                    return tmp;
                }
                if( rhs.isInside(tmp.first) ) {
                    if( rhs.isInside(tmp.last) ) return tmp;
                    tmp.last = rhs.last;
                } else {
                    tmp = 0;
                }
                return tmp;
            }
            template <typename U> RegionType<T> enclosure( RegionType<U> rhs ) const {
                RegionType<T> tmp(*this);
                tmp.normalize();
                rhs.normalize();
                tmp.first = tmp.first.min(rhs.first);
                tmp.last = tmp.last.max(rhs.last);
                return tmp;
            }
            template <typename U> PointType<U> diff( const PointType<U>& rhs ) {
                PointType<U> ret;
                if( isInside(rhs) ) return ret;
                normalize();
                if( rhs.x < (U)first.x ) ret.x = rhs.x-first.x;
                else if( rhs.x > (U)last.x ) ret.x = rhs.x-last.x;
                if( rhs.y < (U)first.y ) ret.y = rhs.y-first.y;
                else if( rhs.y > (U)last.y ) ret.y = rhs.y-last.y;
                return ret;
            }
            template <typename U> void expand( const PointType<U>& rhs ) { normalize(); first-=rhs; last+=rhs; }
            template <typename U> void expand( U rhs ) { expand( PointType<U>(rhs,rhs) ); }
            template <typename U> void shrink( const PointType<U>& rhs ) { normalize(); first+=rhs; last-=rhs; }
            template <typename U> void shrink( U rhs ) { shrink( PointType<U>(rhs,rhs) ); }
            void shrinkSigned( const PointI& rhs ) {
                if( rhs.x > 0 ) last.x -= rhs.x;
                else first.x -= rhs.x;
                if( rhs.y > 0 ) last.y -= rhs.y;
                else first.y -= rhs.y;
            }
            template <typename U> void include( const PointType<U>& rhs ) {
                normalize();
                first = first.min(rhs);
                last = last.max(rhs);
            }
            template <typename U> void include( U y, U x ) { include( PointType<T>(y,x) ); }
            template <typename U> void grow(const RegionType<U>& rhs) { first=first.min(rhs.first); last=last.max(rhs.last); }
            template <typename U> void grow(const PointType<U>& rhs) { first=first.min(rhs); last=last.max(rhs); }
            //void grow(T rhs) { grow( PointType<T>(rhs,rhs) ); }
            template <typename U> void grow(U rhs) { grow( PointType<U>(rhs,rhs) ); }
            void shrink(T rhs) { shrink( PointType<T>(rhs,rhs) ); }
            template <typename U> PointType<U> getAsGlobal( const PointType<U>& rhs ) {
                PointType<U> tmp(rhs);
                if( isNormal() ) return tmp+first;
                if( first.x > last.x ) tmp.x = first.x - tmp.x;
                else tmp.x += first.x;
                if( first.y > last.y ) tmp.y = first.y - tmp.y;
                else tmp.y += first.y;
                return tmp; }
            template <typename U> RegionType<U> getAsGlobal( const RegionType<U>& rhs ) {
                RegionType<U> tmp(rhs);
                tmp.first = getAsGlobal(tmp.first);
                tmp.last = getAsGlobal(tmp.last);
                return tmp; }
            PointType<T> mid( void ) const { return (first+last+1)/2; }
            PointD com( void ) const { return PointD(first+last)/2; }
            Point span( void ) const { return Point(std::abs(last.y-first.y)+1, std::abs(last.x-first.x)+1); }
            template <typename U> RegionType<T> grown(const PointType<U>& rhs) const { RegionType<T> res(*this); res.grow(rhs); return res; }
            RegionType<T> grown(T rhs) const { RegionType<T> res(*this); res.grow(PointType<T>(rhs,rhs)); return res; }
            template <typename U> void restrict(const RegionType<U>& rhs) { first=first.max(rhs.first); last=last.min(rhs.last); }
            template <typename U> PointType<U> restrict( const PointType<U>& rhs ) {
                if( isInside(rhs) ) return rhs;
                PointType<U> ret(rhs);
                ret.min(last);
                ret.max(first);
                return ret;
            }
            operator std::string() const {
                std::string ret = "[" + std::to_string(first.y) + ":" + std::to_string(last.y) +",";
                ret += std::to_string(first.x) + ":" + std::to_string(last.x) + "]";
                return ret;
            }
            template <typename U> RegionType<T>& operator+=(const PointType<U>& rhs) { first += rhs; last += rhs; return *this; }
            template <typename U> RegionType<T>& operator-=(const PointType<U>& rhs) { first -= rhs; last -= rhs; return *this; }
            template <typename U> RegionType<T> operator+(const PointType<U>& rhs) const { RegionType<T> res(*this); return res+=rhs; }
            template <typename U> RegionType<T> operator+(const U& rhs) const { RegionType<T> res(*this); res.first+=rhs; res.last+=rhs; return res; }
            template <typename U> RegionType<T> operator-(const PointType<U>& rhs) const { RegionType<T> res(*this); return res-=rhs; }
            RegionType<T> operator-(void) const { RegionType<T> tmp(first.y,first.x,last.y,last.x); return tmp; }
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
