#ifndef REDUX_UTIL_BOUNDVALUE_HPP
#define REDUX_UTIL_BOUNDVALUE_HPP

#include "redux/exception.hpp"
#include "redux/util/convert.hpp"

#include <iostream>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */

        namespace detail {

            template <typename T>
            static void truncate( T& v, const T& min, const T& max ) {
                if( v < min ) {
                    v = min;
                }
                else if( v > max ) {
                    v = max;
                }
            }

            template <typename T>
            static void wrap( T& v, const T& min, const T& max ) {
                T span( max - min );
                int cnt = static_cast<int>( ( v - min ) / span );
                if( v < min ) cnt--;
                v = v - cnt * span;
            }

            template <typename T>
            static void reflect( T& v, const T& min, const T& max ) {
                T span( max - min );
                int cnt = static_cast<int>( ( v - min ) / span );
                if( v < min ) cnt--;
                bool reverse = ( cnt & 1 );
                if( reverse ) {
                    v = max - v + cnt * span + min;
                }
                else {
                    v -= cnt * span;
                }
            }

            enum TrimType { UNDEFINEDTRIM = 0, TRUNCATE, WRAP, REFLECT };

        }

        /**
         *   @brief     Wrapper around an elementary type, which enforces max/min values.
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         */
        template<class T, detail::TrimType TT = detail::TRUNCATE>
        class BoundValue {
            typedef void ( *trimFunction )( T&, const T&, const T& );

        public:
            BoundValue( T val = 0, T minVal = std::numeric_limits<T>::min(), T maxVal = std::numeric_limits<T>::max() ) : trim_( selectTrimFunction( TT ) ) {
                setLimits( minVal, maxVal );
                assignWithTruncation( val );
            }
            BoundValue( const BoundValue& rhs ) : trim_( selectTrimFunction( TT ) ), val_( rhs.val_ ), minVal_( rhs.minVal_ ), maxVal_( rhs.maxVal_ ) { }
            template<class U, detail::TrimType UU>
            BoundValue( const BoundValue<U, UU>& rhs ) : trim_( selectTrimFunction( UU ) ) {
                T tmpV;
                try {
                    minVal_ = boost::numeric_cast<T>( rhs.minVal_ );
                    maxVal_ = boost::numeric_cast<T>( rhs.maxVal_ );
                    tmpV    = boost::numeric_cast<T>( rhs.val_ );
                }
                catch( const boost::numeric::bad_numeric_cast& ) {
                    minVal_ = std::numeric_limits<T>::min();
                    maxVal_ = std::numeric_limits<T>::max();
                    try {
                        tmpV    = boost::numeric_cast<T>( rhs.val_ );
                    }
                    catch( const boost::numeric::bad_numeric_cast& ) {
                        tmpV = minVal_;
                    }
                }
                assignWithTruncation( tmpV );
            }


            void setLimits( T minVal, T maxVal ) {
                if( minVal > maxVal ) {
                    std::swap( minVal, maxVal );
                }
                minVal_ = minVal;
                maxVal_ = maxVal;
            }

            const T& min( void ) { return minVal_; };
            void setMin( T minVal ) {
                minVal_ = minVal;
                assignWithTruncation( val_ );
            }

            const T& max( void ) { return maxVal_; };
            void setMax( T maxVal ) {
                maxVal_ = maxVal;
                assignWithTruncation( val_ );
            }

            BoundValue& operator= ( const BoundValue& rhs ) {
                trim_ = rhs.trim_;
                val_ = rhs.val_;
                minVal_ = rhs.minVal_;
                maxVal_ = rhs.maxVal_;
                return *this;
            }

            BoundValue& operator= ( T val ) {
                assignWithTruncation( val );
                return *this;
            }

            BoundValue& operator+= ( T val ) {
                assignWithTruncation( val_ + val );
                return *this;
            }

            BoundValue& operator-= ( T val ) {
                assignWithTruncation( val_ - val );
                return *this;
            }

            BoundValue& operator*= ( T val ) {
                assignWithTruncation( val_ * val );
                return *this;
            }

            BoundValue& operator/= ( T val ) {
                assignWithTruncation( val_ / val );
                return *this;
            }

            T& value( void ) { return val_; }

            operator T() {
                return val_;
            }

            template <typename U>
            operator U() {
                return static_cast<U>( val_ );
            }

            bool operator== ( T val ) const {
                return val == val_;
            }

            // Prefix
            BoundValue& operator++() {
                assignWithTruncation( static_cast<T>( val_ + 1 ) );
                return *this;
            }
            BoundValue& operator--() {
                assignWithTruncation( static_cast<T>( val_ - 1 ) );
                return *this;
            }

            // Postfix
            BoundValue operator++ ( int ) {
                BoundValue old( *this );
                assignWithTruncation( static_cast<T>( val_ + 1 ) );
                return old;
            }
            BoundValue operator-- ( int ) {
                BoundValue old( *this );
                assignWithTruncation( static_cast<T>( val_ - 1 ) );
                return old;
            }

        protected:
            trimFunction selectTrimFunction( detail::TrimType tt ) {
                switch( tt ) {
                    case detail::TRUNCATE: return &( detail::truncate<T> );
                    case detail::WRAP: return &( detail::wrap<T> );
                    case detail::REFLECT: return &( detail::reflect<T> );
                    default: throw redux::NotImplemented( "Function not defined for this TrimType" );
                }
            }
            inline T assignWithTruncation( T val ) {
                val_ = val;
                trim_( val_, minVal_, maxVal_ );

                return val_;
            }

            template <class U, detail::TrimType UU> friend std::ostream& operator<< ( std::ostream& os, const BoundValue<U, UU>& bv );
            template <class U, detail::TrimType UU> friend std::istream& operator>> ( std::istream& in, BoundValue<U, UU>& bv );
            template<class U, detail::TrimType UU> friend class BoundValue;

        private:
            trimFunction trim_;
            T val_;
            T minVal_;
            T maxVal_;
        };

        template<class T, detail::TrimType TT>
        std::ostream& operator<< ( std::ostream& os, const BoundValue<T, TT> & bv ) {
            return ( os << bv.val_ );
        }

        template<class T, detail::TrimType TT>
        std::istream& operator>>( std::istream& in, BoundValue<T, TT>& bv ) {
            T tmp;
            in >> tmp;
            bv.assignWithTruncation( tmp );
            return in;
        }

        /*! @} */

    }

}

#endif  // REDUX_UTIL_BOUNDVALUE_HPP
