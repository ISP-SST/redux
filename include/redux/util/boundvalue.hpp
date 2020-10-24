#ifndef REDUX_UTIL_BOUNDVALUE_HPP
#define REDUX_UTIL_BOUNDVALUE_HPP

#include "redux/util/convert.hpp"
#include "redux/util/datautil.hpp"

#include <iostream>
#include <limits>
#include <stdexcept>

#include <boost/numeric/conversion/cast.hpp>

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */

        namespace detail {

            enum RestrictType { UNDEFINEDTRIM=0, TRUNCATE, WRAP, REFLECT };

        }

        /**
         *   @brief     Wrapper around an elementary type, which enforces max/min values.
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         */
        template<class T, detail::RestrictType TT = detail::TRUNCATE>
        class BoundValue {
            typedef T ( *trimFunction ) ( T, T, T );

        public:
            BoundValue( T val = 0, T minVal = std::numeric_limits<T>::lowest(),
                                   T maxVal = std::numeric_limits<T>::max() ) : trim_( setRestrictType(TT) ),
                                   val_(val), minVal_(minVal), maxVal_(maxVal) {
                if( minVal_ > maxVal_ ) std::swap( minVal_, maxVal_ );
                assignWithRestriction( val_ );
            }
            BoundValue( const BoundValue& rhs ) : trim_(rhs.trim_), val_(rhs.val_),
                                                  minVal_(rhs.minVal_), maxVal_(rhs.maxVal_) { }
            template<class U, detail::RestrictType UU>
            BoundValue( const BoundValue<U, UU>& rhs ) : trim_( setRestrictType(UU) ) {
                T tmpV;
                try {
                    minVal_ = boost::numeric_cast<T>( rhs.minVal_ );
                    maxVal_ = boost::numeric_cast<T>( rhs.maxVal_ );
                } catch( const boost::numeric::bad_numeric_cast& ) {
                    minVal_ = std::numeric_limits<T>::lowest();
                    maxVal_ = std::numeric_limits<T>::max();
                }
                try {
                    tmpV = boost::numeric_cast<T>( rhs.val_ );
                } catch ( const boost::numeric::bad_numeric_cast& ) {
                    tmpV = T(0);
                }
                assignWithRestriction( tmpV );
            }


            void setLimits( T minVal, T maxVal ) {
                if( minVal > maxVal ) {
                    std::swap( minVal, maxVal );
                }
                minVal_ = minVal;
                maxVal_ = maxVal;
                assignWithRestriction( val_ );
            }

            const T& min( void ) const { return minVal_; };
            inline void setMin( T minVal ) { setLimits( minVal, maxVal_ ); }
            const T& max( void ) const { return maxVal_; };
            inline void setMax( T maxVal ) { setLimits( minVal_, maxVal ); }

            inline T span(void) const { return (maxVal_ - minVal_); }
            
            T getRelative( float r ) {
                T s = span();
                T ret;
                //try {
                    ret = boost::numeric_cast<T>( minVal_ + s*r );
                //} catch ( const boost::numeric::bad_numeric_cast& ) {
                //    ret = minVal_;
                //}
                return ret;
            }

            template <class U, detail::RestrictType UU>
            BoundValue& operator=( const BoundValue<U,UU>& rhs ) {
                BoundValue<T,TT> tmp(rhs);
                *this = tmp;
                return *this;
            }

            BoundValue& operator=( T val ) {
                assignWithRestriction( val );
                return *this;
            }

            BoundValue& operator+=( T val ) {
                assignWithRestriction( val_ + val );
                return *this;
            }

            BoundValue& operator-=( T val ) {
                assignWithRestriction( val_ - val );
                return *this;
            }

            BoundValue& operator*=( T val ) {
                assignWithRestriction( val_ * val );
                return *this;
            }

            BoundValue& operator/=( T val ) {
                assignWithRestriction( val_ / val );
                return *this;
            }

            T& value( void ) {
                return val_;
            }

            operator T() {
                return val_;
            }

            template <typename U>
            operator U() {
                return static_cast<U> ( val_ );
            }

            bool operator== ( T val ) const {
                return val == val_;
            }

            // Prefix
            BoundValue& operator++() {
                assignWithRestriction ( static_cast<T> ( val_ + 1 ) );
                return *this;
            }
            BoundValue& operator--() {
                assignWithRestriction ( static_cast<T> ( val_ - 1 ) );
                return *this;
            }

            // Postfix
            BoundValue operator++ ( int ) {
                BoundValue old ( *this );
                assignWithRestriction ( static_cast<T> ( val_ + 1 ) );
                return old;
            }
            BoundValue operator-- ( int ) {
                BoundValue old ( *this );
                assignWithRestriction ( static_cast<T> ( val_ - 1 ) );
                return old;
            }

            trimFunction setRestrictType( detail::RestrictType tt ) {
                switch ( tt ) {
                    case detail::TRUNCATE:
                        trim_ = &(redux::util::restrict<T,T>); break;
                    case detail::WRAP:
                        trim_ = &(redux::util::restrict_periodic<T,T>); break;
                    case detail::REFLECT:
                        trim_ = &(redux::util::restrict_reflected<T,T>); break;
                    default:
                        throw std::invalid_argument( "BoundValue: Function not defined for this RestrictType: "+std::to_string(tt) );
                }
                return trim_;
            }
        protected:

            inline void assignWithRestriction( T val ) { val_ = trim_ ( val, minVal_, maxVal_ ); }

            template <class U, detail::RestrictType UU> friend std::ostream& operator<< ( std::ostream& os, const BoundValue<U, UU>& bv );
            template <class U, detail::RestrictType UU> friend std::istream& operator>> ( std::istream& in, BoundValue<U, UU>& bv );
            template<class U, detail::RestrictType UU> friend class BoundValue;

        private:
            trimFunction trim_;
            T val_;
            T minVal_;
            T maxVal_;
        };

        template<class T, detail::RestrictType TT>
        std::ostream& operator<< ( std::ostream& os, const BoundValue<T, TT> & bv ) {
            return ( os << bv.val_ );
        }

        template<class T, detail::RestrictType TT>
        std::istream& operator>> ( std::istream& in, BoundValue<T, TT>& bv ) {
            T tmp;
            in >> tmp;
            bv.assignWithRestriction ( tmp );
            return in;
        }

        /*! @} */

    }

}

#endif  // REDUX_UTIL_BOUNDVALUE_HPP
