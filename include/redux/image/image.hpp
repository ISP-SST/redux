#ifndef REDUX_IMAGE_IMAGE_HPP
#define REDUX_IMAGE_IMAGE_HPP

#include "redux/file/fileinfo.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/convert.hpp"
#include "redux/util/datautil.hpp"

#include <cstring>
#include <memory>

#include "redux/util/stringutil.hpp"
using namespace redux::util;
namespace redux {

    namespace image {

        template <typename T>
        class Image : public redux::util::Array<T> {
        public:
            typedef typename std::shared_ptr<Image> Ptr;

            Image( void ) : nFrames( 1 ) {};
            template <typename U>
            Image( const Image<T>& rhs, const std::vector<U>& indices ) : redux::util::Array<T>( reinterpret_cast<const redux::util::Array<T>&>( rhs ), indices ), nFrames( rhs.nFrames ) {}
            template <typename ...S>
            Image( const Image<T>& rhs, S ...s ) : Image<T>( rhs, std::vector<int64_t>( {static_cast<int64_t>( s )...} ) ) {}
            template <typename ...S>
            Image( S ...s ) : redux::util::Array<T>( s... ), nFrames( 1 ) {}


            template <typename U>
            const Image<T>& operator=( const U& rhs ) {
                redux::util::Array<T>::operator=( rhs );
                return *this;
            }

            template <typename U>
            const Image<T>& operator+=( const Image<U> rhs ) {
                redux::util::Array<T>::operator+=( rhs );
                nFrames += rhs.nFrames;
                return *this;
            }

            template <typename U>
            const Image<T>& operator+=( const U& rhs ) {
                redux::util::Array<T>::operator+=( rhs );
                return *this;
            }

            template <typename U>
            const Image<T>& operator-=( const U& rhs ) {
                redux::util::Array<T>::operator-=( rhs );
                return *this;
            }

            template <typename U>
            const Image<T>& operator*=( const U& rhs ) {
                redux::util::Array<T>::operator*=( rhs );
                return *this;
            }

            template <typename U>
            const Image<T>& operator/=( const U& rhs ) {
                redux::util::Array<T>::operator/=( rhs );
                return *this;
            }

            // operators with scalars
            /*  template <typename U> const Image<T>& operator= ( const U& rhs ) { for( auto & it : *this ) it  = rhs; return *this; }
              template <typename U> const Image<T>& operator+=( const U& rhs ) { for( auto & it : *this ) it += rhs; return *this; }
              template <typename U> const Image<T>& operator-=( const U& rhs ) { for( auto & it : *this ) it -= rhs; return *this; }
              template <typename U> const Image<T>& operator*=( const U& rhs ) { for( auto & it : *this ) it *= rhs; return *this; }
              template <typename U> const Image<T>& operator/=( const U& rhs ) { for( auto & it : *this ) it /= rhs; return *this; }
              */

            void normalize( void ) { if(nFrames) *this /= static_cast<double>(nFrames); nFrames = 1; };
            double mean( void ) const {
                double sum( 0 );
                size_t cnt = 0;
                for( auto & it : *this ) {
                    sum += it;
                    cnt++;
                }
                if( cnt ) {
                    sum /= static_cast<double>( cnt );

                }
                return sum;
            }
//             double mean(const Image<uint8_t>& mask) const {
//                 double sum(0);
//                 size_t cnt=0;
//                 Image<uint8_t>::const_iterator mit = mask.begin();
//                 for( auto & it : *this ) {
//                     if( mit++ ) {
//                         sum += it;
//                         cnt++;
//                     }
//                 }
//                 if(cnt) sum /= static_cast<double>(cnt);
//                 return sum;
//             }
//             template <typename ...S>
//             double mean(S ...s) const {
//                 Image<T> tmp(*this,s...);
//                 return tmp.mean();
//             }

            template <typename U = T>
            Image<U> copy( void ) const {
                std::vector<size_t> newDimSizes;
                for( auto & it : this->dimensions() ) {
                    if( it > 1 ) {
                        newDimSizes.push_back( it );
                    }
                }
                Image<U> tmp( newDimSizes );
                auto cit = this->begin();
                for( auto & it : tmp ) {
                    it = static_cast<U>( *cit );
                    ++cit;
                }
                return tmp;
            }


            template <typename U>
            bool operator==( const redux::util::Array<U>& rhs ) const {
                if( this->sameSize( rhs ) ) {
                    typename Image<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) if( it != *rhsit++ ) return false;
                    return true;
                }
                return false;
            }

            std::shared_ptr<redux::file::FileInfo> hdr;

            uint32_t nFrames;      //! Set this if the image consists of multiple frames (for normalization purposes)

        private:
            template<typename U> friend class Image;
        };


    }   // image

}   // redux


#endif  // REDUX_IMAGE_IMAGE_HPP
