#ifndef REDUX_IMAGE_IMAGE_HPP
#define REDUX_IMAGE_IMAGE_HPP

#include "redux/util/array.hpp"
#include "redux/file/filemeta.hpp"

#include <cstring>
#include <memory>

namespace redux {

    namespace image {

        template <typename T>
        class Image : public redux::util::Array<T> {
        public:
            typedef typename std::shared_ptr<Image> Ptr;

            Image( void ) : nFrames( 1 ) {};
            Image(Image<T>&& rhs) : redux::util::Array<T>(std::move(reinterpret_cast<redux::util::Array<T>&>(rhs))), meta(std::move(rhs.meta)), nFrames(std::move(rhs.nFrames)) { };
            Image(const Image<T>& rhs) : redux::util::Array<T>(reinterpret_cast<const redux::util::Array<T>&>(rhs)), meta(rhs.meta) { };
            template <typename U>
            Image( const Image<T>& rhs, const std::vector<U>& indices ) : redux::util::Array<T>( reinterpret_cast<const redux::util::Array<T>&>( rhs ), indices ), nFrames( rhs.nFrames ) {}
            template <typename ...S>
            Image( const Image<T>& rhs, S ...s ) : Image<T>( rhs, std::vector<int64_t>( {static_cast<int64_t>( s )...} ) ) {}
            template <typename ...S>
            Image( S ...s ) : redux::util::Array<T>( s... ), nFrames( 1 ) {}

            operator const redux::util::Array<T>&() const { return reinterpret_cast<const redux::util::Array<T>&>(*this); }
            
            Image<T>& operator=( const Image<T>& rhs ) {
                redux::util::Array<T>::operator=( reinterpret_cast<const redux::util::Array<T>&>(rhs) );
                meta = rhs.meta;
                nFrames = rhs.nFrames;
                return *this;
            }
            
            template <typename U>
            Image<T>& operator=( const Image<U>& rhs ) {
                redux::util::Array<T>::operator=( reinterpret_cast<const redux::util::Array<U>&>(rhs) );
                meta = rhs.meta;
                nFrames = rhs.nFrames;
                return *this;
            }
            /*template <typename U>
            const Image<T>& operator=( const redux::util::Array<U>& rhs ) {
                redux::util::Array<T>::operator=( rhs );
                return *this;
            }*/

            template <typename U>
            const Image<T>& operator+=( const Image<U>& rhs ) {
                redux::util::Array<T>::operator+=( reinterpret_cast<const redux::util::Array<U>&>(rhs) );
                nFrames += rhs.nFrames;         // TBD: should frames be counted, or assumed to always be normalized before arithmetic ?
                return *this;
            }

            template <typename U>
            const Image<T>& operator-=( const Image<U>& rhs ) {
                redux::util::Array<T>::operator-=( reinterpret_cast<const redux::util::Array<U>&>(rhs) );
                return *this;
            }
            
            template <typename U>
            const Image<T>& operator*=( const Image<U>& rhs ) {
                redux::util::Array<T>::operator*=( reinterpret_cast<const redux::util::Array<U>&>(rhs) );
                return *this;
            }
            
            template <typename U>
            const Image<T>& operator/=( const Image<U>& rhs ) {
                redux::util::Array<T>::operator/=( reinterpret_cast<const redux::util::Array<U>&>(rhs) );
                return *this;
            }

            template <typename U> const Image<T>& operator+=( const U& rhs ) { redux::util::Array<T>::operator+=(rhs); return *this; }
            template <typename U> const Image<T>& operator-=( const U& rhs ) { redux::util::Array<T>::operator-=(rhs); return *this; }
            template <typename U> const Image<T>& operator*=( const U& rhs ) { redux::util::Array<T>::operator*=(rhs); return *this; }
            template <typename U> const Image<T>& operator/=( const U& rhs ) { redux::util::Array<T>::operator/=(rhs); return *this; }


            void normalize( void ) { if(nFrames > 1) redux::util::Array<T>::operator*=(1.0/nFrames); nFrames = 1; };
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

            template <typename U>
            void copy( redux::util::Array<U>& rhs ) const {
                redux::util::Array<T>::copy(rhs);
            }


            template <typename U = T>
            Image<U> copy( bool skipTrivialDims=true ) const {
                Image<U> tmp = redux::util::Array<T>::template copy<U>(skipTrivialDims);
                return std::move(tmp);
//                 std::vector<size_t> newDimSizes = this->dimensions(skipTrivialDims);
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

            std::shared_ptr<redux::file::FileMeta> meta;

            uint32_t nFrames;      //! Set this if the image consists of multiple frames (for normalization purposes)

        private:
            template<typename U> friend class Image;
        };


    }   // image

}   // redux


#endif  // REDUX_IMAGE_IMAGE_HPP
