#ifndef REDUX_IMAGE_IMAGE_HPP
#define REDUX_IMAGE_IMAGE_HPP

#include "redux/file/fileinfo.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/convert.hpp"
#include "redux/util/datautil.hpp"

#include <cstring>
#include <memory>


namespace redux {

    namespace image {

        template <typename T>
        class Image : public redux::util::Array<T> {

        public:
            typedef typename std::shared_ptr<Image> Ptr;

            Image( void ) : weight( 1 ) {};

            template <typename ...S>
            Image( S ...s ) : redux::util::Array<T>(s... ), weight( 1 ) {}
            
            template <typename ...S>
            Image( Image<T>& rhs, S ...s ) : redux::util::Array<T>(reinterpret_cast<redux::util::Array<T>&>(rhs), s... ), weight( rhs.weight ) {}

            // operators with scalars
            const Image<T>& operator= ( const T& rhs ) { for( auto & it : *this ) it  = rhs; return *this; };
            const Image<T>& operator+=( const T& rhs ) { for( auto & it : *this ) it += rhs; return *this; };
            const Image<T>& operator-=( const T& rhs ) { for( auto & it : *this ) it -= rhs; return *this; };
            const Image<T>& operator*=( const T& rhs ) { for( auto & it : *this ) it *= rhs; return *this; };
            const Image<T>& operator/=( const T& rhs ) { for( auto & it : *this ) it /= rhs; return *this; };
            
            void normalize(void) { *this /= weight; weight = 1; };
            double mean(void) const {
                double sum(0);
                size_t cnt=0;
                for( auto & it : *this ) {
                    sum += it;
                    cnt++;
                }
                if(cnt) {
                    sum /= static_cast<double>(cnt);

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

            template <typename U=T>
            Image<U> copy( void ) const {
                std::vector<size_t> newDimSizes;
                for (auto& it: this->dimensions()) {
                    if(it>1) {
                        newDimSizes.push_back(it);
                    }
                }
                Image<U> tmp( newDimSizes );
                auto cit = this->begin();
                for( auto& it : tmp ) {
                    it = static_cast<U>( *cit );
                    ++cit;
                }
                return tmp;
            }
            
            template <typename U>
            void copy( Image<U> img ) const {
                std::vector<size_t> newDimSizes;
                for (auto& it: this->dimensions()) {
                    if(it>1) {
                        newDimSizes.push_back(it);
                    }
                }
                img.resize( newDimSizes );
                auto cit = this->begin();
                for( auto& it : img ) {
                    it = static_cast<U>( *cit );
                    ++cit;
                }

            } 

            template <typename U>
            const Image<T>& operator=( const Image<U>& rhs ) {
                redux::util::Array<T>::operator=(rhs);
                weight = rhs.weight; 
                return *this;
            }

            template <typename U>
            const Image<T>& operator+=( const Image<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    typename Image<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) it += redux::util::bound_cast<T>( *rhsit++ );
                }
                else {
                    throw std::invalid_argument( "image dimensions does not match." );
                }
                return *this;
            }

            template <typename U>
            const Image<T>& operator-=( const Image<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    typename Image<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) it -= redux::util::bound_cast<T>( *rhsit++ );
                }
                else {
                    throw std::invalid_argument( "image dimensions does not match." );
                }
                return *this;
            }

            template <typename U>
            const Image<T>& operator*=( const Image<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    typename Image<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) it *= redux::util::bound_cast<T>( *rhsit++ );
                }
                else {
                    throw std::invalid_argument( "image dimensions does not match." );
                }
                return *this;
            }

            template <typename U>
            const Image<T>& operator/=( const Image<U>& rhs ) {
                if( this->sameSize( rhs ) ) {
                    typename Image<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) it /= redux::util::bound_cast<T>( *rhsit++ );
                }
                else {
                    throw std::invalid_argument( "image dimensions does not match." );
                }
                return *this;
            }


            template <typename U>
            bool operator==( const Image<U>& rhs ) const {
                if( this->sameSize( rhs ) ) {
                    typename Image<U>::const_iterator rhsit = rhs.begin();
                    for( auto & it : *this ) if( it != *rhsit++ ) return false;
                    return true;
                }
                return false;
            }

            std::shared_ptr<redux::file::FileInfo> hdr;

        protected:

            double weight;

        private:
            template<typename U> friend class Image;
        };


    }   // image

}   // redux


#endif  // REDUX_IMAGE_IMAGE_HPP
