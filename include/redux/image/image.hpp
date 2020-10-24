#ifndef REDUX_IMAGE_IMAGE_HPP
#define REDUX_IMAGE_IMAGE_HPP

#include "redux/util/array.hpp"
#include "redux/file/filemeta.hpp"

#include <memory>
#include <mutex>

namespace redux {

    namespace image {

        template <typename T>
        class Image : public redux::util::Array<T> {
        public:
            typedef typename std::shared_ptr<Image> Ptr;

            Image( void ) {};
            Image(Image<T>&& rhs) : redux::util::Array<T>(std::move(reinterpret_cast<redux::util::Array<T>&>(rhs))), meta(std::move(rhs.meta)) { };
            Image(const Image<T>& rhs) : redux::util::Array<T>(reinterpret_cast<const redux::util::Array<T>&>(rhs)), meta(rhs.meta) { };
            template <typename U>
            Image( const Image<T>& rhs, const std::vector<U>& indices ) : redux::util::Array<T>( reinterpret_cast<const redux::util::Array<T>&>( rhs ), indices ) {}
            template <typename ...S>
            Image( const Image<T>& rhs, S ...s ) : Image<T>( rhs, std::vector<int64_t>( {static_cast<int64_t>( s )...} ) ) {}
            template <typename ...S>
            explicit Image( S ...s ) : redux::util::Array<T>( s... ) {}
            
            void reverseX( void ) {
                if( this->nDimensions() != 2 ) return;
                size_t sy = this->dimSize(0);
                size_t sx = this->dimSize(1);
                std::shared_ptr<T*> arrayPtr = util::Array<T>::reshape(sy, sx);
                T** imgPtr = arrayPtr.get();
                redux::util::reverseX( imgPtr, sy, sx );
            }
            void reverseY( void ) {
                if( this->nDimensions() != 2 ) return;
                size_t sy = this->dimSize(0);
                size_t sx = this->dimSize(1);
                std::shared_ptr<T*> arrayPtr = util::Array<T>::reshape(sy, sx);
                T** imgPtr = arrayPtr.get();
                redux::util::reverseY( imgPtr, sy, sx );
            }

            Image<T>& operator=( const Image<T>& rhs ) {
                redux::util::Array<T>::operator=( reinterpret_cast<const redux::util::Array<T>&>(rhs) );
                meta = rhs.meta;
                return *this;
            }
            
            template <typename U>
            Image<T>& operator=( const U& rhs ) {
                redux::util::Array<T>::operator=(rhs);
                return *this;
            }

             template <typename U = T>
            Image<U> copy(void) const {
                Image<U> tmp;
                redux::util::Array<T>::copy(tmp);
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

            std::shared_ptr<redux::file::FileMeta> meta;
            std::mutex imgMutex;

        };

    }   // image
    
}   // redux


#endif  // REDUX_IMAGE_IMAGE_HPP
