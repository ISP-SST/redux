#ifndef REDUX_UTIL_ARRAYUTIL_HPP
#define REDUX_UTIL_ARRAYUTIL_HPP

#include <memory>
#include <iostream>

namespace redux {

    namespace util {

        /*!  @ingroup util
         *  @{
         */

        /*!  @file  arrayutil.hpp
         *   @brief     Set of templates for creating/deleting arrays of arbitrary types and dimensionality.
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         */

        /*! @name delArray
         *  @brief De-allocate pointers and data-block.
         *  @details Will recursively call itself for lower pointer orders and delete the memory block.
         *  @param p Pointer (of any order)
         */
        //@{
        template <class T> void delArray( T*& p ) { delete[] p; }
        template <class T>
        void delArray( T**& p ) {
            delArray( *p );              // Call delArray() for the next pointer-level.
            delete[] p;                  // Delete current pointer-level.
        }
        template <typename T, typename U> void castAndDelete( U* p ) { delete[] reinterpret_cast<T*>(p); }
        //@}

        /*! @name delPointers
         *  @brief De-allocate pointers (but @e NOT memory-block)
         *  @details Will recursively call itself for lower pointer orders.
         *  @param p Array of pointers (of any order)
         */
        //@{
        template <class T>  void delPointers( T* p ) { }
        template <class T>
        void delPointers( T**& p ) {
            delPointers( *p );           // Call delPointers() for the next pointer-level.
            delete[] p;                  // Delete current pointer-level.
        }
        //@}

        namespace detail {

            template <typename T, typename First, typename... Rest>
            struct Dummy {
                typedef typename Dummy<T, Rest...>::dataType* dataType;
                typedef typename Dummy<T, Rest...>::ptrType* ptrType;
            };

            template <typename T, typename First>
            struct Dummy<T, First> {
                typedef T dataType;
                typedef T* ptrType;
            };

        }


        /*! @name newArray
         *  @brief Allocate a multi-dimensional array.
         *  @details Allocate memory-block and generate a pointer structure for it so it can be accessed as a multi-dimensional array.
         *  @param n1 Size of dimension
         *  @param n2 Size of dimension
         *  @param rest Sizes of dimensions
         */
        //@{
        template <class T>
        T* newArray( size_t n1 ) {
            T* ret = new T[ n1 ];
            return ret;
        }
        template <class T, typename... Rest>
        typename detail::Dummy<T, size_t, size_t, Rest...>::ptrType newArray( size_t n1, size_t n2, Rest ...rest ) {
            typedef typename detail::Dummy<T, size_t, Rest...>::ptrType ptrType;
            ptrType* ret = new ptrType[ n1 ];               // Create an array of size n1 with pointers.
            *ret = newArray<T>( n1 * n2, rest... );         // Call newArray() for D = D-1
            for( size_t i = 1; i < n1; ++i ) {
                ret[i] = ret[i - 1] + n2;         // Populate the pointer-array with addresses,
            }                                     // starting from p[1] with step length n2.
            return ret;
        }
        //@}

        
        /*! @name makePointers
         *  @brief Generate a pointer structure for an already allocated memory-block, so it can be accessed as a multi-dimensional array.
         *  @param n1 Size of dimension
         *  @param n2 Size of dimension
         *  @param rest Sizes of dimensions
         */
        //@{
        template <class T>  T* makePointers( T* data, size_t n1 ) { return data; }
        template <class T, typename... Rest>
        typename detail::Dummy<T, size_t, size_t, Rest...>::ptrType makePointers( T* data, size_t n1, size_t n2, Rest ...rest ) {
            typedef typename detail::Dummy<T, size_t, Rest...>::ptrType ptrType;
            ptrType* ret = new ptrType[ n1 ];               // Create an array of size n1 with pointers.
            *ret = makePointers<T>(data, n1 * n2, rest... );         // Call newArray() for D = D-1
            for( size_t i = 1; i < n1; ++i ) {
                ret[i] = ret[i - 1] + n2;         // Populate the pointer-array with addresses,
            }                                     // starting from p[1] with step length n2.
            return ret;
        }
        //@}


        /*! @brief Create a new shared_ptr wrapped array of arbitrary dimensions/type
         *  @details Recursively calls all lower D versions
         */
        template <class T, typename ...S>
        std::shared_ptr<typename detail::Dummy<T,S...>::dataType> sharedArray( S ...sizes ) {
            typedef typename detail::Dummy<T,S...>::dataType dataType;
            return std::shared_ptr<dataType>( newArray<T>( sizes... ), []( dataType* p ) { delArray( p ); } );
        }

        
        /*! @brief Create a new shared_ptr wrapped array of pointers over a data-block. Pointers will be cleaned up on destruction.
         *  @details Recursively calls all lower D versions
         *  @note This wrapper does not store a pointer to the underlying data block, hence it will become invalid if the
         *  block is freed.
         */
        template <class T, typename ...S>
        std::shared_ptr<typename detail::Dummy<T,S...>::dataType> reshapeArray( T* array, S ...sizes ) {
            typedef typename detail::Dummy<T,S...>::dataType dataType;
            return std::shared_ptr<dataType>( makePointers<T>( array, sizes... ), []( dataType* p ) { delPointers( p ); } );
        }


        /*! @} */



    }   // namespace util

} // end namespace redux

#endif // REDUX_UTIL_ARRAYUTIL_HPP
