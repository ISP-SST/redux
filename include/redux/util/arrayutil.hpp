#ifndef REDUX_UTIL_ARRAYUTIL_HPP
#define REDUX_UTIL_ARRAYUTIL_HPP

#include <memory>
#include <iostream>

namespace redux {

    namespace util {

        /*!  @ingroup util
         *  @{
         */

        /*!  @file  array.hpp
         *   @brief     Set of templates for creating/deleting arrays of arbitrary types and dimensionality.
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         */

        template <class T>  void delArray( T*& p ) { delete[] p; }
        template <class T>
        void delArray( T**& p ) {
            delArray( *p );              // Call delArray() for the next pointer-level.
            delete[] p;                  // Delete current pointer-level.
        }


        template <class T>  void delPointers( T* p ) { }
        template <class T>
        void delPointers( T**& p ) {
            delPointers( *p );           // Call delPointers() for the next pointer-level.
            delete[] p;                  // Delete current pointer-level.
        }

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


        template <class T>
        T* newArray( size_t n1 ) {
            T* ret = new T[ n1 ];
            return ret;
        }

        /*! @brief Create an array of arbitrary type and dimensionality
         *  @details Recursively calls all lower D versions to allocate the memoryblock, then asigns pointers
         */
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

        
        template <class T>  T* makePointers( T* data, size_t n1 ) { return data; }

        /*! @brief Create a new pointer-structure over a data-block
         *  @details Recursively calls all lower D versions to allocate/set pointers
         */
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


        /*! @brief Create a new shared_ptr wrapped array of arbitrary dimensions/type
         *  @details Recursively calls all lower D versions
         */
        template <class T, typename ...S>
        std::shared_ptr<typename detail::Dummy<T,S...>::dataType> sharedArray( S ...sizes ) {
            typedef typename detail::Dummy<T,S...>::dataType dataType;
            return std::shared_ptr<dataType>( newArray<T>( sizes... ), []( dataType* p ) { delArray( p ); } );
        }

        
        /*! @brief Create a new shared_ptr wrapped array of pointers over a data-block
         *  @details Recursively calls all lower D versions
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
