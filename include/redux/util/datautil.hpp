#ifndef REDUX_UTIL_DATAUTIL_HPP
#define REDUX_UTIL_DATAUTIL_HPP

#include "redux/util/endian.hpp"

#include <cstring>
#include <string>
#include <vector>

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */

        /*!  @file  datautil.hpp
         *   @brief     Collection of functions for memory/variable manipulation
         *   @details   Functions convenient for manipulation of variables, endianswapping, bit-flipping etc.
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         */

        /*!  Helper for swapping elements between 2 arrays.
         *   @details This is just an iterator which calls std::swap
         */
        template <class T>
        void swap(T* a, T* b, size_t n = 1) {
            if(a == b) {
                return;
            }
            while(n--) {
                std::swap(*a++, *b++);       // uses move semantics since c++11
            }
        }



        /*! @brief Find the new location of element "i" after transposing an MxN matrix.
        *   @param i Original location (as linear offset from the first element in the matrix)
        *   @param sizeY Size of the first dimension ( "row", slow )
        *   @param sizeX Size of the second dimension ( "column", fast )
        *   @returns The new location
        */
        inline size_t nextLocation( int i, size_t sizeY, size_t sizeX ) {
            return ( i % sizeX ) * sizeY + i / sizeX;
        }

        /*! @brief Transpose a 2D array (MxN matrix) of arbitrary type and size.
         *  @details Performs an in-place transpose of the matrix. The Method works by traversing through the permutaion-cycles, swapping the elements pairwise.
         *  @param data Input 2D Array (matrix)
         *  @param sizeY Size of the first dimension ("row",slow)
         *  @param sizeX Size of the second dimension ("column",fast)
         */
        template <class T>
        void transpose( T* data, size_t sizeY, size_t sizeX ) {

            if( sizeY && sizeX && data ) {
                size_t stillToMove = sizeY * sizeX;
                for( size_t i = 0; stillToMove; ++i ) {
                    size_t j, k;
                    for( j = nextLocation( i, sizeY, sizeX ); j > i; j = nextLocation( j, sizeY, sizeX ) ) ; //cycle.push_back(j+1);
                    if( j < i ) continue; // If true, we already traversed this cycle earlier.
                    // Note: j=nextLocation(i,sizeY,sizeX) => i=nextLocation(j,sizeX,sizeY) (cf. transposing MxN vs. NxM)
                    // We need to traverse the cycle backwards to get the elements in the right place by simple swapping, so interchange sizeY & sizeX
                    for( k = i, j = nextLocation( i, sizeX, sizeY ); j != i; k = j, j = nextLocation( j, sizeX, sizeY ) ) {
                        std::swap( data[k], data[j] );
                        --stillToMove;
                    }
                    --stillToMove;
                }
            }
        }

        /*! @name pack
         *  @brief Pack the data into a dense string of characters, suitable for sending across the network or write to file.
         *  @param   ptr Where to store the data
         *  @param   in  Data
         *  @param   count Number of elements to pack
         */
        //@{
        template <typename T>
        inline char* pack(char* ptr, const std::vector<T>& in) {
            size_t sz = in.size()*sizeof(T);
            *reinterpret_cast<size_t*>(ptr) = sz;
            ptr+=sizeof(size_t);
            memcpy(ptr,reinterpret_cast<const char*>(in.data()),sz);
            return ptr+sz;
        }
        
        template <typename T>
        inline char* pack(char* ptr, const T* in, size_t count) {
            size_t sz = count*sizeof(T);
            memcpy(ptr, reinterpret_cast<const char*>(in), sz);
            return ptr+sz;
         }
        
        template <typename T>
        inline char* pack(char* ptr, const T& in) {
            *reinterpret_cast<T*>(ptr) = in;
            return ptr+sizeof(T);
        }
        
        template <>
        inline char* pack<std::string>(char* ptr, const std::string& in) {
            strcpy(ptr,in.c_str());
            return ptr + in.length() + 1;
        }
        //@}
        
        /*! @name unpack
         *  @brief Unpack the data stored in a dense string of characters.
         *  @param   ptr Where to store the data
         *  @param   in  Data
         *  @param   count Number of elements to pack
         *  @param   swap_endian Should the endianess be swapped ?
         */
        //@{
        template <typename T>
        inline const char* unpack(const char* ptr, std::vector<T>& out, bool swap_endian=false) {
            size_t sz = *reinterpret_cast<const size_t*>(ptr);
            ptr+=sizeof(size_t);
            if(swap_endian) swapEndian(sz);
            out.resize(sz/sizeof(T));
            memcpy(reinterpret_cast<char*>(out.data()),ptr,sz);
            if(swap_endian) {
                swapEndian(&out[0],out.size());
            }
            return ptr+sz;
        }
        
        template <typename T>
        inline const char* unpack(const char* ptr, T* out, size_t count, bool swap_endian=false) {
            size_t sz = count*sizeof(T);
            memcpy(reinterpret_cast<char*>(out),ptr,sz);
            if(swap_endian) {
                swapEndian(out,count);
            }
            return ptr+sz;
        }
        
        template <typename T>
        inline const char* unpack(const char* ptr, T& out, bool swap_endian=false) {
            out = *reinterpret_cast<const T*>(ptr);
            if(swap_endian) {
                swapEndian(out);
            }
            return ptr+sizeof(T);
        }
        
        template <>
        inline const char* unpack<std::string>(const char* ptr, std::string& out, bool) {
            out = std::string(ptr);
            return ptr + out.length() + 1;
        }
        //@}
        
        /*! @} */


    }  // namespace util

}  // namespace redux



#endif // REDUX_UTIL_DATAUTIL_HPP
