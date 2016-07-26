#ifndef REDUX_UTIL_DATAUTIL_HPP
#define REDUX_UTIL_DATAUTIL_HPP

#include "redux/util/endian.hpp"

#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <math.h>

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

        
        template <typename T>
        inline bool checkAllZero(const T* ptr, size_t count) {
            for(size_t i=0; i<count; ++i ) if( ptr[i] != T(0) ) return false;
            return true;
        }

        template <typename T>
        inline bool checkAllSmaller(const T* ptr, size_t count, double eps=1E-10 ) {
            for(size_t i=0; i<count; ++i ) if( fabs(ptr[i]) > eps ) return false;
            return true;
        }

        
         /*! @name size
         *  @brief Return the proper size needed for packing/unpacking certain types/structures.
         */
        //@{
        inline uint64_t size( const std::vector<std::string>& data ) {
            uint64_t sz = sizeof(uint64_t);
            for (auto &str: data) {
                sz += str.length() + 1;      // pack as zero-terminated c-style string
            }
            return sz;
        }
       
        template <typename T>
        inline uint64_t size( const std::vector<T>& data ) {
            return data.size()*sizeof(T) + sizeof(uint64_t);
        }
       
        template <class T, class U, class Comp, class Alloc>
        inline uint64_t size(const std::map<T,U,Comp,Alloc>& in) {
            uint64_t sz = sizeof(uint64_t);
	    sz += (sizeof(T) + sizeof(U)) * in.size();
            return sz;
        }
       

        //@}
        
        /*! @name pack
         *  @brief Pack the data into a dense string of characters, suitable for sending across the network or write to file.
         *  @param   ptr Where to store the data
         *  @param   in  Data
         *  @param   count Number of elements to pack
         */
        //@{

        inline uint64_t pack(char* ptr, const std::vector<std::string>& in) {
            uint64_t count = in.size();
            *reinterpret_cast<uint64_t*>(ptr) = count;
            uint64_t totalSize = sizeof(uint64_t);
            for ( auto &str: in ) {
                strncpy( ptr+totalSize, str.c_str(), str.length() + 1 );
                totalSize += str.length() + 1;
            }
            return totalSize;
        }
        
        template <typename T>
        inline uint64_t pack(char* ptr, const std::vector<T>& in) {
            uint64_t sz = in.size();
            *reinterpret_cast<uint64_t*>(ptr) = sz;
            memcpy(ptr+sizeof(uint64_t),reinterpret_cast<const char*>(in.data()),sz*sizeof(T));
            return sz*sizeof(T)+sizeof(uint64_t);
        }
        
        template <class T, class U, class Comp, class Alloc>
        inline uint64_t pack(char* ptr, const std::map<T,U,Comp,Alloc>& in) {
            uint64_t count = in.size();
            *reinterpret_cast<uint64_t*>(ptr) = count;
            uint64_t totalSize = sizeof(uint64_t);
            for (auto &element: in) {
                *reinterpret_cast<T*>(ptr+totalSize) = element.first;
                totalSize += sizeof(T);
                *reinterpret_cast<U*>(ptr+totalSize) = element.second;
                totalSize += sizeof(U);
            }
            return totalSize;
        }
        
        template <class T, class Comp, class Alloc>
        inline uint64_t pack(char* ptr, const std::set<T,Comp,Alloc>& in) {
            uint64_t count = in.size();
            *reinterpret_cast<uint64_t*>(ptr) = count;
            uint64_t totalSize = sizeof(uint64_t);
            for (auto &element: in) {
                *reinterpret_cast<T*>(ptr+totalSize) = element;
                totalSize += sizeof(T);
            }
            return totalSize;
        }
        
        template <typename T>
        inline uint64_t pack(char* ptr, const T* in, size_t count) {
            uint64_t sz = count*sizeof(T);
            memcpy(ptr, reinterpret_cast<const char*>(in), sz);
            return sz;
         }
        
        template <typename T>
        inline uint64_t pack(char* ptr, const T& in) {
            *reinterpret_cast<T*>(ptr) = in;
            return sizeof(T);
        }
        
        template <>
        inline uint64_t pack<std::string>(char* ptr, const std::string& in) {
            strncpy(ptr,in.c_str(), in.length() + 1);
            return in.length() + 1;
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
        inline uint64_t unpack(const char* ptr, std::vector<std::string>& out, bool swap_endian) {
            uint64_t count = *reinterpret_cast<const uint64_t*>(ptr);
            if(swap_endian) swapEndian(count);
            out.resize(count);
            uint64_t totalSize = sizeof(uint64_t);
            for( auto& str: out ) {
                str = ptr+totalSize;
                totalSize += str.length()+1;
            }
            return totalSize;
        }
        
        template <typename T>
        inline uint64_t unpack(const char* ptr, std::vector<T>& out, bool swap_endian=false) {
            uint64_t sz = *reinterpret_cast<const uint64_t*>(ptr);
            if(swap_endian) swapEndian(sz);
            out.resize(sz);
            memcpy(reinterpret_cast<char*>(out.data()),ptr+sizeof(uint64_t),sz*sizeof(T));
            if(swap_endian) {
                swapEndian(&out[0],out.size());
            }
            return sz*sizeof(T)+sizeof(uint64_t);
        }
        
        template <class T, class U, class Comp, class Alloc>
        inline uint64_t unpack(const char* ptr, std::map<T,U,Comp,Alloc>& out, bool swap_endian=false) {
            uint64_t count = *reinterpret_cast<const uint64_t*>(ptr);
            uint64_t totalSize = sizeof(uint64_t);
            if(swap_endian) swapEndian(count);
            for( uint64_t cnt=0; cnt<count; ++cnt) {
                T tmpT = *reinterpret_cast<const T*>(ptr+totalSize);
                totalSize += sizeof(T);
                U tmpU = *reinterpret_cast<const U*>(ptr+totalSize);
                totalSize += sizeof(U);
                if(swap_endian) {
                    swapEndian(tmpT);
                    swapEndian(tmpU);
                }
                out.insert(std::make_pair(tmpT,tmpU));
            }
            return totalSize;
        }
        
        template <class T, class Comp, class Alloc>
        inline uint64_t unpack(const char* ptr, std::set<T,Comp,Alloc>& out, bool swap_endian=false) {
            uint64_t count = *reinterpret_cast<const uint64_t*>(ptr);
            uint64_t totalSize = sizeof(uint64_t);
            if(swap_endian) swapEndian(count);
            for( uint64_t cnt=0; cnt<count; ++cnt) {
                T tmpT = *reinterpret_cast<const T*>(ptr+totalSize);
                totalSize += sizeof(T);
                if(swap_endian) {
                    swapEndian(tmpT);
                }
                out.insert(tmpT);
            }
            return totalSize;
        }
        
        template <typename T>
        inline uint64_t unpack(const char* ptr, T* out, size_t count, bool swap_endian=false) {
            uint64_t sz = count*sizeof(T);
            memcpy(reinterpret_cast<char*>(out),ptr,sz);
            if(swap_endian) {
                swapEndian(out,count);
            }
            return sz;
        }
        
        template <typename T>
        inline uint64_t unpack(const char* ptr, T& out, bool swap_endian=false) {
            out = *reinterpret_cast<const T*>(ptr);
            if(swap_endian) {
                swapEndian(out);
            }
            return sizeof(T);
        }
        
        template <>
        inline uint64_t unpack<std::string>(const char* ptr, std::string& out, bool) {
            out = std::string(ptr);
            return out.length() + 1;
        }
        //@}
        
        
        


        //@{
        /*!
         *   @brief Find the median of **begin, **end and their midpoint, move the median to end.
         *   @details This function is primarily for use with QuickSort/partition.
         */
        template <class T>
        void medianOf3( T** begin, T** end ) {

            if( end < begin + 2 ) {
                return;
            }

            T** mid = begin + ( ( end - begin + 1 ) >> 1 );

            if( **mid > **end ) {
                std::swap( mid, end );  // first sort mid/end
            }

            if( **end > **begin ) {             // if begin >= end, we want end (i.e. we're done.)
                if( **mid > **begin ) {
                    std::swap( mid, end );   // if begin < mid, we want mid
                }
                else {
                    std::swap( begin, end );          // else we want begin
                }
            }

        }

        template <class T>
        void medianOf3( T* begin, T* end ) {

            if( end < begin + 2 ) {
                return;
            }

            T* mid = begin + ( ( end - begin + 1 ) >> 1 );

            if( *mid > *end ) {
                std::swap( mid, end );    // first sort mid/end
            }

            if( *end > *begin ) {           // if begin >= end, we want end (i.e. we're done.)
                if( *mid > *begin ) {
                    std::swap( mid, end ); // if begin < mid, we want mid
                }
                else {
                    std::swap( begin, end );          // else we want begin
                }
            }

        }
        
        template <class T>
        T medianOf3( T* ptr ) {

            if( *ptr > *(ptr+1) ) {
                if( *ptr < *(ptr+2) ) return *ptr;
                else if( *(ptr+1) < *(ptr+2) ) return *(ptr+2);
                else return *(ptr+1);
            } else {
                if( *ptr < *(ptr+2) ) return *(ptr+1);
                else if( *(ptr+1) < *(ptr+2) ) return *(ptr+2);
                else return *(ptr);
            }

        }
        //@}

        //@{
        /*! @fn void partition( T** begin, T** end, bool descending=false )
         *  @details Arrange the array from begin to end in two parts, one with elements > **end, and one with the rest.
         *   I.e. the pivot-element is assumed to be in **end.
         *   This function is for use with QuickSort.
         *  @param begin Pointer to the first element of the array to be sorted.
         *  @param end Pointer to the last element of the array to be sorted.
         *  @param descending Should the part that is larger than the pivot be first ?
         *  @author Tomas Hillberg
         */
        template <class T>
        void partition( T** begin, T** end, bool descending = false ) {

            if( end < begin + 1 ) {
                return;
            }

            T** iter;
            T** store = begin;

            if( descending ) {
                for( iter = begin; iter < end; ++iter ) {
                    if( **iter > **end ) {
                        std::swap( iter, store );
                        store++;
                    }
                }
            }
            else {
                for( iter = begin; iter < end; ++iter ) {
                    if( !( **iter > **end ) ) {
                        std::swap( iter, store );
                        store++;
                    }
                }
            }

            std::swap( store, end );    // Move pivot to it's proper place

            if( store == begin || store >= end ) {
                return;
            }

            if( store > begin ) {
                // Pick a new pivot for the first part, and partition it.
                medianOf3( begin, store );  // stores the median of begin/store/midpoint in store
                partition( begin, store, descending );
            }

            if( ++store < end ) {
                // Pick a new pivot for the second part, and partition it.
                medianOf3( store, end );    // stores the median of store/end/midpoint in end
                partition( store, end, descending );
            }

        }


        template <class T>
        void partition( T* begin, T* end, bool descending = false ) {

            if( end < begin + 1 ) {
                return;
            }

            T* iter;
            T* store = begin;

            if( descending ) {
                for( iter = begin; iter < end; ++iter ) {
                    if( *iter > *end ) {
                        std::swap( iter, store++ );
                    }
                }
            }
            else {
                for( iter = begin; iter < end; ++iter ) {
                    if( !( *iter > *end ) ) {
                        std::swap( iter, store++ );
                    }
                }
            }

            std::swap( store, end );    // Move pivot to it's proper place

            if( store == begin || store >= end ) {
                return;
            }

            if( store > begin ) {
                // Pick a new pivot for the first part, and partition it.
                medianOf3( begin, store );  // stores the median of begin/store/midpoint in store
                partition( begin, store, descending );
            }

            if( ++store < end ) {
                // Pick a new pivot for the second part, and partition it.
                medianOf3( store, end );    // stores the median of store/end/midpoint in end
                partition( store, end, descending );
            }

        }
        //@}


        //@{
        /*!
         *  @brief Sorting algorithm. Version for an array of pointers.
         *  @details Template for sorting arbitrary types/classes, it is written in a way that it only depends
         *   on the ">" operator existing for the class.
         *  @param data First element in the array to be sorted
         *  @param n Length of the array
         *  @param descending Should the largest element be first ?
         *  @param pivotIndex Select an index to be used as pivot-element, else one will be selected using medianOf3()
         *  @author Tomas Hillberg
         *  @todo Optimization: exit partitioning for small arrays (15-20 elements), and use a better sort-method (direct insert?) \n
         *     Check if array is (almost) sorted, if so use another method.
         */
        template <class T>
        void quickSort( T** data, size_t n, bool descending = false, size_t pivotIndex = -1 ) {

            if( n < 2 ) {
                return;
            }

            T** ptr = data;
            T** end = data + n - 1;
            size_t nanCount = 0;

            // If pivotIndex is out of range, call medianOf3, else move pivot to end of array.
            if( pivotIndex < n ) {
                ptr = data + pivotIndex;
                std::swap( ptr, end );
                ptr = data;
            }
            else {
                medianOf3( ptr, end );
            }

            while( ptr < end ) {
                if( !isfinite(**ptr) ) {   // check for uncomparable objects (e.g. NaN's) and place them at the end.
                    std::swap( ptr, --end );
                    nanCount++;
                }
                ++ptr;
            }

            if( nanCount ) {
                ptr = data + n - 1;
                std::swap( end, ptr );
            }

            // recurseviely partition the array
            partition( data, end, descending );

        }


        template <class T>
        void quickSort( T* data, size_t n, bool descending = false, size_t pivotIndex = -1 ) {

            if( n < 2 ) {
                return;
            }

            T* ptr = data;
            T* end = data + n - 1;
            size_t nanCount = 0;

            // If pivotIndex is out of range, call medianOf3, else move pivot to end of array.
            if( pivotIndex < n ) {
                std::swap( data + pivotIndex, end );
            }
            else {
                medianOf3( ptr, end );
            }

            while( ptr < end ) {
                if( !isfinite(*ptr) ) {   // check for uncomparable objects (e.g. NaN's) and place them at the end.
                    std::swap( ptr, --end );
                    nanCount++;
                }
                ++ptr;
            }

            if( nanCount ) {
                std::swap( end, data + n - 1 );
            }

            // recurseviely partition the array
            partition<T>( data, end, descending );

        }

        //@}


        /*! @brief Segment a range from first to last (inklusive) into segments of a given length
        *   @param first First value in range.
        *   @param last Last value in range.
        *   @param segmentLength Size of the segments
        *   @param minimumOverlap Smallest accepted overlap
        *   @returns A vector containing the locations of the segments.
        */
        template <typename T>
        std::vector<T> segment( T first, T last, int segmentLength, int minimumOverlap=0 ) {
            int nSegments=2;
            if( first == last ) return std::move(std::vector<T>());
            if( first > last ) std::swap( first, last );
            double separation = (last-first)/static_cast<double>(nSegments-1);
            double overlap = segmentLength-separation;
            while(overlap < minimumOverlap) {
                ++nSegments;
                separation = (last-first)/static_cast<double>(nSegments-1);
                overlap = segmentLength-separation;
            }
            std::vector<T> ret;
            for(int i = 0; i < nSegments; ++i) {
                ret.push_back(static_cast<T>(i*separation+first) );
            }
            return std::move(ret);
        }

        
        
        /*! @} */


    }  // namespace util

}  // namespace redux


#if defined(__GNUC__)
#  define UNUSED __attribute__ ((unused))
#elif defined(_MSC_VER)
#  define UNUSED __pragma(warning(suppress:4100))
#else
#  define UNUSED
#endif


#endif // REDUX_UTIL_DATAUTIL_HPP
