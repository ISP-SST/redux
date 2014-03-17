#ifndef REDUX_UTIL_BITOPERATIONS_HPP
#define REDUX_UTIL_BITOPERATIONS_HPP

#include <cstdint>
#include <cstdlib>

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */


        /*!  @file      bitoperations.hpp
         *   @brief     Collection of functions for bit-manipulation
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         */


        namespace detail {

            const uint64_t deBruijnMagic32 = 0x077CB531U;
            static const int deBruijnSequence32[32] = {
                0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
                31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
            };
            static const int deBruijnSequence32p1[32] = {   // sequence shifted by 1 to get 1-based indexing
                1, 2, 29, 3, 30, 15, 25, 4, 31, 23, 21, 16, 26, 18, 5, 9,
                32, 28, 14, 24, 22, 20, 17, 8, 27, 13, 19, 7, 12, 6, 11, 10
            };

            const uint64_t deBruijnMagic64 = 0x03f79d71b4cb0a89U;
            static const int deBruijnSequence64[64] = {
                0,  1, 48,  2, 57, 49, 28,  3,
                61, 58, 50, 42, 38, 29, 17,  4,
                62, 55, 59, 36, 53, 51, 43, 22,
                45, 39, 33, 30, 24, 18, 12,  5,
                63, 47, 56, 27, 60, 41, 37, 16,
                54, 35, 52, 21, 44, 32, 23, 11,
                46, 26, 40, 15, 34, 20, 31, 10,
                25, 14, 19,  9, 13,  8,  7,  6
            };
            static const int deBruijnSequence64p1[64] = {   // sequence shifted by 1 to get 1-based indexing
                1,  2, 49,  3, 58, 50, 29,  4,
                62, 59, 51, 43, 39, 30, 18,  5,
                63, 56, 60, 37, 54, 52, 44, 23,
                46, 40, 34, 31, 25, 19, 13,  6,
                64, 48, 57, 28, 61, 42, 38, 17,
                55, 36, 53, 22, 45, 33, 24, 12,
                47, 27, 41, 16, 35, 21, 32, 11,
                26, 15, 20, 10, 14,  9,  8,  7
            };

            static const int binaryMagicNumbers[5] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF, 0x0000FFFF};

        }

        /*! @fn int log2(uint32_t v)
         *  @brief Fast, LUT based, calculation of the base 2 logarithm for unsigned integers.
         *  @param      v value
         *  @returns    log2(v) (rounded down)
         */
        inline int log2( uint32_t v ) {
            v |= v >> 1;
            v |= v >> 2;
            v |= v >> 4;
            v |= v >> 8;
            v |= v >> 16;
            v = ( v >> 1 ) + 1;
            return detail::deBruijnSequence32[( ( uint32_t )( v * detail::deBruijnMagic32 ) ) >> 27];
        }

        /*! @name findLSB  
         *  @brief      Fast, LUT based, calculation of the lowest non-zero bit in an integer
         *  @details    findLSB1 & findLSB64 returns a 1-based bit-index, and findLSB returns 0-based bit-index.
         *  @param      v value
         *  @returns    position of LSB
         *  @code
         *  int a = 17376;      // = 100001111100000
         *  int b = findLSB(a); // = 5
         *  @endcode
         */
        //@{
        inline int findLSB( uint32_t v ) { // returns 0-based index of first non-zero bit.
            return detail::deBruijnSequence32[( ( uint32_t )( ( v & -v ) * detail::deBruijnMagic32 ) ) >> 27];
        }

        inline int findLSB1( uint32_t v ) { // returns 1-based index of first non-zero bit.
            if( !v ) {
                return 0;
            }
            return detail::deBruijnSequence32p1[( ( uint32_t )( ( v & -v ) * detail::deBruijnMagic32 ) ) >> 27];
        }

        inline int findLSB64( uint64_t v ) { // returns 1-based index of first non-zero bit.
            if( !v ) {
                return 0;
            }
            return detail::deBruijnSequence64p1[( ( v & -v ) * detail::deBruijnMagic64 ) >> 58];
        }
        //@}


        // http://burtleburtle.net/bob/hash/doobs.html
        inline uint32_t mix( uint32_t a, uint32_t b, uint32_t c ) {
            a = a - b;  a = a - c;  a = a ^ ( c >> 13 );
            b = b - c;  b = b - a;  b = b ^ ( a << 8 );
            c = c - a;  c = c - b;  c = c ^ ( b >> 13 );
            a = a - b;  a = a - c;  a = a ^ ( c >> 12 );
            b = b - c;  b = b - a;  b = b ^ ( a << 16 );
            c = c - a;  c = c - b;  c = c ^ ( b >> 5 );
            a = a - b;  a = a - c;  a = a ^ ( c >> 3 );
            b = b - c;  b = b - a;  b = b ^ ( a << 10 );
            c = c - a;  c = c - b;  c = c ^ ( b >> 15 );
            return c;
        }


        /*! @fn int nextPowerOfTwo(uint32_t v)
         *  @brief Returns the next power of two larger than (or equal to) v
         *  @param      v value
         *  @returns    next power of two larger than (or equal to) v
         *  @code
         *  int a = 17383;
         *  int b = nextPowerOfTwo(a); // = 32768
         *  @endcode
         */
        inline int nextPowerOfTwo( uint32_t v ) {
            v--;
            v |= v >> 1;
            v |= v >> 2;
            v |= v >> 4;
            v |= v >> 8;
            v |= v >> 16;
            v++;
            return v;
        }

        inline uint32_t nextEven( uint32_t v, int power ) {
            uint32_t tmp = ( 1 << power ) - 1;
            return ( v + tmp ) & ~tmp;
        }

        /*! @name swapBits
         *  @brief      swap the bits set in "mask" between a and b
         *  @param      a Input
         *  @param      b Input
         *  @param      mask a bitmask specifying which bits to swap between a and b
         *  @param      n number of elements in input arrays
         */
        //@{
        template <typename T>
        inline void swapBits( T& a, T& b, const T mask = 0xFFFFFFFFFFFFFFFF ) {
            T masked = ( a ^ b ) & mask;
            a ^= masked;
            b ^= masked;
        }

        template <typename T>
        void swapBits( T* a, T* b, const T mask, size_t n = 1 ) {
            for( size_t i = 0; i < n; ++i ) {
                swapBits( *a++, *b++, mask );
            }
        }
        //@}
        /*! @fn         template <typename T> void flipBits(T* a, const T mask, size_t n = 1)
         *  @brief      flip the bits in a which are set in "mask"
         *  @param      a pointer to the data to be "flipped"
         *  @param      mask a bitmask specifying which bits to flip
         *  @param      n number of elements to flip
         */
        template <typename T>
        void flipBits( T* a, const T mask, size_t n = 1 ) {
            for( size_t i = 0; i < n; ++i ) {
                *a++ ^= mask;
            }
        }

        //@{
        /*!
         *  @brief      Count the number of set bits in v.
         *  @param      v Input
         *  @todo       optimize (arch-specific or asm ??), fix 16/8-bit implementation.
         */
        inline uint64_t countBits( uint64_t v ) {
            v = v - ( ( v >> 1 ) & 0x5555555555555555UL );                // reuse input as temporary
            v = ( v & 0x3333333333333333UL ) + ( ( v >> 2 ) & 0x3333333333333333UL ); // temp
            return ( ( ( v + ( v >> 4 ) ) & 0x0F0F0F0F0F0F0F0FUL ) * 0x0101010101010101UL ) >> 56; // count
        }

        //@}

        /*
        A generalization of the best bit counting method to integers of bit-widths upto 128 (parameterized by type T) is this:
        v = v - ((v >> 1) & (T)~(T)0/3);                           // temp
        v = (v & (T)~(T)0/15*3) + ((v >> 2) & (T)~(T)0/15*3);      // temp
        v = (v + (v >> 4)) & (T)~(T)0/255*15;                      // temp
        c = (T)(v * ((T)~(T)0/255)) >> (sizeof(T) - 1) * CHAR_BIT; // count
        */

        /*! @} */

    }

}

#endif // REDUX_UTIL_BITOPERATIONS_HPP
