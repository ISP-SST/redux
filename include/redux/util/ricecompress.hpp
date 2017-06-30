#ifndef RICECOMPRESS_HPP
#define RICECOMPRESS_HPP

#include <cstddef>
#include <cstdint>


namespace redux {

    namespace util {

        int rice_comp16( const int16_t* in, size_t inSize, uint8_t* out, size_t outSize, size_t blockSize );
        int rice_decomp16( const uint8_t* in, size_t inSize, int16_t* out, size_t outSize, size_t blockSize );

    }   // util
    
}       // redux

#endif      // RICECOMPRESS_HPP
