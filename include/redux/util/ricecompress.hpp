#ifndef RICECOMPRESS_HPP
#define RICECOMPRESS_HPP

#include <cstddef>
#include <cstdint>


namespace redux {

    namespace util {

        int rice_comp16( const int16_t* in, size_t inSize, uint8_t* out, size_t outSize, size_t blockSize );
        
        int rice_decomp16( const uint8_t* in, size_t inSize, int16_t* out, size_t outSize, size_t blockSize );

        uint32_t ana_comp16_old( uint8_t* out, const int16_t* in, int slice, int blockSize, int nBlocks, uint32_t outSize, int t_endian );
        int ana_comp16( const int16_t* in, size_t inSize, uint8_t* out, size_t outSize, int slice, size_t blockSize );
        
        int ana_decomp16( const uint8_t* in, size_t inSize, int16_t* out, size_t outSize, size_t blockSize );

    }   // util
    
}       // redux

#endif      // RICECOMPRESS_HPP
