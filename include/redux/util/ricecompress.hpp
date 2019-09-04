#ifndef RICECOMPRESS_HPP
#define RICECOMPRESS_HPP

#include <cstddef>
#include <cstdint>


namespace redux {

    namespace util {

        int rice_comp16( const int16_t* in, size_t inSize, uint8_t* out, size_t outSize, size_t blockSize );
        int rice_comp16_bad( const int16_t* in, size_t inSize, uint8_t* out, size_t outSize, size_t blockSize );
        
        int rice_decomp16_orig( const uint8_t* in, size_t inSize, int16_t* out, size_t outSize, size_t blockSize );
        // This routine is a hacked version to read corrupted files, written with a buggy version of the camera-code
        // at SST between 2019-08-21 and 2019-09-02
        int rice_decomp16_fix( const uint8_t* in, size_t inSize, int16_t* out, size_t outSize, size_t blockSize );
        int rice_decomp16_dbg( const uint8_t* in, size_t inSize, int16_t* out, size_t outSize, size_t blockSize );
#ifdef FIX_FITS_201908
#define rice_decomp16 rice_decomp16_fix
#else
#define rice_decomp16 rice_decomp16_orig
#endif
        
    }   // util
    
}       // redux

#endif      // RICECOMPRESS_HPP
