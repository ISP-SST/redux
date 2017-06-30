#include "redux/util/ricecompress.hpp"

#include "redux/util/stringutil.hpp"

#include <string.h>             // memset
#include <memory>               // unique_ptr
#include <iostream>
#include <immintrin.h>


using namespace redux::util;
using namespace std;

namespace {
    
    const int nLeadingZeroes[256] = {
        8,                                                          // 0
        7,                                                          // 1
        6, 6,                                                       // 2-3
        5, 5, 5, 5,                                                 // 4-7
        4, 4, 4, 4, 4, 4, 4, 4,                                     // 8-15
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,             // 16-31
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,             // 32-63
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,             // 64-127
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,             // 128-255
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };


    const uint16_t split_map_3[8][9] = { // Lookup-table for 3-bit splits, used for 8-bit input data.  Dimensions: [id][offset]
        // NB. This method assumes little endian ordering.
        // <------------------------------->  for offset < 6 we only need the low byte
        //                                    <-------------------->  for offset >= 6 we use both bytes, but swap order so we can AND directly to the out-stream
        { 0x20, 0x10,  0x8,  0x4, 0x2, 0x1, 0x8000, 0x4000, 0x2000 }, // id =   0    =>   stream: 001    (we store split+1)
        { 0x40, 0x20, 0x10,  0x8, 0x4, 0x2,    0x1, 0x8000, 0x4000 }, // id =   1    =>   stream: 010
        { 0x60, 0x30, 0x18,  0xC, 0x6, 0x3, 0x8001, 0xC000, 0x6000 }, // id =  10    =>   stream: 011
        { 0x80, 0x40, 0x20, 0x10, 0x8, 0x4,    0x2,    0x1, 0x8000 }, // id =  11    =>   stream: 100
        { 0xA0, 0x50, 0x28, 0x14, 0xA, 0x5, 0x8002, 0x4001, 0xA000 }, // id = 100    =>   stream: 101
        { 0xC0, 0x60, 0x30, 0x18, 0xC, 0x6,    0x3, 0x8001, 0xC000 }, // id = 101    =>   stream: 110
        { 0xE0, 0x70, 0x38, 0x1C, 0xE, 0x7, 0x8003, 0xC001, 0xE000 }, // id = 110    =>   stream: 111    (max = 6)
        { 0xE0, 0x70, 0x38, 0x1C, 0xE, 0x7, 0x8003, 0xC001, 0xE000 }  // id = 111    =>   stream: 111
    //    <<5   <<4   <<3   <<2   <<1  aligned   <<7     <<6     <<5
    };

    const uint16_t split_map_4[16][9] = { // Lookup-table for 4-bit splits, used for 16-bit input data.  Dimensions: [id][offset]
        // NB. This method assumes little endian ordering.
        // <------------------------->  for offset < 5 we only need the low byte
        //                              <---------------------------->  for offset >= 5 we use both bytes, but swap order so we can AND directly to the out-stream
        { 0x10,  0x8,  0x4, 0x2,  0x1, 0x8000, 0x4000, 0x2000, 0x1000 }, // id =    0    =>   stream: 0001    (we store split+1)
        { 0x20, 0x10,  0x8, 0x4,  0x2,    0x1, 0x8000, 0x4000, 0x2000 }, // id =    1    =>   stream: 0010
        { 0x30, 0x18,  0xC, 0x6,  0x3, 0x8001, 0xC000, 0x6000, 0x3000 }, // id =   10    =>   stream: 0011
        { 0x40, 0x20, 0x10, 0x8,  0x4,    0x2,    0x1, 0x8000, 0x4000 }, // id =   11    =>   stream: 0100
        { 0x50, 0x28, 0x14, 0xA,  0x5, 0x8002, 0x4001, 0xA000, 0x5000 }, // id =  100    =>   stream: 0101
        { 0x60, 0x30, 0x18, 0xC,  0x6,    0x3, 0x8001, 0xC000, 0x6000 }, // id =  101    =>   stream: 0110
        { 0x70, 0x38, 0x1C, 0xE,  0x7, 0x8003, 0xC001, 0xE000, 0x7000 }, // id =  110    =>   stream: 0111
        { 0x80, 0x40, 0x20, 0x10, 0x8,    0x4,    0x2,    0x1, 0x8000 }, // id =  111    =>   stream: 1000
        { 0x90, 0x48, 0x24, 0x12, 0x9, 0x8004, 0x4002, 0x2001, 0x9000 }, // id = 1000    =>   stream: 1001
        { 0xA0, 0x50, 0x28, 0x14, 0xA,    0x5, 0x8002, 0x4001, 0xA000 }, // id = 1001    =>   stream: 1010
        { 0xB0, 0x58, 0x2C, 0x16, 0xB, 0x8005, 0xC002, 0x6001, 0xB000 }, // id = 1010    =>   stream: 1011
        { 0xC0, 0x60, 0x30, 0x18, 0xC,    0x6,    0x3, 0x8001, 0xC000 }, // id = 1011    =>   stream: 1100
        { 0xD0, 0x68, 0x34, 0x1A, 0xD, 0x8006, 0x4003, 0xA001, 0xD000 }, // id = 1100    =>   stream: 1101
        { 0xE0, 0x70, 0x38, 0x1C, 0xE,    0x7, 0x8003, 0xC001, 0xE000 }, // id = 1101    =>   stream: 1110
        { 0xF0, 0x78, 0x3C, 0x1E, 0xF, 0x8007, 0xC003, 0xE001, 0xF000 }, // id = 1110    =>   stream: 1111    (max = 14)
        { 0xF0, 0x78, 0x3C, 0x1E, 0xF, 0x8007, 0xC003, 0xE001, 0xF000 }  // id = 1111    =>   stream: 1111
    //    <<4   <<3   <<2   <<1  aligned  <<7     <<6     <<5     <<4
    };

    const uint16_t split_map_5[32][9] = { // Lookup-table for 5-bit splits, used for 32-bit input data.  Dimensions: [id][offset]
        // NB. This method assumes little endian ordering.
        // <------------------->  for offset < 4 we only need the low byte
        //                        <------------------------------------>  for offset >= 4 we use both bytes, but swap order so we can AND directly to the out-stream
        {  0x8,  0x4,  0x2,  0x1, 0x8000, 0x4000, 0x2000, 0x1000,  0x800 }, // id =     0    =>   stream: 00001    (we store split+1)
        { 0x10,  0x8,  0x4,  0x2,    0x1, 0x8000, 0x4000, 0x2000, 0x1000 }, // id =     1    =>   stream: 00010
        { 0x18,  0xC,  0x6,  0x3, 0x8001, 0xC000, 0x6000, 0x3000, 0x1800 }, // id =    10    =>   stream: 00011
        { 0x20, 0x10,  0x8,  0x4,    0x2,    0x1, 0x8000, 0x4000, 0x2000 }, // id =    11    =>   stream: 00100
        { 0x28, 0x14,  0xA,  0x5, 0x8002, 0x4001, 0xA000, 0x5000, 0x2800 }, // id =   100    =>   stream: 00101
        { 0x30, 0x18,  0xC,  0x6,    0x3, 0x8001, 0xC000, 0x6000, 0x3000 }, // id =   101    =>   stream: 00110
        { 0x38, 0x1C,  0xE,  0x7, 0x8003, 0xC001, 0xE000, 0x7000, 0x3800 }, // id =   110    =>   stream: 00111
        { 0x40, 0x20, 0x10,  0x8,    0x4,    0x2,    0x1, 0x8000, 0x4000 }, // id =   111    =>   stream: 01000
        { 0x48, 0x24, 0x12,  0x9, 0x8004, 0x4002, 0x2001, 0x9000, 0x4800 }, // id =  1000    =>   stream: 01001
        { 0x50, 0x28, 0x14,  0xA,    0x5, 0x8002, 0x4001, 0xA000, 0x5000 }, // id =  1001    =>   stream: 01010
        { 0x58, 0x2C, 0x16,  0xB, 0x8005, 0xC002, 0x6001, 0xB000, 0x5800 }, // id =  1010    =>   stream: 01011
        { 0x60, 0x30, 0x18,  0xC,    0x6,    0x3, 0x8001, 0xC000, 0x6000 }, // id =  1011    =>   stream: 01100
        { 0x68, 0x34, 0x1A,  0xD, 0x8006, 0x4003, 0xA001, 0xD000, 0x6800 }, // id =  1100    =>   stream: 01101
        { 0x70, 0x38, 0x1C,  0xE,    0x7, 0x8003, 0xC001, 0xE000, 0x7000 }, // id =  1101    =>   stream: 01110
        { 0x78, 0x3C, 0x1E,  0xF, 0x8007, 0xC003, 0xE001, 0xF000, 0x7800 }, // id =  1110    =>   stream: 01111
        { 0x80, 0x40, 0x20, 0x10,    0x8,    0x4,    0x2,    0x1, 0x8000 }, // id =  1111    =>   stream: 10000
        { 0x88, 0x44, 0x22, 0x11, 0x8008, 0x4004, 0x2002, 0x1001, 0x8800 }, // id = 10000    =>   stream: 10001
        { 0x90, 0x48, 0x24, 0x12,    0x9, 0x8004, 0x4002, 0x2001, 0x9000 }, // id = 10001    =>   stream: 10010
        { 0x98, 0x4C, 0x26, 0x13, 0x8009, 0xC004, 0x6002, 0x3001, 0x9800 }, // id = 10010    =>   stream: 10011
        { 0xA0, 0x50, 0x28, 0x14,    0xA,    0x5, 0x8002, 0x4001, 0xA000 }, // id = 10011    =>   stream: 10100
        { 0xA8, 0x54, 0x2A, 0x15, 0x800A, 0x4005, 0xA002, 0x5001, 0xA800 }, // id = 10100    =>   stream: 10101
        { 0xB0, 0x58, 0x2C, 0x16,    0xB, 0x8005, 0xC002, 0x6001, 0xB000 }, // id = 10101    =>   stream: 10110
        { 0xB8, 0x5C, 0x2E, 0x17, 0x800B, 0xC005, 0xE002, 0x7001, 0xB800 }, // id = 10110    =>   stream: 10111
        { 0xC0, 0x60, 0x30, 0x18, 0x800C,    0x6,    0x3, 0x8001, 0xC000 }, // id = 10111    =>   stream: 11000
        { 0xC8, 0x64, 0x32, 0x19,    0xC, 0x4006, 0x2003, 0x9001, 0xC800 }, // id = 11000    =>   stream: 11001
        { 0xD0, 0x68, 0x34, 0x1A, 0x800D, 0x8006, 0x4003, 0xA001, 0xD000 }, // id = 11001    =>   stream: 11010    (max = 25)
        { 0xD8, 0x6C, 0x36, 0x1B,    0xD, 0xC006, 0x6003, 0xB001, 0xD800 }, // id = 11010    =>   stream: 11011
        { 0xE0, 0x70, 0x38, 0x1C, 0x800E,    0x7, 0x8003, 0xC001, 0xE000 }, // id = 11011    =>   stream: 11100
        { 0xE8, 0x74, 0x3A, 0x1D,    0xE, 0x4007, 0xA003, 0xD001, 0xE800 }, // id = 11100    =>   stream: 11101
        { 0xF0, 0x78, 0x3C, 0x1E,    0xF, 0x8007, 0xC003, 0xE001, 0xF000 }, // id = 11101    =>   stream: 11110
        { 0xF8, 0x7C, 0x3E, 0x1F, 0x800F, 0xC007, 0xE003, 0xF001, 0xF800 }, // id = 11110    =>   stream: 11111    (max = 31)
        { 0xF8, 0x7C, 0x3E, 0x1F, 0x800F, 0xC007, 0xE003, 0xF001, 0xF800 }  // id = 11111    =>   stream: 11111
     //   <<3   <<2   <<1 aligned    <<7     <<6     <<5     <<4     <<3
    };
    

    union upack_t {
        upack_t( uint64_t a ) : u64(a) {}
        upack_t( const uint8_t* p, uint8_t cnt, uint8_t rshift ) : u64(0) {
            while(cnt--) u8[cnt] = *p++;
            if( rshift ) u64 >>= rshift;
        }
        uint64_t u64;
        uint32_t u32[2];
        uint16_t u16[4];
        uint8_t u8[8];
    };


    __attribute__ ((target ("default")))
    int preprocess_block( const int16_t* __restrict__ block, size_t blockSize, uint32_t* __restrict__ diff, int16_t refValue ) {
        
        // Naïve implementation
        
        double sum(0.0);
        for( size_t i=0; i<blockSize; i++) {
            const int16_t tmp = block[i] - refValue;
            diff[i] = static_cast<uint32_t>((tmp<0) ? ~(tmp<<1) : (tmp<<1));
            sum += diff[i];
            refValue = block[i];
        }

        if( sum == 0.0 ) return -1;
        sum = (sum - (blockSize/2) - 1)/blockSize;
        if( sum < 0 ) return 0;
        
        uint16_t psum = static_cast<uint64_t>( sum ) >> 1;
        int i(0);
        for(; psum>0; ++i) psum >>= 1;
        
        return i;
        
    }


    __attribute__ ((target ("ssse3")))
    int preprocess_block( const int16_t* block, size_t blockSize, uint32_t* diff, int16_t refValue ) {
        
        // N.B. This function assumes that block and diff are 128-bit aligned.

        const __m128i* __restrict__ blockPtr = (__m128i*)block;
        __m128i* __restrict__ diffPtr = (__m128i*)diff;

        union {
            __m128 v;
            float f[4];
        } acc;
        memset( acc.f, 0, 16 );
        
        static const __m128i tmp0 = _mm_setzero_si128();
        __m128i v32;

        size_t cnt(0);
        while( cnt+8 < blockSize ) {
            const __m128i vdiff = _mm_sub_epi16( *blockPtr, _mm_insert_epi16( _mm_slli_si128( *blockPtr, 2 ), refValue, 0 ) );
            const __m128i fill = _mm_cmplt_epi16( vdiff, tmp0 );            // negative values should be unpacked with leading 1's
            v32 = _mm_slli_epi32( _mm_unpacklo_epi16( vdiff, fill ), 1);    // first 4 values as 32-bit integers
            v32 = _mm_xor_si128( _mm_cmplt_epi32( v32, tmp0 ), v32 );
            _mm_store_si128( diffPtr++, v32 );
            acc.v = _mm_add_ps( _mm_cvtepi32_ps( v32 ), acc.v );
            v32 = _mm_slli_epi32( _mm_unpackhi_epi16( vdiff, fill ), 1);    // last 4 values
            v32 = _mm_xor_si128( _mm_cmplt_epi32( v32, tmp0 ), v32 );
            _mm_store_si128( diffPtr++, v32 );
            acc.v = _mm_add_ps( _mm_cvtepi32_ps( v32 ), acc.v );
            cnt += 8;
            refValue = _mm_extract_epi16( *blockPtr++, 7 );                 // save last value as next reference.
        }
        
        // sum the acc vector
        acc.v = _mm_hadd_ps( acc.v, acc.v );
        acc.v = _mm_hadd_ps( acc.v, acc.v );
        float sum = acc.f[0];   // first value now contains the sum of all elements.

        for( ; cnt<blockSize; cnt++) {  // if blockSize is not a multiple of 8, process the last pixels.
            const int16_t tmp = block[cnt] - refValue;
            diff[cnt] = static_cast<uint32_t>((tmp<0) ? ~(tmp<<1) : (tmp<<1));
            sum += diff[cnt];
            refValue = block[cnt];
        }

        if( sum == 0.0 ) return -1;
        
        // select near-optimal split based on the sum of differences
        sum = (sum - (blockSize/2) - 1)/blockSize;
        if (sum < 0) return 0;
        uint64_t psum = static_cast<uint64_t>( sum ) >> 1;
        int i(0);
        for( ; psum>0; ++i) psum >>= 1;
        return i;
     
    }

    __attribute__ ((target ("default")))
    uint8_t* crunch_block( const uint32_t* __restrict__ block, size_t blockSize, uint8_t* __restrict__ ptr, uint8_t& __restrict__ offset, uint8_t split ) {
        //cout << "." << flush;
        const uint32_t unary1 = (1<<split);
        const uint32_t fsmask = unary1 - 1;
        const uint8_t sh = (31-split);
        for ( size_t j=0; j<blockSize; j++ ) {
            const uint32_t tmp1 = (block[j] >> split)+offset;
            ptr += (tmp1>>3);
            offset = (tmp1&7);
            const upack_t tmp2 = ((block[j]&fsmask)|unary1)<<(sh-offset);
            const uint32_t tmp3 = offset+split+1;       // we are setting the unary "1"-bit together with the entropy bits, hence +1
            //*reinterpret_cast<uint32_t*>(ptr) |= htobe32( tmp2.u32[0] );      // this is slower than the below byte-copy
            *ptr |= tmp2.u8[3];
            *(ptr+1) |= tmp2.u8[2];     // needed when tmp3 > 7
            //*(ptr+2) |= tmp2.u8[1];   // needed when tmp3 > 15
            offset = tmp3&7;
            ptr += tmp3>>3;
        }
        return ptr;
    }

    
    /*__attribute__ ((target ("ssse3")))
    uint8_t* crunch_block( const uint32_t* __restrict__ block, size_t blockSize, uint8_t* __restrict__ ptr, uint8_t& __restrict__ offset, uint8_t split ) {
        //cout << "_" << flush;
        const __m128i* blockPtr = (const __m128i*)block;
        __m128i unary_one = _mm_set1_epi32( 1<<split );
        __m128i lo_mask = _mm_set1_epi32( (1<<(split+1))-1 );   // mask covering entropy bits + "unary_one"
        union data_t {
            __m128i data;
            uint32_t i[4];
            uint8_t b[16];
        } buf;
        uint32_t fsmask = (1<<split) - 1;
        for ( size_t j=0; j<blockSize; j++ ) {
            __m128i a = _mm_loadu_si128( blockPtr++ );             // load (possibly unaligned) data

            // extract entropy bits, add a leading 1 (the ending "1" from the unary part)
            __m128i bits = _mm_and_si128( _mm_or_si128(*blockPtr, unary_one), lo_mask );
            
            // extract hi part, i.e. number of leading zeroes in unary code, add cumulative shifts and offset.
            __m128i offsets= _mm_srli_epi32( *blockPtr++, split );
            
            
            int tmp = (block[j] >> split)+offset;
            ptr += tmp>>3;
            offset = tmp&7;
            *ptr |= 13;
            *(ptr+1) |= 13;
            // *p.p32 |= 13;
            tmp = offset+split+1;       // we are setting the unary "1"-bit together with the entropy bits, hence +1
            ptr += tmp>>3;
            offset = tmp&7;
        }
        return ptr;
    }*/


    uint8_t* copy_block( const uint32_t* __restrict__ in, size_t blockSize, uint8_t* out, uint8_t offset ) {
       
        const uint8_t sh = (8-offset);
        for ( size_t j=0; j<blockSize; j++ ) {
            const upack_t tmp4(in[j]<<sh);
            *(out++) |= tmp4.u8[2];
            *(out++) |= tmp4.u8[1];
            *out |= tmp4.u8[0];
        }
        
        return out;
        
    }
    
}


int redux::util::rice_comp16( const int16_t* in, size_t inSize, uint8_t* out, size_t outSize, size_t blockSize ) {

    // For 16-bit data we allow at most 14 bits of entropy, and encode the number of entropy-bits (split) using 4 bits.
    const int max_entropy_bits = 14;
  
    const int16_t* __restrict__ inPtr = in;
    uint8_t* __restrict__ outPtr = out;

    uint32_t* __restrict__ diffPtr;
    std::unique_ptr<int16_t[]> inTmp;
    std::unique_ptr<uint32_t[]> diff;
    try {
        diffPtr = new uint32_t[ blockSize ];
        diff.reset( diffPtr );
        if( ((unsigned long)inPtr)%16 || blockSize%8 ) {
            inTmp.reset( new int16_t[ blockSize ] );
        }
    } catch ( const std::bad_alloc& ) {
        printf( "rice_comp16: bad_alloc." );
        return(-1);
    }

    memset( outPtr, 0, outSize );
    
    uint8_t offset(0);
    size_t thisBlockSize = blockSize;
    int16_t refValue = inPtr[0];
    
    *outPtr++ = refValue>>8;     // high byte of first value
    *outPtr++ = refValue&0xFF;   // low byte

    for( size_t i=0; i<inSize; i += blockSize ) {

        if( i+thisBlockSize >= inSize ) {
            thisBlockSize = inSize - i;
        }

        int split(0);
        if( ((unsigned long)inPtr)%16 ) {
            memcpy( inTmp.get(), inPtr, thisBlockSize*sizeof(int16_t) );
            split = preprocess_block( inTmp.get(), thisBlockSize, diffPtr, refValue );
        } else {
            split = preprocess_block( inPtr, thisBlockSize, diffPtr, refValue );
        }

        if( split >= max_entropy_bits ) {       // High entropy (hopefully rare): store data as-is without compression.
            *reinterpret_cast<uint16_t*>(outPtr) |= split_map_4[max_entropy_bits][offset];
            offset += 4;
            if( offset > 7 ) {
                outPtr += (offset>>3);
                offset = offset&7;
            }

            outPtr = copy_block( diffPtr, thisBlockSize, outPtr, offset );
        } else if( split == -1 ) {              // All values identical (rare): store split=0 (out is already zeroed, so just increment ptr & offset)
            offset += 4;
            if( offset > 7 ) {
                outPtr += (offset>>3);
                offset = offset&7;
            }

        } else {                                // This is the normal case, i.e. Golomb-Rice compression of the data.

            // first store the split for this block
            *reinterpret_cast<uint16_t*>(outPtr) |= split_map_4[split][offset];
            offset += 4;
            if( offset > 7 ) {
                outPtr += (offset>>3);
                offset = offset&7;
            }

            // do the bit-manipulations and store results in out-array
            outPtr = crunch_block( diffPtr, thisBlockSize, outPtr, offset, split );
        }
        inPtr += thisBlockSize;
        refValue = inPtr[-1];    // reference value for next block is the last value of this block.

    }

    return(outPtr-out+(offset>0));
    
}


int redux::util::rice_decomp16( const uint8_t* in, size_t inSize, int16_t* out, size_t outSize, size_t blockSize ) {
    
    const int max_entropy_bits = 14;
  
    const uint8_t* __restrict__ inPtr = in;
    int16_t* __restrict__ outPtr = out;
    const uint8_t* __restrict__ inEnd = in+inSize;
    const int16_t* __restrict__ outEnd = out+outSize;
    

    uint32_t* __restrict__ diffPtr;

    std::unique_ptr<uint32_t[]> diff;
    try {
        diffPtr = new uint32_t[ blockSize ];
        diff.reset( diffPtr );
    } catch ( const std::bad_alloc& ) {
        printf( "rice_comp16: bad_alloc." );
        return(-1);
    }

    memset( outPtr, 0, outSize*sizeof(int16_t) );
    const uint32_t* diffEnd = diffPtr+blockSize;
    
    int32_t lastValue = inPtr[1] | (inPtr[0]<<8);
    
    inPtr += 2;

    int offset(0);
    
    size_t thisBlockSize = blockSize;
    while( inPtr < inEnd && outPtr < outEnd ) {
        if( outPtr+blockSize >= outEnd ) {
            thisBlockSize = outEnd-outPtr;
        }

        uint16_t split = ((inPtr[1] | (inPtr[0]<<8)) >> (12-offset)) & 0xF;

        offset += 4;
        inPtr += offset>>3;
        offset = offset&7;

        if( split == 0 ) {
            std::fill_n( outPtr, thisBlockSize, lastValue );
            outPtr += thisBlockSize;
            continue;
        } else if( split >= max_entropy_bits ) {
            diffPtr = diff.get();
            if( offset ) {
                const uint8_t sh = (8-offset);
                for ( size_t j=0; j<(thisBlockSize<<1); ) {
                    int nVals = 3;
                    if( diffPtr >= diffEnd ) break;
                    if( diffPtr+nVals > diffEnd ) nVals = diffEnd-diffPtr;
                    if( nVals <= 0 ) break;
                    upack_t tmpU( inPtr+j, 2*nVals+1, sh );
                    j += nVals*sizeof(int16_t);
                    while(nVals--) *diffPtr++ = tmpU.u16[nVals];
                }
            } else {
                const uint16_t* ptr16 = reinterpret_cast<const uint16_t*>(inPtr);
                std::transform( ptr16, ptr16+thisBlockSize, diffPtr, [](const uint16_t&a){
                    return be16toh(a);
                });
            }
            inPtr += thisBlockSize * sizeof(int16_t);
        } else {

            diffPtr = diff.get();
            split--;    // split+1 was stored.
            const uint32_t splitMask = (1<<split) - 1;

            for ( size_t j=0; j<thisBlockSize; j++ ) {

                uint8_t val = *inPtr & (0xFF>>offset);
                int leadingZeroes(0);
                int cnt(0);
                while( !val ) {
                    leadingZeroes += 8;
                    val = *(inPtr + ++cnt);
                }
                if( leadingZeroes ) {
                    leadingZeroes += nLeadingZeroes[val]-offset;
                    offset += leadingZeroes+1;
                } else {
                    leadingZeroes += nLeadingZeroes[val]-offset;
                    offset += leadingZeroes+1;
                }
                inPtr += offset>>3;

                offset = offset&7;
                const int lastBit = offset+split-1;
                const int nBytes = 1+(lastBit>>3);

                const uint8_t rshift = ((nBytes<<3)-lastBit-1);
                const upack_t tmpU( inPtr, nBytes, rshift );

                offset += split;

                *diffPtr++ = (leadingZeroes<<split) | (tmpU.u16[0]&splitMask);
                inPtr += offset>>3;

                offset = offset&7;

            }
        }

        // convert diff back to values.
        diffPtr = diff.get();

        for ( size_t j=0; j<thisBlockSize; j++ ) {
            int32_t d = *diffPtr++;
            if( d & 1 ) {
                d = ~(d>>1);
            } else {
                d >>= 1;
            }
            *outPtr = lastValue + d;
            lastValue = *outPtr++;
        }
    }
    
    return 0;
}

