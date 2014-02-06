#ifndef REDUX_FILE_ANADECOMPRESS_HPP
#define REDUX_FILE_ANADECOMPRESS_HPP

#include <stdint.h>

namespace redux {

    namespace file {

        int anadecrunch32( const unsigned char *x, int32_t *array, int r9, int nx, int ny, int );
        int anadecrunchrun8( const unsigned char *x, int8_t *array, int r9, int nx, int ny, int );
        int anadecrunchrun( const unsigned char *x, int16_t *array, int r9, int nx, int ny, int );
        int anadecrunch8( const unsigned char *x, int8_t *array, int r9, int nx, int ny, int );
        int anadecrunch( const unsigned char *x, int16_t *array, int r9, int nx, int ny, int );

    }
}

#endif  // REDUX_FILE_ANADECOMPRESS_HPP
