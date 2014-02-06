#ifndef REDUX_FILE_ANACOMPRESS_HPP
#define REDUX_FILE_ANACOMPRESS_HPP

#include <cstdint>

namespace redux {

    namespace file {

        uint32_t anacrunchrun8(uint8_t *x, const uint8_t *array, int slice, int nx, int ny, uint32_t limit, int t_endian);
        uint32_t anacrunch8(uint8_t *x, const uint8_t *array, int slice, int nx, int ny, uint32_t limit, int t_endian);
        uint32_t anacrunchrun(uint8_t *x, const int16_t *array, int slice, int nx, int ny, uint32_t limit, int t_endian);
        uint32_t anacrunch(uint8_t *x, const int16_t *array, int slice, int nx, int ny, uint32_t limit, int t_endian);
        uint32_t anacrunch32(uint8_t *x, const int32_t *array, int slice, int nx, int ny, uint32_t limit, int t_endian);

    }
}

#endif  // REDUX_FILE_ANACOMPRESS_HPP
