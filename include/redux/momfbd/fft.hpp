#ifndef REDUX_MOMFBD_FFT_HPP
#define REDUX_MOMFBD_FFT_HPP

#include "redux/types.hpp"
namespace redux {
void fft(redux::complex_t*,int,int);       // 1D
void fft_init(int,int);
void fft_2d(redux::complex_t**,int,int,int);  // 2D
void fft_done(int,int);
void fft_reorder(redux::complex_t **f,int np);
}
#endif                // REDUX_MOMFBD_FFT_HPP
