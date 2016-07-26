#ifndef REDUX_IMAGE_DESCATTER_HPP
#define REDUX_IMAGE_DESCATTER_HPP


#include <redux/types.hpp>
#include <redux/util/array.hpp>

namespace redux {

    namespace image {
        
        /*! Fit and subtract the backscatter patter for semi-transparent CCDs 
         *  @param data Image to be cleaned
         *  @param ccdGain Image containing the backscatter pattern
         *  @param psf PSF for the system
         *  @param maxIterations Exit condition
         *  @param minImprovement Exit if the improvement in one step is less than this factor (1 means no improvement)
         *  @param epsilon Exit if the metric reaches this value.
         */
        void descatter( double* img, size_t imgY, size_t imgX, double* gain, fftw_complex* otf, size_t paddedY, size_t paddedX, int nthreads, int niter, float epsilon );
        void descatter( double* img, size_t szY, size_t szX, double* gain, double* psf, size_t paddedY, size_t paddedX, int nthreads, int niter, float epsilon );

        template <typename T, typename U>
        void descatter(redux::util::Array<T>& data, const redux::util::Array<U>& ccdgain, const redux::util::Array<U>& psf,
                       int maxIterations=50, double minImprovement=1, double epsilon=1E-10 );

    }   // image
    
}   // redux


#endif  // REDUX_IMAGE_DESCATTER_HPP
