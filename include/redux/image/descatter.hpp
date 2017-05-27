#ifndef REDUX_IMAGE_DESCATTER_HPP
#define REDUX_IMAGE_DESCATTER_HPP


#include <redux/types.hpp>
#include <redux/util/array.hpp>

namespace redux {

    namespace image {
        
        /*! Fit and subtract the backscatter patter for semi-transparent CCDs 
         *  @param img Image to be cleaned.
         *  @param nRows Number of pixel rows.
         *  @param nCols Number of pixel columns.
         *  @param nPaddedRows Number of pixel rows, including padding.
         *  @param nPaddedCols Number of pixel columns, including padding.
         *  @param gain Image containing the backscatter pattern.
         *  @param psf Approximated PSF for the backscattering of this CCD.
         *  @param otf OTF for the backscattering of this CCD.
         *  @param nthreads How many threads to use in the multithreaded implementation.
         *  @param maxIterations Exit condition.
         *  @param minImprovement Exit if the improvement in one step is less than this factor (1 means no improvement).
         *  @param epsilon Exit if the metric reaches this value.
         *  @{
         */
        void descatter( double* img, size_t nRows, size_t nCols, double* gain, fftw_complex* otf, size_t nPaddedRows, size_t nPaddedCols, int nthreads, int maxIterations, float epsilon );
        void descatter( double* img, size_t nRows, size_t nCols, double* gain, double* psf, size_t nPaddedRows, size_t nPaddedCols, int nthreads, int maxIterations, float epsilon );

        template <typename T, typename U>
        void descatter(redux::util::Array<T>& img, const redux::util::Array<U>& gain, const redux::util::Array<U>& psf,
                       int maxIterations=50, double minImprovement=1, double epsilon=1E-10 );
        
        /*! @} */

    }   // image
    
}   // redux


#endif  // REDUX_IMAGE_DESCATTER_HPP
