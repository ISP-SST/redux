#include "redux/image/statistics.hpp"

#include "redux/file/fileana.hpp"
#include "redux/image/utils.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/util/bitoperations.hpp"

#include <atomic>
#include <limits>

using namespace redux::image;
using namespace redux::util;
using namespace std;

template <typename T>
void redux::image::Statistics<T>::getStats( const Array<T>& data, int flags ) {

    min = std::numeric_limits<T>::max();
    max = std::numeric_limits<T>::min();
    double sum = 0;
    size_t nElements = data.nElements();
    for( auto & it : data ) {
        if( it > max ) max = it;
        if( it < min ) min = it;
        sum += static_cast<double>( it );
    }
    mean = sum / nElements;
    
    // TODO: median
    
    if( flags & ST_RMS ) {
        rms = stddev = 0;
        for( auto & it : data ) {
            rms += static_cast<double>( it * it );
            double tmp = ( static_cast<double>( it ) - mean );
            stddev += tmp * tmp;
        }
        rms = sqrt( rms / nElements );
        stddev = sqrt( stddev / nElements );
    }
    
    if( flags & ST_NOISE ) {
        Array<T> tmpImage = data.copy();    // make a deep copy because "apodize" is destructive.
        const std::vector<size_t>& dims = tmpImage.dimensions();
        std::vector<int64_t> first = tmpImage.first();
        std::vector<int64_t> last = tmpImage.last();
        size_t nDims = dims.size();
        for( size_t i = 0; i < nDims; ++i ) {   // take a centered subimage that has "next smaller power of two" size
            uint32_t dimSize = tmpImage.dimSize( i );
            uint32_t newSize = nextPowerOfTwo( dimSize ) >> 1;
            if( newSize < dimSize ) {
                int diff = dimSize - newSize;
                first[i] += diff / 2;
                last[i] -= ( diff - diff / 2 );
            }
        }
        tmpImage.setLimits( first, last );
        apodize( tmpImage, 8 );                     // edge smoothing to reduce the FFT boundary-artifacts

        // normalize
        sum = 0;
        for( auto it : tmpImage ) sum += it;
        sum /= tmpImage.nElements();
        tmpImage /= sum;

        FourierTransform ft( tmpImage );

        size_t Ny = ft.dimSize( 0 );
        size_t Nx = ft.dimSize( 1 );
        size_t realNx = tmpImage.dimSize( 1 );  // NB: we are using real-to-complex transform from fftw3, so the transform has dimensions (Ny x Nx/2+1)
        sum = 0;
        size_t count = 0;
        double yxRatio = static_cast<double>( Ny ) / realNx;
        double Nxy4 = realNx * Ny / 4.0;
        for( size_t x = realNx / 6; x < Nx; ++x ) {         // skip the Nx/6 pixels nearest the edge (fft artifacts)
            double xx = x * x * yxRatio;
            for( size_t y = Ny / 6; y < Ny / 2; ++y ) {     // skip the Ny/6 pixels nearest the edge (fft artifacts)
                if( ( xx + y * y ) > Nxy4 ) {               // skip a region corresponding to the cutoff-frequency (diffraction limit)  TODO: verify the numbers used here
                    sum += std::norm( ft( y, x ) );         // norm = abs^2
                    sum += std::norm( ft( Ny-y-1, x ) );
                    count += 2;
                }
            }
        }
        noisePower = sqrt( sum / ( count * Ny * Nx ) );

    }
    //std::cout << "stats0 (" << myID << "):     min = " << min << "  max = " << max << "  mean = " << mean <<  "  rms = " << rms <<  "  stddev = " << stddev << "  pwr = " << noisePower << std::endl;

}
template void redux::image::Statistics<float>::getStats( const Array<float>&, int );
template void redux::image::Statistics<int16_t>::getStats( const Array<int16_t>&, int );

template <typename T>
void redux::image::Statistics<T>::getStats( uint32_t borderClip, const Array<T>& data, int flags ) {

    if( !borderClip ) {
        getStats( data, flags );
        return;
    }

    std::vector<int64_t> first, last;
    const std::vector<size_t>& dims = data.dimensions();
    size_t nDims = dims.size();
    for( size_t i = 0; i < nDims; ++i ) {
        if( dims[i] == 1 ) {            // if dimSize == 1 we don't clip it.
            first.push_back( data.first()[i] );
            last.push_back( data.last()[i] );
        }
        else {
            if( 2 * borderClip > dims[i] ) {    // all data clipped, nothing to do
                return;
            }
            first.push_back( data.first()[i] + borderClip );
            last.push_back( data.last()[i] - borderClip );
        }
    }
    Array<T> clippedImage( data, first, last );
    getStats(clippedImage,flags);
    
}
template void redux::image::Statistics<float>::getStats( uint32_t, const Array<float>&, int );
template void redux::image::Statistics<double>::getStats( uint32_t, const Array<double>&, int );

