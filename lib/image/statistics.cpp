#include "redux/image/statistics.hpp"

#include "redux/file/fileana.hpp"
#include "redux/image/utils.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/util/bitoperations.hpp"

#include <limits>

using namespace redux::image;
using namespace redux::util;

template <typename T>
void redux::image::Statistics<T>::getStats( const Array<T>& data, int flags ) {

    min = std::numeric_limits<T>::max();
    max = std::numeric_limits<T>::min();
    double sum = 0;
    for( auto & it : data ) {
        if( it > max ) max = it;
        if( it < min ) min = it;
        sum += static_cast<double>( it );
    }
    mean = sum / data.nElements();

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
        if( dims[i] == 1 ) {    // if dimSize == 1 we don't clip it.
            first.push_back( data.first()[i] );
            last.push_back( data.last()[i] );
        }
        else {
            if( 2 * borderClip > dims[i] ) {
                return;
            }
            first.push_back( data.first()[i] + borderClip );
            last.push_back( data.last()[i] - borderClip );
        }
    }

    Array<T> clippedImage( data, first, last );
    min = std::numeric_limits<T>::max();
    max = std::numeric_limits<T>::min();
    double sum = 0;
    for( auto & it : clippedImage ) {
        if( it > max ) max = it;
        if( it < min ) min = it;
        sum += static_cast<double>( it );
    }
    // TODO: median
    mean = sum / clippedImage.nElements();

    if( flags & ST_RMS ) {
        rms = stddev = 0;
        for( auto & it : clippedImage ) {
            rms += static_cast<double>(it*it);
            double tmp = ( static_cast<double>( it ) - mean );
            stddev += tmp*tmp;
        }
        rms = sqrt( rms / clippedImage.nElements() );
        stddev = sqrt( stddev / clippedImage.nElements() );
    }
    // FIXME: rms gives wrong value ??!!
    std::cout << "stats2b:   borderClip = " << borderClip << "  min = " << min << "  max = " << max << "  mean = " << mean <<  "  rms = " << rms <<  "  stddev = " << stddev << std::endl;
    static size_t fcnt=0;
    redux::file::Ana::write( "statimage_"+std::to_string(++fcnt)+".f0", clippedImage );

    // TODO verify noise statistics, compare with michiels version
    if( flags & ST_NOISE ) {
        for( size_t i = 0; i < nDims; ++i ) {
           uint32_t dimSize = clippedImage.dimSize(i);
           int diff = dimSize - nextPowerOfTwo(dimSize)/2;
           if( diff && (dims[i] != 1) ) {    // dimensions of size 1 will be mangled in the .copy() below.
                first[i] += diff/2;
                last[i] -= (diff - diff/2);
            }
        }

        clippedImage = data;                // restore the input dimensions (don't use resetLimits since the input data might be a sub-array, so resetting would go back to its "parent")
        clippedImage.setLimits(first,last);
        Array<T> tmpImage = clippedImage.copy();    // make a deep copy
        
        double sum = 0;
        for( auto it: tmpImage ) sum += it;
        sum /= tmpImage.nElements();
        tmpImage /= sum;
        apodize( tmpImage, 8 );         // edge smoothing to reduce the FFT artifacts near the boundary
        
        FourierTransform ft(tmpImage, FT_REORDER);
        
        int Ny = ft.dimSize(0);
        int Nx = ft.dimSize(1);
        double noisePower = 0;
        int N=0, yy, xx;
        double xyRatio = static_cast<double>(Nx)/Ny;
        double Nxy4 = Nx*Ny/4.0;
        for( int y=0; y < Ny; ++y ) {
            if( 6*abs( yy = y-Ny/2 ) > Ny ) {
                double yyyy = yy*yy*xyRatio;
                for( int x=0; x < Nx; ++x ) {
                    if( 6*abs( xx = x - Nx / 2 ) > Nx ) {
                        if( ( yyyy + xx*xx ) > Nxy4 ) {
                            //std::complex<double> elem = ft(y,x);
                            //noisePower += (elem*elem).real();
                            noisePower += std::abs(ft(y,x));
                            ++N;
                        }
                    }
                }
            }
        }

        noisePower = sqrt( noisePower / (N*Ny*Nx) );

 
        std::cout << "stats3b:  noisePower = " << noisePower << std::endl;

    }

}
template void redux::image::Statistics<float>::getStats( uint32_t, const Array<float>&, int );
template void redux::image::Statistics<double>::getStats( uint32_t, const Array<double>&, int );
