#include "redux/image/statistics.hpp"

#include "redux/file/fileana.hpp"
#include "redux/image/utils.hpp"
#include "redux/util/bitoperations.hpp"

#include <limits>

using namespace redux::image;
using namespace redux::util;
using namespace std;

template <typename T>
void redux::image::Statistics::getMinMaxMean(const Array<T>& data) {
    
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::lowest();
    sum = 0;
    size_t nElements = data.nElements();
    for( auto & it : data ) {
        if( it > max ) max = it;
        if( it < min ) min = it;
        sum += static_cast<double>( it );
    }
    mean = sum / nElements;
    
}
template void redux::image::Statistics::getMinMaxMean( const Array<float>& );
template void redux::image::Statistics::getMinMaxMean( const Array<double>& );
template void redux::image::Statistics::getMinMaxMean( const Array<int16_t>& );
template void redux::image::Statistics::getMinMaxMean( const Array<int32_t>& );

           
template <typename T>
void redux::image::Statistics::getRmsStddev(const Array<T>& data, double mean) {
    size_t nElements = data.nElements();
    if(nElements > 1) {
        rms = stddev = 0;
        for( auto & it : data ) {
            rms += static_cast<double>( it * it );
            double tmp = ( static_cast<double>( it ) - mean );
            stddev += tmp * tmp;
        }
        rms = sqrt( rms / nElements );
        stddev = sqrt( stddev / (nElements-1) );
    } else if (nElements == 1) {
        rms = data(0);
        stddev = 0;
    }
}


void redux::image::Statistics::getNoise(const redux::image::FourierTransform& ft) {
    noise = ft.noise(clip, cutoff);
}

       
template <typename T>
void redux::image::Statistics::getNoise(const Array<T>& data) {
    Array<T> tmpImage = data.copy();            // make a deep copy because "apodize" is destructive.
    apodize( tmpImage, 8 );                     // edge smoothing to reduce the FFT boundary-artifacts
    FourierTransform ft( tmpImage, FT_NORMALIZE );
    getNoise(ft);
}


template <typename T>
void redux::image::Statistics::getStats( const Array<T>& data, int flags ) {

    getMinMaxMean(data);
    
    // TODO: median
    
    if( flags & ST_RMS ) {
        getRmsStddev(data,mean);
    }
    
    if( flags & ST_NOISE ) {
        getNoise(data);
    }

}
template void redux::image::Statistics::getStats( const Array<float>&, int );
template void redux::image::Statistics::getStats( const Array<double>&, int );
template void redux::image::Statistics::getStats( const Array<int16_t>&, int );
template void redux::image::Statistics::getStats( const Array<int32_t>&, int );



template <typename T>
void redux::image::Statistics::getStats( uint32_t borderClip, const Array<T>& data, int flags ) {

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
template void redux::image::Statistics::getStats( uint32_t, const Array<float>&, int );
template void redux::image::Statistics::getStats( uint32_t, const Array<double>&, int );
template void redux::image::Statistics::getStats( uint32_t, const Array<int16_t>&, int );
template void redux::image::Statistics::getStats( uint32_t, const Array<int32_t>&, int );

/*
            int clip;
            double cutoff;
            double min, max, median;
            double sum, mean, rms, stddev;
            double noise, noiseRMS;
            uint8_t noiseType;              // flag indicating noise statistics (not used atm.)
*/
size_t redux::image::Statistics::size( void ) {
    return 1 + sizeof(int) + 10*sizeof(double);
}



uint64_t redux::image::Statistics::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack(ptr, clip);
    count += pack(ptr+count, cutoff);
    count += pack(ptr+count, min);
    count += pack(ptr+count, max);
    count += pack(ptr+count, median);
    count += pack(ptr+count, sum);
    count += pack(ptr+count, mean);
    count += pack(ptr+count, rms);
    count += pack(ptr+count, stddev);
    count += pack(ptr+count, noise);
    count += pack(ptr+count, noiseRMS);
    count += pack(ptr+count, noiseType);
    return count;
}


uint64_t redux::image::Statistics::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = unpack(ptr, clip, swap_endian);
    count += unpack(ptr+count, cutoff, swap_endian);
    count += unpack(ptr+count, min, swap_endian);
    count += unpack(ptr+count, max, swap_endian);
    count += unpack(ptr+count, median, swap_endian);
    count += unpack(ptr+count, sum, swap_endian);
    count += unpack(ptr+count, mean, swap_endian);
    count += unpack(ptr+count, rms, swap_endian);
    count += unpack(ptr+count, stddev, swap_endian);
    count += unpack(ptr+count, noise, swap_endian);
    count += unpack(ptr+count, noiseRMS, swap_endian);
    count += unpack(ptr+count, noiseType, swap_endian);
    return count;

}
