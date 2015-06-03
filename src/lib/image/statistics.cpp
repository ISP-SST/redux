#include "redux/image/statistics.hpp"

#include "redux/file/fileana.hpp"
#include "redux/image/utils.hpp"
#include "redux/util/bitoperations.hpp"

#include <limits>

using namespace redux::image;
using namespace redux::util;
using namespace std;

template <typename T>
void redux::image::Statistics::getMinMaxMean(const T* data, size_t n) {
    
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::lowest();
    sum = 0;
    const T* ptr = data + n;
    size_t count(0);
    while( ptr-- > data ) {
        double tmp = static_cast<double>( *ptr );
        if( tmp > max ) {
            max = tmp;
        } else if( tmp < min ) {
            min = tmp;
        }
        if ( tmp == tmp ) { // will exclude NaN
            sum += tmp;
            count++;
        }
    }
    mean = sum;
    if(count) mean /= count;
    
}
template void redux::image::Statistics::getMinMaxMean( const float*, size_t );
template void redux::image::Statistics::getMinMaxMean( const double*, size_t );
template void redux::image::Statistics::getMinMaxMean( const int16_t*, size_t );
template void redux::image::Statistics::getMinMaxMean( const int32_t*, size_t );

           
template <typename T>
void redux::image::Statistics::getRmsStddev(const T* data, size_t n) {
    if(n > 1) {
        size_t count(0);
        rms = stddev = 0;
        const T* ptr = data + n;
        while( ptr-- > data ) {
            long double tmp = static_cast<long double>( *ptr );
            if( tmp == tmp ) {
                rms += static_cast<long double>( tmp * tmp );
                tmp = (tmp - mean);
                stddev += tmp * tmp;
                count++;
            }
        }
        if (count) {
            rms = sqrt( rms / count );
            if (count > 1) stddev = sqrt( stddev / (count-1) );
        }
    } else if (n == 1) {
        rms = (*data == *data)?*data:0;
        stddev = 0;
    }
}
template void redux::image::Statistics::getRmsStddev( const float*, size_t );
template void redux::image::Statistics::getRmsStddev( const double*, size_t );
template void redux::image::Statistics::getRmsStddev( const int16_t*, size_t );
template void redux::image::Statistics::getRmsStddev( const int32_t*, size_t );


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
        getRmsStddev(data);
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
void redux::image::Statistics::getStats( const T* data, size_t count, int flags ) {

    getMinMaxMean(data, count);
    
    // TODO: median
    
    if( flags & ST_RMS ) {
        getRmsStddev(data,count);
    }
/*    
    if( flags & ST_NOISE ) {
        getNoise(data,count);
    }*/

}
template void redux::image::Statistics::getStats( const float*, size_t, int );
template void redux::image::Statistics::getStats( const double*, size_t, int );
template void redux::image::Statistics::getStats( const int16_t*, size_t, int );
template void redux::image::Statistics::getStats( const int32_t*, size_t, int );



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
