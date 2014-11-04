#include "redux/image/utils.hpp"

#include "redux/math/functions.hpp"

#include <map>

using namespace redux::image;
using namespace redux::util;
using namespace std;

namespace {

    const int maxDistance = 32;
    const double deltasqr = 4;
    const double beta = 2;

    float distanceCache[maxDistance*maxDistance];

    bool calculateCache( void ) {
        for( int i = 0; i < maxDistance; ++i ) {
            double i2 = i * i;
            for( int j = 0; j < maxDistance; ++j ) {
                double j2 = j * j;
                distanceCache[i * maxDistance + j] = pow( sqrt( i2 + j2 + deltasqr ), -beta );
            }
        }

        return true;
    }

    map<int, double> getDistanceMap( void ) {
        map<int, double> tmp;
        for( int i = 0; i <= maxDistance; ++i ) {
            int i2 = i * i;
            for( int j = 0; j <= maxDistance; ++j ) {
                int j2 = j * j;
                tmp.insert( pair<int, double>( i2 + j2, pow( i2 + j2 + deltasqr, -beta ) ) );
            }
        }
        return tmp;
    }

}


static __inline int sqr( int x ) {
    return x * x;
}

// nverseDistanceWeight( T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX ) {

//double inv_dist_wght(double **a,int xl,int xh,int yl,int yh,int xb,int yb)
double redux::image::inv_dist_wght( float **a, size_t sizeY, size_t sizeX, size_t posY, size_t posX ) {

    int xl = std::max( 0L, static_cast<int64_t>( posX - maxDistance ) );
    int xh = std::min( sizeX, posX + maxDistance + 1 );
    int yl = std::max( 0L, static_cast<int64_t>( posY - maxDistance ) );
    int yh = std::min( sizeY, posY + maxDistance + 1 );


    double weight = 0.0, res = 0.0;
    for( int y = yl; y < yh; ++y )
        for( int x = xl; x < xh; ++x )
            if( a[y][x] ) {
                double c = pow(  ( double )sqr( x - posX ) + ( double )sqr( y - posY ) + deltasqr, -beta );
                res += c * a[y][x];
                weight += c;
            }
    return res / weight;
}


template <typename T>
double redux::image::inverseDistanceWeight( T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX ) {

    static map<int, double> inverseDistanceSquared = getDistanceMap();
    
    // TODO: verify this function, results look weird

    size_t beginX = std::max( 0L, static_cast<int64_t>( posX - maxDistance ) );
    size_t endX = std::min( sizeX, posX + maxDistance + 1 );
    size_t beginY = std::max( 0L, static_cast<int64_t>( posY - maxDistance ) );
    size_t endY = std::min( sizeY, posY + maxDistance + 1 );

    double normalization = 0.0, weightedSum = 0.0;
    for( int y = beginY; y < endY; ++y ) {
        int y2 = ( y - posY ) * ( y - posY );
        for( int x = beginX; x < endX; ++x ) {
            int x2 = ( x - posX ) * ( x - posX );
            if( array[y][x] > 0 ) {
                double tmp = inverseDistanceSquared.at( y2 + x2 );
                weightedSum += tmp * array[y][x];
                normalization += tmp;
            }
        }
    }
    return weightedSum / normalization;
}
template double redux::image::inverseDistanceWeight( float**, size_t, size_t, size_t, size_t );
template double redux::image::inverseDistanceWeight( double**, size_t, size_t, size_t, size_t );

template <typename T>
double redux::image::horizontalInterpolation( T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX ) {

    T* ptr = array[posY];

    //map the 5 pixel surrounding as bits in a byte
    int val = 0;
    if( ( posX > 1 ) ) val |= ( ( ptr[posX - 2] > 0 ) << 4 );
    if( ( posX > 0 ) )   val |= ( ( ptr[posX - 1] > 0 ) << 3 );
    if( ( posX < sizeX ) )   val |= ( ( ptr[posX + 1] > 0 ) << 1 );
    if( ( posX < sizeX - 1 ) ) val |= ( ( ptr[posX + 2] > 0 ) );
    //now select based on the number
    switch( val ) {
        case( 10 ):     // = 0 1 x 1 0
        case( 11 ):     // = 0 1 x 1 1
        case( 26 ):     // = 1 1 x 1 0
        case( 27 ):     // = 1 1 x 1 1
            return ( ptr[posX - 1] + ptr[posX + 1] ) / 2;
        case( 18 ):     // = 1 0 x 1 0
        case( 19 ):     // = 1 0 x 1 1
            return ( ptr[posX - 2] + 2 * ptr[posX + 1] ) / 3;
        case( 9 ):      // = 0 1 x 0 1
        case( 25 ):     // = 1 1 x 0 1
            return ( 2 * ptr[posX - 1] + ptr[posX + 2] ) / 3;
        default:
            return inverseDistanceWeight<T>( array, sizeY, sizeX, posY, posX );
    }

}
template double redux::image::horizontalInterpolation( float**, size_t, size_t, size_t, size_t );
template double redux::image::horizontalInterpolation( double**, size_t, size_t, size_t, size_t );


template <typename T>
void redux::image::apodize( Array<T>& array, size_t blendRegion ) {
    if(!blendRegion) return;        // nothing to do
    const vector<int64_t>& first = array.first();
    size_t sizeY = array.last()[0]-first[0]+1;
    size_t sizeX = array.last()[1]-first[1]+1;
    blendRegion = std::min(std::min(blendRegion,sizeY),sizeX);
    
    T* tmp = new T[blendRegion+1];
    tmp[0] = 0;
    redux::math::apodize( tmp, blendRegion, T(1) );
//     auto sharedPtrs = array.get(sizeY,sizeX);
//     T** data = sharedPtrs.get();
    for(size_t y=0,yy=sizeY-1; y<sizeY; ++y,--yy ) {
        double yfactor = 1;
        if( y < blendRegion ) yfactor *= tmp[y];
        for(size_t x=0,xx=sizeX-1; x<sizeX; ++x,--xx ) {
            double xfactor = yfactor;
            if( x < blendRegion ) xfactor *= tmp[x];
            if(xfactor < 1) {
                array(y,x)   *= xfactor;
                array(yy,x)  *= xfactor;
                array(y,xx)  *= xfactor;
                array(yy,xx) *= xfactor;
            }
        }
    }
    delete[] tmp;
}
template void redux::image::apodize( Array<int16_t>&, size_t );
template void redux::image::apodize( Array<int32_t>&, size_t );
template void redux::image::apodize( Array<double>&, size_t );
template void redux::image::apodize( Array<float>&, size_t );

