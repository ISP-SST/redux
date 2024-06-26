#include "redux/image/utils.hpp"

#include "redux/image/fouriertransform.hpp"
#include "redux/image/grid.hpp"
#include "redux/image/pupil.hpp"
#include "redux/image/zernike.hpp"

#include "redux/file/fileana.hpp"
#include "redux/math/functions.hpp"
#include "redux/constants.hpp"

#include <functional>
#include <map>
#include <set>
#include <math.h>

#include <gsl/gsl_multifit.h>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#ifdef RDX_WITH_OPENCV
//#       include <opencv2/photo/photo.hpp>
#endif

using namespace redux::image;
//using namespace redux::momfbd;
using namespace redux::util;
using namespace redux;
using namespace std;

namespace {

    const int maxDistance = 64;
    const int maxDistance2 = maxDistance*maxDistance;
    const double deltasqr = 4;
    const double beta = 2;
    
    double distMap[2*maxDistance2+1];
    
    const double* getDistanceMap (void) {
        memset(distMap,0,(2*maxDistance2+1)*sizeof(double));
        for(int i=0; i<=2*maxDistance2; ++i) distMap[i] = pow (i + deltasqr, -beta);
        return distMap;
    }


    template <typename T>
    void maskConnected(T** data, uint8_t** mask, size_t sizeY, size_t sizeX, unsigned int yy, unsigned int xx, unsigned int y_start, unsigned int x_start, T threshold) {
        if (yy >= sizeY || xx >= sizeX) return;
        if ( !mask[yy][xx] && ((yy == y_start && xx == x_start) || data[yy][xx] > threshold)) {       // new point inserted, check neighbours
            mask[yy][xx] = 1;
            maskConnected(data, mask, sizeY, sizeX, yy + 1, xx, y_start, x_start, threshold);
            maskConnected(data, mask, sizeY, sizeX, yy - 1, xx, y_start, x_start, threshold);
            maskConnected(data, mask, sizeY, sizeX, yy, xx + 1, y_start, x_start, threshold);
            maskConnected(data, mask, sizeY, sizeX, yy, xx - 1, y_start, x_start, threshold);
        }
    }

}

template <typename T>
static __inline T sqr (T x) {
    return x * x;
}


double redux::image::makePupil( double** pupil, uint32_t pupilPixels, PointF center, double outerRadius, double innerRadius ) {

    long double area = 0.0;
    size_t pupilPixels2 = pupilPixels * pupilPixels;
    memset( *pupil, 0, pupilPixels2 *sizeof(double) );

    struct Grid grid;
    grid.id.size = pupilPixels;
    grid.id.origin = center;
    grid.init();
    
    auto dist2D = grid.dist2D();
    auto angle2D = grid.angle2D();
    double** distPtr = dist2D.get();      // distance(i,j) is the distance from the centre of the pupil, to the center of pixel (i,j)
    double** anglePtr = angle2D.get();
    bool hasInner = (innerRadius > 0.0);

    for( unsigned int x=0; x< pupilPixels; ++x ) {
        for( unsigned int y=0; y< pupilPixels; ++y ) {
            const long double D = distPtr[y][x];
            const long double absangle = fabs(anglePtr[y][x]);
            long double d;
            if( absangle > M_PI/4 && absangle < 3*M_PI/4 ) {
                d = fabsl(1.0/sinl(absangle));
            } else {
                d = fabsl(1.0/cosl(absangle));
            }
            long double dir = (D-innerRadius);
            if( hasInner && (dir < -d) ) {                  // Entire pixel inside innerRadius
                continue;
            }
            long double dor = (D-outerRadius);
            if( dor > d ) {                                 // Entire pixel outside outerRadius
                continue;
            }
            if( (dor < -d) && (!hasInner || (dir > d)) ) {   // Entire pixel inside outerRadius and outside innerRadius 
                pupil[y][x] = 1;
                area += 1;
                continue;
            }
            long double R = dor;
            d *= 2;
            if( hasInner && (dor < -d) ) {    
                R = dir;
                d = -d;
            }
            long double val = 0.5 - R/d;
            pupil[y][x] = val;
            area += val;

        }
    }
    
    return static_cast<double>(area);


}


template <typename T>
void redux::image::makePlane( T* ptr, size_t ySize, size_t xSize, double* coeffs, uint8_t nC ) {
    if( nC < 1 ) return;
    PointI offset(0);
    double a = coeffs[0];
    double b = (nC>1)?coeffs[1]:0.0;
    double c = 0;
    if( nC > 2 ) {
        c = coeffs[2];
        offset.y = ySize/2;
        offset.x = xSize/2;
    }
    for( size_t y(0); y < ySize; ++y ) {
        double dy = static_cast<double>(y) - offset.y;
        for (size_t x(0); x < xSize; ++x ) {
            double dx = static_cast<double>(x) - offset.x;
            size_t offset = y*xSize+x;
            ptr[offset] = a*dy + b*dx + c;
        }
    }
}


template <typename T>
void redux::image::makePlane( Array<T>& arr, double* coeffs, uint8_t nC ) {
    makePlane( arr.ptr(), arr.dimSize(0), arr.dimSize(1), coeffs, nC );
}
template void redux::image::makePlane(  Array<float>&, double*, uint8_t );
template void redux::image::makePlane(  Array<double>&, double*, uint8_t );


template <typename T>
Array<T> redux::image::makePlane( size_t ySize, size_t xSize, double* coeffs, uint8_t nC ) {
    Array<T> ret( ySize, xSize );
    makePlane( ret.ptr(), ySize, xSize, coeffs, nC );
    return ret;
}
template Array<float> redux::image::makePlane<float>( size_t, size_t, double*, uint8_t nC );
template Array<double> redux::image::makePlane<double>( size_t, size_t, double*, uint8_t nC );

    
template <typename T, typename U>
void redux::image::fitPlane( const T* inPtr, size_t ySize, size_t xSize, const U* maskPtr, double* coeffs, double* chisq ) {
    
    static const U mz(0);
    size_t nPixels = ySize * xSize;
    size_t nSamples = nPixels;
    size_t nParams = 3;                                                        //   fit a plane as:   z = a*x + b*y + c

    if( maskPtr ) {
        nSamples = std::count_if( maskPtr, maskPtr+nPixels, []( U m ){ return ( m > mz ); } );
    }
    
    gsl_vector *data = gsl_vector_alloc( nSamples );
    gsl_vector *coeff = gsl_vector_alloc( nParams );
    gsl_matrix *X = gsl_matrix_alloc( nSamples, nParams );
    gsl_matrix *covar = gsl_matrix_alloc( nParams, nParams );

    size_t yHalf = ySize/2;
    size_t xHalf = xSize/2;
    if( maskPtr ) {
        size_t mCount(0);
        for( size_t i(0); i<nPixels; ++i ) {
            if( maskPtr[i] > mz ) {
                double z_x = static_cast<double>(i%xSize) - xHalf;
                double z_y = static_cast<double>(i/xSize) - yHalf;
                gsl_matrix_set( X, mCount, 0, z_y );
                gsl_matrix_set( X, mCount, 1, z_x );
                gsl_matrix_set( X, mCount, 2, 1 );
                gsl_vector_set( data, mCount, inPtr[i] );
                mCount++;
            }
        }
    } else {
        for( size_t y(0); y <ySize; ++y ) {
            double z_y = static_cast<double>(y) - yHalf;
            for( size_t x(0); x <xSize; ++x ) {
                double z_x = static_cast<double>(x) - xHalf;
                size_t offset = y*xSize + x;
                gsl_matrix_set( X, offset, 0, z_y );
                gsl_matrix_set( X, offset, 1, z_x );
                gsl_matrix_set( X, offset, 2, 1 );
                gsl_vector_set( data, offset, inPtr[offset] );
            }
        }
    }
    
    double this_chisq;
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc( nSamples, nParams );
    gsl_multifit_linear( X, data, coeff, covar, &this_chisq, work );
    gsl_multifit_linear_free( work );

    coeffs[0] = gsl_vector_get( coeff, 0 );
    coeffs[1] = gsl_vector_get( coeff, 1 );
    coeffs[2] = gsl_vector_get( coeff, 2 );

    if( chisq ){
        *chisq = this_chisq/nSamples;
    }
    
    gsl_vector_free( data );
    gsl_vector_free( coeff );
    gsl_matrix_free( X );
    gsl_matrix_free( covar );

}


template <typename T, typename U>
Array<T> redux::image::fitPlane( const Array<T>& in, const Array<U>& mask, bool subtract_mean, double* coeffs, double* chisq ) {

    if( !mask.sameSize(in) ) {
        throw runtime_error("redux::image::fitPlane: image & mask must be the same size" );
    }
    
    size_t ySize = in.dimSize(0);
    size_t xSize = in.dimSize(1);
    
    double tmp_coeffs[4];
    if( !coeffs ) coeffs = &tmp_coeffs[0];
    if( !chisq ) chisq = &tmp_coeffs[3];

    fitPlane( in.ptr(), ySize, xSize, mask.ptr(), coeffs, chisq );

    if( subtract_mean ) {
        coeffs[2] = 0.0;
    }

    Array<T> ret( ySize, xSize );
    makePlane( ret.ptr(), ySize, xSize, coeffs, 3 );
    return ret;
    
}
template Array<float> redux::image::fitPlane( const Array<float>&, const Array<bool>&, bool, double*, double* );
template Array<float> redux::image::fitPlane( const Array<float>&, const Array<uint8_t>&, bool, double*, double* );
template Array<float> redux::image::fitPlane( const Array<float>&, const Array<float>&, bool, double*, double* );
template Array<float> redux::image::fitPlane( const Array<float>&, const Array<double>&, bool, double*, double* );
template Array<double> redux::image::fitPlane( const Array<double>&, const Array<bool>&, bool, double*, double* );
template Array<double> redux::image::fitPlane( const Array<double>&, const Array<uint8_t>&, bool, double*, double* );
template Array<double> redux::image::fitPlane( const Array<double>&, const Array<float>&, bool, double*, double* );
template Array<double> redux::image::fitPlane( const Array<double>&, const Array<double>&, bool, double*, double* );


template <typename T>
Array<T> redux::image::fitPlane( const Array<T>& in, bool subtract_mean, double* coeffs, double* chisq ) {

    size_t ySize = in.dimSize(0);
    size_t xSize = in.dimSize(1);
    
    double tmp_coeffs[4];
    if( !coeffs ) coeffs = &tmp_coeffs[0];
    if( !chisq ) chisq = &tmp_coeffs[3];

    fitPlane( in.ptr(), ySize, xSize, (bool*)nullptr, coeffs, chisq );

    if( subtract_mean ) {
        coeffs[2] = 0.0;
    }

    Array<T> ret( ySize, xSize );
    makePlane( ret.ptr(), ySize, xSize, coeffs, 3 );
    return ret;
    
}
template Array<float> redux::image::fitPlane( const Array<float>&, bool, double*, double* );
template Array<double> redux::image::fitPlane( const Array<double>&, bool, double*, double* );


template <typename T>
void redux::image::connectedRegion(T** data, uint8_t** mask, size_t sizeY, size_t sizeX, unsigned int y_start, unsigned int x_start, T threshold) {
    memset(*mask, 0, sizeY * sizeX);
    maskConnected(data, mask, sizeY, sizeX, y_start, x_start, y_start, x_start, threshold);
}
template void redux::image::connectedRegion(float**, uint8_t**, size_t, size_t, uint, uint, float);
template void redux::image::connectedRegion(double**, uint8_t**, size_t, size_t, uint, uint, double);


template <typename T>
void redux::image::smooth(T** data, size_t sizeY, size_t sizeX, size_t nY, size_t nX) {
    if (nY == 0 || nX == 0) return;
    T** tmp = newArray<T> (sizeY, sizeX);
    size_t n = sizeY*sizeX;
    std::fill_n( *tmp, n, T(0) );
    for (unsigned int y = 0; y < sizeY; ++y) {
        int yl = std::max<int>(y-nY, 0);
        int yh = std::min(y+nY, sizeY);
        for (unsigned int x = 0; x < sizeX; ++x) {
            int xl = std::max<int>(x-nX, 0);
            int xh = std::min(x+nX, sizeX);
            int cnt = 0;
            for (int yy=yl; yy < yh; ++yy) {
                for (int xx=xl; xx < xh; ++xx) {
                    tmp[y][x] += data[yy][xx];
                    cnt++;
                }
            }
            tmp[y][x] /= cnt;
        }
    }
    std::copy_n( *tmp, n, *data );
    delArray (tmp);
}
template void redux::image::smooth (float**, size_t, size_t, size_t, size_t);
template void redux::image::smooth (double**, size_t, size_t, size_t, size_t);
template void redux::image::smooth (complex_t**, size_t, size_t, size_t, size_t);


template <typename T>
void redux::image::ScharmerFilter (T** data, double** q2_inv, size_t sizeY, size_t sizeX, double noise_power, double frequency_cutoff) {

    static const double hi = 1.0;
    static const double lo = 0.2;
    double noiseFactor = noise_power*sizeY*sizeX;
    
    double** pt = newArray<double> (sizeY, sizeX);
    double** rr = newArray<double> (sizeY, sizeX);
    
    std::transform (*data, *data + sizeY * sizeX, *q2_inv, *pt, [](const T& a, const double& b ){ return norm(a)*b; } );
    
    smooth (pt, sizeY, sizeX, 1, 1);
    
    std::transform (*pt, *pt + sizeY * sizeX, *pt, [noiseFactor](const double& a){ return noiseFactor/a;} );
    
    for (unsigned int y=0; y<sizeY; ++y) {
        rr[y][0] = pt[y][0];
        for (unsigned int x=1; x<sizeX; ++x) {
            rr[y][x] = pt[y][sizeX-x];
        }
    }
    
    double* tmp = newArray<double> (sizeY);
    for(unsigned int x=0; x<sizeX; ++x){                        // rr==filter
        rr[0][x] = max((1.0-0.5*(rr[0][x]+pt[0][x])),0.0);
        for(unsigned int y=1; y<sizeY;++y) tmp[y]=rr[y][x];       // temporary storage
        for(unsigned int y=1; y<sizeY;++y) rr[y][x]=max((1.0-0.5*(tmp[sizeY-y]+pt[y][x])), 0.0);
    }
    delArray (tmp);
    
    for(unsigned int y=0; y<sizeY; ++y) {
        for(unsigned int x=0; x<sizeX; ++x){
            if(rr[y][x]<lo) rr[y][x] = 0.0;
            if(rr[y][x]>hi) rr[y][x] = 1.0;
        }
    }
    
    unsigned int yHalf = sizeY/2;
    unsigned int xHalf = sizeX/2;
    rr[yHalf][xHalf] = 1.0;                   // ; DC gain = 1
    
    connectedRegion(rr, sizeY, sizeX, yHalf, xHalf, 0.0);
    smooth(rr, sizeY, sizeX, 4, 4);
    
    frequency_cutoff *= frequency_cutoff;
    for(unsigned int y=0; y<sizeY; ++y) {
        for(unsigned int x=0; x<sizeX; ++x) {
            if( (sqr(y-yHalf)+sqr(x-xHalf)) > frequency_cutoff){
                data[y][x] = 0.0;
            } else{
                data[y][x] *= rr[y][x];
            }
        }
    }
    delArray (pt);
    delArray (rr);
    
}
template void redux::image::ScharmerFilter (float**, double**, size_t, size_t, double, double);
template void redux::image::ScharmerFilter (double**, double**, size_t, size_t, double, double);
template void redux::image::ScharmerFilter (complex_t**, double**, size_t, size_t, double, double);


//double inv_dist_wght(double **a,int xl,int xh,int yl,int yh,int xb,int yb)
double redux::image::inv_dist_wght (float **a, size_t sizeY, size_t sizeX, size_t posY, size_t posX) {

    int xl = std::max (0L, static_cast<int64_t> (posX - maxDistance));
    int xh = std::min (sizeX, posX + maxDistance + 1);
    int yl = std::max (0L, static_cast<int64_t> (posY - maxDistance));
    int yh = std::min (sizeY, posY + maxDistance + 1);


    double weight = 0.0, res = 0.0;
    for (int y = yl; y < yh; ++y)
        for (int x = xl; x < xh; ++x)
            if (a[y][x]) {
                double c = pow ( (double) sqr (x - posX) + (double) sqr (y - posY) + deltasqr, -beta);
                res += c * a[y][x];
                weight += c;
            }
    return res / weight;
}


template <typename T>
double redux::image::inverseDistanceWeight (T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX) {

    static const double* const inverseDistanceSquared = getDistanceMap();     // will only be initialized once.

    // TODO: verify this function, results look weird
    int64_t beginX = std::max (0L, static_cast<int64_t> (posX - maxDistance));
    int64_t endX = std::min (sizeX, posX + maxDistance+1);
    int64_t beginY = std::max (0L, static_cast<int64_t> (posY - maxDistance));
    int64_t endY = std::min (sizeY, posY + maxDistance+1);

    double normalization = 0.0, weightedSum = 0.0;
    for (int y=beginY; y < endY; ++y) {
        int y2 = (y-posY)*(y-posY);
        for (int x = beginX; x < endX; ++x) {
            int x2 = (x-posX)*(x-posX);
            //if( x2+y2 > maxDistance2 ) break;
            if( array[y][x] ) {
                double tmp = inverseDistanceSquared[y2+x2];
                weightedSum += tmp * array[y][x];
                normalization += tmp;
            }
        }
    }
    if( normalization ) {
        return weightedSum / normalization;
    }
    return 0.0;
}
template double redux::image::inverseDistanceWeight (unsigned char**, size_t, size_t, size_t, size_t);
template double redux::image::inverseDistanceWeight (short**, size_t, size_t, size_t, size_t);
template double redux::image::inverseDistanceWeight (int**, size_t, size_t, size_t, size_t);
template double redux::image::inverseDistanceWeight (float**, size_t, size_t, size_t, size_t);
template double redux::image::inverseDistanceWeight (double**, size_t, size_t, size_t, size_t);

template <typename T>
double redux::image::horizontalInterpolation (T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX) {

    T* ptr = array[posY];

    //map the 5 pixel surrounding as bits in a byte
    int val = 0;
    if ( (posX > 1)) val |= ( (ptr[posX - 2] > 0) << 4);
    if ( (posX > 0))   val |= ( (ptr[posX - 1] > 0) << 3);
    if ( (posX+1 < sizeX))   val |= ( (ptr[posX + 1] > 0) << 1);
    if ( (posX+2 < sizeX)) val |= ( (ptr[posX + 2] > 0));
    //now select based on the number
    switch (val) {
        case (10) :     // = 0 1 x 1 0
        case (11) :     // = 0 1 x 1 1
        case (26) :     // = 1 1 x 1 0
        case (27) :     // = 1 1 x 1 1
            return (ptr[posX-1] + ptr[posX+1]) / 2;
        case (18) :     // = 1 0 x 1 0
        case (19) :     // = 1 0 x 1 1
            return (ptr[posX-2] + 2 * ptr[posX+1]) / 3;
        case (9) :      // = 0 1 x 0 1
        case (25) :     // = 1 1 x 0 1
            return (2 * ptr[posX-1] + ptr[posX+2]) / 3;
        default:
            return inverseDistanceWeight<T> (array, sizeY, sizeX, posY, posX);
    }

}
template double redux::image::horizontalInterpolation (uint8_t**, size_t, size_t, size_t, size_t);
template double redux::image::horizontalInterpolation (int16_t**, size_t, size_t, size_t, size_t);
template double redux::image::horizontalInterpolation (uint16_t**, size_t, size_t, size_t, size_t);
template double redux::image::horizontalInterpolation (int32_t**, size_t, size_t, size_t, size_t);
template double redux::image::horizontalInterpolation (uint32_t**, size_t, size_t, size_t, size_t);
template double redux::image::horizontalInterpolation (float**, size_t, size_t, size_t, size_t);
template double redux::image::horizontalInterpolation (double**, size_t, size_t, size_t, size_t);

template <typename T>
void redux::image::apodizeInPlace2( T** data, size_t nRows, size_t nCols, size_t rowBlend, size_t colBlend ) {

    if( rowBlend > nCols/2 ) {
        rowBlend = nCols/2;
    }
    if( colBlend > nRows/2 ) {
        colBlend = nRows/2;
    }
    
    size_t sz = std::max( 2*rowBlend, 2*colBlend )+2;
    double* apod = new double[ sz ];
    T avg1(0);
    Array<T> bla( nCols, 2*rowBlend );
    for( size_t c=0; c<nCols; ++c ) {
        for( size_t r=0; r<rowBlend; ++r ) {
            bla( c, rowBlend+r ) = data[r][c];
            bla( c, rowBlend-1-r ) = data[nRows-1-r][c];
            avg1 += bla( c, rowBlend+r ) + bla( c, rowBlend-1-r );
        }
    }
    file::Ana::write( "aip2_bt.f0", bla );
    avg1 /= (2*nCols*rowBlend);
    
    T avg1b(0);
    memset( apod, 0, sz*sizeof(T) );
    redux::math::apodize( apod, 2*rowBlend+2, 1.0 );
    for( size_t c=0; c<nCols; ++c ) {
        T lastVal = data[nRows-rowBlend-1][c];
        T firstVal = data[rowBlend][c];
        T diff = firstVal-lastVal;
        for( size_t r=0; r<rowBlend; ++r ) {
            double a1 = apod[rowBlend+r+1];
            T db = firstVal-data[2*rowBlend-r][c];
            T de = lastVal-data[nRows-2*rowBlend-1+r][c];
            data[r][c] = a1*(diff+db)+lastVal-0.5*db;
            data[nRows-1-r][c] = a1*(de-diff)+firstVal-0.5*de;
            bla( c, rowBlend+r ) = data[r][c];
            bla( c, rowBlend-1-r ) = data[nRows-1-r][c];
            avg1b += bla( c, rowBlend+r ) + bla( c, rowBlend-1-r );
        }
    }
    file::Ana::write( "aip2_btb.f0", bla );
    avg1b /= (2*nCols*rowBlend);
    
    T avg2(0);
    bla.resize( nRows, 2*colBlend );
    for( size_t r=0; r<nRows; ++r ) {
        for( size_t c=0; c<colBlend; ++c ) {
            bla( r, colBlend-1-c ) = data[r][nCols-1-c];
            bla( r, colBlend+c ) = data[r][c];
            avg2 += bla( r, colBlend-1-c ) + bla( r, colBlend+c );
        }
    }
    file::Ana::write( "aip2_rl.f0", bla );
    avg2 /= (2*nRows*colBlend);
    
    if( colBlend != rowBlend ) {
        memset( apod, 0, sz*sizeof(T) );
        redux::math::apodize( apod, 2*colBlend+2, 1.0 );
    }
    T avg2b(0);
    for( size_t r=0; r<nRows; ++r ) {
        T lastVal = data[r][nCols-colBlend-1];
        T firstVal = data[r][colBlend];
        T diff = firstVal-lastVal;
        for( size_t c=0; c<colBlend; ++c ) {
            double a1 = apod[colBlend+c+1];
            T db = firstVal-data[r][2*colBlend-c];
            T de = lastVal-data[r][nCols-2*colBlend-1+c];
            data[r][c] = a1*(diff+db)+lastVal-0.5*db;
            data[r][nCols-1-c] = a1*(de-diff)+firstVal-0.5*de;
            bla( r, colBlend-1-c ) = data[r][nCols-1-c];
            bla( r, colBlend+c ) = data[r][c];
            avg2b += bla( r, colBlend-1-c ) + bla( r, colBlend+c );
        }
    }
    file::Ana::write( "aip2_rlb.f0", bla );
    avg2b /= (2*nRows*colBlend);
    cout << "Avg1: " << avg1 << "  Avg1b: " << avg1b << endl;
    cout << "Avg2: " << avg2 << "  Avg2b: " << avg2b << endl;
    
    delete[] apod;

}
template void redux::image::apodizeInPlace2(int16_t**, size_t, size_t, size_t, size_t);
template void redux::image::apodizeInPlace2(int32_t**, size_t, size_t, size_t, size_t);
template void redux::image::apodizeInPlace2(float**, size_t, size_t, size_t, size_t);
template void redux::image::apodizeInPlace2(double**, size_t, size_t, size_t, size_t);
template void redux::image::apodizeInPlace2(complex_t**, size_t, size_t, size_t, size_t);


template <typename T>
void redux::image::apodizeInPlace( T** data, size_t nRows, size_t nCols, size_t rowBlend, size_t colBlend, size_t rowMargin, size_t colMargin ) {

    if( rowBlend+rowMargin > nCols/2 ) {
        rowBlend = nCols/2-rowMargin;
    }
    if( colBlend+colMargin > nRows/2 ) {
        colBlend = nRows/2-colMargin;
    }
    
    size_t sz = std::max( rowBlend+rowMargin, colBlend+colMargin ) + 2;
    shared_ptr<double> tmpD = rdx_get_shared<double>( sz );
    double* tmp = tmpD.get();
    
    memset( tmp, 0, sz*sizeof(T) );
    redux::math::apodize( tmp+rowMargin, rowBlend+2, 1.0 );
    for( size_t c=0; c<nCols; ++c ) {
        for( size_t r=0; r<rowBlend+rowMargin; ++r ) {
            data[r][c] *= tmp[r+1];
            data[nRows-r-1][c] *= tmp[r+1];
        }
    }
    memset( tmp, 0, sz*sizeof(T) );
    redux::math::apodize( tmp+colMargin, colBlend+2, 1.0 );
    for( size_t rr=0; rr<nRows; ++rr ) {
        for( size_t c=0; c<colBlend+colMargin; ++c ) {
            data[rr][c] *= tmp[c+1];
            data[rr][nCols-c-1] *= tmp[c+1];
        }
    }

}
template void redux::image::apodizeInPlace(int16_t**, size_t, size_t, size_t, size_t, size_t, size_t);
template void redux::image::apodizeInPlace(int32_t**, size_t, size_t, size_t, size_t, size_t, size_t);
template void redux::image::apodizeInPlace(float**, size_t, size_t, size_t, size_t, size_t, size_t);
template void redux::image::apodizeInPlace(double**, size_t, size_t, size_t, size_t, size_t, size_t);
template void redux::image::apodizeInPlace(complex_t**, size_t, size_t, size_t, size_t, size_t, size_t);


template <typename T>
void redux::image::normalizeIfMultiFrames (redux::image::Image<T>& img) {
    if (img.meta) {
        vector<string> hdrTexts = img.meta->getText();
        string hdr;
        if( !hdrTexts.empty() ) hdr = hdrTexts.front();
        boost::regex re ("(\\d+)[ .]+SUM[= ]+");
        boost::smatch match;
        if (boost::regex_search (hdr, match, re)) {
            int nFrames = boost::lexical_cast<int> (match[1]);
            if (nFrames > 1) {
                img *= (1.0/nFrames);
            }
        }
    }
}
template void redux::image::normalizeIfMultiFrames (redux::image::Image<int16_t>&);
template void redux::image::normalizeIfMultiFrames (redux::image::Image<int32_t>&);
template void redux::image::normalizeIfMultiFrames (redux::image::Image<double>&);
template void redux::image::normalizeIfMultiFrames (redux::image::Image<float>&);


template <typename T, typename U>
void redux::image::inpaint( T* img, U* mask, T* out, size_t ySize, size_t xSize, double radius, int flags ) {
#ifdef RDX_WITH_OPENCV
    cv::Mat imgMat( ySize, xSize, cv::cvType<T>(), img );
    cv::Mat maskMat( ySize, xSize, cv::cvType<U>(), mask );
    cv::Mat outMat( ySize, xSize, cv::cvType<T>(), out );
    std::cerr << "inpaint is not yet implemented for OpenCV builds." << std::endl;
    //cv::inpaint( imgMat, maskMat, outMat, radius, flags );
#else
    std::cerr << "inpaint is not yet implemented for non-OpenCV builds." << std::endl;
#endif            
}
template void redux::image::inpaint( float*, uint8_t*, float*, size_t, size_t, double, int );
template void redux::image::inpaint( double*, uint8_t*, double*, size_t, size_t, double, int );


template <typename T, typename U>
void redux::image::resize( T* in, size_t inSizeY, size_t inSizeX, U* out, size_t outSizeY, size_t outSizeX ) {
#ifdef RDX_WITH_OPENCV
    cv::Mat src( inSizeY, inSizeX, cv::cvType<T>(), in );
    cv::Mat dst( outSizeY, outSizeX, cv::cvType<T>(), out );
    if( (inSizeY > outSizeY) && (inSizeX > outSizeX) ) {
        resize( src, dst, dst.size(), 0, 0, cv::INTER_AREA );
    } else if( (inSizeY < outSizeY) && (inSizeX < outSizeX) ) {
        resize( src, dst, dst.size(), 0, 0, cv::INTER_CUBIC );
    } else {
        resize( src, dst, dst.size(), 0, 0, cv::INTER_CUBIC );
    }
#else
    std::cerr << "redux::image::resize is not yet implemented for non-OpenCV builds." << std::endl;
#endif            
}
template void redux::image::resize( int16_t* in, size_t, size_t, int16_t* out, size_t, size_t );
template void redux::image::resize( int32_t* in, size_t, size_t, int32_t* out, size_t, size_t );
template void redux::image::resize( float* in, size_t, size_t, float* out, size_t, size_t );
template void redux::image::resize( double* in, size_t, size_t, double* out, size_t, size_t );


