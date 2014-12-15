#include "redux/image/utils.hpp"

#include "redux/image/fouriertransform.hpp"

#include "redux/math/functions.hpp"
#include "redux/constants.hpp"
#include "redux/types.hpp"

#include <map>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

namespace {

    const int maxDistance = 32;
    const double deltasqr = 4;
    const double beta = 2;

    map<int, double> getDistanceMap ( void ) {
        map<int, double> tmp;
        for ( int i = 0; i <= maxDistance; ++i ) {
            int i2 = i * i;
            for ( int j = 0; j <= maxDistance; ++j ) {
                int j2 = j * j;
                tmp.insert ( pair<int, double> ( i2 + j2, pow ( i2 + j2 + deltasqr, -beta ) ) );
            }
        }
        return tmp;
    }

}

template <typename T>
static __inline T sqr ( T x ) {
    return x * x;
}


Grid::Grid( uint32_t nPoints ) : size(nPoints,nPoints), origin(nPoints/2,nPoints/2) {
    if( !(nPoints%2) ) {    // for nPoints even, place origin between mid-points
        origin.x += 0.5;
        origin.y += 0.5;
    }
    init();
}

Grid::Grid( uint32_t nPoints, float originY, float originX ) : size(nPoints,nPoints), origin(originY,originX) {
    init();
}

Grid::Grid( uint32_t nPointsY, uint32_t nPointsX, float originY, float originX ) : size(nPointsY, nPointsX), origin(originY,originX) {
    init();
}

void Grid::init(void) {
    
    distance = sharedArray<float>(size.y, size.x);
    angle = sharedArray<float>(size.y, size.x);
    float** distPtr = distance.get();
    float** anglePtr = angle.get();
    for(uint i = 0; i < size.y; ++i) {
        double yDist = i - origin.y;
        double y2 = sqr(yDist);
        for(uint j = 0; j < size.x; ++j) {
            double xDist = j - origin.x;
            if(yDist || xDist) {
                double x2 = sqr(xDist);
                distPtr[i][j] = sqrt(y2 + x2);
                anglePtr[i][j] = atan2(yDist, xDist);       // note: slow index is Y, fast is X
            } else distPtr[i][j] = anglePtr[i][j] = 0;      // this pixel is at the origin -> set the angle to 0. 
        }
    }   

}


double redux::image::makePupil( util::Array<double>& aperture, uint32_t nPoints, float radius, bool normalize ) {
    
    double area = 0.0, origin = 0.0;
    
    uint32_t mid = nPoints/2;
    bool odd = nPoints%2;       // N.B: Don't use odd number of points for pupils!! Implemented just for completion.
    
    odd = true;                 // Pupil should be centered ON pixel (mid,mid), to match the location of the origin for Fourier-transforms. 
                                // Forcing "odd" will center pupil on (mid,mid) rather than making it symmetric around a point between pixels. 
    
    if( odd ) origin = 0.5;
                                
    Grid grid(mid+1,origin,origin);
    float** distPtr = grid.distance.get();
    aperture.resize(nPoints,nPoints);
    aperture.zero();
    
    for(uint x = 0; x < mid; ++x) {
        for(uint y = 0; y <=x; ++y) {
            double val = 0;
            if(distPtr[y+1][x+1] < radius) {
                val = 1;
                //val = distPtr[y+1][x+1];
            } else if(distPtr[y][x] < radius) {
                val = (radius-distPtr[y][x])/(distPtr[y+1][x+1]-distPtr[y][x]); // linear fill-factor from radial ratio
                // TODO: better approximation of pixel fill-factor.
            }
            if(odd) {
                if( x == 0 && y == 0 ) {    // central pixel = 1 for all practical cases
                                            // a pupil of size < sqrt(2) pixel is a bit absurd, this is just for completeness/fun :-)
                    if( radius > sqrt(2)/2.0 ) val = 1;
                    else if( radius < 0.5 ) val = redux::PI*radius*radius;
                    else val = redux::PI*radius*radius + (radius-0.5)/(sqrt(2)/2.0-0.5)*(1-redux::PI*radius*radius);
                    aperture(mid,mid) = val;
                    area += val;
                } else {
                    aperture(mid+y,mid+x) = val;
                    area += val;
                    if ( x != 0 ) {
                        aperture(mid+y,mid-x) = val;
                        area += val;
                        if ( y != 0 ) {
                            aperture(mid-y,mid-x) = val;
                            area += val;
                        }
                    }
                    if ( y != 0 ) {
                        aperture(mid-y,mid+x) = val;
                        area += val;
                    }
                    if( x != y ) {
                        aperture(mid+x,mid+y) = val;
                        area += val;
                        if ( x != 0 ) {
                            aperture(mid-x,mid+y) = val;
                            area += val;
                            if ( y != 0 ) {
                                aperture(mid-x,mid-y) = val;
                                area += val;
                            }
                        }
                        if ( y != 0 ) {
                            aperture(mid+x,mid-y) = val;
                            area += val;
                        }
                    }
                }
            } else {
                aperture(mid+y,mid+x) = val;
                aperture(mid+y,mid-x-1) = aperture(mid-y-1,mid+x) = aperture(mid-y-1,mid-x-1) = val;
                area += 4*val;
                if(x != y) {
                    aperture(mid+x,mid+y) = val;  // Use symmetry to fill quadrants (if off-diagonal)
                    aperture(mid+x,mid-y-1) = aperture(mid-x-1,mid+y) = aperture(mid-x-1,mid-y-1) = val;
                    area += 4*val;
                }
            }
        }
    }   
    
    return area;
    
}


double redux::image::makePupil_mvn( double** pupil, uint32_t nph, float r_c ) {
    
    double area = 0.0;
    memset(*pupil,0,(nph+1)*(nph+1)*sizeof(double));
    int xo = 1 + nph / 2, yo = 1 + nph / 2;
    double dx = 0.5 / r_c, dy = 0.5 / r_c;
    for(uint x = 1; x <= nph; ++x) {
        double xl = fabs((double)(x - xo)) / r_c - dx, xh = fabs((double)(x - xo)) / r_c + dx;
        double xhs = sqr(xh);
        for(uint y = 1; y <= nph; ++y) {
            double yl = fabs((double)(y - yo)) / r_c - dy, yh = fabs((double)(y - yo)) / r_c + dy;
            double yhs = sqr(yh);
            double rsl = sqr(xl) + sqr(yl), rsh = xhs + yhs;
            if(rsl <= 1.0) {   // inside pixel
                if(rsh < 1.0) {   // full pixel
                    pupil[x][y] = 1.0;
                    //pupil[x][y] = sqrt(rsh*r_c*r_c);
                } else {           // partial pixel
                    double x2 = sqrt(max(1.0 - yhs, (double)0.0));
                    double y3 = sqrt(max(1.0 - xhs, (double)0.0));
                    double f = (xh > yh) ? (yh - yl) * (min(xh, max(xl, x2)) - xl) / (4 * dx * dy) :
                             (xh - xl) * (min(yh, max(yl, y3)) - yl) / (4 * dx * dy);
                    pupil[x][y] = f;
                }
                area += pupil[x][y];
            }
            else pupil[x][y] = 0.0; // outside pixel
        }
    }

    return area;
    
}



// nverseDistanceWeight( T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX ) {

//double inv_dist_wght(double **a,int xl,int xh,int yl,int yh,int xb,int yb)
double redux::image::inv_dist_wght ( float **a, size_t sizeY, size_t sizeX, size_t posY, size_t posX ) {

    int xl = std::max ( 0L, static_cast<int64_t> ( posX - maxDistance ) );
    int xh = std::min ( sizeX, posX + maxDistance + 1 );
    int yl = std::max ( 0L, static_cast<int64_t> ( posY - maxDistance ) );
    int yh = std::min ( sizeY, posY + maxDistance + 1 );


    double weight = 0.0, res = 0.0;
    for ( int y = yl; y < yh; ++y )
        for ( int x = xl; x < xh; ++x )
            if ( a[y][x] ) {
                double c = pow ( ( double ) sqr ( x - posX ) + ( double ) sqr ( y - posY ) + deltasqr, -beta );
                res += c * a[y][x];
                weight += c;
            }
    return res / weight;
}


template <typename T>
double redux::image::inverseDistanceWeight ( T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX ) {

    static map<int, double> inverseDistanceSquared = getDistanceMap();

    // TODO: verify this function, results look weird

    int64_t beginX = std::max ( 0L, static_cast<int64_t> ( posX - maxDistance ) );
    int64_t endX = std::min ( sizeX, posX + maxDistance + 1 );
    int64_t beginY = std::max ( 0L, static_cast<int64_t> ( posY - maxDistance ) );
    int64_t endY = std::min ( sizeY, posY + maxDistance + 1 );

    double normalization = 0.0, weightedSum = 0.0;
    for ( int y = beginY; y < endY; ++y ) {
        int y2 = ( y - posY ) * ( y - posY );
        for ( int x = beginX; x < endX; ++x ) {
            int x2 = ( x - posX ) * ( x - posX );
            if ( array[y][x] > 0 ) {
                double tmp = inverseDistanceSquared.at ( y2 + x2 );
                weightedSum += tmp * array[y][x];
                normalization += tmp;
            }
        }
    }
    return weightedSum / normalization;
}
template double redux::image::inverseDistanceWeight ( float**, size_t, size_t, size_t, size_t );
template double redux::image::inverseDistanceWeight ( double**, size_t, size_t, size_t, size_t );

template <typename T>
double redux::image::horizontalInterpolation ( T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX ) {

    T* ptr = array[posY];

    //map the 5 pixel surrounding as bits in a byte
    int val = 0;
    if ( ( posX > 1 ) ) val |= ( ( ptr[posX - 2] > 0 ) << 4 );
    if ( ( posX > 0 ) )   val |= ( ( ptr[posX - 1] > 0 ) << 3 );
    if ( ( posX < sizeX ) )   val |= ( ( ptr[posX + 1] > 0 ) << 1 );
    if ( ( posX < sizeX - 1 ) ) val |= ( ( ptr[posX + 2] > 0 ) );
    //now select based on the number
    switch ( val ) {
        case ( 10 ) :   // = 0 1 x 1 0
        case ( 11 ) :   // = 0 1 x 1 1
        case ( 26 ) :   // = 1 1 x 1 0
        case ( 27 ) :   // = 1 1 x 1 1
            return ( ptr[posX - 1] + ptr[posX + 1] ) / 2;
        case ( 18 ) :   // = 1 0 x 1 0
        case ( 19 ) :   // = 1 0 x 1 1
            return ( ptr[posX - 2] + 2 * ptr[posX + 1] ) / 3;
        case ( 9 ) :    // = 0 1 x 0 1
        case ( 25 ) :   // = 1 1 x 0 1
            return ( 2 * ptr[posX - 1] + ptr[posX + 2] ) / 3;
        default:
            return inverseDistanceWeight<T> ( array, sizeY, sizeX, posY, posX );
    }

}
template double redux::image::horizontalInterpolation ( float**, size_t, size_t, size_t, size_t );
template double redux::image::horizontalInterpolation ( double**, size_t, size_t, size_t, size_t );


template <typename T>
Array<T> redux::image::apodize ( const Array<T>& in, size_t blendRegion ) {
    if ( !blendRegion ) return in;     // nothing to do
    Array<T> array;
    in.copy(array);
    const vector<int64_t>& first = array.first();
    size_t sizeY = array.last() [0]-first[0]+1;
    size_t sizeX = array.last() [1]-first[1]+1;
    blendRegion = std::min ( std::min ( blendRegion,sizeY ),sizeX );

    double* tmp = new double[blendRegion+1];
    tmp[0] = 0;
    redux::math::apodize ( tmp, blendRegion, 1.0 );
//     auto sharedPtrs = array.get(sizeY,sizeX);
//     T** data = sharedPtrs.get();
    for ( size_t y=0,yy=sizeY-1; y<sizeY; ++y,--yy ) {
        double yfactor = 1;
        if ( y < blendRegion ) yfactor *= tmp[y];
        for ( size_t x=0,xx=sizeX-1; x<sizeX; ++x,--xx ) {
            double xfactor = yfactor;
            if ( x < blendRegion ) xfactor *= tmp[x];
            if ( xfactor < 1 ) {
                array ( y,x )   *= xfactor;
                array ( yy,x )  *= xfactor;
                array ( y,xx )  *= xfactor;
                array ( yy,xx ) *= xfactor;
            }
        }
    }
    delete[] tmp;
    return array;
}
template Array<int16_t> redux::image::apodize ( const Array<int16_t>&, size_t );
template Array<int32_t> redux::image::apodize ( const Array<int32_t>&, size_t );
template Array<double> redux::image::apodize ( const Array<double>&, size_t );
template Array<float> redux::image::apodize ( const Array<float>&, size_t );
template Array<complex_t> redux::image::apodize ( const Array<complex_t>&, size_t );

template <typename T>
void redux::image::checkIfMultiFrames( redux::image::Image<T>& img ) {
    if(img.hdr) {
        std::string hdr = img.hdr->getText();
        boost::regex re( "(\\d+)[ .]+SUM[= ]+" );
        boost::smatch match;
        if( boost::regex_search( hdr, match, re ) ) {
            int nFrames = boost::lexical_cast<int>( match[1] );
            if( nFrames > 1) {
                img.nFrames = nFrames;
            } else img.nFrames = 1;
        }
    }
}
template void redux::image::checkIfMultiFrames( redux::image::Image<int16_t>& );
template void redux::image::checkIfMultiFrames( redux::image::Image<int32_t>& );
template void redux::image::checkIfMultiFrames( redux::image::Image<double>& );
template void redux::image::checkIfMultiFrames( redux::image::Image<float>& );


template <typename T, typename U>
void redux::image::descatter ( Array<T>& data, const Array<U>& ccdgain, const Array<U>& psf_in, int maxIterations, double minImprovement, double epsilon ) {

    vector<size_t> dims = data.dimensions ( true );

    if ( dims.size() != 2 || dims != ccdgain.dimensions() || dims != psf_in.dimensions() ) {
        cout << "descatter(): dimensions of gain/psf does not match image." << endl;
        return;
    }

    for ( auto &it: dims ) it *= 2;

    Array<double> img( dims );                 // Twice the size of input
    Array<double> img_center(img, dims[0]/4, 3*dims[0]/4-1, dims[1]/4, 3*dims[1]/4-1 ); // centered subimage of half (i.e. original) size. N.B: shares data with img.
    img.zero();                                 // whole array = 0
    img_center = data;                          // central subimage = img_in

    Array<double> psf( dims );                 // Twice the size of input
    Array<double> psf_center(psf, dims[0]/4, 3*dims[0]/4-1, dims[1]/4, 3*dims[1]/4-1 ); // centered subimage of half (i.e. original) size. N.B: shares data with psf.
    psf.zero();                                 // whole array = 0
    psf_center = psf_in;                        // central subimage = psf_in
    double sum = total ( psf_in );
    if ( sum ) {                                // normalize
        psf_center *= 1.0/sum;
    }

    Array<double> gain;
    ccdgain.copy ( gain );

    redux::image::FourierTransform otf ( psf, FT_REORDER|FT_NORMALIZE );

    Array<double>::const_iterator itg = gain.begin();
    for ( auto& it: img_center ) {
        double g = *itg++;
        it /= ( 1.0 + g*g );
    }

    Array<double> tmp( dims );                 // Twice the size of input
    Array<double> tmp_center(tmp, dims[0]/4, 3*dims[0]/4-1, dims[1]/4, 3*dims[1]/4-1 ); // centered subimage of half (i.e. original) size. N.B: shares data with tmp.

    double metric = std::numeric_limits< double >::max();
    double delta;
    int i = 0;
    do {
        img.copy(tmp);
        tmp_center *= gain;
        otf.convolveInPlace(tmp);
        tmp_center *= gain;
        delta = metric;
        metric = 0.0;
        typename Array<T>::const_iterator in_it = data.begin();
        Array<double>::const_iterator tmp_it = tmp_center.begin();
        for ( auto& img_it: img_center ) {
            double new_img = (*in_it++ - *tmp_it++);
            metric += (img_it - new_img)*(img_it - new_img);
            img_it = new_img;
        }        
        metric /= data.nElements();
        delta /= metric;
        //cout << "descatter:  iter = " << i << ", ChiSq = " << metric <<  ", delta = " << delta << endl;
    } while ( (delta > minImprovement) && ( metric>epsilon ) && ( ++i<maxIterations ) );

    img_center.copy(data);

}
template void redux::image::descatter ( Array<float>&, const Array<float>&, const Array<float>&, int, double, double );
template void redux::image::descatter ( Array<double>&, const Array<float>&, const Array<float>&, int, double, double );
template void redux::image::descatter ( Array<float>&, const Array<double>&, const Array<double>&, int, double, double );


