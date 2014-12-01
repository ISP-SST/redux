#include "redux/image/utils.hpp"

#include "redux/image/fouriertransform.hpp"

#include "redux/math/functions.hpp"

#include <map>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

using namespace redux::image;
using namespace redux::util;
using namespace std;

namespace {

    const int maxDistance = 32;
    const double deltasqr = 4;
    const double beta = 2;

    float distanceCache[maxDistance*maxDistance];

    bool calculateCache ( void ) {
        for ( int i = 0; i < maxDistance; ++i ) {
            double i2 = i * i;
            for ( int j = 0; j < maxDistance; ++j ) {
                double j2 = j * j;
                distanceCache[i * maxDistance + j] = pow ( sqrt ( i2 + j2 + deltasqr ), -beta );
            }
        }

        return true;
    }

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


static __inline int sqr ( int x ) {
    return x * x;
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

    size_t beginX = std::max ( 0L, static_cast<int64_t> ( posX - maxDistance ) );
    size_t endX = std::min ( sizeX, posX + maxDistance + 1 );
    size_t beginY = std::max ( 0L, static_cast<int64_t> ( posY - maxDistance ) );
    size_t endY = std::min ( sizeY, posY + maxDistance + 1 );

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
void redux::image::apodize ( Array<T>& array, size_t blendRegion ) {
    if ( !blendRegion ) return;     // nothing to do
    const vector<int64_t>& first = array.first();
    size_t sizeY = array.last() [0]-first[0]+1;
    size_t sizeX = array.last() [1]-first[1]+1;
    blendRegion = std::min ( std::min ( blendRegion,sizeY ),sizeX );

    T* tmp = new T[blendRegion+1];
    tmp[0] = 0;
    redux::math::apodize ( tmp, blendRegion, T ( 1 ) );
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
}
template void redux::image::apodize ( Array<int16_t>&, size_t );
template void redux::image::apodize ( Array<int32_t>&, size_t );
template void redux::image::apodize ( Array<double>&, size_t );
template void redux::image::apodize ( Array<float>&, size_t );

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

