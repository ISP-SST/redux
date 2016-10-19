#include "redux/math/functions.hpp"

#include "redux/constants.hpp"

#include <cmath>
#include <complex>
#include <vector>

using namespace redux::math;
using namespace redux;

namespace {
    
    const double log_16 = log( 16 );
    
    /* coefficients of the rational approximation formula
       to the complementary error function */
    const double A0 = 122.607931777104326;
    const double A1 = 214.382388694706425;
    const double A2 = 181.928533092181549;
    const double A3 = 93.155580458138441;
    const double A4 = 30.180142196210589;
    const double A5 = 5.912626209773153;
    const double A6 = 0.564189583562615;
    const double B0 = 122.60793177387535;
    const double B1 = 352.730625110963558;
    const double B2 = 457.334478783897737;
    const double B3 = 348.703917719495792;
    const double B4 = 170.354001821091472;
    const double B5 = 53.992906912940207;
    const double B6 = 10.479857114260399;

}


template <typename T>
void redux::math::apodize( T* data, size_t n, const T& target ) {

    double mean = 0.5 * static_cast<double>(*data + target);
    double amplitude = 0.5 * static_cast<double>(*data - target);
    double step = PI / (n-1);     // FIXME: this is the correct apodization step, using MvN's version below for now...
    //double step = 2*PI / (2*n-1);

    for( size_t i = 1; i < n; ++i ) {
        *(data + i) = mean + amplitude * cos( i * step );
    }

}
template void redux::math::apodize( int16_t*, size_t, const int16_t&);
template void redux::math::apodize( int32_t*, size_t, const int32_t&);
template void redux::math::apodize( float*, size_t, const float&);
template void redux::math::apodize( double*, size_t, const double&);


void redux::math::cauchylorentz( double* data, size_t n, double fwhm, float ncav ) {    // lorentzian: ncav = 1.  for prefilters ncav is number of cavities.
    
    if( n < 3 ) {
        return;
    }
    
    double* ptr1 = data + (n>>1);
    double* ptr2 = ptr1 - 1;
    double a = 4 / ( fwhm * fwhm );
    double i(0);

    if( n & 1 ) {
        *ptr1 = 1;
        ptr1++;
        i = 0.5;
    }
    
    if( ncav == 1.0 ) {
        while( ptr2 >= data ) {
            *ptr2 = *ptr1 = 1.0 / ( i * i * a + 1 );
            ptr1++;
            ptr2--;
            i += 1;
        }
    } else {
        while( ptr2 >= data ) {
            *ptr2 = *ptr1 = 1.0 / ( pow( i * i * a, ncav ) + 1 );
            ptr1++;
            ptr2--;
            i += 1;
        }

    }

}


void redux::math::gauss( double* data, size_t n, double dispersion, double fwhm ) {

    if( n < 3 ) {
        return;
    }
    
    double* ptr1 = data + (n>>1);
    double* ptr2 = ptr1 - 1;
    double a = log_16 / (fwhm*fwhm);

    if( n & 1 ) {
        *ptr1 = 1;
        ptr1++;
    } else {
        *ptr2 = *ptr1 = exp( -a * dispersion*dispersion / 4 );
        ptr1++;
        ptr2--;
    }
    //return;
    a = exp( -a * dispersion*dispersion );
    double b = a;

    while( ptr2 >= data ) {
        //*ptr1++ = exp( -a * i * i );
        *ptr2 = *ptr1 = *(ptr1-1) * b;
        ptr1++;
        ptr2--;
        b *= a * a;

    }


}

template <typename T>
void redux::math::gauss( T* data, size_t sizeY, size_t sizeX, double fwhmY, double fwhmX, double centerY, double centerX) {
    // TODO: this implementation is naĩve, optimize if necessary...
    std::vector<double> ydiffsquared(sizeY);
    std::vector<double> xdiffsquared(sizeX);
    double ay = log_16 / (fwhmY*fwhmY);
    double ax = log_16 / (fwhmX*fwhmX);
    for(size_t y=0; y<sizeY; ++y) {
        double tmp = (y-centerY);
        ydiffsquared[y] = tmp*tmp;
    }
    for(size_t x=0; x<sizeX; ++x) {
        double tmp = (x-centerX);
        xdiffsquared[x] = tmp*tmp;
    }
    for(size_t y=0; y<sizeY; ++y) {
        for(size_t x=0; x<sizeX; ++x) {
            *(data+y*sizeX+x) = exp( -(xdiffsquared[x]*ax + ydiffsquared[y]*ay) );
        }
    }
}
template void redux::math::gauss( float*, size_t, size_t, double, double, double, double);
template void redux::math::gauss( double*, size_t, size_t, double, double, double, double);


double redux::math::hann( double x ) {
    return 0.5 * ( 1.0 - cos( PI * x ) );
}

void redux::math::hann( double* data, size_t n ) {

    if( !data || n < 3 ) {
        return;
    }

    double arg = ( PI / ( n - 1 ) );
    if( n & 1 ) {
        data[n >> 1] = 0.5;
    }
    n--;
    for( size_t i=0; i<n; ++i,--n ) {
        double cai =  0.5 * cos( arg * i );
        data[i] = 0.5 - cai;
        data[n] = 0.5 + cai;
    }

}


template <class T>
void redux::math::hann( T& data, size_t n ) {
    if( n < 3 ) {
        return;
    }
    data.resize( n );
    double arg = PI / (n-1);
    if( n & 1 ) {
        data[n >> 1] = 0.5;
    }
    n--;
    for( size_t i=0; i<n; ++i,--n ) {
        double cai =  0.5 * cos( arg * i );
        data[i] = static_cast<typename T::value_type>(0.5 - cai);
        data[n] = static_cast<typename T::value_type>(0.5 + cai);
    }
}
template void redux::math::hann( std::vector<double>&, size_t n );
template void redux::math::hann( std::vector<float>&, size_t n );


double lorentz( double x, double center, double fwhm, double ncav ) {    // lorentzian: ncav = 1.  for prefilters ncav is number of cavities.
    double a = x - center;
    double b = 2*a/fwhm;
    b *= b;
    if( ncav != 1.0 ) {
        b = pow( b, ncav );
    }
    
    return 1.0 / ( b + 1 );
}


void redux::math::faraday_voigt( double* vgt, double* far, size_t n, double dispersion, double damping ) {

    if( !vgt || n < 2 ) {
        return;
    }

    double* vgt1 = vgt + (n>>1);
    double* vgt2 = vgt1 - 1;
    
    double* far1 = far + (n>>1);
    double* far2 = far1 - 1;
    
    double offset(0);
    
    if( n & 1 ) {
        *vgt1 = 1;
        vgt1++;
        *far1 = 0;
        far1++;
        offset += dispersion;
    } else {
        offset += 0.5*dispersion;
    }
    
    
    if( damping == 0 ) {
        while( vgt2 >= vgt ) {
            *vgt1++ = *vgt2-- = exp( -offset*offset );
            offset += dispersion;
            *far1 = offset * ( 1. - .66666667 * ( offset * offset ) ) * 5.641895836E-1;
            *far2-- = - *far1++;
       }
    } else {
        std::complex<double> Z(damping),ZZ;
        while( vgt2 >= vgt ) {
            Z.imag(-offset);
            ZZ = ( ( ( ( ( ( A6 * Z + A5 ) * Z + A4 ) * Z + A3 ) * Z + A2 ) * Z + A1 ) * Z + A0 ) /
                ( ( ( ( ( ( ( Z + B6 ) * Z + B5 ) * Z + B4 ) * Z + B3 ) * Z + B2 ) * Z + B1 ) * Z + B0 );
            *vgt1++ = *vgt2-- = ZZ.real();
            offset += dispersion;
            *far1 = .5 * ZZ.imag();
            *far2-- = - *far1++;
        }

    }

}

template <class T>
void redux::math::faraday_voigt( T& vgt, T& far, double dispersion, double damping ) {

    size_t n = vgt.size();
    if( n != far.size() || n < 2 ) {
        return;
    }

    size_t i(n>>1);
    size_t ii(i-1);
    
    double offset(0);
    if( n & 1 ) {
        vgt[i] = 1;
        far[i] = 0;
        i++;
        offset += dispersion;
    } else {
        offset += 0.5*dispersion;
    }
    
    if( damping == 0 ) {
        do {
            vgt[i] = vgt[ii] = static_cast<typename T::value_type>(exp( -offset*offset ));
            offset += dispersion;
            far[i] = static_cast<typename T::value_type>(offset * ( 1. - .66666667 * ( offset * offset ) ) * 5.641895836E-1);
            far[ii--] = - far[i++];
       } while( ii > 0 );
    } else {
        std::complex<double> Z(damping),ZZ;
         do {
            Z.imag(-offset);
            ZZ = ( ( ( ( ( ( A6 * Z + A5 ) * Z + A4 ) * Z + A3 ) * Z + A2 ) * Z + A1 ) * Z + A0 ) /
                ( ( ( ( ( ( ( Z + B6 ) * Z + B5 ) * Z + B4 ) * Z + B3 ) * Z + B2 ) * Z + B1 ) * Z + B0 );
            vgt[i] = vgt[ii] = static_cast<typename T::value_type>(ZZ.real());
            offset += dispersion;
            far[i] = static_cast<typename T::value_type>(.5 * ZZ.imag());
            far[ii--] = - far[i++];
        } while( ii > 0 );

    }

}
template void redux::math::faraday_voigt( std::vector<double>&, std::vector<double>&, double dispersion, double damping );
template void redux::math::faraday_voigt( std::vector<float>&, std::vector<float>&, double dispersion, double damping );

template <class T>
void redux::math::faraday_voigt( T& vgt, T& far, T& lambda, double dispersion, double damping ) {

    size_t n = vgt.size();
    if( n != far.size() || n < 2 ) {
        return;
    }

    size_t i(n>>1);
    size_t ii(i-1);
    
    double offset(0);
    if( n & 1 ) {
        vgt[i] = 1;
        far[i] = 0;
        lambda[i] = 0;
        i++;
        offset += dispersion;
    } else {
        offset += 0.5*dispersion;
    }
    
    if( damping == 0 ) {
        do {
            lambda[i] = static_cast<typename T::value_type>(offset);
            lambda[ii] = -lambda[i];
            vgt[i] = vgt[ii] = static_cast<typename T::value_type>(exp( -offset*offset ));
            offset += dispersion;
            far[i] = static_cast<typename T::value_type>(offset * ( 1. - .66666667 * ( offset * offset ) ) * 5.641895836E-1);
            far[ii--] = - far[i++];
       } while( ii > 0 );
    } else {
        std::complex<double> Z(damping),ZZ;
        do {
            Z.imag(-offset);
            ZZ = ( ( ( ( ( ( A6 * Z + A5 ) * Z + A4 ) * Z + A3 ) * Z + A2 ) * Z + A1 ) * Z + A0 ) /
                ( ( ( ( ( ( ( Z + B6 ) * Z + B5 ) * Z + B4 ) * Z + B3 ) * Z + B2 ) * Z + B1 ) * Z + B0 );
            lambda[i] = static_cast<typename T::value_type>(offset);
            lambda[ii] = -lambda[i];
            vgt[i] = vgt[ii] = static_cast<typename T::value_type>(ZZ.real());
            offset += dispersion;
            far[i] = static_cast<typename T::value_type>(.5 * ZZ.imag());
            far[ii--] = -far[i++];
        } while( ii > 0 );

    }

}
template void redux::math::faraday_voigt( std::vector<double>&, std::vector<double>&, std::vector<double>&, double dispersion, double damping );
template void redux::math::faraday_voigt( std::vector<float>&, std::vector<float>&, std::vector<float>&, double dispersion, double damping );

void redux::math::faraday_voigt(double* vgt, double* far, double offset, double damping ) {

    int offsetsign = 1;
    std::complex<double> Z;
    
    if( offset < 0 ) {
        offset = -offset;
        offsetsign = -1;
    }
    
    // If the damping is zero, then the convolution is dominated by the gaussian
    if( damping == 0 ) {
        *vgt = exp( - ( offset * offset ) );
        *far = offset * ( 1. - .66666667 * ( offset * offset ) ) * offsetsign * 5.641895836E-1;
    }
    else {
        Z = std::complex<double> ( damping, -offset ) ;
        Z = ( ( ( ( ( ( A6 * Z + A5 ) * Z + A4 ) * Z + A3 ) * Z + A2 ) * Z + A1 ) * Z + A0 ) /
            ( ( ( ( ( ( ( Z + B6 ) * Z + B5 ) * Z + B4 ) * Z + B3 ) * Z + B2 ) * Z + B1 ) * Z + B0 );
        *vgt = Z.real();
        *far = .5 * offsetsign * Z.imag();
    }
}


unsigned long redux::math::n_choose_k( int n, int k ) {

    if( k > n ) return 0;
    if( k == n ) return 1;

    unsigned long temp = 1;
    int imax, base;
    if( n < (k<<1) ) {
        imax = k + 1;
        base = n - k;
    } else {
        imax = n - k + 1;
        base = k;
    }
    for( int i = 1; i < imax; ++i ) {
        temp *= ( base + i );
        temp /= i;
    }

    return temp;
}

