#include "zernike.hpp"

#include "defs.hpp"

#include <cstdint>
#include <cstdlib>
#include <vector>

#include <gsl/gsl_math.h>   // for M_PI

using namespace redux::speckle;


#define RADIAL_SAMPLES     200001                // radial sampling points for zernikes
#define ANGULAR_SAMPLES    400001                // angular sampling points for zernikes

/*
    Functions for the numerical integration of SLE and STF using the
    Wang & Markey model for partially compensated wavefront errors

    Helper library to compute Zernike polynomials
 */

namespace {     // Anonymuous namespace -> only visible from this compilation unit (.cpp file)
    
    const double pi_15(1.5 * M_PI);
    const double pi_2(2.0 * M_PI);
    
    const double stepr(1.0 / (RADIAL_SAMPLES-1));
    const double stepphi(2*M_PI / (ANGULAR_SAMPLES-1));
    
    double rz[SPECKLE_IMAX+1][RADIAL_SAMPLES];
    double cs[ANGULAR_SAMPLES];
    
    bool doInit(true);

}



/*
 *      convert Noll/Wang index i to Born/Wolf n and m with proper signs
 */
void redux::speckle::noll_to_nm(int i, int& n, int& m) {
    n = 0;
    int len = 1;
    for(int j = 1; len < i; ++j) {
        len += (n = j) + 1;
    }
    int dl = n + 1 - len + i;
    m = 2 * ((dl + (n % 2)) / 2) + !(n % 2) - 1;
    m *= ((i % 2) ? -1 : 1);  // sign of m
}


double redux::speckle::calcZernikeCovariance(int i, int j) {

    if((i < 2) || (j < 2)) return 0.0;
    int m, n, o, p;
    noll_to_nm(i, n, m);
    noll_to_nm(j, p, o);
    if(m != o) return 0.0;
    if(m && ((i + j) % 2)) return 0.0;

    //  ; Now deal with the numerical terms: Dai
    int isign;
    static double kk = pow(4.8 * exp(lgamma_r(6.0 / 5.0, &isign)), 5.0 / 6.0) *
                       exp(lgamma_r(14.0 / 3.0, &isign) + 2.0 * lgamma_r(11.0 / 6.0, &isign)) / (pow(2.0, (8.0 / 3.0)) * M_PI);
    double k = kk*sqrt((double)((n + 1) * (p + 1)));

    double g1 = lgamma_r(((double)(n + p) -  5.0 / 3.0) / 2.0, &isign);
    double g2 = lgamma_r(((double)(n - p) + 17.0 / 3.0) / 2.0, &isign);
    double g3 = lgamma_r(((double)(p - n) + 17.0 / 3.0) / 2.0, &isign);
    double g4 = lgamma_r(((double)(n + p) + 23.0 / 3.0) / 2.0, &isign);
    
    isign = (((n + p - 2 * m)/2) % 2) ? -1 : +1;

    return k * exp(g1 - g2 - g3 - g4) * isign;

}


void redux::speckle::init_zernike(void) {
    
    if (doInit) {
        for (int i = 0; i < ANGULAR_SAMPLES; ++i) {
            cs[i] = cos (i * stepphi);
        }
        int n, l;
        for (int j = 1; j <= SPECKLE_IMAX; j++) {
            noll_to_nm (j, n, l);
            for (int i = 0; i < RADIAL_SAMPLES; i++) {
                rz[j][i] = zernrad( i*stepr, n, l )*sqrt((double)(n+1));
            }
        }
        doInit = false;
    }
    
}

/*
 *      fak: compute k! recursively
 */
long double fak(uint32_t k) {

    static std::vector<long double> factorials;
    
    uint32_t fsz = factorials.size();
    if(fsz <= k) {
        factorials.resize(k + 1, 1);
        for(uint32_t i = std::max(fsz, 1U); i <= k; ++i) factorials[i] = i * factorials[i-1];
    }
    return factorials[k];

}

/*
 *      zernrad: compute radial Zernike polynomial for degree (n,l)
 *               at radius r
 */
double redux::speckle::zernrad (double r, int n, int l) {
    
    int sign = 1;
    int m = std::abs (l);
    long double result = 0;

    for (int s = 0; s <= (n - m) / 2; s++)  {
        long double denominator = fak (s) * fak ( (n + m) / 2 - s) * fak ( (n - m) / 2 - s);
        long double ilg = fak (n - s) / denominator;
        result  += sign * ilg * pow (r, (double) (n - 2 * s));
        sign *= -1;
    }

    return (double) result;
}

/*
 *      zernike:        Return value of Zernike polynominal of degree
 *                      (dn,dl) with magnitude mc (cosine term) and ms
 *                      (sine term) for a Cartesian position x and y.
 *
 *      i:              Zernike index (Noll)
 *      x, y:           Co-ordinates in pupil with radius 1.
 */
double redux::speckle::zernike(int i, int m, double x, double y) {

    double rho = sqrt(x*x + y*y);

    if ( rho > 1.0 || i < 1 || i > SPECKLE_IMAX ) {
        return 0.0;
    }
    
    int idxr = static_cast<int>(rho / stepr);
    double z = rz[i][idxr];
    if (idxr < (RADIAL_SAMPLES-1)) {          // linear interpolation
        double remainder = (rho / stepr) - idxr;
        z = (1.0 - remainder) * z + remainder * rz[i][idxr+1];
    }
    
    if (m == 0) {
        return z;
    } else {
        z *= M_SQRT2;
    }
    
    if (rho > SPECKLE_EPS)  {
        double phi = m*(atan2(y, x) + 2*M_PI);
        if( i % 2 ) {  // odd mode -> convert cosine to sine by adding 3/2*\pi 
            phi += 1.5 * M_PI; // sine
        }
        int factor = static_cast<int>(phi/pi_2);
        if(phi < 0) factor--;
        phi -= factor*pi_2;
        int idxphi = static_cast<int>(phi/stepphi+0.5);
        z *= cs[idxphi];
    }

    return z;

}


double redux::speckle::zernikePolar(int i, int m, double rho, double phi) {
    
    if ( rho > 1.0 || i < 1 || i > SPECKLE_IMAX ) {
        return 0.0;
    }

    int idxr = static_cast<int>(rho / stepr);
    double z = rz[i][idxr];
    if ( idxr < (RADIAL_SAMPLES-1)) {          // linear interpolation
        double remainder = (rho / stepr) - idxr;
        z = (1.0 - remainder) * z + remainder * rz[i][idxr+1];
    }
    
    if (m == 0) {
        return z;
    } else {
        z *= M_SQRT2;
    }
    
    if (rho > SPECKLE_EPS)  {
        phi *= m;
        if( i % 2 ) {  // odd mode -> convert cosine to sine by adding 3/2*\pi 
            phi += pi_15;
        }
        int factor = static_cast<int>(phi/pi_2);
        if(phi < 0) factor--;
        phi -= factor*pi_2;
        int idxphi = static_cast<int>(phi/stepphi+0.5);
        z *= cs[idxphi];
    }

    return z;
}


double redux::speckle::zernikeDiffPolar(int i, int m, double rho1, double phi1, double rho2, double phi2) {
    
    if ( i < 1 || i > SPECKLE_IMAX ) {
        return 0.0;
    }

    double z1(0);
    double z2(0);
    if ( rho1 <= 1.0 ) {
        int idxr1 = static_cast<int>(rho1 / stepr);
        z1 = rz[i][idxr1];
        if (idxr1 < (RADIAL_SAMPLES-1)) {          // linear interpolation
            double remainder = (rho1 / stepr) - idxr1;
            z1 = (1.0 - remainder)*z1 + remainder * rz[i][idxr1+1];
        }
    }
    
    if ( rho2 <= 1.0 ) {
        int idxr2 = static_cast<int>(rho2 / stepr);
        z2 = rz[i][idxr2];
        if (idxr2 < (RADIAL_SAMPLES-1)) {          // linear interpolation
            double remainder = (rho2 / stepr) - idxr2;
            z2 = (1.0 - remainder)*z2 + remainder * rz[i][idxr2+1];
        }
    }
    
    if (m == 0) {
        return (z1-z2);
    }
    
    z1 *= M_SQRT2;
    z2 *= M_SQRT2;
    
    if (rho1 > SPECKLE_EPS)  {
        phi1 *= m;
        if( i % 2 ) {  // odd mode -> convert cosine to sine by adding 3/2*\pi 
            phi1 += pi_15; // sine
        }
        int factor = static_cast<int>(phi1/pi_2);
        if(phi1 < 0) factor--;
        phi1 -= factor*pi_2;
        int idxphi = static_cast<int>(phi1/stepphi+0.5); // (phi + 0.5*step)/step
        z1 *= cs[idxphi];
    }
    
    if (rho2 > SPECKLE_EPS)  {
        phi2 *= m;
        if( i % 2 ) {  // odd mode -> convert cosine to sine by adding 3/2*\pi 
            phi2 += pi_15; // sine
        }
        int factor = static_cast<int>(phi2/pi_2);
        if(phi2 < 0) factor--;
        phi2 -= factor*pi_2;
        int idxphi = static_cast<int>(phi2/stepphi+0.5);
        z2 *= cs[idxphi];
    }

    return (z1-z2);
}




