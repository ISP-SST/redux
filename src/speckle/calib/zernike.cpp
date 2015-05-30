#include "zernike.hpp"

#include "defs.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include <gsl/gsl_math.h>   // for M_PI


#include <iostream>
#include "redux/util/array.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/momfbd/cache.hpp"
#include "redux/file/fileana.hpp"
using namespace redux::file;
using namespace redux::momfbd;
using namespace redux::util;

using namespace redux::speckle;
using namespace std;

#define RADIAL_SAMPLES     200001                // radial sampling points for zernikes
#define ANGULAR_SAMPLES    400001                // angular sampling points for zernikes

/*
    Functions for the numerical integration of SLE and STF using the
    Wang & Markey model for partially compensated wavefront errors

    Helper library to compute Zernike polynomials
 */

namespace {     // Anonymuous namespace -> only visible from this compilation unit (.cpp file)

    const double pi_15 (1.5 * M_PI);
    const double pi_2 (2.0 * M_PI);
    const double pi_2_inv (1.0 / pi_2);

    const double stepr_inv (RADIAL_SAMPLES - 1);
    const double stepr (1.0 / stepr_inv);
    const double stepphi (2 * M_PI / (ANGULAR_SAMPLES - 1));
    const double stepphi_inv (1.0 / stepphi);

    double zernikeRadialValues[ (SPECKLE_IMAX + 1)][RADIAL_SAMPLES];
    double cosineValues[ANGULAR_SAMPLES];

    bool doInit (true);


    /*
     *      factorial: compute k! recursively
     */
    long double factorial(uint32_t k) {
        static vector<long double> factorials;
        uint32_t fsz = factorials.size();
        if (fsz <= k) {
            factorials.resize (k + 1, 1);
            for (uint32_t i = std::max (fsz, 1U); i <= k; ++i) factorials[i] = i * factorials[i - 1];
        }
        return factorials[k];
    }

    int initAngular (void) {
        for (int i = 0; i < ANGULAR_SAMPLES; ++i) {
            cosineValues[i] = cos (i * stepphi);
        }
        return 0;
    }

}



/*
 *      convert Noll/Wang index i to Born/Wolf n and m with proper signs
 */
void redux::speckle::noll_to_nm (int i, int& n, int& m) {
    n = 0;
    int len = 1;

    for (int j = 1; len < i; ++j) {
        len += (n = j) + 1;
    }

    int dl = n + 1 - len + i;
    m = 2 * ( (dl + (n % 2)) / 2) + ! (n % 2) - 1;
    m *= ( (i % 2) ? -1 : 1);    // sign of m
}



double redux::speckle::calcZernikeCovariance (int i, int j) {

    if ( (i < 1) || (j < 1)) return 0.0;

    int m, n, o, p;
    noll_to_nm (i, n, m);
    noll_to_nm (j, p, o);

    if (m != o) return 0.0;                         // only the same azimuthal order is non-zero.

    if (m && (i + j) % 2) return 0.0;               // only sine-sine or cosine-cosine is non-zero (or m=0 => only radial component)

    //  ; Now deal with the numerical terms: Dai
    int isign;
    static double kk = pow (4.8 * exp (lgamma_r (6.0 / 5.0, &isign)), 5.0 / 6.0) *
                       exp (lgamma_r (14.0 / 3.0, &isign) + 2.0 * lgamma_r (11.0 / 6.0, &isign)) / (pow (2.0, (8.0 / 3.0)) * M_PI);
    double k = kk * sqrt ( (double) ( (n + 1) * (p + 1)));

    double g1 = lgamma_r ( ( (double) (n + p) -  5.0 / 3.0) / 2.0, &isign);
    double g2 = lgamma_r ( ( (double) (n - p) + 17.0 / 3.0) / 2.0, &isign);
    double g3 = lgamma_r ( ( (double) (p - n) + 17.0 / 3.0) / 2.0, &isign);
    double g4 = lgamma_r ( ( (double) (n + p) + 23.0 / 3.0) / 2.0, &isign);

    isign = ( ( (n + p - 2 * m) / 2) % 2) ? -1 : +1;

    return k * exp (g1 - g2 - g3 - g4) * isign;

}


void redux::speckle::calcRadialZernike (double* out, uint32_t nPoints, uint16_t n, uint16_t abs_m) {

    long double step = 1.0L / static_cast<long double> (nPoints - 1);

    if (n == 0) {
        for (uint32_t i = 0; i < nPoints; ++i) out[i] = 1;
        return;
    } else if (n == 1) {
        for (uint32_t i = 0; i < nPoints; ++i) out[i] = i * step;
    }

    long double* rn = new long double[nPoints];
    long double* r_2 = new long double[nPoints];

    for (uint32_t i = 0; i < nPoints; ++i) {
        long double tmp = i * step;
        r_2[i] = 1.0L / (tmp * tmp);
        rn[i] = pow (tmp, static_cast<long double> (n));                  // initialize to r^n
    }

    int32_t nmm = (n - abs_m) / 2;
    int32_t npm = (n + abs_m) / 2;
    vector<long double> coeff (nmm + 1);

    for (int32_t s = 0, pm = -1; s <= nmm; ++s) {
        coeff[s] = (double) ( (pm *= -1) * factorial (n - s)) / (factorial (s) * factorial (npm - s) * factorial (nmm - s))*sqrt ( (double) (n + 1));
    }

    memset (out, 0, nPoints * sizeof (double));

    for (auto & c : coeff) {
        for (uint32_t i = 0; i < nPoints; ++i) {
            out[i] += c * rn[i];
            rn[i] *= r_2[i];                                      //  next term r^{n-2} etc.
        }
    }
    if(abs_m == 0) {
                Array<double> wrapper (out, nPoints);
                Ana::write ("newr_"+to_string(n)+".f0", wrapper);
    }
    delete[] rn;
    delete[] r_2;

}


// TODO: test accuracy of samplings and pick some smart sizes...
uint32_t selectSampling (uint32_t n) {
//     if(n<2) return 2;
//     if(n<4) return 2000;
//     if(n<8) return 20000;
//     if(n<16) return 100000;
//     return 200000;
    return RADIAL_SAMPLES;
}


ZernikeData::ZernikeData::cov_t::cov_t (uint32_t n1_in, uint32_t n2_in, int32_t m_in, double cov_in) :
    n1 (max (n1_in, n2_in)), n2 (min (n1_in, n2_in)), m (m_in), cov (cov_in) {
    sampling1 = selectSampling (n1);
    sampling2 = selectSampling (n2);
}

void ZernikeData::ZernikeData::cov_t::swap (void) {
    std::swap(n1,n2);
    std::swap(sampling1,sampling1);
}

bool ZernikeData::ZernikeData::cov_t::operator< (const cov_t& rhs) const {
    if (m == rhs.m) {
        if (n1 == rhs.n1) {
            return (n2 < rhs.n2);
        }
        return (n1 < rhs.n1);
    }
    return (m < rhs.m);
}


bool ZernikeData::ZernikeData::data_index::operator< (const data_index& rhs) const {
    if (n == rhs.n) {
        if (abs_m == rhs.abs_m) {
            return (sampling < rhs.sampling);
        }
        return (abs_m < rhs.abs_m);
    }
    return (n < rhs.n);
}


ZernikeData& ZernikeData::get (void) {
    static ZernikeData zd;
    return zd;
}


void ZernikeData::init (const vector<float>& eff, uint32_t Imax, float cov_cutoff) {
    {
        unique_lock<mutex> lock (mtx);
        static int dummy = initAngular();                   // static variable will only initialize once, i.e. only call initAngular once.
        int32_t n1,n2,nEff = eff.size();
        int32_t m;                                          // m is the same, or cov=0.0
        covariances.clear();

        for (uint16_t i1 = 2; i1 <= nEff; ++i1) {                   // i1,i2 matches the Noll-indices
            for (uint16_t i2 = 2; i2 <= Imax; ++i2) {
                double a = calcZernikeCovariance (i1, i2);
                if (fabs (a) > cov_cutoff) {
                    noll_to_nm (i1, n1, m);
                    noll_to_nm (i2, n2, m);
                    cov_t tmpc(n1, n2, m, a);
                    if (i2 > nEff) {
                        tmpc.cov *= eff[i1 - 1];                    // ...but the list of efficiencies start with Z_2, so subtract 1.
                        tmpc.swap();                                // Needed to use the right datapointer for the QPcorr part.
                    } else {
                        tmpc.cov *= eff[i2 - 1] * (1.0 - 0.5 * eff[i1 - 1]);
                    }

                    auto ret = covariances.emplace (tmpc);

                    if (!ret.second) {  // already existing => this is the same Zernikes, but swapped  => add the covariances
                        ret.first->cov += tmpc.cov;
                    }
                }
            }
        }
    }   // end of locked scope

    for (auto & it : covariances) { // get pointers to data
        it.data1 = radialData (it.n1, abs (it.m), it.sampling1);
        it.data2 = radialData (it.n2, abs (it.m), it.sampling2);
    }

}




std::vector<float> ZernikeData::init (const string& filename, uint32_t Imax, float cov_cutoff) {
    
    ifstream fp( filename );
    string tmp;
    float val;

    vector<float> eff;
    while( fp.good() ) {
        std::getline( fp, tmp );
        stringstream ss( tmp );
        ss >> val;

        if( !ss.fail() ) {
            eff.push_back( val );
        }
    }
    
    init(eff,Imax,cov_cutoff);
    
    return std::move(eff);
    
}


std::shared_ptr<double> ZernikeData::radialData (uint32_t n, uint32_t abs_m, uint32_t nPoints) {
    unique_lock<mutex> lock (mtx);
    data_index index = {n, abs_m, nPoints};
    auto it = radialPolynomials.find (index);
    if (it != radialPolynomials.end()) {
        return it->second;
    }
    shared_ptr<double> data (new double[nPoints], [](double* p) { delete[] p; } );
    calcRadialZernike (data.get(), nPoints, n, abs_m);
    return radialPolynomials.emplace(index, data).first->second;
}


double ZernikeData::QRsumQPcorr (const double& plus1_rho, const double& plus1_phi, const double& minus1_rho, const double& minus1_phi,
                                 const double& plus2_rho, const double& plus2_phi, const double& minus2_rho, const double& minus2_phi) const {

    double val = 0.0;
    for (const auto & it : covariances) {
        if (it.data1 && it.data2) {
            double tmp = zernikePolar(it.data1.get(), it.sampling1, it.m, plus1_rho, plus1_phi)
                       - zernikePolar(it.data1.get(), it.sampling1, it.m, minus1_rho, minus1_phi);
            tmp *= zernikePolar(it.data2.get(), it.sampling2, it.m, plus2_rho, plus2_phi)
                 - zernikePolar(it.data2.get(), it.sampling2, it.m, minus2_rho, minus2_phi);
            val += it.cov * tmp;
        }
    }
    return val;

}


void redux::speckle::init_zernike (void) {

    if (doInit) {
        for (int i = 0; i < ANGULAR_SAMPLES; ++i) {
            cosineValues[i] = cos (i * stepphi);
        }

        int32_t n;
        int32_t m;
        for (int j = 1; j <= SPECKLE_IMAX; j++) {
            noll_to_nm (j, n, m);
            for (int i = 0; i < RADIAL_SAMPLES; i++) {
                zernikeRadialValues[j][i] = zernrad ((i * stepr), n, m) * sqrt ( (double) (n + 1));
            }
            if( m == 0 ) {
                Array<double> wrapper (zernikeRadialValues[j], RADIAL_SAMPLES);
                Ana::write ("oldr_"+to_string(n)+".f0", wrapper);
            }
        }
        doInit = false;
    }

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
        long double denominator = factorial (s) * factorial ( (n + m) / 2 - s) * factorial ( (n - m) / 2 - s);
        long double ilg = factorial (n - s) / denominator;
        result  += sign * ilg * pow (r, (double) (n - 2 * s));
        sign *= -1;
    }

    return (double) result;
}


const double* redux::speckle::getZRV(uint32_t i) {
    if (i > SPECKLE_IMAX) return nullptr;
    return zernikeRadialValues[i];
}

uint32_t redux::speckle::getSamplesR (void) {
    return RADIAL_SAMPLES;
}

/*
 *      zernike:        Return value of Zernike polynominal of degree
 *                      (dn,dl) with magnitude mc (cosine term) and ms
 *                      (sine term) for a Cartesian position x and y.
 *
 *      i:              Zernike index (Noll)
 *      x, y:           Co-ordinates in pupil with radius 1.
 */
double redux::speckle::zernike (int i, int m, double x, double y) {
    
    if (i < 1 || i > SPECKLE_IMAX) {
        return 0.0;
    }

    return zernike(zernikeRadialValues[i], RADIAL_SAMPLES, m, x, y);

}


double redux::speckle::zernikePolar (int i, int m, double rho, double phi) {

    if (i < 1 || i > SPECKLE_IMAX) {
        return 0.0;
    }

    return zernikePolar (zernikeRadialValues[i], RADIAL_SAMPLES, m, rho, phi);

}


double redux::speckle::zernikeDiffPolar (int i, int m, const double& rho1, const double& phi1, const double& rho2, const double& phi2) {

    if (i < 1 || i > SPECKLE_IMAX) {
        return 0.0;
    }

    return zernikeDiffPolar (zernikeRadialValues[i], RADIAL_SAMPLES, m, rho1, phi1, rho2, phi2);

}


double redux::speckle::zernike(const double* rdata, uint32_t sampling, int32_t m, double x, double y) {

    double rho = sqrt (x * x + y * y);
    double phi = 0;
    if (rho > 1.0 ) {
        return 0.0;
    } else if(rho > SPECKLE_EPS) {
        phi = atan2 (y, x);
    }
    return zernikePolar(rdata,sampling,m,rho,phi);
}


double redux::speckle::zernikePolar (const double* rdata, uint32_t sampling, int32_t m, double rho, double phi) {

    if (rho > 1.0) {
        return 0.0;
    }

    sampling--;     // we use (sampling-1) below, so make it easier...

    double tmp = rho * sampling;
    uint32_t idx = static_cast<uint32_t> (tmp);
    double z = rdata[idx];

    if (idx < sampling) {         // linear interpolation
        tmp -= idx;
        z = (1.0 - tmp) * z + tmp * rdata[idx + 1];
    }

    if (m == 0) {
        return z;
    }

    z *= M_SQRT2;
    if (rho > SPECKLE_EPS)  {
        phi *= abs (m);
        if (m < 0) {   // odd mode -> convert cosine to sine by adding 3/2*\pi
            phi += pi_15;
        }
        int factor = static_cast<int> (phi * pi_2_inv);
        if (phi < 0) factor--;
        phi -= factor * pi_2;
        idx = static_cast<uint32_t> (phi * stepphi_inv + 0.5);
        z *= cosineValues[idx];
    }

    return z;
}


// TODO: This is basically 2 calls to the above function, see if it can be done smarter...
double redux::speckle::zernikeDiffPolar (const double* rdata, uint32_t sampling, int32_t m, const double& rho1, const double& phi1, const double& rho2, const double& phi2) {

    double z1 (0);
    double z2 (0);

    sampling--;     // we use (sampling-1) below, so make it easier...

    double tmp = rho1 * sampling;
    uint32_t idx = static_cast<uint32_t> (tmp);
    z1 = rdata[idx];

    if (idx < sampling) {        // linear interpolation
        tmp -= idx;
        z1 = (1.0 - tmp) * z1 + tmp * rdata[idx + 1];
    }

    tmp = rho2 * sampling;
    idx = static_cast<uint32_t> (tmp);
    z2 = rdata[idx];

    if (idx < sampling) {        // linear interpolation
        tmp -= idx;
        z2 = (1.0 - tmp) * z2 + tmp * rdata[idx + 1];
    }

    if (m == 0) {
        return (z1 - z2);
    }

    z1 *= M_SQRT2;
    z2 *= M_SQRT2;

    if (rho1 > SPECKLE_EPS) {
        tmp = abs (m) * phi1;

        if (m < 0) {   // odd mode -> convert cosine to sine by adding 3/2*\pi
            tmp += pi_15;
        }

        int factor = static_cast<int> (tmp * pi_2_inv);

        if (tmp < 0) factor--;

        tmp -= factor * pi_2;
        idx = static_cast<uint32_t> (tmp * stepphi_inv + 0.5);
        z1 *= cosineValues[idx];
    }

    if (rho2 > SPECKLE_EPS) {
        tmp = abs (m) * phi2;

        if (m < 0) {   // odd mode -> convert cosine to sine by adding 3/2*\pi
            tmp += pi_15;
        }

        int factor = static_cast<int> (tmp * pi_2_inv);

        if (tmp < 0) factor--;

        tmp -= factor * pi_2;
        idx = static_cast<uint32_t> (tmp * stepphi_inv + 0.5);
        z2 *= cosineValues[idx];
    }

    return (z1 - z2);
}
