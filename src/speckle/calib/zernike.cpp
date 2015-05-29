#include "zernike.hpp"

#include "defs.hpp"
#include "atmcov.hpp"

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

    double rz[ (SPECKLE_IMAX + 1)][RADIAL_SAMPLES];
    double cs[ANGULAR_SAMPLES];

    bool doInit (true);


    /*
     *      fak: compute k! recursively
     */
    long double fak (uint32_t k) {

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
            cs[i] = cos (i * stepphi);
        }

        return 0;
    }


    size_t dataCount (0);
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

    // cout << "calcRadialZernike()   nPoints = " << nPoints << "  n = " << n << "  abs_m = " << abs_m << endl;

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
        coeff[s] = (double) ( (pm *= -1) * fak (n - s)) / (fak (s) * fak (npm - s) * fak (nmm - s))*sqrt ( (double) (n + 1));
    }
// * sqrt ( (double) (n + 1));
//    std::reverse (coeff.begin(), coeff.end());      // reverse so that coeff[0] is the coefficient for r^n

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


// TODO: test accuracy of samplings and pick some smart size...
uint32_t selectSampling (uint16_t n) {
//     if(n<2) return 2;
//     if(n<4) return 2000;
//     if(n<8) return 20000;
//     if(n<16) return 100000;
//     return 200000;
    return RADIAL_SAMPLES;
}


ZernikeData::ZernikeData::cov_t::cov_t (uint16_t n1_in, uint16_t n2_in, int32_t m_in, double cov_in) :
    n1 (max (n1_in, n2_in)), n2 (min (n1_in, n2_in)), m (m_in), cov (cov_in) {
    sampling1 = selectSampling (n1);
    sampling2 = selectSampling (n2);
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


void ZernikeData::init (const vector<float>& eff, uint16_t Imax, float cov_cutoff) {
   // cout << " ZernikeData::init()   eff.sz = " << eff.size() << endl;
    {
        unique_lock<mutex> lock (mtx);
        static int dummy = initAngular();                   // static variable will only initialize once, i.e. only call initAngular once.
        uint16_t nEff = eff.size();
        covariances.clear();

        for (uint16_t i1 = 1; i1 <= nEff; ++i1) {                   // i1,i2 matches the Noll-indices
            for (uint16_t i2 = 1; i2 <= Imax; ++i2) {
                //double a = calcZernikeCovariance (i1, i2);
                double a = redux::speckle::atmcov (i1, i2);

                if (fabs (a) > cov_cutoff) {
                    int m, n1, n2;                                  // m is the same, or cov=0.0
                    noll_to_nm (i1, n1, m);
                    noll_to_nm (i2, n2, m);
                    cov_t tmpc (n1, n2, m, a);

                    if (i2 > nEff) {
                        tmpc.cov *= eff[i1 - 1];                  // ...but the list of efficiencies start with Z_2, so subtract 1.
                    } else {
                        tmpc.cov *= eff[i2 - 1] * (1.0 - 0.5 * eff[i1 - 1]);
                    }

                    auto ret = covariances.emplace (tmpc);

                    if (!ret.second) {  // already existing => this is the same Zernikes, but swapped  => add the covariances
                     //   cout << "DUP:   n1 = " << tmpc.n1 << "  n2 = " << tmpc.n2 << "  m = " << tmpc.m << "  cov = " << tmpc.cov << "  cov2 = " << ret.first->cov << flush;
                        ret.first->cov += tmpc.cov;
                     //   cout << "  newcov2 = " << ret.first->cov << endl;
                    }
                }
            }
        }
    }   // end of locked scope
 //   cout << " ZernikeData::init()   cov.sz = " << covariances.size() << endl;

    for (auto & it : covariances) { // get pointers to data
        //cout << "COV:   n1 = " << it.n1 << "  n2 = " << it.n2 << "  m = " << it.m << "  cov = " << it.cov << "  s1 = " << it.sampling1 << "  s2 = " << it.sampling2 << endl;
        it.data1 = radialData (it.n1, abs (it.m), it.sampling1);
        it.data2 = radialData (it.n2, abs (it.m), it.sampling2);
    }

 //   cout << " ZernikeData::init()   DataSize = " << dataCount << "   orig = " << ( (SPECKLE_IMAX + 1) *RADIAL_SAMPLES)
 //        << "   (" << (100 * dataCount / ( (double) ( (SPECKLE_IMAX + 1) *RADIAL_SAMPLES))) << "%)" << endl;

}


std::shared_ptr<double>& ZernikeData::radialData (uint16_t n, uint16_t abs_m, uint32_t nPoints) {
    unique_lock<mutex> lock (mtx);
    data_index index = {n, abs_m, nPoints};
    auto it = radialPolynomials.find (index);

    if (it != radialPolynomials.end()) {
        //cout << "radialData() EXISTING  nPoints = " << nPoints << "  n = " << n << "  abs_m = " << abs_m << endl;
        return it->second;
    }

    dataCount += nPoints;
    shared_ptr<double> data (new double[nPoints], [](double* p) {  /*cout << "FREEING DATA @ "  << hexString(p) << endl;*/  delete[] p; } );
    //cout << "radialData()    NEW DATA @ "  << hexString(data.get())<< "  nPoints = " << nPoints << "  n = " << n << "  abs_m = " << abs_m << endl;
    calcRadialZernike (data.get(), nPoints, n, abs_m);
   // radialPolynomials.emplace (index, data);
   // return data;
    return radialPolynomials.emplace(index, data).first->second;

}


double ZernikeData::QRsumQPcorr (const double& plus1_rho, const double& plus1_phi, const double& minus1_rho, const double& minus1_phi,
                                 const double& plus2_rho, const double& plus2_phi, const double& minus2_rho, const double& minus2_phi) const {

    double val = 0.0;
    for (const auto & it : covariances) {
        if (it.data1 && it.data2) {
//        cout << "QRsumQPcorr:                  data1 = " << hexString(it.data1.get()) << "  data2 = " << hexString(it.data2.get())
//              << "  m = " << it.m << "  sampling1 = " << it.sampling1 << "  sampling2 = " << it.sampling2 << endl;
            //double tmp =  zernikeDiffPolar (it.data1.get(), it.sampling1, it.m, plus1_rho, plus1_phi, minus1_rho, minus1_phi);
            //tmp *= zernikeDiffPolar (it.data2.get(), it.sampling2, it.m, plus2_rho, plus2_phi, minus2_rho, minus2_phi);
            double tmp =  zernikePolar(it.data1.get(), it.sampling1, it.m, plus1_rho, plus1_phi) - zernikePolar(it.data1.get(), it.sampling1, it.m, minus1_rho, minus1_phi);
            tmp *= zernikePolar(it.data2.get(), it.sampling2, it.m, plus2_rho, plus2_phi) - zernikePolar(it.data2.get(), it.sampling2, it.m, minus2_rho, minus2_phi);
            val += it.cov * tmp;
        } else {
        cout << "QRsumQPcorr:NULL at  n1 = " << it.n1 << "  n2 = " << it.n2
              << "  m = " << it.m << "  sampling1 = " << it.sampling1 << "  sampling2 = " << it.sampling2 << endl;
            
        }
    }

    return val;

}


void redux::speckle::init_zernike (void) {

 //   cout << " init_zernike()   doInit = " << doInit << endl;

    if (doInit) {
        for (int i = 0; i < ANGULAR_SAMPLES; ++i) {
            cs[i] = cos (i * stepphi);
        }

        int n, l;

//         vector<double> r(RADIAL_SAMPLES);
//         int cnt(0);
//         std::generate (r.begin(), r.end(), [&cnt](void){ return cnt++*stepr; } );
//cout << printArray(r,"RRRRRR") << endl;

        for (int j = 1; j <= SPECKLE_IMAX; j++) {
            noll_to_nm (j, n, l);

//             zernrad( r.data(), rz[j], RADIAL_SAMPLES, n, l );
            for (int i = 0; i < RADIAL_SAMPLES; i++) {
                rz[j][i] = zernrad (/*sqrt*/ (i * stepr), n, l) * sqrt ( (double) (n + 1));
            }
            if( l == 0 ) {
                Array<double> wrapper (rz[j], RADIAL_SAMPLES);
                Ana::write ("oldr_"+to_string(n)+".f0", wrapper);
            }
        }

  //      cout << endl << " init_zernike()   doInit = done.  " << rz[17][100] << "  " << rz[17][100] << endl;


//         {
//             Array<double> awrapper (cs, ANGULAR_SAMPLES);
//             Ana::write ("ang_old.f0", awrapper);
//         }
// 
//         {
//             Array<double> wrapper (rz[0], SPECKLE_IMAX, RADIAL_SAMPLES);
//             Ana::write ("rad_new.f0", wrapper);
//         }
// 
        /*        Cache& cache = Cache::getCache();

                for (int j = 1; j <= SPECKLE_IMAX; j++) {
                    noll_to_nm (j, n, l);
                    auto coeff = cache.zernikeRadialPolynomial (n, l);
                    cout << "j = " << j << "  nC = " << coeff.size() << printArray (coeff, "  coeff") << endl;

                    if (true || coeff.empty()) continue;

                    double r2;
                    int m = abs (l);

                    for (int i = 0; i < RADIAL_SAMPLES; i++) {
                        double r = i * stepr;
                        double r2 = r * r;
                        double val = 1;

                        if (m == 0) val = 1;
                        else if (m == 1) val = r;
                        else val = pow (r, m);                  // leading term ~ r^{m}

                        val *= coeff[0];

        //                 for (auto it = coeff.begin() + 1; it < coeff.end(); it++) {
        //                     rPtr[y][x] *= r2Ptr[y][x];                          //  next term ~ r^{m-2}
        //                     modePtr[y][x] += rPtr[y][x] * *it;
        //
        //                 }







                        rz[j][i] = zernrad (i * stepr, n, l) * sqrt ( (double) (n + 1));
                    }
                }
        */
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
        long double denominator = fak (s) * fak ( (n + m) / 2 - s) * fak ( (n - m) / 2 - s);
        long double ilg = fak (n - s) / denominator;
        result  += sign * ilg * pow (r, (double) (n - 2 * s));
        sign *= -1;
    }

    return (double) result;
}

void redux::speckle::zernrad (const double* r, double* out, int count, int n, int l) {

    int sign = 1;
    int m = std::abs (l);
    int nCoeffs = (n - m) / 2 + 1;
    double* coeffs = new double[nCoeffs];

    for (int s = 0; s < nCoeffs; ++s)  {
        long double denominator = fak (s) * fak ( (n + m) / 2 - s) * fak ( (n - m) / 2 - s);
        coeffs[s] = sign * fak (n - s) / denominator * sqrt ( (double) (n + 1));
        sign *= -1;
    }

    memset (out, 0, count * sizeof (double));

    for (int s = 0; s < nCoeffs; ++s)  {
        for (int i = 0; i < count; ++i)  {
            out[i] += coeffs[s] * pow (r[i], (double) (n - 2 * s));
        }

//         transform (r, r+count, out,
//                       [&coeffs,s,n](double rr){
//                          return coeffs[s] * pow( rr, ( double )( n - 2 * s ) ) * sqrt( ( double )( n + 1 ) );
//                       } );
    }

    delete[] coeffs;
}
void redux::speckle::zernrad2 (const double* r, double* out, int count, int n, int m) {

    int32_t nmm = (n - m) / 2;
    int32_t npm = (n + m) / 2;

    vector<double> coeff (nmm + 1);

    for (int32_t s = 0, pm = -1; s <= nmm; ++s) {
        coeff[s] = (double) ( (pm *= -1) * fak (n - s)) / (fak (s) * fak (npm - s) * fak (nmm - s));
    }

    std::reverse (coeff.begin(), coeff.end());

    /*        for(int y = 0; y < nPoints; ++y) {
                for(int x = 0; x < nPoints; ++x) {
                    //if(pupPtr[y][x]>0) {
                        double tmp = distPtr[y][x] / r_c;                       // normalize radial distance
                        //if(tmp>1) continue;
                        r2Ptr[y][x] = tmp * tmp;
                        if(m == 0) rPtr[y][x] = 1;
                        else if(m == 1) rPtr[y][x] = tmp;
                        else rPtr[y][x] = pow(tmp, m);                          // leading term ~ r^{m}
                        modePtr[y][x] = rPtr[y][x] * coeff[0];                  // add leading order to mode
                    //}
                }
            }

            // generate polynomial part
            for(auto it = coeff.begin() + 1; it < coeff.end(); it++) {
                for(int y = 0; y < nPoints; ++y) {
                    for(int x = 0; x < nPoints; ++x) {
                        //if(pupPtr[y][x]>0) {
                            rPtr[y][x] *= r2Ptr[y][x];                          //  next term ~ r^{m-2}
                            modePtr[y][x] += rPtr[y][x] * *it;
                        //}
                    }
                }
            }
    */

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

    double rho = sqrt (x * x + y * y);

    if (rho > 1.0 || i < 1 || i > SPECKLE_IMAX) {
        return 0.0;
    }

    double remainder = rho * stepr_inv;
    int idxr = static_cast<int> (remainder);
    double z = rz[i][idxr];

    if (idxr < stepr_inv) {        // linear interpolation
        remainder -= idxr;
        z = (1.0 - remainder) * z + remainder * rz[i][idxr + 1];
    }

    if (m == 0) {
        return z;
    } else {
        z *= M_SQRT2;
    }

    if (rho > SPECKLE_EPS)  {
        double phi = abs (m) * (atan2 (y, x) + 2 * M_PI);

        if (i % 2) {   // odd mode -> convert cosine to sine by adding 3/2*\pi
            phi += 1.5 * M_PI; // sine
        }

        int factor = static_cast<int> (phi * pi_2_inv);

        if (phi < 0) factor--;

        phi -= factor * pi_2;
        int idxphi = static_cast<int> (phi * stepphi_inv + 0.5);
        z *= cs[idxphi];
    }

    return z;

}


double redux::speckle::zernikePolar (int i, int m, double rho, double phi) {

    if (i < 1 || i > SPECKLE_IMAX) {
        return 0.0;
    }

    return zernikePolar (rz[i], RADIAL_SAMPLES, m, rho, phi);

}


double redux::speckle::zernikeDiffPolar (int i, int m, const double& rho1, const double& phi1, const double& rho2, const double& phi2) {

    if (i < 1 || i > SPECKLE_IMAX) {
        return 0.0;
    }

    return zernikeDiffPolar (rz[i], RADIAL_SAMPLES, m, rho1, phi1, rho2, phi2);

}


double redux::speckle::zernikePolar (const double* rdata, uint32_t sampling, int32_t m, double rho, double phi) {

    if (rho > 1.0) {
        return 0.0;
    }

    sampling--;     // we use (sampling-1) below, so make it easier...

    double tmp = rho * sampling;
    uint32_t idxr = static_cast<uint32_t> (tmp);
    double z = rdata[idxr];

    if (idxr < sampling) {         // linear interpolation
        tmp -= idxr;
        z = (1.0 - tmp) * z + tmp * rdata[idxr + 1];
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
        uint32_t idxphi = static_cast<uint32_t> (phi * stepphi_inv + 0.5);
        z *= cs[idxphi];
    }

    return z;
}


double redux::speckle::zernikeDiffPolar (const double* rdata, uint32_t sampling, int32_t m, const double& rho1, const double& phi1, const double& rho2, const double& phi2) {

    double z1 (0);
    double z2 (0);

    sampling--;     // we use (sampling-1) below, so make it easier...

    double tmp = rho1 * sampling;
    uint32_t idx = static_cast<uint32_t> (tmp);
//    cout << "zernikeDiffPolar()   data = " << hexString(rdata) << "  sampling = " << sampling << "   m = " << m << "   tmp = " << tmp << "   idx = " <<  idx << endl;
    z1 = rdata[idx];

    if (idx < sampling) {        // linear interpolation
        tmp -= idx;
        z1 = (1.0 - tmp) * z1 + tmp * rdata[idx + 1];
    }

    tmp = rho2 * sampling;
    idx = static_cast<uint32_t> (tmp);
//    cout << "zernikeDiffPolar()                                               tmp2 = " << tmp << "  idx2 = " <<  idx << endl;
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
        z1 *= cs[idx];
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
        z2 *= cs[idx];
    }

    return (z1 - z2);
}
