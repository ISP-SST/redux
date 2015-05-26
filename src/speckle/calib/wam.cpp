#include "wam.hpp"

#include "defs.hpp"
#include "zernike.hpp"

#include <fstream>
#include <sstream>

#include <gsl/gsl_math.h>   // for M_PI and other constants

using namespace std;

/*
    Functions for the numerical integration of SLE and STF using the
    Wang & Markey model for partially compensated wavefront errors

    Helper library to compute mean Wang & Markey SE transfer function
 */


namespace {     // Anonymuous namespace -> only visible from this compilation unit (.cpp file)
    
    const double five_thirds( 5.0 / 3.0 );

    vector<float> eff;                      // these are the efficency factors, first value corresponds to Z_2 (tilt)
    struct weight_t {
        int i1,i2,m1,m2;
        double val;
    };
    // vector<weight_t> Kweights;              // eff[i2]*(1.0-0.5*eff[i1])*a           used in QRsum
    // vector<weight_t> Lweights;              // eff[i1]*a                             used in QPcorr
    vector<weight_t> KLweights;             // eff[i2]*(1.0-0.5*eff[i1])*a, OR eff[i1]*a if i2 > nEfficiencies      used in QRsumQPcorr

}


/*
 *   init_wam_master: read in the efficiency factors
 *   ===============================================
 */
vector<float>& redux::speckle::init_wam_master( char *filename ) {
    eff.clear();
    ifstream fp(filename);
    string tmp;
    float val;
    while (fp.good()) {
        std::getline(fp,tmp);
        stringstream ss(tmp);
        ss >> val;
        if (!ss.fail()) {
            eff.push_back(val);
        }
    }
    return eff;
}


/*
 *   init_wam_slave: set the efficiency factors
 *   ==========================================
 */
void redux::speckle::init_wam_slave( const vector<float>& data ) {
    
    eff = data;
    
    uint16_t nEff = eff.size();
    for(uint16_t i1=2; i1<=nEff; ++i1) {                        // i1,i2 matches the Noll-indices
        for(uint16_t i2=2; i2<=nEff; ++i2) {
            double a = calcZernikeCovariance( i1, i2 );
            if (fabs(a) > SPECKLE_EPS2) {
                int n,m1,m2;
                noll_to_nm(i1,n,m1);
                noll_to_nm(i2,n,m2);
                a *= eff[i2-1]*(1.0-0.5*eff[i1-1]);             // ...but the list of efficiencies start with Z_2, so subtract 1.
                //Kweights.push_back( {i1,i2,abs(m1),abs(m2),a});
                KLweights.push_back( {i1,i2,abs(m1),abs(m2),a});
            }
        }
    }

    for(uint16_t i1=2; i1<=nEff; ++i1) {
        for(uint16_t i2=nEff+1; i2<=SPECKLE_IMAX; ++i2) {
            double a = calcZernikeCovariance( i1, i2 );
            if (fabs(a) > SPECKLE_EPS2) {
                int n,m1,m2;
                noll_to_nm(i1,n,m1);
                noll_to_nm(i2,n,m2);
                a *= eff[i1-1];
                //Lweights.push_back( {i1,i2,abs(m1),abs(m2),a});
                KLweights.push_back( {i1,i2,abs(m1),abs(m2),a});
            }
        }
    }

}


/*
 *   QRsum: sum over compensated terms
 *          this is implemented with efficency factor, so if
 *          efficency factor is 1 everywhere the value is 1/2
 *          of the Wang & Markey Q
 *   ========================================================
 *   We're using:  plus1 = r+s  minus1 = r-s  plus2 = r'+s  minus2 = r'-s
 */
// double redux::speckle::QRsum( double plus1_rho, double plus1_phi, double minus1_rho, double minus1_phi,
//                                double plus2_rho, double plus2_phi, double minus2_rho, double minus2_phi ) {
// 
//     double val = 0.0;
//     for( const auto& it: Kweights ) {
//         double diff1 = zernikeDiffPolar( it.i1, it.m1, plus1_rho, plus1_phi, minus1_rho, minus1_phi );
//         double diff2 = zernikeDiffPolar( it.i2, it.m2, plus2_rho, plus2_phi, minus2_rho, minus2_phi );
//         val += it.val * diff1 * diff2;
// 
//     }
//     return val;
// 
// }


/*
 *   QPcorr:  sum over correlated terms
 *   ==================================
 *   We're using:  plus1 = r+s  minus1 = r-s  plus2 = r'+s  minus2 = r'-s
 */
// double redux::speckle::QPcorr( double plus1_rho, double plus1_phi, double minus1_rho, double minus1_phi,
//                                 double plus2_rho, double plus2_phi, double minus2_rho, double minus2_phi ) {
// 
//     double val = 0.0;
//     for( const auto& it: Lweights ) {
//         double diff1 = zernikeDiffPolar( it.i1, it.m1, plus1_rho, plus1_phi, minus1_rho, minus1_phi );
//         double diff2 = zernikeDiffPolar( it.i2, it.m2, plus2_rho, plus2_phi, minus2_rho, minus2_phi );
//         val += it.val * diff1 * diff2;
// 
//     }
//     return val;
// 
// }

/*
 *   QRsumQPcorr:  QRsum + QPcorr
 *   ==================================
 *   We're using:  plus1 = r+s  minus1 = r-s  plus2 = r'+s  minus2 = r'-s
 */
double redux::speckle::QRsumQPcorr( double plus1_rho, double plus1_phi, double minus1_rho, double minus1_phi,
                                    double plus2_rho, double plus2_phi, double minus2_rho, double minus2_phi ) {

    double val = 0.0;
    for( auto& it: KLweights ) {
        double diff1 = zernikeDiffPolar( it.i1, it.m1, plus1_rho, plus1_phi, minus1_rho, minus1_phi );
        double diff2 = zernikeDiffPolar( it.i2, it.m2, plus2_rho, plus2_phi, minus2_rho, minus2_phi );
        val += it.val * diff1 * diff2;

    }
    return val;

}


/*
 *   The integrand of the Wang & Markey SE transfer function
 *   =======================================================
 *   TODO: optimize. e.g. do the vector addition directly in polar coordinates to avoid transforming to/from cartesian.
 */
double redux::speckle::wam( double *k, size_t dim, void *params ) {

    if( dim < 2 || k[1] > 1.0 ) return 0.0;

    parameters *p = reinterpret_cast<parameters*>(params);

    complex_t r1 = polar<double>(k[1],k[0]);

    complex_t tmp = r1+p->delta;
    double plus1_rho = abs(tmp);
    if( plus1_rho>1 ) return 0.0;
    double plus1_phi = arg(tmp);

    tmp = r1-p->delta;
    double minus1_rho = abs(tmp);
    if( minus1_rho>1 ) return 0.0;
    double minus1_phi = arg(tmp);

    double val = -3.44 * p->q_five_thirds;

    if( !KLweights.empty() ) {
        val += QRsumQPcorr( plus1_rho, plus1_phi, minus1_rho, minus1_phi,
                            plus1_rho, plus1_phi, minus1_rho, minus1_phi );
    }
    val = exp( val*p->alpha_five_thirds );

    return k[1] * val * M_1_PI;

}


/*
 *   The integrand of the Wang & Markey SE transfer function
 *   =======================================================
 *   TODO: optimize. e.g. do the vector addition directly in polar coordinates to avoid transforming to/from cartesian.
 */
double redux::speckle::wam2( double *k, size_t dim, void *params ) {

    if( dim < 2 || k[1] > 1.0 || k[3] > 1.0 ) return 0.0;

    parameters *p = reinterpret_cast<parameters*>(params);
    const complex_t& delta = p->delta;

    complex_t r1 = polar<double>(k[1],k[0]);        // quickfix: using complex type to get abs/arg-functions as well as cartesian vector addition.
    complex_t r2 = polar<double>(k[3],k[2]);

    complex_t tmp = r1+delta;
    double plus1_rho = abs(tmp);
    if( plus1_rho>1 ) return 0.0;
    double plus1_phi = arg(tmp);

    tmp = r1-delta;
    double minus1_rho = abs(tmp);
    if( minus1_rho>1 ) return 0.0;
    double minus1_phi = arg(tmp);

    tmp = r2+delta;
    double plus2_rho = abs(tmp);
    if( plus2_rho>1 ) return 0.0;
    double plus2_phi = arg(tmp);

    tmp = r2-delta;
    double minus2_rho = abs(tmp);
    if( minus2_rho>1 ) return 0.0;
    double minus2_phi = arg(tmp);

    tmp = (r1-r2)/2.0;

    double val = 3.44*(pow( abs(tmp+delta), five_thirds )+pow( abs(tmp-delta), five_thirds ));
    val -= 6.88*(p->q_five_thirds + pow( abs(tmp), five_thirds ));

    if( !KLweights.empty() ) {
        val += QRsumQPcorr( plus1_rho, plus1_phi, minus1_rho, minus1_phi,
                            plus1_rho, plus1_phi, minus1_rho, minus1_phi );
        val += QRsumQPcorr( plus2_rho, plus2_phi, minus2_rho, minus2_phi,
                            plus2_rho, plus2_phi, minus2_rho, minus2_phi );
        val -= QRsumQPcorr( plus1_rho, plus1_phi, minus1_rho, minus1_phi,
                            plus2_rho, plus2_phi, minus2_rho, minus2_phi );
        val -= QRsumQPcorr( plus2_rho, plus2_phi, minus2_rho, minus2_phi,
                            plus1_rho, plus1_phi, minus1_rho, minus1_phi );
    }

    val *= p->alpha_five_thirds;

    return k[1] * k[3] * exp( val ) * M_1_PI * M_1_PI;

}
