#include "wam.hpp"

#include "zernike.hpp"

#include <gsl/gsl_math.h>   // for M_PI and other constants

using namespace std;


namespace {     // Anonymuous namespace -> only visible from this compilation unit (.cpp file)

    const double five_thirds( 5.0 / 3.0 );

}


double redux::speckle::QRsumQPcorr( const double& plus1_rho, const double& plus1_phi, const double& minus1_rho, const double& minus1_phi,
                                    const double& plus2_rho, const double& plus2_phi, const double& minus2_rho, const double& minus2_phi ) {
    
    static ZernikeData& zd = ZernikeData::get();
    double val = 0.0;
    auto it = zd.getCovariances().begin();
    auto end = zd.getCovariances().end();
    for (; it != end; ++it) {
        if (it->data1 && it->data2) {
            double tmp = zernikePolar(it->data1.get(), it->sampling1, it->m, plus1_rho, plus1_phi)
                       - zernikePolar(it->data1.get(), it->sampling1, it->m, minus1_rho, minus1_phi);
            tmp *= zernikePolar(it->data2.get(), it->sampling2, it->m, plus2_rho, plus2_phi)
                 - zernikePolar(it->data2.get(), it->sampling2, it->m, minus2_rho, minus2_phi);
            val += it->cov * tmp;
        }
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

    static ZernikeData& zd = ZernikeData::get();
    parameters *p = reinterpret_cast<parameters*>( params );

    complex_t r1 = polar<double>( k[1], k[0] );

    complex_t tmp = r1 + p->delta;
    double plus1_rho = abs( tmp );
    if( plus1_rho > 1 ) return 0.0;
    double plus1_phi = arg( tmp );

    tmp = r1 - p->delta;
    double minus1_rho = abs( tmp );
    if( minus1_rho > 1 ) return 0.0;
    double minus1_phi = arg( tmp );

    double val = -3.44 * p->q_five_thirds;

    if( !zd.empty() ) {
        val += QRsumQPcorr( plus1_rho, plus1_phi, minus1_rho, minus1_phi,
                               plus1_rho, plus1_phi, minus1_rho, minus1_phi );
    }

    val = exp( val * p->alpha_five_thirds );

    return k[1] * val * M_1_PI;

}


/*
 *   The integrand of the Wang & Markey SE transfer function
 *   =======================================================
 *   TODO: optimize. e.g. do the vector addition directly in polar coordinates to avoid transforming to/from cartesian.
 */
double redux::speckle::wam2( double *k, size_t dim, void *params ) {

    if( dim < 2 || k[1] > 1.0 || k[3] > 1.0 ) return 0.0;

    static ZernikeData& zd = ZernikeData::get();
    parameters *p = reinterpret_cast<parameters*>( params );
    const complex_t& delta = p->delta;

    complex_t r1 = polar<double>( k[1], k[0] );     // quickfix: using complex type to get abs/arg-functions as well as cartesian vector addition.
    complex_t r2 = polar<double>( k[3], k[2] );

    complex_t tmp = r1 + delta;
    double plus1_rho = abs( tmp );
    if( plus1_rho > 1 ) return 0.0;
    double plus1_phi = arg( tmp );

    tmp = r1 - delta;
    double minus1_rho = abs( tmp );
    if( minus1_rho > 1 ) return 0.0;
    double minus1_phi = arg( tmp );

    tmp = r2 + delta;
    double plus2_rho = abs( tmp );
    if( plus2_rho > 1 ) return 0.0;
    double plus2_phi = arg( tmp );

    tmp = r2 - delta;
    double minus2_rho = abs( tmp );

    if( minus2_rho > 1 ) return 0.0;

    double minus2_phi = arg( tmp );

    tmp = ( r1 - r2 ) / 2.0;

    double val = 3.44 * ( pow( abs( tmp + delta ), five_thirds ) + pow( abs( tmp - delta ), five_thirds ) );
    val -= 6.88 * ( p->q_five_thirds + pow( abs( tmp ), five_thirds ) );

    if( !zd.empty() ) {
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
