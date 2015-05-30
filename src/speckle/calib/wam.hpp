#ifndef REDUX_SPECKLE_CALIB_WAM_HPP
#define REDUX_SPECKLE_CALIB_WAM_HPP

#include <complex>

namespace redux {
    
    namespace speckle {
        
        typedef std::complex<double> complex_t;
        struct parameters { 
            double alpha_five_thirds;
            double q_five_thirds;
            complex_t delta;
        };


        double QRsumQPcorr(const double& plus1_rho, const double& plus1_phi, const double& minus1_rho, const double& minus1_phi,
                           const double& plus2_rho, const double& plus2_phi, const double& minus2_rho, const double& minus2_phi);
        
        double wam(double *k, size_t dim, void *para);
        double wam2(double *k, size_t dim, void *params);


    } // speckle
    
} // redux

#endif  // REDUX_SPECKLE_CALIB_WAM_HPP
