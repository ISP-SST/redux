#ifndef REDUX_SPECKLE_CALIB_WAM_HPP
#define REDUX_SPECKLE_CALIB_WAM_HPP

#include <complex>
#include <vector>

namespace redux {
    
    namespace speckle {
        
        typedef std::complex<double> complex_t;
        struct parameters { 
            double alpha_five_thirds;
            double q_five_thirds;
            complex_t delta;
        };

        std::vector<float>& init_wam_master(char *file);
        void init_wam_slave(const std::vector<float>& data);

        double QRsum(double plus1_rho, double plus1_phi, double minus1_rho, double minus1_phi,
                      double plus2_rho, double plus2_phi, double minus2_rho, double minus2_phi);
        double QPcorr(double plus1_rho, double plus1_phi, double minus1_rho, double minus1_phi,
                       double plus2_rho, double plus2_phi, double minus2_rho, double minus2_phi);
        double QRsumQPcorr(const double& plus1_rho, const double& plus1_phi, const double& minus1_rho, const double& minus1_phi,
                           const double& plus2_rho, const double& plus2_phi, const double& minus2_rho, const double& minus2_phi);
        
        double wam(double *k, size_t dim, void *para);
        double wam2(double *k, size_t dim, void *params);


    } // speckle
    
} // redux

#endif  // REDUX_SPECKLE_CALIB_WAM_HPP
