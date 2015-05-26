#ifndef REDUX_SPECKLE_CALIB_ZERNIKE_HPP
#define REDUX_SPECKLE_CALIB_ZERNIKE_HPP

namespace redux {
    
    namespace speckle {
        
        void noll_to_nm(int i, int& n, int& m);
        double calcZernikeCovariance(int i, int j);
        void init_zernike(void);
        double zernrad(double r, int n, int l);
        double zernike(int i, int m, double x, double y);
        double zernikePolar(int i, int m, double rho, double phi);
        double zernikeDiffPolar(int i, int m, double rho1, double phi1, double rho2, double phi2);

    } // speckle
    
} // redux

#endif  // REDUX_SPECKLE_CALIB_ZERNIKE_HPP
