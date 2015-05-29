#ifndef REDUX_SPECKLE_CALIB_ZERNIKE_HPP
#define REDUX_SPECKLE_CALIB_ZERNIKE_HPP

#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <vector>

namespace redux {

    namespace speckle {

        class ZernikeData {
            
            struct cov_t {
                cov_t(uint16_t n1, uint16_t n2, int32_t m, double cov);
                bool operator<(const cov_t&) const;
                uint16_t n1;                            // radial order of first Zernike
                uint16_t n2;                            // radial order of second Zernike
                int32_t m;                              // azimuthal order of Zernikes (must be same or covariance is zero)
                mutable double cov;                     // "effective covariance", i.e. weighted with the efficiencies passed to "init()".  

                // for convenience in the loops, we keep pointers to the data here too.
                mutable std::shared_ptr<double> data1,data2;    // pointers to radial polynomials for the zernikes.
                uint32_t sampling1, sampling2;           // allow for different samplings for the different polynomials (not used yet)
            };

        public:
            static ZernikeData& get(void);
            void init(const std::vector<float>& eff, uint16_t Imax=231, float cov_cutoff=0.0);
            std::shared_ptr<double>& radialData(uint16_t n, uint16_t abs_m, uint32_t nPoints);
            
            double QRsumQPcorr( const double& plus1_rho, const double& plus1_phi, const double& minus1_rho, const double& minus1_phi,
                                const double& plus2_rho, const double& plus2_phi, const double& minus2_rho, const double& minus2_phi ) const;

        private:
            ZernikeData (void) {};      // private constructor, only 1 instance possible (get reference by calling ZernikeData::get())
            ZernikeData (const ZernikeData&) {};
            
            
            struct data_index {
                bool operator<(const data_index&) const;
                uint16_t n, abs_m;
                uint32_t sampling;
            };
            
            std::mutex mtx;
            std::set<cov_t> covariances;         
            std::map<data_index, std::shared_ptr<double>> radialPolynomials;

        };

        void noll_to_nm (int i, int& n, int& m);
        double calcZernikeCovariance (int i, int j);
        void calcRadialZernike(double* out, uint32_t nPoints, uint16_t n, uint16_t abs_m);
        void init_zernike (void);
        double zernrad (double r, int n, int l);
        void zernrad (const double* r, double* out, int count, int n, int l);
        void zernrad2 (const double* r, double* out, int count, int n, int l);
        double zernike (int i, int m, double x, double y);
        double zernikePolar (int i, int m, double rho, double phi);
        double zernikeDiffPolar (int i, int m, const double& rho1, const double& phi1, const double& rho2, const double& phi2);
        double zernikePolar (const double* rdata, uint32_t sampling, int32_t m, double rho, double phi);
        double zernikeDiffPolar (const double* rdata, uint32_t sampling, int32_t m, const double& rho1, const double& phi1, const double& rho2, const double& phi2);

    } // speckle

} // redux

#endif  // REDUX_SPECKLE_CALIB_ZERNIKE_HPP
