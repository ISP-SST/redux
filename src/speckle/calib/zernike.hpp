#ifndef REDUX_SPECKLE_CALIB_ZERNIKE_HPP
#define REDUX_SPECKLE_CALIB_ZERNIKE_HPP

#include "defs.hpp"

#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace redux {

    namespace speckle {

        class ZernikeData {
            
            struct cov_t {
                cov_t(uint32_t n1, uint32_t n2, int32_t m, double cov);
                void swap(void);
                bool operator<(const cov_t&) const;
                uint32_t n1;                            // radial order of first Zernike
                uint32_t n2;                            // radial order of second Zernike
                int32_t m;                              // azimuthal order of Zernikes (must be same or covariance is zero)
                mutable double cov;                     // "effective covariance", i.e. weighted with the efficiencies passed to "init()".  

                // for convenience in the loops, we keep pointers to the data here too.
                mutable std::shared_ptr<double> data1, data2;    // pointers to radial polynomials for the zernikes.
                uint32_t sampling1, sampling2;           // allow for different samplings for the different polynomials (not used yet)
            };

        public:
            bool empty(void) const { return covariances.empty(); }
            static ZernikeData& get(void);
            void init(const std::vector<float>& eff, uint32_t Imax=SPECKLE_IMAX, float cov_cutoff=SPECKLE_COV_CUTOFF);
            std::vector<float> init(const std::string& filename, uint32_t Imax=SPECKLE_IMAX, float cov_cutoff=SPECKLE_COV_CUTOFF);
            std::shared_ptr<double> getRadialData(uint32_t n, int32_t abs_m, uint32_t nPoints=0);
            std::shared_ptr<double> getRadialData(uint32_t i);
            uint32_t getSampling(uint32_t n, int32_t m);
            uint32_t getSampling(uint32_t i);
            const std::set<cov_t>& getCovariances(void) const { return covariances; };

            double zernike (int i, double x, double y);
            double zernikePolar (int i, double rho, double phi);
            double zernikeDiffPolar(int i, const double& rho1, const double& phi1, const double& rho2, const double& phi2);
            
        private:
            ZernikeData (void) {};      // private constructor, only 1 instance possible (get reference by calling ZernikeData::get())
            ZernikeData (const ZernikeData&) {};
            
            
            struct data_index {
                data_index(uint32_t n, uint32_t m, uint32_t sampling=0) : n(n), abs_m(m), sampling(sampling) {};
                bool operator<(const data_index&) const;
                uint32_t n, abs_m, sampling;
            };
            
            std::set<cov_t> covariances;         
            std::map<data_index, std::shared_ptr<double>> radialPolynomials;

        };

        void noll_to_nm (int i, int& n, int& m);
        double calcZernikeCovariance (uint32_t i, uint32_t j);
        void calcRadialZernike(double* out, uint32_t nPoints, uint16_t n, uint16_t abs_m);
        
        double zernike (const double* rdata, uint32_t sampling, int32_t m, double x, double y);
        double zernikePolar (const double* rdata, uint32_t sampling, int32_t m, double rho, double phi);
        double zernikeDiffPolar (const double* rdata, uint32_t sampling, int32_t m, const double& rho1, const double& phi1, const double& rho2, const double& phi2);

    } // speckle

} // redux

#endif  // REDUX_SPECKLE_CALIB_ZERNIKE_HPP
