#ifndef REDUX_MOMFBD_MODES_HPP
#define REDUX_MOMFBD_MODES_HPP

#include "redux/util/array.hpp"

#define MT_NONE   0
#define MT_PUPIL  1
#define MT_OTF_BF 2

namespace redux {

    namespace momfbd {

        struct PupilMode : public redux::util::Array<double> {
            
            // Karhunen-Loeve expansion coefficients
            struct KL {
                std::vector< std::pair<uint16_t, double> > zernikeWeights; //< zernike mode numbers and corresponding weights
                double covariance;
            };
            typedef std::shared_ptr<PupilMode> Ptr;
            typedef std::shared_ptr<KL> KLPtr;

            PupilMode() : atm_rms(0) {};
            PupilMode ( uint16_t modeNumber, uint16_t nPoints, double r_c = 1.0, double angle = 0.0 ); // Zernike
            PupilMode ( uint16_t firstMode, uint16_t lastMode, uint16_t klModeNumber, uint16_t nPoints, double r_c = 1.0, double angle = 0.0, double cutoff=0.0 ); // KL

            operator const redux::util::Array<double>&() const { return reinterpret_cast<const redux::util::Array<double>&>(*this); }
            double atm_rms;                         //!< = sqrt(covariance), used in metric computations.
        };

    }

}

#endif      // REDUX_MOMFBD_MODES_HPP
