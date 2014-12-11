#ifndef REDUX_MOMFBD_MODES_HPP
#define REDUX_MOMFBD_MODES_HPP

#include "redux/util/array.hpp"

#define MT_NONE   0
#define MT_PUPIL  1
#define MT_OTF_BF 2

namespace redux {

    namespace momfbd {

/*        class Modes {

        public:
            


*/
            struct PupilMode : public redux::util::Array<double> {

                typedef std::shared_ptr<PupilMode> Ptr;
                
                enum Basis { MB_NONE=0, MB_ZERNIKE, MB_KARHUNEN_LOEVE };

                struct KL_cfg {
                    // zernike mode numbers and corresponding weights
                    std::vector< std::pair<int, double> > zernikeWeights; 
                    double covariance;
                };

                PupilMode() {};
                PupilMode ( int modeNumber, int nPoints, double r_c = 1.0, double lambda = 1.0, double angle = 0.0 ); // Zernike
                PupilMode ( int firstMode, int lastMode, int klModeNumber, int nPoints, double r_c = 1.0, double lambda = 1.0, double angle = 0.0, double cutoff=0.0 ); // KL

                Basis base;
            };

/*
            Modes();
            Modes ( KL_cfg*, double, double, int, int, int, int*, int, int*, int**, int**, int, int, double, double, double** );
            virtual ~Modes ( void );

            void init ( KL_cfg*, double, double, int, int, int, int*, int, int*, int**, int**, int, int, double, double, double** );

        private:
            PupilMode*** mode;

            std::map<int,PupilMode::Ptr> zernikeModes, klModes, defaultModes;

            double **pupil, area;
            double **covar;
            int32_t nph, zmin, zmax, kmin, kmax, first_mode, last_mode;
        };
*/
    }

}

#endif      // REDUX_MOMFBD_MODES_HPP
