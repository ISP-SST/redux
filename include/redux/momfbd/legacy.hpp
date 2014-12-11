#ifndef REDUX_MOMFBD_LEGACY_HPP
#define REDUX_MOMFBD_LEGACY_HPP

#include "redux/momfbd/modes.hpp"

#include <cmath>
#include <map>

namespace redux {

    namespace momfbd {
        
        namespace legacy {

            // limbo for old functions that are kept until replaced


            double gammln( const double& xx, double& sign );

            std::map<int, int> mappingVector( double func( int, int ), int, int );
            std::map<int, int> reverseMappingVector( std::map<int, int>& );

            double*** reorderedMatrix( double func( int, int ), const std::map<int, int>&, int, int, int &, int *&, int *& );
            void blockwiseSVD( double ***blockMatrix, double *singular_values, int nBlocks, int *blockFirst, int *blockLast );

            PupilMode::KL_cfg* klCoefficients(double ***blockMatrix,double *values,int *first,int *last,const std::map<int, int>& map,const std::map<int, int>& rmap,int nl,int nh);
            PupilMode::KL_cfg* klConfig(int kl_min_mode,int kl_max_mode);
            
            void svdcmp( double * const * const, int, int, double *, double * const * const );

            double **qr(double **c, int m, int n, uint8_t fast_QR);
 

        }

    }

}

#endif // REDUX_MOMFBD_LEGACY_HPP
