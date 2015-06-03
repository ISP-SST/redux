#ifndef REDUX_SPECKLE_CALIB_MPIWRAPPER_HPP
#define REDUX_SPECKLE_CALIB_MPIWRAPPER_HPP

#include "zernike.hpp"

#include <memory>
#include <mutex>
#include <ostream>
#include <vector>

#include <mpi.h>

namespace redux {

    namespace speckle {

        class MpiWrapper {
            
            public:
                MpiWrapper (int argc, char *argv[]);
                ~MpiWrapper();
                
                bool next(const double*, double *&);
                void run(void);
                
                void writeResults(std::ostream&) const;

            private:
                
                void makeCube(void);

                int myRank;                    // Rank of process
                int worldSize;                 // Number of available processes
                MPI_Status status;             // status variable for communication
                
                size_t nCalls;
                size_t nIterations;
                bool running;
                
                std::vector<int> nThreads;
                std::vector<int> blockSizes;
                std::vector<int> blockOffsets;
                int nTotalThreads;

                std::vector<float> efficiencies;

                std::vector<double> alpha, freq;
                std::shared_ptr<double> data;  // [alpha][q][4]  input:  [0]=pow( freq, 5.0/3.0 ) ), [1]=pow( 1.0 / alpha, 5.0/3.0 )
                                                              //         [2]=delta_x                 [3]= delta_y
                                                              // result: [0]=ltf                     [1]=ltf_err
                                                              //         [2]=stf                     [3]=stf_Err
                std::mutex mtx;
                size_t currentOffset;
                ZernikeData& zd;
        };


    } // speckle

} // redux

#endif  // REDUX_SPECKLE_CALIB_MPIWRAPPER_HPP
