#ifndef REDUX_MOMFBD_SOLVER_HPP
#define REDUX_MOMFBD_SOLVER_HPP


#include "redux/momfbd/data.hpp"
#include "redux/momfbd/subimage.hpp"

#include "redux/util/gsl.hpp"
#include "redux/util/progresswatch.hpp"
#include "redux/util/stopwatch.hpp"

#include <memory>


namespace redux {
    
    namespace logging {
        class Logger;
    }

    namespace network {
        class Host;
    }

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
        
        namespace thread {
        
            struct TmpStorage {
                TmpStorage() : patchSize(0), pupilSize(0) {}
                void resize( uint16_t patchSz, uint16_t pupSz ) {
                    if( patchSz != patchSize ) {
                        patchSize = patchSz;
                        if( patchSize ) {
                            D.resize( patchSize, patchSize );
                            FT.resize( patchSize, patchSize );
                        } else {
                            D.clear();
                            FT.clear();
                        }
                    }
                    if( pupSz != pupilSize ) {
                        pupilSize = pupSz;
                        size_t otfSize = 2*pupilSize;
                        if( otfSize ) {
                            D2.resize( otfSize, otfSize );
                            C.resize( otfSize, otfSize );
                            C2.resize( otfSize, otfSize );
                            OTF.resize( otfSize, otfSize, redux::image::FT_FULLCOMPLEX );
                        } else {
                            D2.clear();
                            C.clear();
                            C2.clear();
                            OTF.clear();
                        }
                    }
                }
                uint16_t patchSize, pupilSize;
                redux::util::Array<double> D,D2;
                redux::util::Array<complex_t> C,C2;
                redux::image::FourierTransform FT,OTF;
            };
    
        }

        class MomfbdJob;
        class SubImage;
         /*! Container used during processing. Basically temporary arrays and reorganized references to original data.
         */
        struct Solver {
            
            typedef std::shared_ptr<Solver> Ptr;
            
            Solver(redux::momfbd::MomfbdJob&, boost::asio::io_service&, uint16_t nThreads);
            ~Solver();
            
            void init(void);
            
            void getMetric(boost::asio::io_service&, uint8_t nThreads);
            void reset(void);
            void dumpImages( boost::asio::io_service&, std::string );
            
            double my_f( const gsl_vector*, void* );
            void my_df( const gsl_vector*, void*, gsl_vector* );
            void my_fdf( const gsl_vector*, void*, double*, gsl_vector* );
            void my_precalc( const gsl_vector*, const gsl_vector* );
            
            void run(PatchData::Ptr);
            
            void shiftAndInit( double* a, bool doReset=false );
            void applyAlpha( double* a );
            inline void applyAlpha(void) { applyAlpha( alpha ); } ;
            void applyBeta( const gsl_vector* beta );
            void applyBeta( const gsl_vector* beta, double scale );
            void applyBeta( double scale );
            
            void applyConstraints( const double* alpha, double* beta );
            void reverseConstraints( const double* beta, double* alpha );

            void loadInit( const PatchData::Ptr pd, double* a) const;
            void initImages( double* a );
            
            double metric(void);
            double metric2(void);
            double metricAt(double step);       // evaluate metric at alpha + step*grad
            void calcPQ(void);
            void calcPQ2(void);
            void gradient(void);
            void gradient2(void);
            void gradient(gsl_vector* out);
            void gradient2(gsl_vector* out);
          
            void clear(void);
            
            void dump( std::string tag );
            
            MomfbdJob& job;
            network::Host& myInfo;
            logging::Logger& logger;
            const std::vector<std::shared_ptr<Object>>& objects;
            boost::asio::io_service& service;
            
            redux::util::Array<double> window, noiseWindow;
            redux::util::Array<double> tmpPhi, tmpPhiGrad;
            redux::util::Array<complex_t> tmpOTF;
            
            uint16_t patchSize;
            uint16_t pupilSize;
            uint16_t nModes;
            uint16_t nThreads;
            uint32_t nParameters;
            uint32_t nFreeParameters;
            uint32_t nTotalImages;
            
            uint16_t *modeNumbers;
            uint16_t *enabledModes;

            double *alpha, *alpha_offset, *grad_alpha, *tmp_alpha;
            double *beta, *grad_beta, *search_dir, *tmp_beta;
            double grad_beta_norm;
            double *regAlphaWeights;
            
            double max_wavelength;
            
            grad_t gradientMethod;

            redux::util::StopWatch timer;
            redux::util::ProgressWatch progWatch;
            
            static thread_local thread::TmpStorage tmp;
            std::mutex mtx;
            
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_SOLVER_HPP
