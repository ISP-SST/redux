#ifndef REDUX_MOMFBD_SOLVER_HPP
#define REDUX_MOMFBD_SOLVER_HPP


#include "redux/momfbd/data.hpp"
#include "redux/momfbd/subimage.hpp"
#include "redux/momfbd/wavefront.hpp"

#include "redux/util/gsl.hpp"
#include "redux/util/progresswatch.hpp"
#include "redux/util/stopwatch.hpp"

#include <memory>


namespace redux {
    
    namespace logging {
        class Logger;
    }

    namespace network {
        struct Host;
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
                        } else {
                            D.clear();
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
                            FT.resize( otfSize, otfSize );
                        } else {
                            D2.clear();
                            C.clear();
                            C2.clear();
                            OTF.clear();
                            FT.clear();
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
        struct SubImage;
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
            
            template <typename T>
            void shiftAndInit( const T* a, bool doReset=false );
            void shiftAndInit( bool doReset=false ){ shiftAndInit( alpha.get(), doReset ); }
            
            void alignWavefronts( void );
            void zeroAlphas( void );
            
            template <typename T> void applyAlpha( T* a );
            inline void applyAlpha(void) { applyAlpha( alpha.get() ); } ;
            void applyBeta( const gsl_vector* beta );
            void applyBeta( const gsl_vector* beta, double scale );
            void applyBeta( double scale );
            
            void applyConstraints( const double* a, double* b );
            void reverseConstraints( const double* b, double* a );
            
            void zeroAvgTilt( double* a, int m );
            void zeroAvgTilts( double* a, int m1, int m2 );

            void loadInit( const PatchData::Ptr pd, double* a) const;
            void initImages( double* a );
            
            double metric(void);
            double metricAt(double step);       // evaluate metric at alpha + step*grad
            void calcPQ(void);
            void gradient(void);
            void gradient(gsl_vector* out);
          
            void clear(void);
            
            void dump( std::string tag );
            
            MomfbdJob& job;
            network::Host& myInfo;
            logging::Logger& logger;
            std::map<uint32_t,std::shared_ptr<WaveFront>> wavefronts;
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
            
            std::shared_ptr<bool> enabledModes;
            std::shared_ptr<double> alpha, alpha_offset, grad_alpha;

            double *tmp_alpha;
            double *beta, *grad_beta, *search_dir, *tmp_beta;
            double grad_beta_norm;
            double *regAlphaWeights;
            
            double max_wavelength;
            size_t nTotalPixels;
            
            grad_t gradientMethod;

            redux::util::StopWatch timer;
            redux::util::ProgressWatch progWatch;
            
            std::vector< thread::TmpStorage > tmps;
            static thread_local thread::TmpStorage* tmp;
            std::mutex mtx;
            
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_SOLVER_HPP
