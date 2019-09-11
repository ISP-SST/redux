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
                TmpStorage() : m_patchSize(0), m_pupilSize(0) {}
                TmpStorage( const TmpStorage& ) = delete;
                TmpStorage( TmpStorage&& ) = delete;
                ~TmpStorage() {}
                static void setSize( uint16_t patchSz, uint16_t pupSz ) { patchSize=patchSz; pupilSize=pupSz; }
                void resize( void ) {
                    using namespace redux::util;
                    using namespace std;
                    size_t otfSize = 2*pupilSize;
                    size_t tmpSize = max<size_t>(patchSize,otfSize);
                    size_t tmpSize2 = tmpSize*tmpSize;
                    if( m_patchSize != patchSize ) {
                        m_patchSize = patchSize;
                        if( patchSize ) {
                            D.reset( new double[tmpSize2] );
                        } else {
                            D.reset();
                        }
                    }
                    if( m_pupilSize != pupilSize ) {
                        m_pupilSize = pupilSize;
                        if( tmpSize ) {
                            D2.reset( new double[tmpSize2] );
                            C.resize( tmpSize, tmpSize );
                            C2.resize( tmpSize, tmpSize );
                            OTF.resize( tmpSize, tmpSize, redux::image::FT_FULLCOMPLEX );
                            FT.resize( tmpSize, tmpSize );
                        } else {
                            D2.reset();
                            C.clear();
                            C2.clear();
                            OTF.clear();
                            FT.clear();
                        }
                    }
                }
                uint16_t m_patchSize, m_pupilSize;
                static uint16_t patchSize, pupilSize;
                std::unique_ptr<double[]> D,D2;
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
            uint16_t maxThreads;
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
            
            static thread::TmpStorage* tmp(void);
            
            std::mutex mtx;
            
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_SOLVER_HPP
