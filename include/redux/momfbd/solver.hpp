#ifndef REDUX_MOMFBD_SOLVER_HPP
#define REDUX_MOMFBD_SOLVER_HPP


#include "redux/momfbd/data.hpp"
#include "redux/momfbd/subimage.hpp"

#include "redux/util/gsl.hpp"
#include "redux/util/stopwatch.hpp"

#include <memory>


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
         
        class MomfbdJob;
        class SubImage;
         /*! Container used during processing. Basically temporary arrays and reorganized references to original data.
         */
        struct Solver {
            
            typedef std::shared_ptr<Solver> Ptr;
            
            Solver(const redux::momfbd::MomfbdJob&, boost::asio::io_service&, uint16_t nThreads);
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
            
            void applyAlpha(void);
            void applyBeta( const gsl_vector* beta );
            void applyBeta( const gsl_vector* beta, double scale );
            void applyBeta( double scale );
            
            double metric(void);
            double metricAt(double step);       // evaluate metric at alpha + step*grad
            void calcPQ(void);
            void gradient(void);
            void gradient(gsl_vector* out);
          
            void clear(void);
            
            void dump( std::string tag );
            
            const MomfbdJob& job;
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
            
            double max_wavelength;
            
            grad_t gradientMethod;

            redux::util::StopWatch timer;
            
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_SOLVER_HPP
