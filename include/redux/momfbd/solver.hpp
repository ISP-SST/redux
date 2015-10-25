#ifndef REDUX_MOMFBD_SOLVER_HPP
#define REDUX_MOMFBD_SOLVER_HPP


#include "redux/momfbd/data.hpp"
#include "redux/momfbd/subimage.hpp"

#include "redux/util/gsl.hpp"

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

            const std::vector<std::shared_ptr<Object>>& objects;
            
            Solver(const redux::momfbd::MomfbdJob&);
            ~Solver();
            
            void init(void);
            
            void getMetric(boost::asio::io_service&, uint8_t nThreads);
            void resetAll( boost::asio::io_service& );
            void dumpImages( boost::asio::io_service&, std::string );
            double my_test( boost::asio::io_service&, const gsl_vector*, void*, gsl_vector*, std::string );
            double my_f( boost::asio::io_service&, const gsl_vector*, void* );
            void my_df( boost::asio::io_service&, const gsl_vector*, void*, gsl_vector* );
            void my_fdf( boost::asio::io_service&, const gsl_vector*, void*, double*, gsl_vector* );
            void run(PatchData::Ptr, boost::asio::io_service&, uint16_t nThreads);
            

            double objectMetric(boost::asio::io_service&);
            
            void clear(void);
            void getAlpha(void);
            
            void dump( std::string tag );
            PatchData::Ptr data;
            
            const MomfbdJob& job;
            
            redux::util::Array<double> window, noiseWindow;
            
            uint16_t nModes;
            uint16_t nThreads;
            uint32_t nParameters;
            uint32_t nFreeParameters;
            uint32_t nTotalImages;
            
            uint16_t *modeNumbers;
            uint16_t *enabledModes;
            double *alpha,*grad_alpha;
            
            grad_t gradientMethod;
            
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_SOLVER_HPP
