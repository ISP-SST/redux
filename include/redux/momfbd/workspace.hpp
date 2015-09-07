#ifndef REDUX_MOMFBD_WORKSPACE_HPP
#define REDUX_MOMFBD_WORKSPACE_HPP


#include "redux/momfbd/data.hpp"
#include "redux/momfbd/tilts.hpp"
#include "redux/momfbd/wavefront.hpp"

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
        struct WorkSpace : public std::enable_shared_from_this<WorkSpace> {
            
            typedef std::shared_ptr<WorkSpace> Ptr;

            const std::vector<std::shared_ptr<Object>>& objects;
            std::map<uint32_t, std::shared_ptr<WaveFront>> wavefronts;              //!< Constrained groups, using imageNumber/wf_num as identifier.
            std::map<uint16_t, std::shared_ptr<Tilts>> tilts;                        //!< Constrained tilts.
            
            WorkSpace(const redux::momfbd::MomfbdJob&);
            ~WorkSpace();
            
            void init(void);
            
            void getMetric(boost::asio::io_service&, uint8_t nThreads);
            void resetAll( boost::asio::io_service& );
            void dumpImages( boost::asio::io_service&, std::string );
            double my_test( boost::asio::io_service&, const gsl_vector*, void*, gsl_vector*, std::string );
            double my_f( boost::asio::io_service&, const gsl_vector*, void* );
            void my_df( boost::asio::io_service&, const gsl_vector*, void*, gsl_vector* );
            void my_fdf( boost::asio::io_service&, const gsl_vector*, void*, double*, gsl_vector* );
            void run(PatchData::Ptr, boost::asio::io_service&, uint8_t nThreads);
            
            double coefficientMetric(boost::asio::io_service&);
            void calcOTFs(boost::asio::io_service&);
            double objectMetric(boost::asio::io_service&);
           
            
            void clear(void);
            void resetAlpha(void);
           // PatchResult::Ptr& getResult(void);
            
            void dump( std::string tag );
            PatchData::Ptr data;
           // PatchResult::Ptr result;
            
            const MomfbdJob& job;
            
            redux::util::Array<double> window, noiseWindow;
            size_t nFreeParameters;
            std::shared_ptr<double> alpha, grad_alpha,saved_alpha;
            gsl_vector *beta, *grad_beta, *beta_init;
            gsl_vector alpha_init;
            
            grad_t gradient;
            
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_WORKSPACE_HPP
