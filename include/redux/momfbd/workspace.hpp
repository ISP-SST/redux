#ifndef REDUX_MOMFBD_WORKSPACE_HPP
#define REDUX_MOMFBD_WORKSPACE_HPP


#include "redux/momfbd/data.hpp"
#include "redux/momfbd/result.hpp"
#include "redux/momfbd/wavefront.hpp"

#include <memory>


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
         
        class MomfbdJob;
        
         /*! Container used during processing. Basically temporary arrays and reorganized references to original data.
         */
        struct WorkSpace : public std::enable_shared_from_this<WorkSpace> {
            
            typedef std::shared_ptr<WorkSpace> Ptr;
            
            const std::vector<std::shared_ptr<Object>>& objects;
            std::map<uint32_t, std::shared_ptr<WaveFront>> wavefronts;              //!< Constrained groups, using imageNumber/wf_num as identifier.
            
            WorkSpace(const redux::momfbd::MomfbdJob&);
            ~WorkSpace();
            
            void init(void);
            
            void run(PatchData::Ptr, boost::asio::io_service&, uint8_t nThreads);
            
            double coefficientMetric(boost::asio::io_service&);
            double objectMetric(boost::asio::io_service&);
           
            
            void clear(void);
            void resetAlpha(void);
            PatchResult::Ptr& getResult(void);
            
            PatchData::Ptr data;
            PatchResult::Ptr result;
            
            const MomfbdJob& job;
            
            redux::util::Array<double> window, noiseWindow;
            size_t nFreeParameters;

            std::shared_ptr<double> alpha,grad_alpha,saved_alpha;
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_WORKSPACE_HPP
