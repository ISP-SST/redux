#ifndef REDUX_MOMFBD_WORKSPACE_HPP
#define REDUX_MOMFBD_WORKSPACE_HPP


#include "redux/momfbd/data.hpp"
#include "redux/momfbd/patch.hpp"

#include "redux/momfbd/momfbdjob.hpp"

#include <memory>


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
            
         /*! Container used during processing. Basically temporary arrays and reorganized references to original data.
         */
        struct WorkSpace {
            
            //typedef std::shared_ptr<WaveFrontData> WfPtr;
            //typedef std::shared_ptr<ChannelData> ChPtr;
            //typedef std::shared_ptr<ObjectData> ObjPtr;

            
            std::vector<ObjPtr> objects;
            std::map<uint32_t, WfPtr> wavefronts;           //! Constrained groups, using imageNumber/wf_num as identifier.
            
            WorkSpace(const redux::momfbd::MomfbdJob&, PatchData::Ptr);
            ~WorkSpace();
            
            void init(boost::asio::io_service&);
            void clear(void);
            void collectResults(void);
            
            PatchData::Ptr data;
            const MomfbdJob& cfg;
            redux::util::Array<PupilMode::Ptr> modes;
            redux::util::Array<double> window,noiseWindow;
            redux::util::Array<complex_t> imgFTs;               //! Stack of Fourier-transforms of input-images ([nTotalImages][pupilSize][pupilSize])
            redux::util::Array<complex_t> objFT;                //! Fourier-transforms of the approximate object ([pupilSize][pupilSize])
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_WORKSPACE_HPP
