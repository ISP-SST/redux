#ifndef REDUX_MOMFBD_WORKSPACE_HPP
#define REDUX_MOMFBD_WORKSPACE_HPP

#include "redux/types.hpp"
#include "redux/momfbd/patch.hpp"
#include "redux/momfbd/modecache.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/statistics.hpp"

#include <memory>
#include <mutex>

#include <boost/asio.hpp>

namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */

        class ObjectCfg;
        class ChannelCfg;
            
         /*! Container used during processing. Basically temporary arrays and reorganized references to original data.
         */
        struct WorkSpace {
            
            struct Image;
            struct WaveFront;
            struct Channel;
            struct Object;
            typedef std::shared_ptr<Image> ImPtr;
            typedef std::shared_ptr<WaveFront> WfPtr;
            typedef std::shared_ptr<Channel> ChPtr;
            typedef std::shared_ptr<Object> ObjPtr;
            
            struct Image : public redux::util::Array<float>, public std::enable_shared_from_this<Image> {
                
                Image(const WorkSpace&, const ObjPtr&, const ChPtr&, const redux::util::Array<float>& stack, uint32_t index, int firstY, int lastY, int firstX, int lastX);
                ~Image(void);
                
                void init(void);
                void setWaveFront( WfPtr wf ) { wfg = wf; };
                ImPtr shared(void) { return shared_from_this(); };
                void clear(void);
                
                Point offset;
                // alpha (sptr)
                // OTF
                // imgFT
                uint32_t index;
                
                const WorkSpace& ws;
                ObjPtr object;
                ChPtr channel;
                WfPtr wfg;
                
                redux::util::Array<double> img;
                redux::util::Array<complex_t> SJ;
                redux::image::FourierTransform ft;
                redux::image::Statistics stats;
                
            };


            struct WaveFront : public std::enable_shared_from_this<WaveFront> {
                
                WaveFront(ImPtr im) { if(im) { images.push_back( im->shared() ); } };
                
                void addImage(ImPtr im) { images.push_back(im); im->setWaveFront( shared() ); };
                WfPtr shared(void) { return shared_from_this(); };
                
                std::vector<double> alpha;
                std::vector<ImPtr> images;          // list of images sampling this wavefront

            };
            

            struct Object : public std::enable_shared_from_this<Object> {
                
                Object(const WorkSpace&, const std::shared_ptr<ObjectCfg>&);
                ~Object();
                
                void prepareData(const PatchData::Ptr&, std::map<uint32_t, WfPtr>&, boost::asio::io_service&);
                ObjPtr shared(void) { return shared_from_this(); };
                void clear(void);
                
                void addToFT(const redux::image::FourierTransform&);
                void addToPQ(const redux::image::FourierTransform&, const redux::util::Array<complex_t>);
                
                // modes
                // approxObj
                const WorkSpace& ws;
                std::vector<ChPtr> channels;
                std::shared_ptr<ObjectCfg> cfg;
                
                std::mutex mtx;
                redux::image::FourierTransform ftSum;
                redux::util::Array<double> Q;
                redux::util::Array<complex_t> P;
            };

            
            struct Channel : public std::enable_shared_from_this<Channel> {
                
                Channel(const WorkSpace&, const ObjPtr&, const std::shared_ptr<ChannelCfg>&);
                ~Channel(void);
                
                void prepareData(const PatchData::Ptr&, std::map<uint32_t, WfPtr>&, boost::asio::io_service&);
                ChPtr shared(void) { return shared_from_this(); };
                void clear(void);

                const WorkSpace& ws;
                std::vector<ImPtr> images;
                std::shared_ptr<ChannelCfg> cfg;
                ObjPtr obj;
                // phi_fixed
            };

            
            std::vector<ObjPtr> objects;
            std::map<uint32_t, WfPtr> wavefronts;           //! Constrained groups, using imageNumber/wf_num as identifier.
            
            WorkSpace(const redux::momfbd::MomfbdJob&, PatchData::Ptr);
            ~WorkSpace();
            
            void init(boost::asio::io_service&);
            void clear(void);
            
            PatchData::Ptr data;
            const MomfbdJob& cfg;
            redux::util::Array<ModeCache::ModePtr> modes;
            redux::util::Array<double> window,noiseWindow;
            redux::util::Array<complex_t> imgFTs;               //! Stack of Fourier-transforms of input-images ([nTotalImages][pupilSize][pupilSize])
            redux::util::Array<complex_t> objFT;                //! Fourier-transforms of the approximate object ([pupilSize][pupilSize])
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_WORKSPACE_HPP
