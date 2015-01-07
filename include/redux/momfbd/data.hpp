#ifndef REDUX_MOMFBD_DATA_HPP
#define REDUX_MOMFBD_DATA_HPP

#include "redux/momfbd/patch.hpp"

#include "redux/types.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/statistics.hpp"
#include "redux/util/array.hpp"

#include <memory>


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
        struct WorkSpace;
        struct Object;
        struct Channel;
        struct ImageData;
        struct WaveFrontData;
        struct ChannelData;
        struct ObjectData;
        typedef std::shared_ptr<ImageData> ImPtr;
        typedef std::shared_ptr<WaveFrontData> WfPtr;
        typedef std::shared_ptr<ChannelData> ChPtr;
        typedef std::shared_ptr<ObjectData> ObjPtr;
        
        struct ImageData : public redux::util::Array<float>, public std::enable_shared_from_this<ImageData> {
            
            ImageData(const WorkSpace&, const ObjPtr&, const ChPtr&, const redux::util::Array<float>& stack, uint32_t index, int firstY, int lastY, int firstX, int lastX);
            ~ImageData(void);
            
            void init(void);
            void setWaveFront( WfPtr wf ) { wfg = wf; };
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


        struct WaveFrontData : public std::enable_shared_from_this<WaveFrontData> {
            
            WaveFrontData(ImPtr im) { if(im) { images.push_back(im); } };
            
            void addImage(ImPtr im) { images.push_back(im); im->setWaveFront( shared_from_this() ); };
            
            std::vector<double> alpha;
            std::vector<ImPtr> images;          // list of images sampling this wavefront

        };
        

        struct ObjectData : public std::enable_shared_from_this<ObjectData> {
            
            ObjectData(const WorkSpace&, const std::shared_ptr<Object>&);
            ~ObjectData();
            
            void prepareData(const PatchData::Ptr&, std::map<uint32_t, WfPtr>&, boost::asio::io_service&);
            void clear(void);
            
            void addToFT(const redux::image::FourierTransform&);
            void addToPQ(const redux::image::FourierTransform&, const redux::util::Array<complex_t>);
            
            // modes
            // approxObj
            const WorkSpace& ws;
            std::vector<ChPtr> channels;
            std::shared_ptr<Object> cfg;
            
            std::mutex mtx;
            redux::image::FourierTransform ftSum;
            redux::util::Array<double> Q;
            redux::util::Array<complex_t> P;
        };

        
        struct ChannelData : public std::enable_shared_from_this<ChannelData> {
            
            ChannelData(const WorkSpace&, const ObjPtr&, const std::shared_ptr<Channel>&);
            ~ChannelData(void);
            
            void prepareData(const PatchData::Ptr&, std::map<uint32_t, WfPtr>&, boost::asio::io_service&);
            void clear(void);

            const WorkSpace& ws;
            std::vector<ImPtr> images;
            std::shared_ptr<Channel> cfg;
            ObjPtr obj;
            // phi_fixed
        };


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_DATA_HPP
