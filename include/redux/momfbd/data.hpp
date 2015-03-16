#ifndef REDUX_MOMFBD_DATA_HPP
#define REDUX_MOMFBD_DATA_HPP

#include "redux/momfbd/cache.hpp"

#include "redux/types.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/statistics.hpp"
#include "redux/util/array.hpp"
#include "redux/work.hpp"

#include <memory>
#include <mutex>


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
        
        
                
        /*! Datastructures to contain the data/info sent from the master to the slaves.
         *  Note: the configuration/settings are mostly sent via the MomfbdJob class, these structures
         *  just holds the patch-specific information and the input data.
         */
                
        struct ChannelData {
            typedef std::shared_ptr<ChannelData> Ptr;
            
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);

            /*! Remainders of the x/y offsets after aligning (on master) to nearest pixel,
             *  this will be added to the fixed modes, as tilts, before processing (on slaves).
             */
            PointF residualOffset;
            PointI offset;              //!< Local offset of this channel relative to the anchor
            
            //! Image data for this channel
            redux::util::Array<float> images;
        };

        struct ObjectData {
            
            typedef std::shared_ptr<ObjectData> Ptr;
            
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            std::vector<ChannelData> channels;
            
        };

        enum PartType { PT_DEFAULT=0, PT_GLOBAL };
        
        class MomfbdJob;
        struct PatchData : public Part {
            typedef std::shared_ptr<PatchData> Ptr;
            const MomfbdJob& myJob;
            std::vector<ObjectData> objects;
            Point16 index;                      //! Patch-index in mozaic
            Point16 pos;                        //! Position of patch, coordinates in the "anchor channel"
            PatchData( const MomfbdJob& j, uint16_t yid=0, uint16_t xid=0);

            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            bool operator==(const PatchData&);
        };
       
        struct GlobalData : public Part {   // Used for sending e.g. modes/pupils to the slaves.
            typedef std::shared_ptr<GlobalData> Ptr;
            std::map<std::pair<uint32_t, float>, const std::pair<redux::util::Array<double>, double>> pupils;
            std::map<Cache::ModeID, const PupilMode::Ptr> modes;
            GlobalData(void){ partType = PT_GLOBAL; };
            void fetch(const Cache::ModeID&);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
        };


        

        struct WaveFrontData;
        struct ChannelResult;
        struct ObjectResult;
        typedef std::shared_ptr<WaveFrontData> WfPtr;
        typedef std::shared_ptr<ChannelResult> ChPtr;
        typedef std::shared_ptr<ObjectResult> ObjPtr;
        
        struct ImageData : public redux::util::Array<float>, public std::enable_shared_from_this<ImageData> {
            typedef std::shared_ptr<ImageData> Ptr;
            
            ImageData(const ObjPtr&, const ChPtr&, const redux::util::Array<float>& stack, uint32_t index, int firstY, int lastY, int firstX, int lastX);
            ~ImageData(void);
            
            void init(void);
            void setWaveFront( WfPtr wf ) { wfg = wf; };
            void clear(void);
            
            Point offset;
            // alpha (sptr)
            // OTF
            // imgFT
            uint32_t index;
            
            ObjPtr object;
            ChPtr channel;
            WfPtr wfg;
            
            redux::util::Array<double> img;
            redux::util::Array<complex_t> SJ;
            redux::image::FourierTransform ft;
            redux::image::Statistics stats;
            
        };


        struct WaveFrontData : public std::enable_shared_from_this<WaveFrontData> {
            
            WaveFrontData(ImageData::Ptr im) { if(im) { images.push_back(im); } };
            
            void addImage(ImageData::Ptr im) { images.push_back(im); im->setWaveFront( shared_from_this() ); };
            
            std::vector<double> alpha;
            std::vector<ImageData::Ptr> images;          // list of images sampling this wavefront

        };
        

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_DATA_HPP
