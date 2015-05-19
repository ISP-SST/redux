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
                
        class Channel;
        struct ChannelData {
            
            typedef std::shared_ptr<ChannelData> Ptr;
            std::shared_ptr<Channel> myChannel;
            ChannelData( std::shared_ptr<Channel> c );
            void getPatchData(Point16&);
            void initPatch(void);
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

        class Object;
        struct ObjectData {
            
            typedef std::shared_ptr<ObjectData> Ptr;
            std::shared_ptr<Object> myObject;
            ObjectData( std::shared_ptr<Object> o );
            void getPatchData(Point16&);
            void initPatch(void);
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
            void getData(void);
            void initPatch(void);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            bool operator==(const PatchData&);
        };
       
        struct GlobalData : public Part {   // Used for sending e.g. modes/pupils to the slaves.
            typedef std::shared_ptr<GlobalData> Ptr;
            std::mutex mtx;
            std::map<std::pair<uint32_t, float>, const std::pair<redux::util::Array<double>, double>&> pupils;
            std::map<Cache::ModeID, const PupilMode::Ptr> modes;
            GlobalData(void){ partType = PT_GLOBAL; };
            const std::pair<redux::util::Array<double>, double>& fetch(uint16_t,double);
            const PupilMode::Ptr fetch(const Cache::ModeID&);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
        };


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_DATA_HPP
