#ifndef REDUX_MOMFBD_PATCH_HPP
#define REDUX_MOMFBD_PATCH_HPP

#include "redux/momfbd/cache.hpp"
#include "redux/momfbd/modes.hpp"

#include "redux/image/fouriertransform.hpp"
#include "redux/types.hpp"
#include "redux/work.hpp"
#include "redux/util/array.hpp"

#include <memory>

namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
        class MomfbdJob;
        class Object;
        class Channel;
        struct ObjectData;
        /*! Datastructures to contain the data/info sent from the master to the slaves.
         *  Note: the configuration/settings are mostly sent via the MomfbdJob class, these structures
         *  just holds the patch-specific information and the input data.
         */
                
        struct ChannelData : public std::enable_shared_from_this<ChannelData> {
            
            typedef std::shared_ptr<ChannelData> Ptr;
            
            ChannelData(const std::shared_ptr<ObjectData>&, const std::shared_ptr<Channel>&);
            ~ChannelData(void);
            
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);

            void init(uint16_t yid, uint16_t xid);
            void clear(void);

            PointI offset;                      //! Local offset of this channel relative to the anchor
            PointF residualOffset;              //! Remaining part of the x/y offsets after aligning to nearest pixel
            
            redux::util::Array<float> images;

            std::shared_ptr<Channel> channel;
            std::shared_ptr<ObjectData> object;

        };


        struct ObjectData : public std::enable_shared_from_this<ObjectData> {
            
            typedef std::shared_ptr<ObjectData> Ptr;
        
            ObjectData(const std::shared_ptr<Object>&);
            ~ObjectData();
            
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);

            void init(uint16_t yid, uint16_t xid);
            void clear(void);
            
            void addToFT(const redux::image::FourierTransform&);
            void addToPQ(const redux::image::FourierTransform&, const redux::util::Array<complex_t>);
            
            std::mutex mtx;
            redux::image::FourierTransform ftSum;
            redux::util::Array<double> Q;
            redux::util::Array<complex_t> P;
            
            std::vector<std::shared_ptr<ChannelData>> channels;
            std::shared_ptr<Object> object;

            
        };

        
        struct PatchData : public Part {
            typedef std::shared_ptr<PatchData> Ptr;
            Point16 index;                      //! Patch-index in mozaic
            Point16 pos;                        //! Position of patch, coordinates in the "anchor channel"
            std::vector<ObjectData::Ptr> objects;
            const MomfbdJob& myJob;
            PatchData( const MomfbdJob& j ) : myJob(j) {};
            void setIndex(uint16_t yid, uint16_t xid);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            bool operator==(const PatchData&);
        };

        
        struct GlobalData : public Part {   // Used for sending e.g. modes/pupils to the slaves.
            typedef std::shared_ptr<GlobalData> Ptr;
            std::map<std::pair<uint32_t, float>, const std::pair<redux::util::Array<double>, double>> pupils;
            std::map<Cache::ModeID, const PupilMode::Ptr> modes;
            void fetch(const Cache::ModeID&);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
        };


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_PATCH_HPP
