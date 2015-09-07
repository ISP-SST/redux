#ifndef REDUX_MOMFBD_DATA_HPP
#define REDUX_MOMFBD_DATA_HPP

#include "redux/momfbd/cache.hpp"
#include "redux/momfbd/constraints.hpp"

#include "redux/types.hpp"
#include "redux/image/fouriertransform.hpp"
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
        struct PatchData;
        struct ChannelData {
            
            ChannelData( std::shared_ptr<Channel> c );
            
            /***** PreProcessing (manager) *****/
            void getPatchData(const PatchData&);
            /***********************************/
            
            /**** Local processing (slaves) ****/
            void initPatch(void);
            void collectResults(void);
            /***********************************/
            
            /********* Input *********/
            redux::util::Array<float> images;                   //!< Image stack for this channel.
            PointI shift;                                       //!< The displacement of this channel. Initially just from the x/y-offset files, but during the fitting the cutout might get shifted.
            PointI offset;                                      //!< Local offset, i.e. where in the "block" we cut out the subimage for processing.
            PointF residualOffset;                              //!< Remainders of the x/y offsets after aligning (on master) to nearest pixel.
            /*************************/
            
            /********* Output ********/
            redux::util::Array<float> alpha;                    //!< Mode coefficients  (nImages x nModes)
            /*************************/
            
            
            uint64_t size(void) const;
            uint64_t pack(char*) const;                         //!< Pack channel data to char-array (for sending/storing)
            uint64_t unpack(const char*, bool);
            void cclear(void);

            std::shared_ptr<Channel> myChannel;
            
        };

        
        class Object;
        struct ObjectData {
            
            typedef std::shared_ptr<ObjectData> Ptr;
            std::shared_ptr<Object> myObject;
            ObjectData( std::shared_ptr<Object> o );
            void getPatchData(const PatchData& patch);
            void initPatch(void);
            void collectResults(void);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            void cclear(void);
            
            std::vector<ChannelData> channels;
            
            redux::util::Array<float> results;      //! stack of images to be returned to the master
            
        };

        enum PartType { PT_DEFAULT=0, PT_GLOBAL };
        
        class MomfbdJob;
        struct PatchData : public Part, public redux::util::CacheItem {
            typedef std::shared_ptr<PatchData> Ptr;
            const MomfbdJob& myJob;
            std::vector<ObjectData> objects;
            Point16 index;                      //! Patch-index in mozaic
            Region16 roi;                       //! Region/position of this patch in the full image
            PatchData( const MomfbdJob& j, uint16_t yid=0, uint16_t xid=0);
            void getData(void);
            void initPatch(void);
            void collectResults(void);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            size_t csize(void) const { return size(); };
            uint64_t cpack(char* p) const { return pack(p); };
            uint64_t cunpack(const char* p, bool e) { return unpack(p,e); };
            void cclear(void);
            bool operator==(const PatchData&);
        };
       
        struct GlobalData : public Part {   // Used for sending e.g. modes/pupils to the slaves.
            typedef std::shared_ptr<GlobalData> Ptr;
            std::mutex mtx;
            std::map<std::pair<uint32_t, double>, const std::pair<redux::util::Array<double>, double>&> pupils;
            std::map<Cache::ModeID, const PupilMode::Ptr> modes;
            Constraints constraints;
            GlobalData(const MomfbdJob& j ) : constraints(j) { partType = PT_GLOBAL; };
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
