#ifndef REDUX_MOMFBD_DATA_HPP
#define REDUX_MOMFBD_DATA_HPP

#include "redux/momfbd/constraints.hpp"
#include "redux/momfbd/modes.hpp"

#include "redux/types.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/pupil.hpp"
#include "redux/util/array.hpp"
#include "redux/work.hpp"
#include "redux/util/compress.hpp"

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
        struct ChannelData : public redux::util::CacheItem {
            
            explicit ChannelData( std::shared_ptr<Channel> c );
            ~ChannelData();
            
            void setPath(const std::string& path);
            
            /**** Local processing (slaves) ****/
            void initPatch(void);
            /***********************************/
            
            /********* Input *********/
            redux::util::Array<float> images;                   //!< Image stack for this channel.
            RegionI cutout;                                     //!< Coordinates of the cutout region.
            PointI shift;                                       //!< The displacement of this channel. Initially just from the x/y-offset files, but during the fitting the cutout might get shifted.
            PointI offset;                                      //!< Local offset, i.e. where in the "block" we cut out the subimage for processing.
            PointF residualOffset;                              //!< Remainders of the x/y offsets after aligning (on master) to nearest pixel.
            /*************************/
            
            virtual uint64_t size(void) const;
            virtual uint64_t pack(char*) const;                         //!< Pack channel data to char-array (for sending/storing)
            virtual uint64_t unpack(const char*, bool);
            size_t csize(void) const { return size(); };
            uint64_t cpack(char* p) const { return pack(p); };
            uint64_t cunpack(const char* p, bool e) { return unpack(p,e); };
            void cclear(void);

            std::shared_ptr<Channel> myChannel;
            
        };

        
        class Object;
        struct ObjectData : public redux::util::CacheItem {
            
            typedef std::shared_ptr<ObjectData> Ptr;
            std::shared_ptr<Object> myObject;
            ObjectData( std::shared_ptr<Object> o );
            ~ObjectData();
            
            void setPath(const std::string& path);
            void initPatch(void);
            virtual uint64_t size(void) const;
            virtual uint64_t pack(char*) const;
            virtual uint64_t unpack(const char*, bool);
            size_t csize(void) const { return size(); };
            uint64_t cpack(char* p) const { return pack(p); };
            uint64_t cunpack(const char* p, bool e) { return unpack(p,e); };
            bool cacheLoad(bool removeAfterLoad=false);
            bool cacheStore(bool clearAfterStore=false);
            void cclear(void);
            
            std::vector<redux::util::Compressed<ChannelData,5>> channels;
            
            /********* Results ********/
            redux::util::Array<float> img;              //!< Restored image.
            redux::util::Array<float> psf;              //!< PSFs.               ((1|nImages) x nPixels x nPixels)
            redux::util::Array<float> cobj;             //!< Convolved objects.  (nImages x nPixels x nPixels)
            redux::util::Array<float> res;              //!< Residuals.          (nImages x nPixels x nPixels)
            redux::util::Array<float> alpha;            //!< Mode coefficients   (nImages x nModes)
            redux::util::Array<float> div;              //!< Diversity           (nCh x pupilPixels x pupilPixels)
            /**************************/
            
        };

        enum PartType { PT_DEFAULT=0, PT_GLOBAL };
        
        class MomfbdJob;
        struct PatchData : public Part {
            typedef std::shared_ptr<PatchData> Ptr;
            const MomfbdJob& myJob;
            std::vector<redux::util::Compressed<ObjectData,5>> objects;
            Point16 index;                      //! Patch-index in mozaic
            Region16 roi;                       //! Region/position of this patch in the full image
            PatchData( const MomfbdJob& j, uint16_t yid=0, uint16_t xid=0);
            PatchData( const PatchData& ) = delete;
            ~PatchData();

            void setPath(const std::string& path);
            void initPatch(void);
            virtual uint64_t size(void) const;
            virtual uint64_t pack(char*) const;
            virtual uint64_t unpack(const char*, bool);
            size_t csize(void) const { return size(); };
            uint64_t cpack(char* p) const { return pack(p); };
            uint64_t cunpack(const char* p, bool e) { return unpack(p,e); };
            void cclear(void);
            bool cacheLoad(bool removeAfterLoad=false);
            bool cacheStore(bool clearAfterStore=false);

            bool operator==(const PatchData&);
        };
       
        struct GlobalData : public Part {   // Used for sending e.g. modes/pupils to the slaves.
            typedef std::shared_ptr<GlobalData> Ptr;
            std::mutex mtx;
            std::map<ModeInfo, ModeSet&> modes;
            std::map<redux::image::PupilInfo, redux::image::Pupil&> pupils;
            Constraints constraints;
            GlobalData(const MomfbdJob& j ) : constraints(j) { partType = PT_GLOBAL; };
            ModeSet& get(const ModeInfo&, const ModeSet& ms=ModeSet());
            redux::image::Pupil& get(const redux::image::PupilInfo&, const redux::image::Pupil& ms=redux::image::Pupil());
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);

        };


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_DATA_HPP
