#ifndef REDUX_MOMFBD_DATA_HPP
#define REDUX_MOMFBD_DATA_HPP

#include "redux/momfbd/constraints.hpp"
#include "redux/momfbd/modes.hpp"

#include "redux/util/region.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/pupil.hpp"
#include "redux/util/array.hpp"
#include "redux/work.hpp"
#include "redux/util/compress.hpp"

#ifdef RDX_TRACE_PARTS
#   include "redux/util/trace.hpp"
#endif

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
        struct ChannelData
#ifdef RDX_TRACE_PARTS
            : public redux::util::TraceObject<ChannelData>
#endif
        {
            
            typedef std::shared_ptr<ChannelData> Ptr;
            
            explicit ChannelData( std::shared_ptr<Channel> c );
            ChannelData( const ChannelData& ) = delete;
            ChannelData( ChannelData&& ) = delete;
            ~ChannelData();
            
            /**** Local processing (slaves) ****/
            void initPatch(void);
            /***********************************/
            
            /********* Input *********/
            redux::util::Array<float> images;               //!< Image stack for this channel.
            redux::util::PointF exactPatchPosition;         //!< True coordinates for the patch in the (un-clipped) camera reference frame.
            redux::util::PointI cutoutPosition;             //!< The integer position of the patch, restricted by edges and rounded
            redux::util::RegionI cutoutRegion;              //!< Coordinates of the cutout region. (size is usually patchSize+2*maxLocalShift, but might be constrained by edges)
            redux::util::PointI channelOffset;              //!< The displacement of this channel relative to the reference channel.
            redux::util::PointI patchStart;                 //!< Local offset, i.e. where in the cutout we cut out the subimage for processing.
            redux::util::PointF residualOffset;             //!< Remainders of the x/y offsets after aligning (on master) to nearest pixel.
             /*************************/
            
            void clear(void);
            void load(void);
            uint64_t size(void) const;
            uint64_t pack(char*) const;                         //!< Pack channel data to char-array (for sending/storing)
            uint64_t unpack(const char*, bool);
            ChannelData& operator=(const ChannelData&);

            void copyResults( const ChannelData& rhs );
            void dump( std::string ) const;
            
            std::shared_ptr<Channel> myChannel;
            
        };

        
        class Object;
        struct ObjectData
#ifdef RDX_TRACE_PARTS
            : public redux::util::TraceObject<ObjectData>
#endif
         {
            
            typedef std::shared_ptr<ObjectData> Ptr;
            std::shared_ptr<Object> myObject;
            
            ObjectData(void);
            explicit ObjectData( std::shared_ptr<Object> o );
            ~ObjectData();
            
            void initPatch(void);
            void clear(void);
            void load(void);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            ObjectData& operator=(const ObjectData&);
            
            void copyResults( const ObjectData& rhs );
            void dump( std::string ) const;
            
            std::vector<ChannelData::Ptr> channels;
            
            /********* Results ********/
            redux::util::Array<float> img;              //!< Restored image.
            redux::util::Array<float> psf;              //!< PSFs.               ((1|nImages) x nPixels x nPixels)
            redux::util::Array<float> cobj;             //!< Convolved objects.  (nImages x nPixels x nPixels)
            redux::util::Array<float> res;              //!< Residuals.          (nImages x nPixels x nPixels)
            redux::util::Array<float> alpha;            //!< Mode coefficients   (nImages x nModes)
            redux::util::Array<float> div;              //!< Diversity           (nCh x pupilPixels x pupilPixels)
            /**************************/
            
        };

        struct WavefrontData
#ifdef RDX_TRACE_PARTS
            : public redux::util::TraceObject<WavefrontData>
#endif
        {

            WavefrontData();
            ~WavefrontData();
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            WavefrontData& operator=(const WavefrontData&);
            void copyResults( const WavefrontData& rhs );
            void dump( const std::string& ) const;
            
            /********* Results ********/
            std::vector<uint32_t> ids;                  //!< Wavefront indices   (nWf)
            redux::util::Array<float> alpha;            //!< Mode coefficients   (nWf x nModes)
            /**************************/
            
        };

        enum PartType { PT_DEFAULT=0, PT_GLOBAL };
        
        class MomfbdJob;
        struct PatchData : public Part {
            typedef std::shared_ptr<PatchData> Ptr;
            MomfbdJob& myJob;
            std::vector<ObjectData::Ptr> objects;
            std::vector<ObjectData::Ptr> trace_data;
            WavefrontData waveFronts;
            redux::util::Point16 index;                      //! Patch-index in mozaic
            redux::util::Point16 position;                   //! Patch position (centre coordinates in the reference channel)
            redux::util::Region16 roi;                       //! Region/position of this patch in the full image
            float finalMetric;
            std::vector<float> metrics;
            explicit PatchData( MomfbdJob& j, uint16_t yid=0, uint16_t xid=0);
            PatchData( const PatchData& ) = delete;
            ~PatchData();

            ObjectData::Ptr getObjectData( uint16_t id ) const;
            ChannelData::Ptr getChannelData( uint16_t oid, uint16_t cid ) const;
            void setPath(const std::string& path) override;
            void initPatch(void);
            void load(void) override;
            void unload(void) override;
            void prePack( bool force=false ) override;
            void clear(void);
            uint64_t size(void) const override;
            uint64_t pack(char*) const override;
            uint64_t unpack(const char*, bool) override;
            size_t csize(void) const override { return size(); };
            uint64_t cpack(char* p) const override { return pack(p); };
            uint64_t cunpack(const char* p, bool e) override { return unpack(p,e); };
            void cclear(void) override;
            bool cacheLoad(bool removeAfterLoad=false) override;
            bool cacheStore(bool clearAfterStore=false) override;

            bool operator==(const PatchData&);
            PatchData& operator=(const PatchData&);
            
            template <typename T>
            void loadAlpha( T* a ) const {
                if( waveFronts.alpha.nElements() ) {
                    waveFronts.alpha.copyTo<T>(a);
                }
            }
            template <typename T>
            void storeAlpha( T* a ) {
                if( waveFronts.alpha.nElements() ) {
                    waveFronts.alpha.copyFrom<T>(a);
                }
            }

            void copyResults( const PatchData& rhs );
            void dump( std::string ) const;
            
        };
       
        struct GlobalData : public Part {   // Used for sending e.g. modes/pupils to the slaves.
            typedef std::shared_ptr<GlobalData> Ptr;
            typedef std::pair<ModeInfo, double> ScaledModeInfo;
            mutable std::mutex mtx;
            std::map<ModeInfo, std::shared_ptr<ModeSet>> modes;
            std::map<ScaledModeInfo, std::shared_ptr<ModeSet>> scaled_modes;
            std::map<redux::image::PupilInfo, std::shared_ptr<redux::image::Pupil>> pupils;
            Constraints constraints;
            std::shared_ptr<ModeSet> get(const ModeInfo&, double scale, const std::shared_ptr<redux::image::Pupil>& pupil, const std::shared_ptr<ModeSet>& ms=nullptr);
            explicit GlobalData( MomfbdJob& j );
            ~GlobalData();
            std::shared_ptr<ModeSet> get(const ModeInfo&, const std::shared_ptr<ModeSet>& ms=nullptr);
            std::shared_ptr<redux::image::Pupil> get(const redux::image::PupilInfo&,
                                                     const std::shared_ptr<redux::image::Pupil>& ms=nullptr);
            bool verify(void) const;
            uint64_t size(void) const override;
            uint64_t pack(char*) const override;
            uint64_t unpack(const char*, bool) override;
            void unload(void) override;
            void prePack( bool force=false ) override;
            void dump( const std::string& tag="gd" ) const;
            
        };


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_DATA_HPP
