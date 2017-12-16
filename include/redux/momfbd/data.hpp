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
            
            explicit ChannelData( std::shared_ptr<Channel> c );
            ~ChannelData();
            
            /**** Local processing (slaves) ****/
            void initPatch(void);
            /***********************************/
            
            /********* Input *********/
            redux::util::Array<float> images;                           //!< Image stack for this channel.
            redux::util::RegionI cutout;                                //!< Coordinates of the cutout region. (size is usually patchSize+2*maxLocalShift, but might be constrained by edges)
            redux::util::PointI channelOffset;                          //!< The displacement of this channel relative to the reference channel.
            redux::util::PointI patchStart;                             //!< Local offset, i.e. where in the cutout we cut out the subimage for processing.
            redux::util::PointF residualOffset;                         //!< Remainders of the x/y offsets after aligning (on master) to nearest pixel.
            /*************************/
            
            void clear(void);
            void load(void);
            uint64_t size(void) const;
            uint64_t pack(char*) const;                         //!< Pack channel data to char-array (for sending/storing)
            uint64_t unpack(const char*, bool);
            const ChannelData& operator=(const ChannelData&);

            void copyResults( const ChannelData& rhs );
            void dump( std::string ) const;
            
            std::shared_ptr<Channel> myChannel;
            
        };

        
        class Object;
        struct ObjectData {
            
            typedef std::shared_ptr<ObjectData> Ptr;
            std::shared_ptr<Object> myObject;
            ObjectData(void) {};
            explicit ObjectData( std::shared_ptr<Object> o );
            ~ObjectData();
            
            void initPatch(void);
            void clear(void);
            void load(void);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            const ObjectData& operator=(const ObjectData&);
            
            void copyResults( const ObjectData& rhs );
            void dump( std::string ) const;
            
            std::vector<std::shared_ptr<ChannelData>> channels;
            
            /********* Results ********/
            redux::util::Array<float> img;              //!< Restored image.
            redux::util::Array<float> psf;              //!< PSFs.               ((1|nImages) x nPixels x nPixels)
            redux::util::Array<float> cobj;             //!< Convolved objects.  (nImages x nPixels x nPixels)
            redux::util::Array<float> res;              //!< Residuals.          (nImages x nPixels x nPixels)
            redux::util::Array<float> alpha;            //!< Mode coefficients   (nImages x nModes)
            redux::util::Array<float> div;              //!< Diversity           (nCh x pupilPixels x pupilPixels)
            /**************************/
            
        };

        struct WavefrontData {

            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            const WavefrontData& operator=(const WavefrontData&);
            void copyResults( const WavefrontData& rhs );
            void dump( std::string ) const;
            
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
            std::vector<std::shared_ptr<ObjectData>> objects;
            std::vector<std::shared_ptr<ObjectData>> trace_data;
            WavefrontData waveFronts;
            redux::util::Point16 index;                      //! Patch-index in mozaic
            redux::util::Point16 position;                   //! Patch position (centre coordinates in the reference channel)
            redux::util::Region16 roi;                       //! Region/position of this patch in the full image
            float finalMetric;
            std::vector<float> metrics;
            PatchData( MomfbdJob& j, uint16_t yid=0, uint16_t xid=0);
            PatchData( const PatchData& ) = delete;
            ~PatchData();

            void setPath(const std::string& path);
            void initPatch(void);
            void clear(void);
            void load(void);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            size_t csize(void) const { return size(); };
            uint64_t cpack(char* p) const { return pack(p); };
            uint64_t cunpack(const char* p, bool e) { return unpack(p,e); };
            void cclear(void);
            bool cacheLoad(bool removeAfterLoad=false);
            bool cacheStore(bool clearAfterStore=false);

            bool operator==(const PatchData&);
            const PatchData& operator=(const PatchData&);
            
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
            std::mutex mtx;
            std::map<ModeInfo, std::shared_ptr<ModeSet>> modes;
            std::map<redux::image::PupilInfo, std::shared_ptr<redux::image::Pupil>> pupils;
            Constraints constraints;
            explicit GlobalData(MomfbdJob& j ) : constraints(j) { partType = PT_GLOBAL; };
            std::shared_ptr<ModeSet> get(const ModeInfo&, const std::shared_ptr<ModeSet>& ms=nullptr);
            std::shared_ptr<redux::image::Pupil> get(const redux::image::PupilInfo&,
                                                     const std::shared_ptr<redux::image::Pupil>& ms=nullptr);
            bool verify(void) const;
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            void dump( std::string tag="gd" ) const;

        };


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_DATA_HPP
