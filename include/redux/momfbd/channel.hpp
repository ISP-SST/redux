#ifndef REDUX_MOMFBD_CHANNEL_HPP
#define REDUX_MOMFBD_CHANNEL_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/solver.hpp"

#include "redux/image/image.hpp"
#include "redux/util/arraystats.hpp"
#include "redux/util/progresswatch.hpp"
#include "redux/util/region.hpp"

#ifdef RDX_TRACE_JOB
#   include "redux/util/trace.hpp"
#endif

#include <future>

#include <boost/asio.hpp>
#include <boost/property_tree/ptree.hpp>
namespace bpt = boost::property_tree;


namespace redux {

    namespace logging {
        class Logger;
    }

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
        

        class MomfbdJob;
        class Object;
        struct Solver;
        struct ChannelData;
        struct SubImage;
        struct Tilts;
        struct WaveFront;
        /*! @brief Class containing the channel-specific configuration for a MomfbdJob/Object
         * 
         */
        class Channel : public ChannelCfg
#ifdef RDX_TRACE_JOB
#ifndef RDX_TRACE_PROC
            , public redux::util::TraceObject<Channel>
#endif
#endif
        {

        public:

            typedef std::shared_ptr<Channel> Ptr;
            Channel( Object&, MomfbdJob&, uint16_t id=0 );
            ~Channel();

            void parsePropertyTree( bpt::ptree&, redux::logging::Logger& );
            bpt::ptree getPropertyTree( bpt::ptree&, bool showAll=false );

            uint16_t id(void) const { return ID; };
            void cleanup(void);
            size_t size(void) const override;
            uint64_t pack(char*) const override;
            uint64_t unpack(const char*, bool) override;
            double getMaxMean(void);
            double getMaxMedian(void);
            double getMedianMedian(void);
            void getFileNames(std::vector<std::string>&, const std::vector<uint32_t>&) const;
            void getFileNames(std::vector<std::string>& list) const { getFileNames( list, waveFrontList ); };
            uint32_t nImages( const std::vector<uint32_t>& );
            uint32_t nImages(void);
            redux::util::PointF getOffsetAt( const redux::util::Point16& pos, size_t sz=0 ) const;
            void adjustCutout(ChannelData&, const PatchData::Ptr&) const;
            void adjustCutouts(redux::util::Array<PatchData::Ptr>&);
            void setClip( const std::vector<int16_t>& clip ) { alignClip = clip; checkFlip(); };
            void setMap( const std::vector<float>& map ) { alignMap = map; checkFlip(); };
            void setMaps( const std::vector<float>& xmap, const std::vector<float>& ymap ) { alignMapX = xmap; alignMapY=ymap; checkFlip(); };
            void setImageSize( const redux::util::Point16& sz ) { imgSize = sz; };
            void setImageSize( uint16_t sz ) { imgSize = sz; };
            void clearOffsets( void ) { xOffset.clear(); yOffset.clear(); }
            void setOffsets( image::Image<int16_t>& xoff, image::Image<int16_t>& yoff ) {
                xOffset = xoff; // N.B. will use shared data.
                yOffset = yoff;
            }
            
            /*************   Processing on slave   ***************/
            /*************         Methods         ***************/
            void initProcessing( Solver& );
            void initPatch(ChannelData&);
            const std::vector<std::shared_ptr<SubImage>>& getSubImages(void) const { return subImages; };
            void initPhiFixed(void);
            void addAllFT(redux::util::Array<double>&);
            //void addAllPQ(void) const;
            /*****************************************************/

            std::string idString( void ) const;
            void dump( std::string tag );
        private:
            
            void checkFlip( void );
            bool checkCfg(void);
            bool checkData( bool verbose=false );
            
            void initChannel(void);
            
            void loadCalib(boost::asio::io_service&);
            void loadData(boost::asio::io_service&, redux::util::Array<PatchData::Ptr>&);
            void unloadCalib(void);

            void addTimeStamps( const bpx::ptime& newStart, const bpx::ptime& newEnd );
            void loadFile( size_t fileIndex, size_t offset );
            void preprocessImage( size_t index );
            void maybeLoadImages( void );          
            void getStorage(ChannelData&);          
            
            redux::util::Point16 getImageSize( bool force=false );
            void logAndThrow( std::string );


            /*************   Local variables for   ***************/
            /************* PreProcessing on master ***************/
            std::vector<redux::util::ArrayStats::Ptr> imageStats;
            redux::image::Image<float> images;
            redux::image::Image<float> dark, gain;
            redux::image::Image<float> ccdResponse, ccdScattering;
            redux::image::Image<float> psf, modulationMatrix;
            redux::image::Image<int16_t> xOffset, yOffset;
            std::shared_ptr<uint8_t> gainMask;
            bpx::ptime startT, endT;
            std::future<bool> patchWriteFail;
            std::vector<size_t> nFrames;                                    //!< Number of frames in each file
            std::vector<std::vector<size_t>> frameNumbersPerFile;            //!< for storing the frameNumbers obtained from each file.
            std::string cacheFile;
            /*****************************************************/
            
            /*************   Local variables for   ***************/
            /*************   Processing on slave   ***************/
            std::mutex mtx;
            redux::util::ProgressWatch progWatch;
            std::vector<std::shared_ptr<SubImage>> subImages;
            redux::util::Array<double> phi_fixed;                           //!< The fixed aberrations for this channel (i.e. phase diversity)
            redux::util::Array<double> phi_channel;                         //!< The fixed part + tilt-corrections for this channel
            double otfNormalization;
            /*****************************************************/
            
            uint16_t ID;
            redux::util::Point16 imgSize;
            bool flipX, flipY;
            uint32_t nTotalFrames;
            Object& myObject;
            MomfbdJob& myJob;
            logging::Logger& logger;

            friend struct PatchData;
            friend struct ObjectData;
            friend struct ChannelData;
            friend class Object;
            friend class MomfbdJob;
            friend struct ModeSet;
            friend struct Solver;
            friend struct SubImage;
            
        };

        /*! @} */
                
    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_CHANNEL_HPP
