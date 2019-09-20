#ifndef REDUX_MOMFBD_CHANNEL_HPP
#define REDUX_MOMFBD_CHANNEL_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/solver.hpp"

#include <redux/image/image.hpp>
#include <redux/util/arraystats.hpp>
#include <redux/util/progresswatch.hpp>
#include <redux/util/region.hpp>

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
        class Solver;
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

            Channel( Object&, MomfbdJob&, uint16_t id=0 );
            ~Channel();

            void parsePropertyTree( bpt::ptree&, redux::logging::Logger& );
            bpt::ptree getPropertyTree( bpt::ptree&);

            uint16_t id(void) const { return ID; };
            void cleanup(void);
            size_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            double getMaxMean(void);
            void getFileNames(std::vector<std::string>&, const std::vector<uint32_t>&) const;
            void getFileNames(std::vector<std::string>& list) const { getFileNames( list, waveFrontList ); };
            uint32_t nImages( const std::vector<uint32_t>& );
            uint32_t nImages(void);
            void adjustCutout(ChannelData&, const PatchData::Ptr&) const;
            void adjustCutouts(redux::util::Array<PatchData::Ptr>&);          
            
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
            
            bool checkCfg(void);
            bool checkData( bool verbose=false );
            
            void initCache(void);
            
            void loadCalib(boost::asio::io_service&);
            void loadData(boost::asio::io_service&, redux::util::Array<PatchData::Ptr>&);
            void unloadCalib(void);

            void addTimeStamps( const bpx::ptime& newStart, const bpx::ptime& newEnd );
            void loadFile( size_t fileIndex, size_t offset );
            void preprocessImage( size_t index );
            void maybeLoadImages( void );          
            void getStorage(ChannelData&);          
            
            redux::util::Point16 getImageSize(void);
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

            friend struct ChannelData;
            friend class Object;
            friend class MomfbdJob;
            friend struct ModeSet;
            friend class Solver;
            friend struct SubImage;
            
        };

        /*! @} */
                
    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_CHANNEL_HPP
