#ifndef REDUX_MOMFBD_CHANNEL_HPP
#define REDUX_MOMFBD_CHANNEL_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/solver.hpp"

#include <redux/image/image.hpp>
#include <redux/util/arraystats.hpp>

#include <future>

#include <boost/asio.hpp>
#include <boost/property_tree/ptree.hpp>
namespace bpt = boost::property_tree;


namespace redux {

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
        class Channel : public ChannelCfg {

        public:

            Channel( Object&, MomfbdJob&, uint16_t id=0 );
            ~Channel();

            void parsePropertyTree( bpt::ptree& );
            bpt::ptree getPropertyTree( bpt::ptree&);

            uint16_t id(void) const { return ID; };
            size_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            double getMaxMean(void);
            void getFileNames(std::vector<std::string>&) const;
            uint32_t nImages(void) const { return imageNumbers.size(); } 
            void adjustCutout(ChannelData&, const Region16&) const;
            void adjustCutouts(redux::util::Array<PatchData::Ptr>&);          
            void storePatchData(boost::asio::io_service& service, redux::util::Array<PatchData::Ptr>&);          
            
            /*************   Processing on slave   ***************/
            /*************         Methods         ***************/
            void initProcessing( const Solver& );
            void initPatch(ChannelData&);
            const std::vector<std::shared_ptr<SubImage>>& getSubImages(void) const { return subImages; };
            void initPhiFixed(void);
            void computePhi(void);
            //void addMode(redux::util::Array<double>&, uint16_t, double) const;
            //void getPhi(redux::util::Array<double>&, const WaveFront&) const;
            //void addAllFT(redux::image::FourierTransform&);
            void addAllFT(redux::util::Array<double>&);
            void addAllPQ(void) const;
            double metric(void);
            /*****************************************************/

            void dump( std::string tag );
        private:
            
            bool checkCfg(void);
            bool checkData(void);
            
            void initCache(void);
            
            void loadCalib(boost::asio::io_service&);
            void loadData(boost::asio::io_service&, redux::util::Array<PatchData::Ptr>&);
            void unloadCalib(void);

            void addTimeStamps( const bpx::ptime& newStart, const bpx::ptime& newEnd );
            void loadImage(size_t index);
            void preprocessImage(size_t index, redux::image::Image<float>& img);
            void copyImagesToPatch(ChannelData&);          
            
            Point16 getImageSize(void);
            void logAndThrow( std::string );


            /*************   Local variables for   ***************/
            /************* PreProcessing on master ***************/
            std::vector<redux::util::ArrayStats::Ptr> imageStats;
            redux::image::Image<float> images;
            redux::image::Image<float> dark, gain;
            redux::image::Image<float> ccdResponse, ccdScattering;
            redux::image::Image<float> psf, modulationMatrix;
            redux::image::Image<int16_t> xOffset, yOffset;
            bpx::ptime startT, endT;
            std::future<bool> patchWriteFail;
            /*****************************************************/
            
            /*************   Local variables for   ***************/
            /*************   Processing on slave   ***************/
            std::mutex mtx;
            std::vector<std::shared_ptr<SubImage>> subImages;
            redux::util::Array<double> phi_fixed;                           //!< The fixed aberrations for this channel (i.e. phase diversity)
            redux::util::Array<double> phi_channel;                         //!< The fixed part + tilt-corrections for this channel
            /*****************************************************/
            
            uint16_t ID;
            Point16 imgSize;
            Object& myObject;
            MomfbdJob& myJob;

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
