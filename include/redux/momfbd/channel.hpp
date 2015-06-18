#ifndef REDUX_MOMFBD_CHANNEL_HPP
#define REDUX_MOMFBD_CHANNEL_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/workspace.hpp"

#include <redux/image/image.hpp>
#include <redux/image/statistics.hpp>

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
        class WorkSpace;
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

            size_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            double getMaxMean(void) const;
            void collectImages(redux::util::Array<float>&) const;
            void getPatchData(ChannelData&, Point16 patchID) const;
            void calcPatchPositions(const std::vector<uint16_t>&, const std::vector<uint16_t>&);
            
            /*************   Processing on slave   ***************/
            /*************         Methods         ***************/
            void initProcessing(WorkSpace::Ptr);
            void initPatch(ChannelData&);
            void initPhiFixed(void);
            void computePhi(void);
            void addMode(redux::util::Array<double>&, uint16_t, double) const;
            void getPhi(redux::util::Array<double>&, const WaveFront&) const;
            //void addAllFT(redux::image::FourierTransform&);
            void addAllFT(redux::util::Array<double>&);
            double metric(void);
            /*****************************************************/
            
            uint32_t dataOffset;

        private:
            
            bool checkCfg(void);
            bool checkData(void);
            
            void init(void);
            void initCache(void);
            void cleanup(void);
            
            size_t nImages(size_t offset=0) { dataOffset=offset; return images.dimSize(0); } 

            void loadData(boost::asio::io_service&);
            void preprocessData(boost::asio::io_service&);
            void normalizeData(boost::asio::io_service&, double value);

            void loadImage(size_t index);
            void preprocessImage(size_t index, double avgMean);
            void normalizeImage(size_t index, double value);
            
            size_t sizeOfPatch(uint32_t) const;
            
            Point16 clipImages(void);


            /*************   Local variables for   ***************/
            /************* PreProcessing on master ***************/
            std::vector<redux::image::Statistics::Ptr> imageStats;
            redux::image::Image<float> images, dark, gain;
            redux::image::Image<float> ccdResponse, ccdScattering;
            redux::image::Image<float> psf, modulationMatrix;
            redux::image::Image<int16_t> xOffset, yOffset;
            /*****************************************************/
            
            /*************   Local variables for   ***************/
            /*************   Processing on slave   ***************/
            std::mutex mtx;
            WorkSpace::Ptr workspace;
            std::map<uint16_t, const PupilMode::Ptr> modes;                 //!< modes used in this channel
            std::pair<redux::util::Array<double>, double> pupil;            //!< pupil & area of pupil
            std::vector<std::shared_ptr<SubImage>> subImages;
            redux::image::Image<float> fittedPlane;
            redux::util::Array<double> phi_fixed;                           //!< The fixed aberrations for this channel (i.e. phase diversity)
            redux::util::Array<double> phi_channel;                         //!< The fixed part + tilt-corrections for this channel
            std::set<size_t> pupilIndices, otfIndices;                      //!< Arrays with the offsets where the pupil/otf are greater than some threshold
            std::set<std::pair<size_t,size_t>> pupilInOTF;                  //!< Arrays with the offsets where the pupil/img are located/centered in the OTF grid
            /*****************************************************/
            
            double frequencyCutoff, pupilRadiusInPixels;

            uint16_t ID;
            Object& myObject;
            MomfbdJob& myJob;

            friend class Object;
            friend struct Tilts;
            friend struct SubImage;
            
        };

        /*! @} */
                
    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_CHANNEL_HPP
