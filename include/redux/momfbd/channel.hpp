#ifndef REDUX_MOMFBD_CHANNEL_HPP
#define REDUX_MOMFBD_CHANNEL_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/workspace.hpp"

#include <redux/image/image.hpp>
#include <redux/util/arraystats.hpp>

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

            uint16_t id(void) const { return ID; };
            size_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            double getMaxMean(void) const;
            void getFileNames(std::vector<std::string>&) const;
            size_t nImages(void) const { return images.size(); } 
            void getPatchData(ChannelData&, const PatchData&) const;
            
            
            /*************   Processing on slave   ***************/
            /*************         Methods         ***************/
            void initProcessing(WorkSpace::Ptr);
            void initPatch(ChannelData&);
            void getResults(ChannelData&);
            const std::vector<std::shared_ptr<SubImage>>& getSubImages(void) const { return subImages; };
            void writeAna(const redux::util::Array<PatchData::Ptr>&);
            void initPhiFixed(void);
            void computePhi(void);
            void addMode(redux::util::Array<double>&, uint16_t, double) const;
            void getPhi(redux::util::Array<double>&, const WaveFront&) const;
            //void addAllFT(redux::image::FourierTransform&);
            void addAllFT(redux::util::Array<double>&);
            double metric(void);
            /*****************************************************/

            void dump( std::string tag );
        private:
            
            bool checkCfg(void);
            bool checkData(void);
            
            void init(void);
            void initCache(void);
            void cleanup(void);
            
            void loadData(boost::asio::io_service&);
            void unloadData(void);
            void unloadCalib(void);

            void preprocessData(boost::asio::io_service&);
            void normalizeData(boost::asio::io_service&, double value);

            void loadImage(size_t index);
            void preprocessImage(size_t index, double avgMean);
            void normalizeImage(size_t index, double value);
            
            size_t sizeOfPatch(uint32_t) const;
            
            Point16 getImageSize(void);


            /*************   Local variables for   ***************/
            /************* PreProcessing on master ***************/
            std::vector<redux::util::ArrayStats::Ptr> imageStats;
            std::vector<redux::image::Image<float>> images;
            redux::image::Image<float> dark, gain;
            redux::image::Image<float> ccdResponse, ccdScattering;
            redux::image::Image<float> psf, modulationMatrix;
            redux::image::Image<int16_t> xOffset, yOffset;
            bpx::ptime startT, endT;
            /*****************************************************/
            
            /*************   Local variables for   ***************/
            /*************   Processing on slave   ***************/
            std::mutex mtx;
            WorkSpace::Ptr workspace;
            std::map<uint16_t, const PupilMode::Ptr> modes;                 //!< modes used in this channel
            std::pair<redux::util::Array<double>, double> pupil;            //!< pupil & area of pupil
            std::vector<std::shared_ptr<SubImage>> subImages;
            redux::util::Array<float> fittedPlane;
            redux::util::Array<double> phi_fixed;                           //!< The fixed aberrations for this channel (i.e. phase diversity)
            redux::util::Array<double> phi_channel;                         //!< The fixed part + tilt-corrections for this channel
            std::set<size_t> pupilIndices, otfIndices;                      //!< Arrays with the offsets where the pupil/otf are greater than some threshold
            std::set<std::pair<size_t,size_t>> pupilInOTF;                  //!< Arrays with the offsets where the pupil/img are located/centered in the OTF grid
            /*****************************************************/
            
            double frequencyCutoff, pupilRadiusInPixels;

            uint16_t ID;
            Point16 imgSize;
            Object& myObject;
            MomfbdJob& myJob;

            friend class Object;
            friend struct Tilts;
            friend class WorkSpace;
            friend struct SubImage;
            
        };

        /*! @} */
                
    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_CHANNEL_HPP
