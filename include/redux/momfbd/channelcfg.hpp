#ifndef REDUX_MOMFBD_CHANNELCFG_HPP
#define REDUX_MOMFBD_CHANNELCFG_HPP

#include "redux/momfbd/patch.hpp"

#include <redux/image/image.hpp>
#include <redux/image/statistics.hpp>
#include <redux/types.hpp>

#include <boost/asio.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/thread.hpp>

namespace po = boost::program_options;
namespace bpt = boost::property_tree;


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
        
        class MomfbdJob;
        class ObjectCfg;
        class WorkSpace;
        /*! @brief Class containing the channel-specific configuration for a MomfbdJob/Object
         * 
         */
        class ChannelCfg {

        public:

            ChannelCfg( const ObjectCfg&, const MomfbdJob& );
            ~ChannelCfg();

            void parseProperties( bpt::ptree& tree );
            bpt::ptree getPropertyTree( bpt::ptree* root=nullptr );

            size_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            double getMaxMean(void) const;
            void collectImages(redux::util::Array<float>&) const;
            
            void initWorkSpace( WorkSpace& ws );

            std::vector<uint32_t> imageNumbers;
            std::vector<uint32_t> wf_num;
            std::vector<double> stokesWeights;
            std::vector<double> diversity;
            std::vector<uint32_t> diversityOrders;
            std::vector<uint32_t> diversityTypes;

            uint8_t fillpix_method, mmRow, mmWidth, incomplete;
            uint32_t flags, image_num_offs, sequenceNumber, imageOffset;
            double noiseFudge;
            
        private:
            
            bool checkCfg(void);
            bool checkData(void);
            size_t nImages(size_t offset=0) { imageOffset=offset; return images.dimSize(0); } 

            void loadData(boost::asio::io_service&);
            void preprocessData(boost::asio::io_service&);
            void normalizeData(boost::asio::io_service&, double value);

            void loadImage(size_t index);
            void preprocessImage(size_t index, double avgMean);
            void normalizeImage(size_t index, double value);
            
            size_t sizeOfPatch(uint32_t) const;
            void applyLocalOffsets(PatchData::Ptr) const;
            
            Point clipImages(void);

            std::vector<uint32_t> darkNumbers;
            std::vector<int16_t> alignClip;     // {firstX,lastX,firstY,lastY}
            std::string imageDataDir, imageTemplate, darkTemplate, gainFile;
            std::string responseFile, backgainFile, psfFile, mmFile;
            std::string offxFile, offyFile;
            std::vector<redux::image::Statistics::Ptr> imageStats;

            redux::image::Image<float> images, dark, gain;
            redux::image::Image<float> ccdResponse, ccdScattering;
            redux::image::Image<float> psf, modulationMatrix;
            redux::image::Image<int16_t> xOffset, yOffset;
            
            const ObjectCfg& myObject;
            const MomfbdJob& myJob;

            friend class ObjectCfg;
            
        };

        /*! @} */
                
    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_CHANNELCFG_HPP
