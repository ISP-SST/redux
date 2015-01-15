#ifndef REDUX_MOMFBD_CHANNEL_HPP
#define REDUX_MOMFBD_CHANNEL_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/patch.hpp"

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
        /*! @brief Class containing the channel-specific configuration for a MomfbdJob/Object
         * 
         */
        class Channel : public ChannelCfg {

        public:

            Channel( Object&, MomfbdJob& );
            ~Channel();

            void parsePropertyTree( bpt::ptree& );
            bpt::ptree getPropertyTree( bpt::ptree&);

            size_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            double getMaxMean(void) const;
            void collectImages(redux::util::Array<float>&) const;
            void getPatchData(ChannelData::Ptr, uint16_t yid, uint16_t xid) const;
            void calcPatchPositions(const std::vector<uint16_t>&, const std::vector<uint16_t>&);
            
            void initWorkSpace( WorkSpace& ws );
            
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

            std::vector<redux::image::Statistics::Ptr> imageStats;

            redux::image::Image<float> images, dark, gain;
            redux::image::Image<float> ccdResponse, ccdScattering;
            redux::image::Image<float> psf, modulationMatrix;
            redux::image::Image<int16_t> xOffset, yOffset;
            
            /*
            Object& myObject;
            MomfbdJob& myJob;

            friend class Object;
            
        };

        /*! @} */
                
    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_CHANNEL_HPP
