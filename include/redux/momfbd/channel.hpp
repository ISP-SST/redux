#ifndef REDUX_MOMFBD_CHANNEL_HPP
#define REDUX_MOMFBD_CHANNEL_HPP

#include <redux/image/image.hpp>
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
        
        class Object;
        class MomfbdJob;
        /*! @brief Class containing the channel-specific configuration for a MomfbdJob/Object
         * 
         */
        class Channel {

        public:

            Channel( const Object&, const MomfbdJob& );
            ~Channel();

            void parseProperties( bpt::ptree& tree );
            bpt::ptree getPropertyTree( bpt::ptree* root=nullptr );

            size_t size(void) const;
            char* pack(char*) const;
            const char* unpack(const char*, bool);
        
        private:
            
            bool isValid(void);
            void loadData(boost::asio::io_service&, boost::thread_group&);
            void preprocessData(boost::asio::io_service&, boost::thread_group&);

            void loadImage(size_t index);
            void preprocessImage(size_t index, double avgMean);

            std::vector<uint32_t> imageNumbers, darkNumbers;
            std::vector<int16_t> alignClip;     // {xl,xh,yl,yh}
            std::vector<uint32_t> wf_num;
            std::string imageDataDir, imageTemplate, darkTemplate, gainFile;
            std::string responseFile, backgainFile, psfFile, mmFile;
            std::string offxFile, offyFile;
            std::vector<double> imageMeans;
            std::vector<double> stokesWeights;
            std::vector<double> diversity;
            std::vector<uint32_t> diversityOrders;
            std::vector<uint32_t> diversityTypes;

            uint32_t flags;
            uint8_t mmRow, mmWidth;
            uint8_t fillpix_method;
            uint32_t image_num_offs, sequenceNumber;
            double nf;

            redux::image::Image<float> images, dark, gain;
            redux::image::Image<float> ccdResponse, ccdScattering;
            redux::image::Image<float> psf, modulationMatrix;
            redux::image::Image<int16_t> xOffset, yOffset;
            
            const Object& myObject;
            const MomfbdJob& myJob;

            friend class Object;
        };

        /*! @} */
                
    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_CHANNEL_HPP
