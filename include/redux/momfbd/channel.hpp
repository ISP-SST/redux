#ifndef REDUX_MOMFBD_CHANNEL_HPP
#define REDUX_MOMFBD_CHANNEL_HPP

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

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

            Channel( Object&, MomfbdJob& );
            ~Channel();

            void parseProperties( bpt::ptree& tree );
            bpt::ptree getPropertyTree( bpt::ptree* root=nullptr );

            size_t size(void) const;
            char* pack(char*) const;
            const char* unpack(const char*, bool);
        
        private:

            std::vector<uint32_t> imageNumbers, darkNumbers;
            std::vector<int16_t> alignClip;     // {xl,xh,yl,yh}
            std::vector<uint32_t> wf_num;
            std::string imageDataDir, filenameTemplate, darkTemplate, gainFile;
            std::string responseFile, backgainFile, psfFile, mmFile;
            std::string offxFile, offyFile;
            std::vector<double> stokesWeights;
            std::vector<double> diversity;
            std::vector<uint32_t> diversityOrders;
            std::vector<uint32_t> diversityTypes;

            uint32_t flags;
            uint8_t mmRow, mmWidth;
            uint8_t fillpix_method;
            uint32_t image_num_offs, sequenceNumber;
            double nf;

            Object& myObject;
            MomfbdJob& myJob;

            friend class Object;
        };

        /*! @} */
                
    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_CHANNEL_HPP
