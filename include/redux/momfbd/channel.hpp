#ifndef REDUX_MOMFBD_CHANNEL_HPP
#define REDUX_MOMFBD_CHANNEL_HPP

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

namespace po = boost::program_options;
namespace bpt = boost::property_tree;


namespace redux {

    namespace momfbd {

        class Object;
        class MomfbdJob;
        class Channel {

        public:

            Channel( Object&, MomfbdJob& );
            ~Channel();

            void parseProperties( bpt::ptree& tree );
            bpt::ptree getPropertyTree( bpt::ptree* root=nullptr );

        private:

            std::vector<uint32_t> imageNumbers, darkNumbers;
            std::vector<int16_t> alignClip;     // {xl,xh,yl,yh}
            std::vector<uint32_t> wf_num;
            std::string imageDataDir, filenameTemplate, darkTemplate, gainFile;
            std::string responseFile, backgainFile, psfFile, mmFile;
            std::string offxFile, offyFile;
            std::vector<double> stokesWeights;
            std::vector<double> diversity;
            std::vector<int> diversityOrders;
            std::vector<int> diversityTypes;

            uint32_t flags;
            uint8_t mmRow, mmWidth;
            uint8_t fillpix_method;
            int image_num_offs, sequenceNumber;
            double nf;

            Object& myObject;
            MomfbdJob& myJob;

            template <typename Archive>
            void serialize( Archive& ar, const unsigned int version ) {

                ar & imageNumbers & darkNumbers;
                ar & alignClip;     // {xl,xh,yl,yh}
                ar & wf_num;
                ar & imageDataDir & filenameTemplate & darkTemplate & gainFile;
                ar & responseFile & backgainFile & psfFile & mmFile;
                ar & offxFile & offyFile;
                ar & stokesWeights;
                ar & diversity;
                ar & diversityOrders;
                ar & diversityTypes;

                ar & flags;
                ar & mmRow & mmWidth;
                ar & fillpix_method;
                ar & image_num_offs & sequenceNumber;
                ar & nf;

            }

            friend class boost::serialization::access;
            friend class Object;
        };

    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_CHANNEL_HPP
