#ifndef REDUX_MOMFBD_OBJECT_HPP
#define REDUX_MOMFBD_OBJECT_HPP

#include "redux/momfbd/channel.hpp"
#include "redux/serialization.hpp"
#include "redux/util/array.hpp"

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

namespace po = boost::program_options;
namespace bpt = boost::property_tree;


namespace redux {

    namespace momfbd {

        class MomfbdJob;
        class Object {

        public:

            Object( MomfbdJob& );
            ~Object();

            void parseProperties( bpt::ptree& tree, const std::string& fn );
            bpt::ptree getProperties( bpt::ptree* );

            std::vector<uint32_t> imageNumbers, sequenceNumbers, darkNumbers;
            double reg_gamma, weight, angle, lambda;
            int nPoints, sequenceNumber, nph;
            std::vector<double> stokesWeights;

        private:

            uint32_t flags;

            std::vector<std::shared_ptr<Channel>> channels;
            std::vector<uint32_t> wf_num;
            std::string imageDataDir, outputFileName;
            uint8_t fillpix_method;
            double lim_freq, r_c;
            MomfbdJob& myJob;
            redux::util::Array<double> pupil;
            uint8_t output_data_type;

            template <typename Archive>
            void serialize( Archive& ar, const unsigned int version ) {

                std::string tmp = "Object";
                ar & tmp;
                ar & imageNumbers & sequenceNumbers & darkNumbers;
                ar & reg_gamma & weight & angle & lambda;
                ar & nPoints & sequenceNumber & nph;
                ar & stokesWeights;
                ar & flags;

                ar & channels;
                ar & wf_num;
                ar & imageDataDir & outputFileName;
                ar & fillpix_method;
                ar & lim_freq & r_c;
                ar & output_data_type;

                // TODO
                //redux::util::Array<double> pupil;

            }

            friend class boost::serialization::access;
            friend class Channel;

        };

    }

}

#endif  // REDUX_MOMFBD_OBJECT_HPP
