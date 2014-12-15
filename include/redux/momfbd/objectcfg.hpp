#ifndef REDUX_MOMFBD_OBJECTCFG_HPP
#define REDUX_MOMFBD_OBJECTCFG_HPP

#include "redux/momfbd/channelcfg.hpp"
#include "redux/momfbd/patch.hpp"
#include "redux/momfbd/workspace.hpp"
#include "redux/util/array.hpp"
#include "redux/types.hpp"

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
        
        /*! @brief Class containing the object-specific configuration a MOMFBD job
         * 
         */
        class ObjectCfg {

        public:

            ObjectCfg( const MomfbdJob& );
            ~ObjectCfg();

            void parseProperties( bpt::ptree& tree, const std::string& fn );
            bpt::ptree getPropertyTree( bpt::ptree* root=nullptr );

            size_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            size_t nImages(size_t offset=0);
            void collectImages(redux::util::Array<float>&) const;
            
            const std::vector<std::shared_ptr<ChannelCfg>>& getChannels(void) const { return channels; };
            const MomfbdJob& getJob(void) const { return myJob; };
            
            void initWorkSpace( WorkSpace& ws );
        
            uint8_t fillpix_method, output_data_type;
            uint32_t objectSize, sequenceNumber, objectPupilSize;
            uint32_t flags;
            double reg_gamma, weight, angle, wavelength;
            double lim_freq, r_c, cf2pix, pix2cf, cf2def;

            std::vector<uint32_t> imageNumbers, sequenceNumbers, darkNumbers;
            std::vector<uint32_t> wf_num;
            std::vector<double> stokesWeights;
            
            std::string imageDataDir, outputFileName, pupilFile;

            redux::util::Array<double> pupil;
            std::map<uint32_t, const ModeCache::ModePtr> modes;

            const MomfbdJob& myJob;
            std::vector<std::shared_ptr<ChannelCfg>> channels;
            
        private:

            
            bool checkCfg(void);
            bool checkData(void);

            void init(void);
            void cleanup(void);
            void loadData(boost::asio::io_service&);
            void preprocessData(boost::asio::io_service&);
            void normalize(boost::asio::io_service&);
            void prepareStorage(void);
            
            size_t sizeOfPatch(uint32_t) const;
            void applyLocalOffsets(PatchData::Ptr) const;
            
            Point clipImages(void);

            friend class MomfbdJob;
            friend class ChannelCfg;

        };


        /*! @} */
        
    }

}

#endif  // REDUX_MOMFBD_OBJECTCFG_HPP
