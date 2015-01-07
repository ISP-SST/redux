#ifndef REDUX_MOMFBD_OBJECT_HPP
#define REDUX_MOMFBD_OBJECT_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/modecache.hpp"
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
        class Object : public ObjectCfg {

        public:

            Object( const MomfbdJob& );
            ~Object();

            void parsePropertyTree( bpt::ptree& );
            bpt::ptree getPropertyTree( bpt::ptree& );

            size_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            size_t nImages(size_t offset=0);
            void collectImages(redux::util::Array<float>&) const;
            
            const std::vector<std::shared_ptr<Channel>>& getChannels(void) const { return channels; };
            const MomfbdJob& getJob(void) const { return myJob; };
            
            void initWorkSpace( WorkSpace& ws );
        
            double lim_freq, r_c, cf2pix, pix2cf, cf2def;
            redux::util::Array<double> pupil;
            std::map<uint32_t, const ModeCache::ModePtr> modes;

            const MomfbdJob& myJob;
            std::vector<std::shared_ptr<Channel>> channels;
            
        private:

            
            bool checkCfg(void);
            bool checkData(void);

            void init(void);
            void cleanup(void);
            void loadData(boost::asio::io_service&);
            void preprocessData(boost::asio::io_service&);
            void normalize(boost::asio::io_service&);
            void prepareStorage(void);
            void storePatches( WorkInProgress&, boost::asio::io_service&, uint8_t );
            
            size_t sizeOfPatch(uint32_t) const;
            void applyLocalOffsets(PatchData::Ptr) const;
            
            Point16 clipImages(void);

            friend class MomfbdJob;
            friend class Channel;

        };


        /*! @} */
        
    }

}

#endif  // REDUX_MOMFBD_OBJECT_HPP
