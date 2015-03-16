#ifndef REDUX_MOMFBD_OBJECT_HPP
#define REDUX_MOMFBD_OBJECT_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/cache.hpp"

#include "redux/util/array.hpp"
#include "redux/work.hpp"

#include <boost/asio.hpp>
#include <boost/property_tree/ptree.hpp>
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

            Object( MomfbdJob& );
            ~Object();

            void parsePropertyTree( bpt::ptree& );
            bpt::ptree getPropertyTree( bpt::ptree& );

            size_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            size_t nImages(size_t offset=0);
            void collectImages(redux::util::Array<float>&) const;
            void calcPatchPositions(const std::vector<uint16_t>&, const std::vector<uint16_t>&);
            
            const std::vector<std::shared_ptr<Channel>>& getChannels(void) const { return channels; };
            const MomfbdJob& getJob(void) const { return myJob; };
            
            void initWorkSpace( WorkSpace& ws );
        
            MomfbdJob& myJob;
            std::vector<std::shared_ptr<Channel>> channels;
            
        private:

            
            bool checkCfg(void);
            bool checkData(void);

            void init(void);
            void initCache(void);
            void cleanup(void);
            void loadData(boost::asio::io_service&);
            void preprocessData(boost::asio::io_service&);
            void normalize(boost::asio::io_service&);
            void prepareStorage(void);
            void storePatches( WorkInProgress&, boost::asio::io_service&, uint8_t );
            
            size_t sizeOfPatch(uint32_t) const;
            
            Point16 clipImages(void);

            friend class MomfbdJob;
            friend class Channel;

        };


        /*! @} */
        
    }

}

#endif  // REDUX_MOMFBD_OBJECT_HPP
