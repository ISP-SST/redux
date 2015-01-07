#ifndef REDUX_MOMFBD_MOMFBDJOB_HPP
#define REDUX_MOMFBD_MOMFBDJOB_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/object.hpp"
#include "redux/momfbd/patch.hpp"
#include "redux/job.hpp"
#include "redux/util/array.hpp"

#include <map>

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

namespace redux {

    namespace momfbd {

        /*! @defgroup momfbd MOMFBD
         *  @{
         */

        class Channel;
        /*! @brief Class containing the configuration settings for a MOMFBD job.
         *
         */
        class MomfbdJob : public Job, public GlobalCfg {

        enum Step { JSTEP_SUBMIT = 1,
                    JSTEP_PREPROCESS = 2,
                    JSTEP_QUEUED = 4,
                    JSTEP_RUNNING = 8,
                    JSTEP_POSTPROCESS = 16,
                    JSTEP_COMPLETED = 32,
                    JSTEP_ERR = 128 };
        
        public:
            static size_t jobType;

            MomfbdJob( void );
            ~MomfbdJob( void );

            uint64_t unpackParts(const char* ptr, std::vector<Part::Ptr>&, bool);
            
            void parsePropertyTree( bpo::variables_map& vm, bpt::ptree& tree );
            bpt::ptree getPropertyTree( bpt::ptree* root = nullptr );

            uint64_t size( void ) const;
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );

            bool getWork( WorkInProgress&, uint8_t );
            void ungetWork( WorkInProgress& );
            void returnResults( WorkInProgress& );

            void init(void);
            void cleanup(void);
            bool run( WorkInProgress&, boost::asio::io_service&, uint8_t );
            
            bool check(void);
            bool checkCfg(void);
            bool checkData(void);
            const std::vector<std::shared_ptr<Object>>& getObjects(void) const { return objects; };

            const MomfbdJob& operator=(const GlobalCfg&);

        private:

            void checkParts( void );
            void preProcess( boost::asio::io_service& );
            void applyLocalOffsets( PatchData::Ptr );
            void runMain( WorkSpace& );
            void storePatches( WorkInProgress&, boost::asio::io_service&, uint8_t );
            void postProcess( boost::asio::io_service& );

            std::vector<std::string> outputFiles;

            std::vector<std::shared_ptr<Object>> objects;

            redux::util::Array<double> pupil;
            redux::util::Array<float> imageStack;
            std::map<uint32_t, const ModeCache::ModePtr> modes;

            std::map<size_t,PatchData::Ptr> patches;
        
            friend class Object;
            friend class Channel;


        };

        const size_t dummy = MomfbdJob::jobType;       // this will trigger the registration of MomfbdJob in Job::jobMap

        /*! @} */

    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_MOMFBDJOB_HPP
