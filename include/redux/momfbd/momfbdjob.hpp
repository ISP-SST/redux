#ifndef REDUX_MOMFBD_MOMFBDJOB_HPP
#define REDUX_MOMFBD_MOMFBDJOB_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/object.hpp"

#include "redux/momfbd/solver.hpp"

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

            static Job* create(void) { return new MomfbdJob(); }

            uint64_t unpackParts(const char* ptr, WorkInProgress&, bool);
            
            void parsePropertyTree( bpo::variables_map& vm, bpt::ptree& tree );
            bpt::ptree getPropertyTree( bpt::ptree* root = nullptr );

            uint64_t size( void ) const;
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );

            size_t nImages(void) const;
            
            bool getWork( WorkInProgress&, uint16_t );
            void ungetWork( WorkInProgress& );
            void returnResults( WorkInProgress& );

            bool run( WorkInProgress&, boost::asio::io_service&, uint16_t );
            
            bool check(void);
            bool checkCfg(void);
            bool checkData(void);
            const std::vector<std::shared_ptr<Object>>& getObjects(void) const { return objects; };

            const MomfbdJob& operator=(const GlobalCfg&);

        private:

            void checkParts( void );
            void preProcess( boost::asio::io_service&, uint16_t nThreads );
            void initCache( void );
            void storePatches( WorkInProgress&, boost::asio::io_service&, uint8_t );
            void postProcess( boost::asio::io_service&, uint16_t nThreads );

            std::vector<std::shared_ptr<Object>> objects;

            redux::util::Array<PatchData::Ptr> patches;
            
            GlobalData::Ptr globalData;
            Solver::Ptr solver;
        
            friend class Object;
            friend class Channel;
            friend class Solver;
            friend struct ModeSet;
            friend struct redux::image::Pupil;


        };

        const MomfbdJob MomfbdJobDummy UNUSED;       // this will trigger the registration of MomfbdJob in Job::jobMap

        /*! @} */

    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_MOMFBDJOB_HPP
