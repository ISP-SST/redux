#ifndef REDUX_MOMFBD_MOMFBDJOB_HPP
#define REDUX_MOMFBD_MOMFBDJOB_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/object.hpp"

#include "redux/momfbd/solver.hpp"

#include "redux/job.hpp"
#include "redux/util/array.hpp"
#include "redux/util/progresswatch.hpp"
#include "redux/util/region.hpp"

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

        enum Step : uint16_t { 
                    JSTEP_CHECKING = 1,
                    JSTEP_CHECKED,
                    JSTEP_PREPROCESS,
                    JSTEP_QUEUED,
                    JSTEP_RUNNING,
                    JSTEP_DONE,
                    JSTEP_VERIFY,
                    JSTEP_VERIFIED,
                    JSTEP_POSTPROCESS,
                    JSTEP_WRITING
                 // JSTEP_SUBMIT = 0x100,              Defined in job.hpp
                 // JSTEP_COMPLETED = 0x4000,          Defined in job.hpp
                 // JSTEP_ERR = 0x8000                 Defined in job.hpp
                  };
        
        public:
            static size_t jobType;

            static int staticInit(void);
            
            MomfbdJob( void );
            ~MomfbdJob( void );

            static Job* create(void) { return new MomfbdJob(); }

            size_t getTypeID(void) { return Job::MOMFBDJOB; }
            uint64_t unpackParts(const char* ptr, WorkInProgress::Ptr, bool);
            
            void parsePropertyTree( bpo::variables_map& vm, bpt::ptree& tree, redux::logging::Logger& );
            bpt::ptree getPropertyTree( bpt::ptree* root = nullptr );

            uint64_t size( void ) const;
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );

            size_t nImages(void) const;
            
            bool getWork( WorkInProgress::Ptr, uint16_t, const std::map<Job::StepID,Job::CountT>& );
            void ungetWork( WorkInProgress::Ptr );
            void failWork( WorkInProgress::Ptr );
            void returnResults( WorkInProgress::Ptr );

            void cleanup(void);
            bool run( WorkInProgress::Ptr, boost::asio::io_service&, uint16_t );
            
            void setLogChannel(std::string channel);
            
            bool mayBeDeleted(void);
            
            size_t memUsage(void);       //!< Approximate current memory usage of this job
            size_t diskUsage(void);      //!< Approximate current disk usage of this job
            size_t procUsage(void);      //!< Approximate memory usage for processing 1 part
        
            bool active(void);
            bool check(void);
            bool checkCacheUsage(void);
            bool checkOutputUsage(void);
            bool checkCfg(void);
            bool checkData(bool verbose=false);
            bool checkPre(void);
            bool checkPost(void);
            bool checkWriting(void);
            uint16_t getNextStep( uint16_t s=JSTEP_NONE ) const;
            const std::vector<std::shared_ptr<Object>>& getObjects(void) const { return objects; };

            const MomfbdJob& operator=(const GlobalCfg&);

        private:

            void checkParts( void );
            bool checkPatchPositions(void);
            void unloadCalib( boost::asio::io_service& );
            void preProcess( boost::asio::io_service&, uint16_t nThreads );
            void initCache( void );
            void clearPatches( void );
            void verifyPatches( void );
            void writeOutput( boost::asio::io_service& );
            void loadPatchResults( boost::asio::io_service&, uint16_t nThreads );
            void postProcess( boost::asio::io_service&, uint16_t nThreads );
            
            void updateProgressString(void);

            std::vector<std::shared_ptr<Object>> objects;

            redux::util::Array<PatchData::Ptr> patches;
            redux::util::Region16 roi;            // patch locations etc. are specified relative to the reference align-clip, this variable stores that region.
            
            GlobalData::Ptr globalData;
            Solver::Ptr solver;
            
            redux::util::ProgressWatch progWatch;
            
            bool cfgChecked;
            bool dataChecked;
        
            friend class Constraints;
            friend class Object;
            friend class Channel;
            friend class Solver;
            friend class SubImage;
            friend struct ModeSet;
            friend struct redux::image::Pupil;


        };

        /*! @} */

    }   // momfbd

    /*! @} */


}   // redux

#endif  // REDUX_MOMFBD_MOMFBDJOB_HPP
