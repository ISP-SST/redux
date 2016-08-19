#ifndef REDUX_MOMFBD_MOMFBDJOB_HPP
#define REDUX_MOMFBD_MOMFBDJOB_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/object.hpp"

#include "redux/momfbd/solver.hpp"

#include "redux/job.hpp"
#include "redux/util/array.hpp"

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

namespace redux {

    namespace momfbd {
        class MomfbdJob;
    }
    
    /*!
    *  This is just a trick to trigger the registration of the JobType MomfbdJob by only including the header file.
    *  See below for the actual triggering.
    */
    template<typename Dummy>
    class StaticJobTypeInit<momfbd::MomfbdJob, Dummy> {
        public:
            static size_t const jobType;
    };
    
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
                    JSTEP_CHECKED = 1,
                    JSTEP_PREPROCESS = 2,
                    JSTEP_QUEUED = 4,
                    JSTEP_RUNNING = 8,
                    JSTEP_POSTPROCESS = 0x40,
                    JSTEP_WRITING = 0x80
                 // JSTEP_SUBMIT = 0x100,              Defined in job.hpp
                 // JSTEP_COMPLETED = 0x4000,          Defined in job.hpp
                 // JSTEP_ERR = 0x8000                 Defined in job.hpp
                  };
        
        public:
            static size_t jobType;

            MomfbdJob( void );
            ~MomfbdJob( void );

            static Job* create(void) { return new MomfbdJob(); }

            void setProgressString(void);
            uint64_t unpackParts(const char* ptr, WorkInProgress&, bool);
            
            void parsePropertyTree( bpo::variables_map& vm, bpt::ptree& tree );
            bpt::ptree getPropertyTree( bpt::ptree* root = nullptr );

            uint64_t size( void ) const;
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );

            size_t nImages(void) const;
            
            bool getWork( WorkInProgress&, uint16_t, bool );
            void ungetWork( WorkInProgress& );
            void failWork( WorkInProgress& );
            void returnResults( WorkInProgress& );

            void cleanup(void);
            bool run( WorkInProgress&, boost::asio::io_service&, uint16_t );
            
            void setLogChannel(std::string channel);
            
            bool active(void);
            bool check(void);
            bool checkCfg(void);
            bool checkData(void);
            const std::vector<std::shared_ptr<Object>>& getObjects(void) const { return objects; };

            const MomfbdJob& operator=(const GlobalCfg&);

        private:

            uint8_t checkParts( void );
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
            friend class SubImage;
            friend struct ModeSet;
            friend struct redux::image::Pupil;


        };

        /*! @} */

    }   // momfbd

    template<typename Dummy>
    size_t const StaticJobTypeInit<momfbd::MomfbdJob, Dummy>::jobType = Job::registerJob( "momfbd", momfbd::MomfbdJob::create );

    /*! @} */
    
    namespace momfbd {
        // this will trigger the registration of MomfbdJob in Job::jobMap
        const size_t MomfbdJobDummy RDX_UNUSED = StaticJobTypeInit<MomfbdJob, void>::jobType;
    }

}   // redux

#endif  // REDUX_MOMFBD_MOMFBDJOB_HPP
