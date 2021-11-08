#ifndef REDUX_MOMFBD_MOMFBDJOB_HPP
#define REDUX_MOMFBD_MOMFBDJOB_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/object.hpp"
#include "redux/momfbd/wavefront.hpp"

#include "redux/momfbd/solver.hpp"

#include "redux/job.hpp"
#include "redux/util/array.hpp"
#include "redux/util/progresswatch.hpp"
#include "redux/util/region.hpp"

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

// define RDX_DO_TRANSPOSE to enable old way of transposing input
//#define RDX_DO_TRANSPOSE

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

            size_t getTypeID(void) override { return Job::MOMFBDJOB; }
            uint64_t packParts(char* ptr, WorkInProgress::Ptr) const override;
            uint64_t unpackParts(const char* ptr, WorkInProgress::Ptr, bool) override;
            
            void parsePropertyTree( bpo::variables_map& vm, bpt::ptree& tree, redux::logging::Logger& ) override;
            bpt::ptree getPropertyTree( bpt::ptree* root = nullptr, bool showAll=false ) override;

            uint64_t size( void ) const override;
            uint64_t pack( char* ) const override;
            uint64_t unpack( const char*, bool ) override;
            void prePack( bool force=false ) override;

            size_t nImages(void) const;
            
            bool getWork( WorkInProgress::Ptr, uint16_t, const std::map<Job::StepID,Job::CountT>& ) override;
            void ungetWork( WorkInProgress::Ptr ) override;
            void failWork( WorkInProgress::Ptr ) override;
            void returnResults( WorkInProgress::Ptr ) override;

            void cleanup(void) override;
            bool run( WorkInProgress::Ptr, uint16_t ) override;
            static void setRunningCount( const CountT& c ) { Job::setStepCount(StepID(Job::MOMFBDJOB,JSTEP_RUNNING), c); }
            
            bool mayBeDeleted(void) override;
            
            size_t memUsage(void) override;       //!< Approximate current memory usage of this job
            size_t diskUsage(void) override;      //!< Approximate current disk usage of this job
            size_t procUsage(void) override;      //!< Approximate memory usage for processing 1 part
        
            bool active(void) override;
            bool check(void) override;
            bool checkCacheUsage(void);
            bool checkOutputUsage(void);
            bool checkCfg(void);
            bool checkData(bool verbose=false);
            bool checkPre(void);
            bool checkPost(void);
            bool checkWriting(void);
            uint16_t getNextStep( uint16_t s=JSTEP_NONE ) const override;
            const std::vector<Object::Ptr>& getObjects(void) const { return objects; };
            const Object::Ptr getObject( uint16_t id ) const;
            void addObject( Object::Ptr obj ) { objects.push_back(obj); }
            Object::Ptr addObject( void ) { Object::Ptr obj( new Object( *this, objects.size() ) ); addObject(obj); return obj; }
            void addTraceObject( Object::Ptr obj ) { trace_objects.push_back(obj); }
            const std::vector<Channel::Ptr>& getChannels(uint16_t objID) const;
            redux::util::Point16 getSmallestImageSize( void );

            const MomfbdJob& operator=(const GlobalCfg&);

        private:

            void checkParts( void );
            bool checkPatchPositions(void);
            void unloadCalib( void );
            void preProcess( void );
            void initCache( void );
            void clearPatches( void );
            void verifyPatches( void );
            void writeOutput( void );
            void loadPatchResults( void );
            int getReferenceObject( void );
            void generateTraceObjects( void );
            void generateTraceData(PatchData::Ptr);
            void postProcess( void );
            
            void dumpConfig(void);
            void updateProgressString(void);
            void updateStatus(void);

            std::vector<Object::Ptr> objects;
            std::vector<Object::Ptr> trace_objects;
            WaveFronts waveFronts;

            redux::util::Array<PatchData::Ptr> patches;
            redux::util::Region16 roi;            // patch locations etc. are specified relative to the reference align-clip, this variable stores that region.
            
            GlobalData::Ptr globalData;
            Solver::Ptr solver;
            
            redux::util::ProgressWatch progWatch;
            
            bool cfgChecked;
            bool dataChecked;
        
            friend struct Constraints;
            friend class WaveFronts;
            friend class Object;
            friend class Channel;
            friend struct Solver;
            friend struct SubImage;
            friend struct PatchData;
            friend struct ModeSet;
            friend struct redux::image::Pupil;


        };

        /*! @} */

    }   // momfbd

    /*! @} */


}   // redux

#endif  // REDUX_MOMFBD_MOMFBDJOB_HPP
