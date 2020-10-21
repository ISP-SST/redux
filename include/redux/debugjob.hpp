#ifndef REDUX_DEBUGJOB_HPP
#define REDUX_DEBUGJOB_HPP

#include "redux/job.hpp"
#include "redux/work.hpp"
#include "redux/util/array.hpp"

#include <complex>

namespace redux {

    /*! @ingroup redux
     *  @{
     */

    
    /*!
     *  Dummy job for testing the framework. Also usable as a tutorial for how to implement a new/custom job.
     */
    class DebugJob : public Job {
        
        enum Step { JSTEP_SUBMIT=1, JSTEP_PREPROCESS=2, JSTEP_QUEUED=4, JSTEP_RUNNING=8, JSTEP_POSTPROCESS=16, JSTEP_COMPLETED=32, JSTEP_ERR=255 };
        
        struct DebugPart : public Part {
            uint32_t xPixelL, xPixelH, yPixelL, yPixelH;
            double beginX, endX, beginY, endY;
            size_t sortedID;
            redux::util::Array<int64_t> result;      // use int64_t as temporary storage, cast to int16_t in post-processing
            size_t size(void) const override;
            uint64_t pack(char*) const override;
            uint64_t unpack(const char*, bool) override;
        };
        typedef std::shared_ptr<DebugPart> PartPtr;
        
        uint64_t unpackParts(const char*, WorkInProgress::Ptr, bool) override;

    public:

        DebugJob(void);
        DebugJob( const DebugJob& ) = delete;
        ~DebugJob(void);
        
        static Job* create(void) { return new DebugJob(); }

        void parsePropertyTree( bpo::variables_map& vm, bpt::ptree& tree, redux::logging::Logger& ) override;
        bpt::ptree getPropertyTree( bpt::ptree* root=nullptr, bool showAll=false ) override;
        
        size_t size(void) const override;
        uint64_t pack(char*) const override;
        uint64_t unpack(const char*, bool) override;
        
        bool check(void) override;
        
        size_t getTypeID(void) override { return Job::DEBUGJOB; }
        bool getWork( WorkInProgress::Ptr, uint16_t, const std::map<Job::StepID,Job::CountT>& ) override;
        void ungetWork(WorkInProgress::Ptr) override;
        void returnResults(WorkInProgress::Ptr) override;
        
        bool run( WorkInProgress::Ptr, uint16_t ) override;

    private:

        void preProcess(void);
        void runMain(Part::Ptr&);
        void postProcess(void);
        
        void checkParts(void);
        uint32_t maxIterations, patchSize;
        double gamma;
        uint32_t xSize, ySize;
        double coordinates[4];
        std::map<size_t,PartPtr> jobParts;

    };
   
    /*! @} */


}

#endif  // REDUX_DEBUGJOB_HPP
