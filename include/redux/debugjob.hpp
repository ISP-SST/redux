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

    /*! Dummy job for testing the framework. Also usable as a tutorial for how to implement a new/custom job.
     * 
     */
    class DebugJob : public Job {
        
        struct DebugPart : public Part {
            uint32_t xPixelL, xPixelH, yPixelL, yPixelH;
            double beginX, endX, beginY, endY;
            size_t sortedID;
            redux::util::Array<int64_t> result;      // use int64_t as temporary storage, cast to int16_t in post-processing
            size_t size(void) const;
            char* pack(char*) const;
            const char* unpack(const char*, bool);
        } dPart;
        typedef std::shared_ptr<DebugPart> PartPtr;
        
        const char* unpackParts(const char* ptr, std::vector<Part::Ptr>&, bool);

    public:

        static size_t jobType;

        DebugJob(void);
        ~DebugJob(void);

        void parseProperties( po::variables_map& vm, bpt::ptree& tree );
        bpt::ptree getPropertyTree( bpt::ptree* root=nullptr );
        
        size_t size(void) const;
        char* pack(char*) const;
        const char* unpack(const char*, bool);
        
        size_t getParts(WorkInProgress&);
        void ungetParts(WorkInProgress&);
        void returnParts(WorkInProgress&);
        
        bool run(WorkInProgress&);

    private:

        void preProcess(void);
        void runMain(Part::Ptr&,int);
        void postProcess(void);
        
        void checkParts(void);

        uint32_t maxIterations, patchSize;
        double gamma;
        uint32_t xSize, ySize;
        double coordinates[4];
        std::map<size_t,PartPtr> jobParts;

    };
   
    /*! @} */
    
    namespace debugjob {
        const size_t dummy = redux::DebugJob::jobType;       // this will trigger the registration of DebugJob in Job::jobMap
    }

}

#endif  // REDUX_DEBUGJOB_HPP
