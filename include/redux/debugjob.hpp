#ifndef REDUX_DEBUGJOB_HPP
#define REDUX_DEBUGJOB_HPP

#include "redux/job.hpp"


namespace redux {

    /*! @ingroup redux
     *  @{
     */

    /*! "Dummy" job for testing the framework.
     * 
     */
    class DebugJob : public Job {

    public:

        static size_t jobType;

        DebugJob(void);
        ~DebugJob(void);

        void parseProperties( po::variables_map& vm, bpt::ptree& tree );
        bpt::ptree getPropertyTree( bpt::ptree* root=nullptr ) { return bpt::ptree(); };
        
        size_t size(void) const;
        char* pack(char*) const;
        const char* unpack(const char*, bool);
        
    private:
        void* prePart( void );
        void postPart( void* );
        void* runPreJob( void );
        void runPostJob( void* );

        uint32_t preProcess( void );
        uint32_t postProcess( void );
        uint32_t runJob( void );

    };
   
    /*! @} */
    
    namespace debugjob {
        const size_t dummy = redux::DebugJob::jobType;       // this will trigger the registration of DebugJob in Job::jobMap
    }

}

#endif  // REDUX_DEBUGJOB_HPP
