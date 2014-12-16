#ifndef REDUX_JOB_HPP
#define REDUX_JOB_HPP

#include "redux/work.hpp"

#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <string>

#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/thread/thread.hpp>

namespace po = boost::program_options;
namespace bpt = boost::property_tree;


namespace redux {
        
    /*! @ingroup redux
     *  @{
     */

    struct WorkInProgress;
    
    /*! Base class for a "job" to be processed by the redux framework.
     * 
     */
    class Job {
        
    public:

        enum State : uint8_t { JSTATE_IDLE      = 1,
                               JSTATE_ACTIVE    = 2,    
                               JSTATE_PAUSED    = 4,
                               JSTATE_CANCELLED = 64,
                               JSTATE_ERR       = 128           // Error flag. Should not be used, dynamic error handling is better...
                             };
        typedef std::shared_ptr<Job> JobPtr;
        struct JobCompare {
            bool operator()( const JobPtr &a, const JobPtr &b ) const { return ( *a < *b ); }
        };
        typedef std::set < JobPtr, JobCompare> JobSet;
        typedef Job* ( *JobCreator )( void );
        typedef std::map<std::string, std::pair<size_t, JobCreator>> MapT;

        static MapT& getMap(void) { static MapT m; return m; };
        static size_t registerJob( const std::string&, JobCreator f );
        static std::vector<JobPtr> parseTree( po::variables_map& vm, bpt::ptree& tree );
        static JobPtr newJob( const std::string& );
        
        virtual std::string stepString(uint8_t) { return "-"; };
        static std::string stateString(uint8_t);
        static std::string stateTag(uint8_t);
        
        struct Info {
            size_t id;
            uint8_t priority, verbosity, maxThreads, maxPartRetries;
            std::atomic<uint8_t> step;
            std::atomic<uint8_t> state;
            std::string typeString, name, user, host;
            std::string logFile;
            boost::posix_time::ptime submitTime;
            Info();
            size_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            static std::string printHeader(void);
            std::string print(void);
        } info;
        
        virtual uint64_t unpackParts(const char* ptr, std::vector<Part::Ptr>& p, bool) { p.clear(); return 0; };
        
        virtual void parseProperties( po::variables_map&, bpt::ptree& );
        /*! @brief Returns a boost::property_tree containing the settings for this job.
         *  @details If a ptree pointer is provided, the tree will also be appended to it,
         *  inside a named block
         *  @code
         *  "jobtag"
         *  {
                // settings //
            }
         *  @endcode
         */
        virtual bpt::ptree getPropertyTree( bpt::ptree* root=nullptr );

        Job(void);
        Job(const Job&) = delete;
        virtual ~Job(void);

        virtual size_t size(void) const;
        virtual uint64_t pack(char*) const;
        virtual uint64_t unpack(const char*, bool);

        virtual bool check(void) { return false; };         //! will be called several times during processing, should return true if all is ok.
        
        virtual bool getWork(WorkInProgress&, uint8_t) { return false; };
        virtual void ungetWork(WorkInProgress&) { };
        virtual void returnResults(WorkInProgress&) { };

        virtual void init(void) {};
        virtual void cleanup(void) {};
        virtual bool run(WorkInProgress&, boost::asio::io_service&, uint8_t) = 0;
        
        bool operator<(const Job& rhs);
        bool operator!=(const Job& rhs);
    
    protected:
        
        std::mutex jobMutex;

    };


    /*! @} */
        
}  // redux


#endif // REDUX_JOB_HPP
