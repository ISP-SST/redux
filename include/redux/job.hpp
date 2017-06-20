#ifndef REDUX_JOB_HPP
#define REDUX_JOB_HPP

#include "redux/logging/logger.hpp"
#include "redux/work.hpp"

#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include <boost/asio.hpp>
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/thread/thread.hpp>

namespace bpo = boost::program_options;
namespace bpt = boost::property_tree;


namespace redux {
        
    /*! @ingroup redux
     *  @{
     */
    
    void runThreadsAndWait(boost::asio::io_service& service, uint16_t nThreads );
    
    class job_error: public std::exception {
        
    public:
        explicit job_error(const char* message) : msg_(message) { }
        explicit job_error(const std::string& message) : msg_(message) { }
        virtual ~job_error() throw () {}

        virtual const char* what() const throw () { return msg_.c_str(); }

    protected:
        std::string msg_;

    };
        
    class job_fatal: public job_error {


    };
        

    
    /*! Base class for a "job" to be processed by the redux framework.
     * 
     */
    class Job {
        
    public:

        enum State : uint8_t { JSTATE_NONE      = 0,
                               JSTATE_IDLE      = 1,
                               JSTATE_ACTIVE    = 2,    
                               JSTATE_PAUSED    = 4,
                               JSTATE_CANCELLED = 64,
                               JSTATE_ERR       = 128           // Error flag. Should not be used, dynamic error handling is better...
                             };
        enum Step : uint16_t { JSTEP_SUBMIT = 0,
                               // Intermediate steps can be defined per job. The ones listed here are needed by the Deamon class.
                               JSTEP_NONE = 0x100,
                               JSTEP_COMPLETED = 0x4000,
                               JSTEP_ERR = 0x8000
                            };
        constexpr static uint16_t StepUserMask = 0xFF;      // the job-specific steps

        typedef std::shared_ptr<Job> JobPtr;
        struct JobCompare {
            bool operator()( const JobPtr &a, const JobPtr &b ) const { return ( *a < *b ); }
        };
        typedef std::set < JobPtr, JobCompare> JobSet;
        typedef Job* ( *JobCreator )( void );
        typedef std::map<std::string, std::pair<size_t, JobCreator>> MapT;

        static MapT& getMap(void) { static MapT m; return m; };
        static size_t registerJob( const std::string&, JobCreator f );
        static std::vector<JobPtr> parseTree( bpo::variables_map& vm, bpt::ptree& tree, redux::logging::Logger&, bool check=false );
        static JobPtr newJob( const std::string& );
        
        static std::string stateString(uint8_t);
        static std::string stateTag(uint8_t);
        
        struct Info {
            uint32_t id, timeout, maxProcessingTime;
            uint8_t priority, verbosity, maxPartRetries;
            uint16_t maxThreads;
            std::atomic<uint16_t> step;
            std::atomic<uint8_t> state;
            std::string typeString, name, user, host;
            std::string logFile;
            std::string outputDir;                  //!< Where the output goes (defaults to current directory of rsub)
            std::string progressString;
            boost::posix_time::ptime submitTime;
            boost::posix_time::ptime startedTime;
            boost::posix_time::ptime completedTime;
            Info();
            Info(const Info&);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            static std::string printHeader(void);
            std::string print(void);
        } info;
        
        virtual size_t getTypeID(void) { return 0; };
        virtual uint64_t unpackParts(const char* ptr, WorkInProgress::Ptr wip, bool) { wip->parts.clear(); return 0; };
        
        virtual void parsePropertyTree( bpo::variables_map&, bpt::ptree&, redux::logging::Logger&);
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

        virtual uint64_t size(void) const;
        virtual uint64_t pack(char*) const;
        virtual uint64_t unpack(const char*, bool);

        virtual bool active(void) { return false; };
        virtual bool check(void) { return false; };         //! will be called several times during processing, should return true if all is ok.
        
        virtual bool getWork(WorkInProgress::Ptr, uint16_t, const std::map<uint16_t,uint16_t>&) { return false; };
        virtual void ungetWork(WorkInProgress::Ptr) { };
        virtual void failWork(WorkInProgress::Ptr) { };
        virtual void returnResults(WorkInProgress::Ptr) { };

        virtual void init(void) {};
        virtual void cleanup(void) {};
        virtual bool run(WorkInProgress::Ptr, boost::asio::io_service&, uint16_t) = 0;
        
        std::string cfg(void);
        
        virtual void setLogChannel(std::string channel) { jobLogChannel = channel; };
        std::string getLogChannel(void) { return jobLogChannel; }
        redux::logging::Logger& getLogger(void) { return logger; }
        void startLog(bool overwrite=false);
        void printJobInfo(void);
        void stopLog(void);
        
        std::unique_lock<std::mutex> getLock(bool trylock=false) {
            if(trylock) return std::move( std::unique_lock<std::mutex>(jobMutex,std::try_to_lock) );
            return std::move( std::unique_lock<std::mutex>(jobMutex) );
        }
        
        virtual bool mayBeDeleted(void) { return true; }
        
        bool operator<(const Job& rhs);
        bool operator!=(const Job& rhs);
    
    protected:
       
        std::mutex jobMutex;
        std::string cachePath;
        std::string jobLogChannel;
        redux::logging::Logger logger;
        
        friend class Daemon;
        friend class Worker;

    };

    template<class JobType, typename Dummy>
    class StaticJobTypeInit { };

    /*! @} */
        
}  // redux


#endif // REDUX_JOB_HPP
