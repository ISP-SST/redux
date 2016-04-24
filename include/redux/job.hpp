#ifndef REDUX_JOB_HPP
#define REDUX_JOB_HPP

#include "redux/logger.hpp"
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

namespace bpo = boost::program_options;
namespace bpt = boost::property_tree;


namespace redux {
        
    /*! @ingroup redux
     *  @{
     */
    
    void runThreadsAndWait(boost::asio::io_service& service, uint16_t nThreads );

    struct WorkInProgress;
    
    class job_check_failed: public std::exception {
        
    public:
        explicit job_check_failed(const char* message) : msg_(message) { }
        explicit job_check_failed(const std::string& message) : msg_(message) { }
        virtual ~job_check_failed() throw () {}

        virtual const char* what() const throw () { return msg_.c_str(); }

    protected:
        std::string msg_;

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
        enum Step : uint8_t { JSTEP_NONE = 0,
                              JSTEP_SUBMIT = 1,
                              // Intermediate steps can be defined per job. The ones listed here are needed by the Deamon class.
                              JSTEP_COMPLETED = 64,
                              JSTEP_ERR = 128
                            };
        constexpr static uint8_t StepUserMask = 0x3E;      // the job-specific steps

        typedef std::shared_ptr<Job> JobPtr;
        struct JobCompare {
            bool operator()( const JobPtr &a, const JobPtr &b ) const { return ( *a < *b ); }
        };
        typedef std::set < JobPtr, JobCompare> JobSet;
        typedef Job* ( *JobCreator )( void );
        typedef std::map<std::string, std::pair<size_t, JobCreator>> MapT;

        static MapT& getMap(void) { static MapT m; return m; };
        static size_t registerJob( const std::string&, JobCreator f );
        static std::vector<JobPtr> parseTree( bpo::variables_map& vm, bpt::ptree& tree, bool check=false );
        static JobPtr newJob( const std::string& );
        
        virtual void setProgressString(void) { info.progressString = ""; };
        static std::string stateString(uint8_t);
        static std::string stateTag(uint8_t);
        
        struct Info {
            uint32_t id, timeout, maxProcessingTime;
            uint8_t priority, verbosity, maxPartRetries;
            uint16_t maxThreads;
            std::atomic<uint8_t> step;
            std::atomic<uint8_t> state;
            std::atomic<uint32_t> progress[2];
            std::string typeString, name, user, host;
            std::string progressString, logFile;
            std::string outputDir;                  //!< Where the output goes (defaults to current directory of rsub)
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
        
        virtual uint64_t unpackParts(const char* ptr, WorkInProgress& wip, bool) { wip.parts.clear(); return 0; };
        
        virtual void parsePropertyTree( bpo::variables_map&, bpt::ptree& );
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
        
        virtual bool getWork(WorkInProgress&, uint16_t, bool) { return false; };
        virtual void ungetWork(WorkInProgress&) { };
        virtual void failWork(WorkInProgress&) { };
        virtual void returnResults(WorkInProgress&) { };

        virtual void init(void) {};
        virtual void cleanup(void) {};
        virtual bool run(WorkInProgress&, boost::asio::io_service&, uint16_t) = 0;
        
        virtual void setLogChannel(std::string channel) { jobLogChannel = channel; };
        std::string getLogChannel(void) { return jobLogChannel; }
        void startLog(bool overwrite=false);
        void stopLog(void);
        
        bool operator<(const Job& rhs);
        bool operator!=(const Job& rhs);
    
    protected:
        
        std::mutex jobMutex;
        std::string cachePath;
        std::string jobLogChannel;
        std::shared_ptr<LogSink> jlog;

    };


    /*! @} */
        
}  // redux


#endif // REDUX_JOB_HPP
