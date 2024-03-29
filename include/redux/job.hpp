#ifndef REDUX_JOB_HPP
#define REDUX_JOB_HPP

#include "redux/logging/logger.hpp"
#include "redux/work.hpp"

#ifdef RDX_TRACE_JOB
#   include "redux/util/trace.hpp"
#endif

#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#ifndef Q_MOC_RUN
#include <boost/asio.hpp>
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/thread/thread.hpp>
#endif

namespace bpo = boost::program_options;
namespace bpt = boost::property_tree;

#define RDX_JOB_PROGSTRING_LENGTH 30

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
    class Job
#ifdef RDX_TRACE_JOB
            : public redux::util::TraceObject<Job>
#endif
    {
    public:

        enum JobTypes : size_t { DEBUGJOB = 1,
                                 MOMFBDJOB
                               };
        
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
        enum Flags : uint16_t { CHECKED = 1,
                                NOCHECK = 2
                            };
        constexpr static uint16_t StepUserMask = 0xFF;      // the job-specific steps

        typedef std::shared_ptr<Job> JobPtr;
        struct JobCompare {
            bool operator()( const JobPtr &a, const JobPtr &b ) const { return ( *a < *b ); }
        };
        typedef std::set < JobPtr, JobCompare> JobSet;
        typedef Job* ( *JobCreator )( void );
        typedef std::map<std::string, std::pair<size_t, JobCreator>> MapT;

        struct StepID {
            StepID( size_t i, uint16_t s ) : id(i), step(s) {}
            bool operator<( const StepID& rhs ) const {
                if(!id || !rhs.id || (id==rhs.id)) return (step<rhs.step);
                return (id<rhs.id);
            }
            const size_t id;
            const uint16_t step;
        };
        struct CountT {
            CountT( void ) : min(0), max(std::numeric_limits<int64_t>::max()), active(0) {}
            explicit CountT( int64_t a ) : min(0), max(a), active(0) {}
            CountT( int64_t a, int64_t b ) : min(a), max(b), active(0) {}
            int64_t min, max;
            int64_t active;
        };
        static std::map<StepID,CountT> counts;
        static MapT jobCreators;
        static MapT& getMap(void) { static MapT m; return m; };
        static size_t registerJob( const std::string&, JobCreator f );
        static std::vector<JobPtr> parseTree( bpo::variables_map& vm, bpt::ptree& tree, redux::logging::Logger&, bool check=false );
        static JobPtr newJob( const std::string& );
        
        static std::string stateString(uint8_t);
        static std::string stateTag(uint8_t);
        
        static void setStepCount( const StepID& s, const CountT& c ) { counts[ s ] = c; };

        
        struct Info {
            uint32_t id, timeout, maxProcessingTime;
            uint8_t priority, verbosity, maxPartRetries;
            uint16_t maxThreads;
            uint16_t flags;
            std::atomic<uint16_t> step;
            std::atomic<uint8_t> state;
            std::string typeString, name, user, host;
            std::string logFile;
            std::string outputDir;                  //!< Where the output goes (defaults to current directory of rsub)
            char progressString[RDX_JOB_PROGSTRING_LENGTH];     // fixed size to avoid crashed when it's updated while packing/sending.
            boost::posix_time::ptime submitTime;
            boost::posix_time::ptime startedTime;
            boost::posix_time::ptime completedTime;
            std::map<uint16_t,boost::posix_time::ptime> times;
            Info();
            Info(const Info&);
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            static std::string printHeader(void);
            std::string print(void);
        } info;
        
        virtual size_t getTypeID(void) { return 0; };
        virtual uint64_t packParts( char* , WorkInProgress::Ptr ) const { return 0; };
        virtual uint64_t unpackParts( const char* , WorkInProgress::Ptr wip, bool ) { wip->parts.clear(); return 0; };
        
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
        virtual bpt::ptree getPropertyTree( bpt::ptree* root=nullptr, bool showAll=false );

        Job(void);
        Job(const Job&) = delete;
        virtual ~Job(void);

        virtual uint64_t size(void) const;
        virtual uint64_t pack(char*) const;
        virtual uint64_t unpack(const char*, bool);
        virtual void prePack( bool force=false ) {};

        virtual bool active(void) { return false; };
        virtual bool check(void) { return false; };         //! will be called several times during processing, should return true if all is ok.
        virtual uint16_t getNextStep( uint16_t s=JSTEP_NONE ) const;
        boost::posix_time::time_duration getElapsed( uint16_t from, uint16_t to ) { return (info.times[to] - info.times[from]); }
        
        virtual bool getWork(WorkInProgress::Ptr, uint16_t, const std::map<StepID,CountT>& ) { return false; };
        virtual void ungetWork(WorkInProgress::Ptr) { };
        virtual void failWork(WorkInProgress::Ptr) { };
        virtual void returnResults(WorkInProgress::Ptr) { };

        void setFailed(void);
        bool isOK(void);
        
        virtual void init(void) {};
        virtual void cleanup(void) {};
        virtual bool run(WorkInProgress::Ptr, uint16_t) = 0;
        
        std::string cfg(void);
        
        redux::logging::Logger& getLogger(void) { return logger; }
        void startLog(bool overwrite=false);
        void printJobInfo(void);
        void stopLog(void);
        
        std::unique_lock<std::mutex> getLock(bool trylock=false) {
            if(trylock) return std::unique_lock<std::mutex>(jobMutex,std::try_to_lock);
            return std::unique_lock<std::mutex>(jobMutex);
        }
        
        static void moveTo( Job* job, uint16_t to );
        static std::unique_lock<std::mutex> getGlobalLock(bool trylock=false) {
            if(trylock) return std::unique_lock<std::mutex>(globalMutex,std::try_to_lock);
            return std::unique_lock<std::mutex>(globalMutex);
        }
        
        virtual bool mayBeDeleted(void) { return true; }
        
        virtual size_t memUsage(void) { return 0; }       //!< Approximate current memory usage of this job
        virtual size_t diskUsage(void) { return 0; }      //!< Approximate current disk usage of this job
        virtual size_t procUsage(void) { return 0; }      //!< Approximate memory usage for processing 1 part
        virtual void updateStatus(void) {};
        
        void addThread( uint16_t n );
        void delThread( uint16_t n );
        void cleanupThreads( void );
        void threadLoop( void );
        
        bool operator<(const Job& rhs);
        bool operator!=(const Job& rhs);
    
    protected:
       
        mutable PackedData packed;
        std::mutex jobMutex;
        static std::mutex globalMutex;
        std::string cachePath;
        redux::logging::Logger logger;
        
        boost::asio::io_service ioService;
        std::shared_ptr<boost::asio::io_service::work> workLoop;
        boost::thread_group pool;
        
        friend class Daemon;
        friend class Worker;

    };

    template<class JobType, typename Dummy>
    class StaticJobTypeInit { };

    /*! @} */
        
}  // redux


#endif // REDUX_JOB_HPP
