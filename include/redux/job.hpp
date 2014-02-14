#ifndef REDUX_JOB_HPP
#define REDUX_JOB_HPP

#include <map>
#include <memory>
#include <string>

#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/archive/text_oarchive.hpp>

namespace po = boost::program_options;
namespace bpt = boost::property_tree;


namespace redux {
        
    /*! @ingroup redux
     *  @{
     */

    /*! Base class for a "job" to be processed by the redux framework.
     * 
     */
    class Job {
        
    public:
        
        enum State : uint8_t { JST_UNDEFINED = 0,       // NB: specific-type enums are supported by gcc/clang, but may not behave well with other compilers
                               JST_QUEUED,    
                               JST_IDLE,
                               JST_ACTIVE,
                               JST_PAUSED,
                               JST_COMPLETED,    
                               JST_ERR = 255
                             };
                             
        typedef std::shared_ptr<Job> JobPtr;
        typedef Job* ( *JobCreator )( void );
        typedef std::map<std::string, std::pair<size_t, JobCreator>> MapT;

        static MapT& getMap(void) { static MapT m = []( void ) {
            MapT tmp;
            std::cout << "NEW" << std::endl;
            tmp.clear();
            return tmp; }();
            return m;
        };
        static size_t registerJob( const std::string&, JobCreator f );
        static std::vector<JobPtr> parseTree( po::variables_map& vm, bpt::ptree& tree );
        static JobPtr newJob( const std::string& );
        
        struct Info {
            size_t id;
            uint8_t priority, verbosity;
            State state;
            std::string typeString, name, user, host;
            std::string logFile;
            boost::posix_time::ptime submitTime;
            Info();
            size_t size(void) const;
            char* pack(char*) const;
            const char* unpack(const char*, bool);
            static std::string printHeader(void);
            std::string print(void);
        } info;
        
        virtual void parseProperties( po::variables_map&, bpt::ptree& ) {};
        
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
        virtual bpt::ptree getPropertyTree( bpt::ptree* root=nullptr ) { return bpt::ptree(); };

        Job(void);
        virtual ~Job(void);

        virtual size_t size(void) const;
        virtual char* pack(char*) const;
        virtual const char* unpack(const char*, bool);
        
        const size_t& id(void) { return info.id; };
        void setID(size_t id) { info.id = id; };
         
    private:
        //static std::map<std::string, std::pair<size_t, JobCreator>> jobMap;
        static size_t nJobTypes;
    
    protected:
        
        template <typename Archive>
        void serialize( Archive& ar, const unsigned int version ) {
            std::string tmp = "Job";
            ar & tmp;
            // no members yet...
        }
        
        friend class boost::serialization::access;
    };


    /*! @} */
        
}  // redux


#endif // REDUX_JOB_HPP
