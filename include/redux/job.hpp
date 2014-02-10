#ifndef REDUX_JOB_HPP
#define REDUX_JOB_HPP

#include <map>
#include <memory>
#include <string>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/archive/text_oarchive.hpp>

namespace po = boost::program_options;
namespace bpt = boost::property_tree;


namespace redux {

    class Job {
        
    public:
        
        typedef std::shared_ptr<Job> JobPtr;
        typedef Job* ( *JobCreator )( void );
       
        static size_t registerJob( const std::string, JobCreator f );

        static std::vector<JobPtr> parseTree(po::variables_map& vm, bpt::ptree& tree);
        virtual void parseProperties( po::variables_map&, bpt::ptree& ) {};
        
        /*! @brief Returns a boost::property_tree containing the settings for this job.
         *  @details If 
         * 
         */
        virtual bpt::ptree getPropertyTree( bpt::ptree* root=nullptr ) { return bpt::ptree(); };

        Job(void);
        virtual ~Job(void);

    private:
        static std::map<std::string, std::pair<size_t, JobCreator>> jobMap;
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


}  // redux


#endif // REDUX_JOB_HPP
