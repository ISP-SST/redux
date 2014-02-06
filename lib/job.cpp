#include "redux/job.hpp"

#include "redux/logger.hpp"
#include "redux/util/stringutil.hpp"

using namespace redux;
using namespace std;

map<string, pair<size_t, Job::JobCreator>> Job::jobMap;
size_t Job::nJobTypes( Job::jobMap.size() );

#define lg Logger::lg
namespace {
    const string thisChannel = "job";
}

size_t Job::registerJob( const string name, JobCreator f ) {
    // BUG in GCC (still in 4.8.1), inserting into a map immediately after initialization causes a segfault.
    // A workaround is to call clear() before inserting.
    static bool dummy = []( void ) { jobMap.clear(); return false; }();
    auto ret = jobMap.insert( { redux::util::uppercase( name ), {nJobTypes + 1, f}} );
    if( ret.second == dummy ) {
        return ret.first->second.first;
    }
    return nJobTypes++;
}

vector<Job::JobPtr> Job::parseTree( po::variables_map& vm, bpt::ptree& tree ) {
    vector<JobPtr> tmp;
    for( auto & it : tree ) {
        string nm = it.first;
        auto it2 = jobMap.find( redux::util::uppercase( nm ) );
        if( it2 != jobMap.end() ) {
            LOG_DEBUG << "Parsing configuration \"" << nm << "\".";
            Job* tmpJob = it2->second.second();
            tmpJob->parseProperties( vm, it.second );
            tmp.push_back( shared_ptr<Job>( tmpJob ) );
        }
    }
    return tmp;
}



Job::Job( void ) {

}

Job::~Job( void ) {

}

