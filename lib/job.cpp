#include "redux/job.hpp"

#include "redux/logger.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stringutil.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/date_time/posix_time/time_formatters.hpp>

using namespace redux::util;
using namespace redux;
using namespace std;

size_t Job::nJobTypes( Job::getMap().size() );

#define lg Logger::lg
namespace {

    const string thisChannel = "job";

    const std::string StateNames[7] = { "undefined", "queued", "idle", "active", "paused", "completed", "error" };
    const std::string StateTags[7] = { "-", "Q", "I", "A", "P", "C", "E" };

}

size_t Job::registerJob( const string& name, JobCreator f ) {
    auto ret = getMap().insert( { boost::to_upper_copy( name ), {nJobTypes + 1, f}} );
    if( ret.second ) {
        return ret.first->second.first;
    }
    return nJobTypes++;
}

vector<Job::JobPtr> Job::parseTree( po::variables_map& vm, bpt::ptree& tree ) {
    vector<JobPtr> tmp;
    for( auto & it : tree ) {
        string nm = it.first;
        auto it2 = getMap().find( boost::to_upper_copy( nm ) );       // check if the current tag matches a registered (Job-derived) class.
        if( it2 != getMap().end() ) {
            Job* tmpJob = it2->second.second();
            tmpJob->parseProperties( vm, it.second );
            tmp.push_back( shared_ptr<Job>( tmpJob ) );
        }
    }
    return tmp;
}

Job::JobPtr Job::newJob( const string& name ) {
    JobPtr tmp;
    auto it = getMap().find( boost::to_upper_copy( name ) );
    if( it != getMap().end() ) {
        tmp.reset( it->second.second() );
    } else LOG_WARN << "No job with tag \"" << name << "\" registered.";

    return tmp;
}


Job::Info::Info(void) : id(0), priority(0), verbosity(0), state(JST_UNDEFINED) { }

size_t Job::Info::size(void) const {
    size_t sz = sizeof(size_t) + sizeof(time_t) + sizeof(State) + 2;
    sz += typeString.length() + name.length() + user.length() + host.length() + logFile.length() + 5 ;
    return sz;
}

char* Job::Info::pack(char* ptr) const {

    using redux::util::pack;
    
    ptr = pack(ptr,typeString);         // NB: the type-name has to be first in the packed block,
    ptr = pack(ptr,name);               //   it is used to identify which job-class to instantiate on the receiving side.
    ptr = pack(ptr,user);  
    ptr = pack(ptr,host);  
    ptr = pack(ptr,logFile);  
    ptr = pack(ptr,id);  
    ptr = pack(ptr,priority);  
    ptr = pack(ptr,verbosity);  
    ptr = pack(ptr,state);  
    ptr = pack(ptr,to_time_t( submitTime ));  

    return ptr;
   
}

const char* Job::Info::unpack(const char* ptr, bool doSwap) {
    
    using redux::util::unpack;
    
    ptr = unpack(ptr,typeString);
    ptr = unpack(ptr,name);  
    ptr = unpack(ptr,user);  
    ptr = unpack(ptr,host);  
    ptr = unpack(ptr,logFile);  
    ptr = unpack(ptr,id);
    ptr = unpack(ptr,priority);
    ptr = unpack(ptr,verbosity);
    ptr = unpack(ptr,state);
    time_t timestamp;
    ptr = unpack(ptr,timestamp);
    if( doSwap ) {
        swapEndian(id);
        swapEndian(timestamp);
    }
    submitTime = boost::posix_time::from_time_t( timestamp );
    
    return ptr;
}

std::string Job::Info::printHeader(void) {
    string hdr = alignRight("ID",5) + alignCenter("type",10) + alignCenter("submitted",20);
    hdr += alignCenter("name",15) + alignLeft("user",15) + alignCenter("priority",8) + alignCenter("state",8);
    return hdr;
}

std::string Job::Info::print(void) {
    string info = alignRight(std::to_string(id),5) + alignCenter(typeString,10);
    info += alignCenter(to_iso_extended_string( submitTime ),20);
    info += alignCenter(name,15) + alignLeft(user+"@"+host,15) + alignCenter(std::to_string(priority),8) + alignLeft(StateTags[state],8);
    return info;
}

Job::Job( void ) {
    info.user = getUname();
    info.host = boost::asio::ip::host_name();
}

Job::~Job( void ) {

}

size_t Job::size(void) const {
    size_t sz = info.size();
    return sz;
}

char* Job::pack(char* ptr) const {
    ptr = info.pack(ptr);
    return ptr;
}

const char* Job::unpack(const char* ptr, bool doSwap) {
    ptr = info.unpack(ptr,doSwap);
    return ptr;
}
