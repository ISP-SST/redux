#include "redux/job.hpp"

#include "redux/logger.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stringutil.hpp"

#include <mutex>

#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/date_time/posix_time/time_formatters.hpp>

using namespace redux::util;
using namespace redux;
using namespace std;


#ifdef DEBUG_
//#define DBG_JOB_
#endif

#define lg Logger::lg
namespace {

    const string thisChannel = "job";

    const std::string StateNames[10] = { "undefined", "preprocessed", "queued", "done", "completed", "postprocess", "idle", "active", "paused", "error" };
    const std::string StateTags[10] = { "-", "Pre", "Q", "D", "C", "Post", "I", "A", "P", "E" };

    mutex globalJobMutex;
    /*const*/ Job::Info globalDefaults;
    
#ifdef DBG_JOB_
    static atomic<int> jobCounter(0);
#endif

}


void redux::runThreadsAndWait(boost::asio::io_service& service, uint16_t nThreads) {
    boost::thread_group pool;
    for(uint16_t t = 0; t < nThreads; ++t) {
        pool.create_thread(boost::bind(&boost::asio::io_service::run, &service));
    }
    pool.join_all();
    service.reset();        // reset service so that next call will not fail
}


size_t Job::registerJob(const string& name, JobCreator f) {
    static size_t nJobTypes(0);
    std::unique_lock<mutex> lock(globalJobMutex);
    auto ret = getMap().insert({ boost::to_upper_copy(name), {nJobTypes + 1, f}});
    if(ret.second) {
        return ret.first->second.first;
    }
    return nJobTypes++;
}


vector<Job::JobPtr> Job::parseTree(bpo::variables_map& vm, bpt::ptree& tree, bool check) {
    vector<JobPtr> tmp;
    std::unique_lock<mutex> lock(globalJobMutex);
    for(auto & property : tree) {
        string name = property.first;
        auto it2 = getMap().find(boost::to_upper_copy(name));       // check if the current tag matches a registered (Job-derived) class.
        if(it2 != getMap().end()) {
            Job* tmpJob = it2->second.second();
            tmpJob->parsePropertyTree(vm, property.second);
            if(!check || tmpJob->check()) {
                tmp.push_back(shared_ptr<Job>(tmpJob));
            } else LOG_WARN << "Job \"" << tmpJob->info.name << "\" of type " << tmpJob->info.typeString << " failed cfgCheck, skipping.";
        } else {
#ifdef DEBUG_
            LOG_WARN << "No job-class with tag \"" << name << "\" registered.";
#endif
        }
    }
    return tmp;
}


Job::JobPtr Job::newJob(const string& name) {
    JobPtr tmp;
    std::unique_lock<mutex> lock(globalJobMutex);
    auto it = getMap().find(boost::to_upper_copy(name));
    if(it != getMap().end()) {
        tmp.reset(it->second.second());
    } else {
#ifdef DEBUG__
        LOG_WARN << "No job with tag \"" << name << "\" registered.";
#endif
    }
    return tmp;
}


string Job::stateString(uint8_t state) {

    switch(state) {
        case JSTATE_IDLE: return "idle";
        case JSTATE_ACTIVE: return "active";
        case JSTATE_PAUSED: return "paused";
        case JSTATE_CANCELLED: return "cancelled";
        case JSTATE_ERR: return "error";
        default: return "-";
    }

}


string Job::stateTag(uint8_t state) {

    switch(state) {
        case JSTATE_IDLE: return "I";
        case JSTATE_ACTIVE: return "A";
        case JSTATE_PAUSED: return "P";
        case JSTATE_CANCELLED: return "C";
        case JSTATE_ERR: return "E";
        default: return "";
    }

}


Job::Info::Info(void) : id(0), timeout(1000), priority(10), verbosity(0), maxPartRetries(1), maxThreads(255), 
         step(0), state(JSTATE_IDLE) {

}


uint64_t Job::Info::size(void) const {
    uint64_t sz = 2*sizeof(uint32_t) + sizeof(uint16_t) + 5;
    sz += typeString.length() + name.length() + user.length() + host.length() + logFile.length() + 5;
    sz += 3*sizeof(time_t);
    return sz;
}


uint64_t Job::Info::pack(char* ptr) const {
    using redux::util::pack;
    uint64_t count = pack(ptr, typeString);    // NB: the type-string has to be first in the packed block,
    count += pack(ptr+count, id);              //   it is used to identify which job-class to instantiate on the receiving side.
    count += pack(ptr+count, timeout);
    count += pack(ptr+count, priority);
    count += pack(ptr+count, verbosity);
    count += pack(ptr+count, maxThreads);
    count += pack(ptr+count, maxPartRetries);
    count += pack(ptr+count, step.load());
    count += pack(ptr+count, state.load());
    count += pack(ptr+count, name);
    count += pack(ptr+count, user);
    count += pack(ptr+count, host);
    count += pack(ptr+count, logFile);
    count += pack(ptr+count, redux::util::to_time_t(submitTime));
    count += pack(ptr+count, redux::util::to_time_t(startedTime));
    count += pack(ptr+count, redux::util::to_time_t(completedTime));
    return count;
}


uint64_t Job::Info::unpack(const char* ptr, bool swap_endian) {
    using redux::util::unpack;
    uint64_t count = unpack(ptr, typeString);
    count += unpack(ptr+count, id, swap_endian);
    count += unpack(ptr+count, timeout);
    count += unpack(ptr+count, priority);
    count += unpack(ptr+count, verbosity);
    count += unpack(ptr+count, maxThreads);
    count += unpack(ptr+count, maxPartRetries);
    uint8_t tmp;
    count += unpack(ptr+count, tmp);
    step.store(tmp);
    count += unpack(ptr+count, tmp);
    state.store(tmp);
    count += unpack(ptr+count, name);
    count += unpack(ptr+count, user);
    count += unpack(ptr+count, host);
    count += unpack(ptr+count, logFile);
    time_t timestamp;
    count += unpack(ptr+count, timestamp, swap_endian);
    submitTime = boost::posix_time::from_time_t(timestamp);
    count += unpack(ptr+count, timestamp, swap_endian);
    startedTime = boost::posix_time::from_time_t(timestamp);
    count += unpack(ptr+count, timestamp, swap_endian);
    completedTime = boost::posix_time::from_time_t(timestamp);
    return count;
}


std::string Job::Info::printHeader(void) {
    string hdr = alignRight("ID", 5) + alignCenter("type", 10) + alignCenter("submitted", 20);
    hdr += alignCenter("name", 15) + alignLeft("user", 15) + alignCenter("priority", 8) + alignCenter("state", 8);
    return hdr;
}


std::string Job::Info::print(void) {
    string info = alignRight(std::to_string(id), 5) + alignCenter(typeString, 10);
    info += alignCenter(to_iso_extended_string(submitTime), 20);
    info += alignCenter(name, 15) + alignLeft(user + "@" + host, 15) + alignCenter(std::to_string(priority), 8)
    + alignCenter(stateTag(state), 3); // + alignLeft(stepString(step), 15);
    return info;
}


void Job::parsePropertyTree(bpo::variables_map&, bpt::ptree& tree) {
    
    info.timeout = tree.get<uint32_t>("TIMEOUT", globalDefaults.timeout);
    info.priority = tree.get<uint8_t>("PRIORITY", globalDefaults.priority);
    info.verbosity = tree.get<uint8_t>("VERBOSITY", globalDefaults.verbosity);
    info.maxThreads = tree.get<uint16_t>("MAX_THREADS", globalDefaults.maxThreads);
    info.maxPartRetries = tree.get<uint8_t>("MAX_PART_RETRIES", globalDefaults.maxPartRetries);
    info.logFile = tree.get<string>("LOGFILE", globalDefaults.logFile);
    
}


bpt::ptree Job::getPropertyTree(bpt::ptree* root) {
    bpt::ptree tree;
    if(info.timeout != globalDefaults.timeout) tree.put("TIMEOUT", info.timeout);
    if(info.priority != globalDefaults.priority) tree.put("PRIORITY", info.priority);
    if(info.verbosity != globalDefaults.verbosity) tree.put("VERBOSITY", info.verbosity);
    if(info.maxThreads != globalDefaults.maxThreads) tree.put("MAX_THREADS", info.maxThreads);
    if(info.maxPartRetries != globalDefaults.maxPartRetries) tree.put("MAX_PART_RETRIES", info.maxPartRetries);
    if(info.logFile != globalDefaults.logFile) tree.put("LOGFILE", info.logFile);
    if(root) {
        root->push_back(bpt::ptree::value_type("job", tree));
    }
    return tree;
}


Job::Job(void) {
    info.user = getUname();
    info.host = boost::asio::ip::host_name();
#ifdef DBG_JOB_
    LOG_DEBUG << "Constructing Job: (" << hexString(this) << ") new instance count = " << (jobCounter.fetch_add(1)+1);
#endif
}


Job::~Job(void) {
#ifdef DBG_JOB_
    LOG_DEBUG << "Destructing Job: (" << hexString(this) << ") new instance count = " << (jobCounter.fetch_sub(1)-1);
#endif
}


uint64_t Job::size(void) const {
    return info.size();
}


uint64_t Job::pack(char* ptr) const {
    return info.pack(ptr);
}


uint64_t Job::unpack(const char* ptr, bool swap_endian) {
    return info.unpack(ptr, swap_endian);
}


bool Job::operator<(const Job& rhs) {
    return (info.id < rhs.info.id);
}


bool Job::operator!=(const Job& rhs) {
    return (info.id != rhs.info.id);
}
