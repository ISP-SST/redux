#include "redux/job.hpp"

#include "redux/util/cache.hpp"
#include "redux/util/convert.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/version.hpp"

#include <mutex>

#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/date_time/posix_time/time_formatters.hpp>
#include <boost/property_tree/info_parser.hpp>

using namespace redux::logging;
using namespace redux::util;
using namespace redux;
using namespace std;

std::mutex Job::globalMutex;
std::map<Job::StepID, Job::CountT> Job::counts = { { {0,JSTEP_SUBMIT}, {1,5} },
                                                   { {0,JSTEP_NONE}, {1,5} },
                                                   { {0,JSTEP_COMPLETED}, {} },
                                                   { {0,JSTEP_ERR}, {} } };

#ifdef DEBUG_
//#define DBG_JOB_
#endif


namespace {

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


size_t Job::registerJob( const string& name, JobCreator f ) {
    
    static size_t nJobTypes(0);
    std::unique_lock<mutex> lock(globalJobMutex);
    auto ret = getMap().insert({ boost::to_upper_copy(name), {nJobTypes + 1, f}});
    if( !ret.second ) {
        return ret.first->second.first;
    }
    return nJobTypes++;
}


vector<Job::JobPtr> Job::parseTree(bpo::variables_map& vm, bpt::ptree& tree, redux::logging::Logger& logger, bool check) {

    vector<JobPtr> tmp;
    std::unique_lock<mutex> lock(globalJobMutex);
    for(auto & property : tree) {
        string name = property.first;
        auto it2 = getMap().find(boost::to_upper_copy(name));       // check if the current tag matches a registered (Job-derived) class.
        if(it2 != getMap().end()) {
            Job* tmpJob = it2->second.second();
            tmpJob->parsePropertyTree( vm, property.second, logger );
            tmpJob->getLogger().addLogger( logger );
            bool checkResult(false);
            if(!check || (checkResult=tmpJob->check()) ) {
                if(checkResult) tmpJob->info.flags |= Job::CHECKED;
                tmp.push_back(shared_ptr<Job>(tmpJob));
            } else LOG_ERR << "Job \"" << tmpJob->info.name << "\" of type " << tmpJob->info.typeString << " failed cfgCheck, skipping." << ende;
        } else {
#ifdef DEBUG_
            LOG_ERR << "No job-type with tag \"" << name << "\" registered." << ende;
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
        cerr << "No job with tag \"" << name << "\" registered." << endl;
#endif
    }
    return tmp;
}


string Job::stateString(uint8_t state) {

    switch(state) {
        case JSTATE_NONE: return "";
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
        case JSTATE_NONE: return "";
        case JSTATE_IDLE: return "I";
        case JSTATE_ACTIVE: return "A";
        case JSTATE_PAUSED: return "P";
        case JSTATE_CANCELLED: return "C";
        case JSTATE_ERR: return "E";
        default: return "";
    }

}


Job::Info::Info(void) : id(0), timeout(36000), maxProcessingTime(0), priority(10), verbosity(0), maxPartRetries(5),
         maxThreads(255), flags(0), step(JSTEP_NONE), state(JSTATE_NONE) {

}


Job::Info::Info(const Info& rhs) {

    id = rhs.id;
    timeout = rhs.timeout;
    maxProcessingTime = rhs.maxProcessingTime;
    priority = rhs.priority;
    verbosity = rhs.verbosity;
    maxPartRetries = rhs.maxPartRetries;
    maxThreads = rhs.maxThreads;
    step = rhs.step.load();
    state = rhs.state.load();
    typeString = rhs.typeString;
    name = rhs.name;
    user = rhs.user;
    host = rhs.host;
    logFile = rhs.logFile;
    outputDir = rhs.outputDir;
    progressString = "";
    submitTime = rhs.submitTime;
    startedTime = rhs.startedTime;
    completedTime = rhs.completedTime;

}

uint64_t Job::Info::size(void) const {
    
    static uint64_t fixed_sz = 2*sizeof(uint32_t)       // id, timeout, maxProcessingTime
                             + 3*sizeof(uint16_t)       // maxThreads, step, flags
                             + 4                        // priority, verbosity, maxPartRetries, state
                             + 3*sizeof(time_t)         // submitTime, startedTime, completedTime
                             + 20 + 6;                  // progressString + \0 for typeString/name/user/host/logFile/outputDir
                                                        // fixed size for progressString, so it doesn't matter if it's updated between size() & pack()

    uint64_t sz = fixed_sz + typeString.length() + name.length() + user.length() + host.length();
    sz += logFile.length() + outputDir.length();
    
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
    count += pack(ptr+count, flags);
    count += pack(ptr+count, maxPartRetries);
    count += pack(ptr+count, step.load());
    count += pack(ptr+count, state.load());
    count += pack(ptr+count, name);
    count += pack(ptr+count, user);
    count += pack(ptr+count, host);
    count += pack(ptr+count, progressString);
    count += pack(ptr+count, logFile);
    count += pack(ptr+count, outputDir);
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
    count += unpack(ptr+count, flags);
    count += unpack(ptr+count, maxPartRetries);
    uint16_t tmp16(0);
    count += unpack(ptr+count, tmp16);
    step.store(tmp16);
    uint8_t tmp8(0);
    count += unpack(ptr+count, tmp8);
    state.store(tmp8);
    count += unpack(ptr+count, name);
    count += unpack(ptr+count, user);
    count += unpack(ptr+count, host);
    count += unpack(ptr+count, progressString);
    count += unpack(ptr+count, logFile);
    count += unpack(ptr+count, outputDir);
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
    string hdr = alignRight("#", 4) + alignRight("ID", 5) + alignCenter("type", 10) + alignCenter("submitted", 20); // + alignCenter("started", 20);
    hdr += alignCenter("name", 15) + alignLeft("user", 15) + alignCenter("priority", 8) + alignCenter("state", 8);
    return hdr;
}


std::string Job::Info::print(void) {
    string info = alignRight(std::to_string(id), 5) + alignCenter(typeString, 10);
    string startedString = "";
    if ( !startedTime.is_not_a_date_time() ) startedString = to_iso_extended_string(startedTime);
    info += alignCenter(to_iso_extended_string(submitTime), 20); // + alignCenter(startedString, 20);
    info += alignCenter(name, 15) + alignLeft(user + "@" + host, 15) + alignCenter(std::to_string(priority), 8)
    + alignCenter(stateTag(state), 3) + alignLeft(progressString, 15);
    return info;
}


void Job::parsePropertyTree(bpo::variables_map& vm, bpt::ptree& tree, redux::logging::Logger& logger) {
    
    Info defaults = globalDefaults;

    if( vm.count( "verbosity" ) > 0 ) {         // if --verbosity N is specified, use it.
        defaults.verbosity = vm["verbosity"].as<int>();
    } else defaults.verbosity = Logger::getDefaultLevel();
    
    if( vm.count( "name" ) > 0 ) {
        defaults.name = vm["name"].as<string>();
    }
    if( vm.count( "log-file" ) ) {
        vector<string> logfiles = vm["log-file"].as<vector<string>>();
        for( auto & filename : logfiles ) { // get first non-empty logfile
            if( filename.empty() )  continue;
            defaults.logFile = filename;
            break;
        }
    } else if( !defaults.name.empty() ) {
        defaults.logFile = defaults.name + ".log";
    }
    
    info.timeout = tree.get<uint32_t>("TIMEOUT", defaults.timeout);
    info.priority = tree.get<uint8_t>("PRIORITY", defaults.priority);
    info.verbosity = tree.get<uint8_t>("VERBOSITY", defaults.verbosity);
    info.maxThreads = tree.get<uint16_t>("MAX_THREADS", defaults.maxThreads);
    info.maxPartRetries = tree.get<uint8_t>("MAX_PART_RETRIES", defaults.maxPartRetries);
    info.logFile = tree.get<string>("LOGFILE", defaults.logFile);
    info.outputDir = tree.get<string>( "OUTPUT_DIR", defaults.outputDir );
    if( info.outputDir.empty() ) info.outputDir = bfs::current_path().string();
    info.name = tree.get<string>("NAME", defaults.name);

}


bpt::ptree Job::getPropertyTree(bpt::ptree* root) {
    bpt::ptree tree;
    if(info.timeout != globalDefaults.timeout) tree.put("TIMEOUT", info.timeout);
    if(info.priority != globalDefaults.priority) tree.put("PRIORITY", info.priority);
    if(info.verbosity != globalDefaults.verbosity) tree.put("VERBOSITY", info.verbosity);
    if(info.maxThreads != globalDefaults.maxThreads) tree.put("MAX_THREADS", info.maxThreads);
    if(info.maxPartRetries != globalDefaults.maxPartRetries) tree.put("MAX_PART_RETRIES", info.maxPartRetries);
    if(info.logFile != globalDefaults.logFile) tree.put("LOGFILE", info.logFile);
    if(info.outputDir != globalDefaults.outputDir) tree.put("OUTPUT_DIR", info.outputDir);
    if(info.name != globalDefaults.name) tree.put("NAME", info.name);
    if(root) {
        root->push_back(bpt::ptree::value_type("job", tree));
    }
    return tree;
}


Job::Job(void) : cachePath("") {
    info.user = getUname();
    info.host = boost::asio::ip::host_name();
#ifdef DBG_JOB_
    LOG_DEBUG << "Constructing Job: (" << hexString(this) << ") new instance count = " << (jobCounter.fetch_add(1)+1) << ende;
#endif

}


Job::~Job(void) {
#ifdef DBG_JOB_
    LOG_DEBUG << "Destructing Job#" << info.id << ": (" << hexString(this) << ") new instance count = " << (jobCounter.fetch_sub(1)-1) << ende;
#endif

    if( !cachePath.empty() ) {
        try {
            bfs::remove_all( bfs::path(Cache::get().path()) / bfs::path(cachePath) );
        } catch( const exception& e) {
            cerr << "Failed to remove path: " << bfs::path(Cache::get().path()) / bfs::path(cachePath) << endl
                 << "  reason: " << e.what() << endl;
        }
    }
    
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


uint16_t Job::getNextStep( uint16_t step ) const {

    switch( step ) {
        //
        default: return JSTEP_ERR;
    }
    
}

        
void Job::setFailed(void) {
    
    unique_lock<mutex>(jobMutex);
    info.state |= JSTATE_ERR;
    
}


bool Job::isOK(void) {
    
    unique_lock<mutex>(jobMutex);
    return !(info.state&JSTATE_ERR);
    
}

string Job::cfg(void) {

    bpt::ptree pt;
    this->getPropertyTree( &pt );
    stringstream ss;
    bpt::write_info( ss, pt );
    return ss.str();
    
}


void Job::startLog( bool overwrite ) {
    
    bfs::path logFilePath = bfs::path( info.logFile );

    string tmpChan = "job "+to_string( info.id );
    if( isRelative( logFilePath ) ) {
        bfs::path outDir( info.outputDir );
        if( !outDir.empty() && !bfs::exists(outDir) ) {
            if( !bfs::create_directories(outDir) ) {
                throw job_error( info.name + ": Failed to create directory for output: " + outDir.string() );
            }
        }
        logFilePath = bfs::path(info.outputDir) / logFilePath;
    }

    logger.setLevel( info.verbosity );
    logger.setContext( "job "+to_string(info.id) );
    try {
        logger.addFile( logFilePath.string(), 0, overwrite );
    } catch( exception& e ) {
        throw job_error( info.name + ": " + string(e.what()) );
    }
    

}


void Job::printJobInfo(void) {
    
    LOG << "Redux version = " << getLongVersionString();
    LOG << "\nJob configuration:\n" << cfg() << ende;
    
}


void Job::stopLog(void) {
    logger.flushAll();
    logger.removeAllOutputs();
    //jlog.reset();
}


void Job::moveTo( Job* job, uint16_t to ) {
    if( !job ) return;
    auto glock = getGlobalLock();
    CountT& c_old = counts[StepID(job->getTypeID(),job->info.step)];
    CountT& c_new = counts[StepID(job->getTypeID(),to)];
    c_old.active--;
    c_new.active++;
    glock.unlock();
    job->info.step = to;
}


bool Job::operator<(const Job& rhs) {
    return (info.id < rhs.info.id);
}


bool Job::operator!=(const Job& rhs) {
    return (info.id != rhs.info.id);
}
