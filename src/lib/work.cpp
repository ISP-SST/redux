#include "redux/work.hpp"

#ifdef DEBUG_
#   define TRACE_THREADS
#endif

#include "redux/job.hpp"
#include "redux/logging/logger.hpp"
#include "redux/util/convert.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/trace.hpp"

using namespace redux::util;
using namespace redux;
using namespace std;

#ifdef DEBUG_
//#define DBG_PART_
//#define DBG_WIP_
#endif

namespace {

#ifdef DBG_PART_
    static atomic<int> partCounter(0);
#endif
#ifdef DBG_WIP_
    static atomic<int> wipCounter(0);
#endif
    
}
Part::Part() : id( 0 ), partStarted(boost::posix_time::not_a_date_time), step( 0 ), nRetries( 0 ), partType(0),
        nThreads(0), runtime_wall(0.0), runtime_cpu(0.0) {
#ifdef DBG_PART_
    LOG_DEBUG << "Constructing Part: (" << hexString(this) << ") new instance count = " << (partCounter.fetch_add(1)+1);
#endif
}

Part::~Part() {
#ifdef DBG_PART_
    LOG_DEBUG << "Destructing Part: (" << hexString(this) << ") new instance count = " << (partCounter.fetch_sub(1)-1);
#endif
}


uint64_t Part::pack( char* ptr ) const {

    using redux::util::pack;
    
    uint64_t count = pack( ptr, id );
    count += pack( ptr+count, step );
    count += pack( ptr+count, nRetries );
    count += pack( ptr+count, partType );
    count += pack( ptr+count, nThreads );
    count += pack( ptr+count, runtime_wall );
    count += pack( ptr+count, runtime_cpu );

    return count;

}


uint64_t Part::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;
    
    uint64_t count = unpack( ptr, id, swap_endian );
    count += unpack( ptr+count, step, swap_endian );
    count += unpack( ptr+count, nRetries, swap_endian );
    count += unpack( ptr+count, partType, swap_endian );
    count += unpack( ptr+count, nThreads, swap_endian );
    count += unpack( ptr+count, runtime_wall, swap_endian );
    count += unpack( ptr+count, runtime_cpu, swap_endian );
    
    return count;
}


void Part::load(void) {
    
    CacheItem::cacheLoad();
    
}


void Part::unload(void) {
    
    CacheItem::cacheClear();
    
}


WorkInProgress::WorkInProgress(void) : job(), jobID(0),
    workStarted(boost::posix_time::not_a_date_time), isRemote(false), hasResults(false), nParts(0), nCompleted(0) {
#ifdef DBG_WIP_
    LOG_DEBUG << "Constructing WIP: (" << hexString(this) << ") new instance count = " << (wipCounter.fetch_add(1)+1);
#endif
}

WorkInProgress::WorkInProgress(const WorkInProgress& rhs) :  job(), jobID(rhs.jobID), parts(rhs.parts), workStarted(rhs.workStarted),
    isRemote(rhs.isRemote), hasResults(false), nParts(rhs.nParts), nCompleted(rhs.nCompleted) {
#ifdef DBG_WIP_
    LOG_DEBUG << "Constructing WIP: (" << hexString(this) << ") new instance count = " << (wipCounter.fetch_add(1)+1);
#endif
}

WorkInProgress::~WorkInProgress() {
#ifdef DBG_WIP_
    LOG_DEBUG << "Destructing WIP: (" << hexString(this) << ") new instance count = " << (wipCounter.fetch_sub(1)-1);
#endif
}


uint64_t WorkInProgress::size(void) const {
    static uint64_t sz = 2*sizeof(uint16_t)+sizeof(uint32_t) // nParts + nCompleted + jobID
                       + sizeof(time_t);                     // workStarted is converted and transferred as time_t
    return sz;
}


uint64_t WorkInProgress::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack( ptr, jobID );
    count += pack( ptr+count, nParts );
    count += pack( ptr+count, nCompleted );
    count += pack(ptr+count,redux::util::to_time_t( workStarted ));
    return count;
}


uint64_t WorkInProgress::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = unpack( ptr, jobID, swap_endian );
    count += unpack( ptr+count, nParts, swap_endian );
    count += unpack( ptr+count, nCompleted, swap_endian );
    time_t timestamp;
    count += unpack(ptr+count,timestamp,swap_endian);
    workStarted = boost::posix_time::from_time_t( timestamp );
    return count;
}


void WorkInProgress::reset( void ) {
    
    job.reset();
    jobID = 0;
    resetParts();
    workStarted = boost::posix_time::not_a_date_time;
    
}


void WorkInProgress::resetParts( void ) {
    
    parts.clear();
    hasResults = false;
    nParts = nCompleted = 0;
    
}


uint64_t WorkInProgress::workSize(void) {
    uint64_t sz = this->size() + 1; // + newJob
    Job::JobPtr thisJob = job.lock();
    if( thisJob && (jobID != thisJob->info.id)) {
        sz += thisJob->size();
    }
    nParts = 0;
    for( const auto& part: parts ) {
        if( part ) {
            nParts++;
            sz += part->size();
        }
    }
    return sz;
}


uint64_t WorkInProgress::packWork( char* ptr ) {
    
    Job::JobPtr thisJob = job.lock();
    if( !thisJob ) {
        throw invalid_argument( "Can't pack WIP without a job instance..." );
    }
    THREAD_MARK
    using redux::util::pack;
    uint64_t count = this->size();                                  // skip WIP info, it's packed below.
    bool newJob = (jobID != thisJob->info.id);
    count += pack( ptr+count, newJob );
    THREAD_MARK
    if( newJob ) {
        count += thisJob->pack( ptr+count );                            // pack Job-info
    }
    count += thisJob->packParts( ptr+count, shared_from_this() );       // pack parts
    this->pack( ptr );                                              // write WIP info, (done last, in case nParts changed).
    THREAD_MARK
    return count;
}


uint64_t WorkInProgress::unpackWork( const char* ptr, std::shared_ptr<Job>& tmpJob, bool swap_endian ) {

    using redux::util::unpack;
    uint64_t count = this->unpack( ptr, swap_endian );
    bool newJob(false);
    count += unpack( ptr+count, newJob );
    if( newJob ) {
        string tmpS = string( ptr+count );
        tmpJob = Job::newJob( tmpS );
        job = tmpJob;
        if( tmpJob ) {
            count += tmpJob->unpack( ptr+count, swap_endian );
        } else throw invalid_argument( "Unrecognized Job tag: \"" + tmpS + "\"" );
    } else {
        tmpJob = job.lock();
    }
    
    if( tmpJob ) {
        count += tmpJob->unpackParts( ptr+count, shared_from_this(), swap_endian );
    } else throw invalid_argument( "Can't unpack parts without a job instance..." );
    
    return count;
}


void WorkInProgress::returnResults(void) {
    Job::JobPtr thisJob = job.lock();
    if( thisJob ) {
        thisJob->returnResults( shared_from_this() );
        resetParts();
    }
}


bool WorkInProgress::operator<( const WorkInProgress& rhs ) const {
    Job::JobPtr thisJob = job.lock();
    Job::JobPtr rhsJob = rhs.job.lock();
    if( thisJob != rhsJob ) return thisJob < rhsJob;
    size_t sz = std::min( parts.size(), rhs.parts.size() );
    for( size_t i(0); i<sz; ++i ) {
        if( parts[i] != rhs.parts[i] ) return thisJob < rhsJob;
    }
    //return (parts.size() < rhs.parts.size());
    return false;
    
}


std::string WorkInProgress::print( void ) {
    Job::JobPtr thisJob = job.lock();
    string ret = "\"" + ( thisJob ? thisJob->info.name : string( "undefined" ) ) + "\"";
    if( parts.size() ) {
        ret += " (job: #" + (thisJob ? to_string(thisJob->info.id) : "0" ) + " part(s):";
        for( auto & part : parts ) {
            if( part ) {
                ret += " " + to_string( part->id );
            }
        }
        ret += ")";
    }
    
    return ret;
    
}
