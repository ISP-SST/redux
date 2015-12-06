#include "redux/work.hpp"

#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/job.hpp"
#include "redux/logger.hpp"

using namespace redux::util;
using namespace redux;
using namespace std;

#ifdef DEBUG_
//#define DBG_PART_
//#define DBG_WIP_
#define lg Logger::mlg
#endif

namespace {

    const string thisChannel = "work";
    
#ifdef DBG_PART_
    static atomic<int> partCounter(0);
#endif
#ifdef DBG_WIP_
    static atomic<int> wipCounter(0);
#endif
    
}
Part::Part() : id( 0 ), step( 0 ), nRetries( 0 ), partType(0) {
#ifdef DBG_PART_
    LOG_DEBUG << "Constructing Part: (" << hexString(this) << ") new instance count = " << (partCounter.fetch_add(1)+1);
#endif
}

Part::~Part() {
#ifdef DBG_PART_
    LOG_DEBUG << "Destructing Part: (" << hexString(this) << ") new instance count = " << (partCounter.fetch_sub(1)-1);
#endif
}

uint64_t Part::size( void ) const {
    return sizeof( id ) + sizeof( step ) + sizeof( nRetries )  + sizeof( partType ) ;
}


uint64_t Part::pack( char* ptr ) const {

    using redux::util::pack;
    
    uint64_t count = pack( ptr, id );
    count += pack( ptr+count, step );
    count += pack( ptr+count, nRetries );
    count += pack( ptr+count, partType );

    return count;

}


uint64_t Part::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;
    
    uint64_t count = unpack( ptr, id, swap_endian );
    count += unpack( ptr+count, step, swap_endian );
    count += unpack( ptr+count, nRetries, swap_endian );
    count += unpack( ptr+count, partType, swap_endian );
    
    return count;
}


WorkInProgress::WorkInProgress( network::TcpConnection::Ptr c ) : isRemote( false ), hasResults(false), nParts(0), nCompleted(0) {
#ifdef DBG_WIP_
    LOG_DEBUG << "Constructing WIP: (" << hexString(this) << ") new instance count = " << (wipCounter.fetch_add(1)+1);
#endif
}

WorkInProgress::WorkInProgress(const WorkInProgress& rhs) :  job(rhs.job), parts(rhs.parts), isRemote(rhs.isRemote),
    hasResults(false), nParts(rhs.nParts), nCompleted(rhs.nCompleted) {
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
    static uint64_t sz = 2*sizeof(uint16_t); // nParts + nCompleted
    return sz;
}


uint64_t WorkInProgress::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack( ptr, nParts );
    count += pack( ptr+count, nCompleted );
    return count;
}


uint64_t WorkInProgress::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = unpack( ptr, nParts, swap_endian );
    count += unpack( ptr+count, nCompleted, swap_endian );
    return count;
}


uint64_t WorkInProgress::workSize(void) {
    uint64_t sz = this->size() + 1; // + newJob
    if(job && job != previousJob) {
        sz += job->size();
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


uint64_t WorkInProgress::packWork( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = this->pack( ptr );
    bool newJob = (job && job != previousJob);
    count += pack( ptr+count, newJob );
    if(newJob) {
        count += job->pack(ptr+count);
    }
    for( auto& part: parts ) {
        if( part ) {
            count += part->pack(ptr+count);
        }
    }
    return count;
}


uint64_t WorkInProgress::unpackWork( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = this->unpack( ptr, swap_endian );
    bool newJob;
    count += unpack( ptr+count, newJob );
    if( newJob ) {
        swap(job,previousJob);
        string tmpS = string( ptr+count );
        job = Job::newJob( tmpS );
        if( job ) {
            count += job->unpack( ptr+count, swap_endian );
        }
        else throw invalid_argument( "Unrecognized Job tag: \"" + tmpS + "\"" );
    } else job = previousJob;
    if( job ) {
        count += job->unpackParts( ptr+count, *this, swap_endian );
    } else throw invalid_argument( "Can't unpack parts without a job instance..." );
    
    return count;
}


std::string WorkInProgress::print( void ) {
    
    string ret = "\"" + ( job ? job->info.name : string( "undefined" ) ) + "\"";
    if( parts.size() ) {
        ret += " (job: #" + (job ? to_string(job->info.id) : "0" ) + " part(s):";
        for( auto & part : parts ) {
            if( part ) {
                ret += " " + to_string( part->id );
            }
        }
        ret += ")";
    }
    
    return ret;
    
}
