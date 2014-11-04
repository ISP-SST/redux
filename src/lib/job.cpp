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


#define lg Logger::lg
namespace {

    const string thisChannel = "job";

    const std::string StateNames[10] = { "undefined", "preprocessed", "queued", "done", "completed", "postprocess", "idle", "active", "paused", "error" };
    const std::string StateTags[10] = { "-", "Pre", "Q", "D", "C", "Post", "I", "A", "P", "E" };

    mutex globalJobMutex;

}


size_t Job::registerJob( const string& name, JobCreator f ) {
    static size_t nJobTypes( 0 );
    std::unique_lock<mutex> lock( globalJobMutex );
    auto ret = getMap().insert( { boost::to_upper_copy( name ), {nJobTypes + 1, f}} );
    if( ret.second ) {
        //LOG_DEBUG << "Job with tag \"" << name << "\" successfully registered.";
        return ret.first->second.first;
    } //else LOG_WARN << "Failed to register job with tag \"" << name << "\".";
    return nJobTypes++;
}


vector<Job::JobPtr> Job::parseTree( po::variables_map& vm, bpt::ptree& tree ) {
    vector<JobPtr> tmp;
    std::unique_lock<mutex> lock( globalJobMutex );
    for( auto & it : tree ) {
        string name = it.first;
        auto it2 = getMap().find( boost::to_upper_copy( name ) );       // check if the current tag matches a registered (Job-derived) class.
        if( it2 != getMap().end() ) {
            Job* tmpJob = it2->second.second();
            tmpJob->parseProperties( vm, it.second );
            tmp.push_back( shared_ptr<Job>( tmpJob ) );
        }
        //else LOG_WARN << "No job with tag \"" << name << "\" registered.";
    }
    return tmp;
}


Job::JobPtr Job::newJob( const string& name ) {
    JobPtr tmp;
    std::unique_lock<mutex> lock( globalJobMutex );
    auto it = getMap().find( boost::to_upper_copy( name ) );
    if( it != getMap().end() ) {
        tmp.reset( it->second.second() );
    }
    else LOG_WARN << "No job with tag \"" << name << "\" registered.";

    return tmp;
}

string Job::stepString( uint8_t step ) {
    switch( step ) {
        case JSTEP_SUBMIT: return "submit";
        case JSTEP_RECEIVED: return "received";
        case JSTEP_QUEUED: return "queued";
        case JSTEP_RUNNING: return "running";
        case JSTEP_POSTPROCESS: return "postprocess";
        case JSTEP_COMPLETED: return "completed";
        case JSTEP_ERR: return "error";
        default: return "-";
    }
}


string Job::stateString( uint8_t state ) {

    switch( state ) {
        case JSTATE_IDLE: return "idle";
        case JSTATE_ACTIVE: return "active";
        case JSTATE_PAUSED: return "paused";
        case JSTATE_CANCELLED: return "cancelled";
        case JSTATE_ERR: return "error";
        default: return "-";
    }

}

string Job::stateTag( uint8_t state ) {

    switch( state ) {
        case JSTATE_IDLE: return "I";
        case JSTATE_ACTIVE: return "A";
        case JSTATE_PAUSED: return "P";
        case JSTATE_CANCELLED: return "C";
        case JSTATE_ERR: return "E";
        default: return "";
    }

}


Job::Info::Info( void ) : id( 0 ), priority( 0 ), verbosity( 0 ), step( JSTEP_PRE_SUBMIT ), state( JSTATE_IDLE ) {

}


size_t Job::Info::size( void ) const {
    size_t sz = sizeof( size_t ) + sizeof( time_t ) + 6;
    sz += typeString.length() + name.length() + user.length() + host.length() + logFile.length() + 5 ;
    return sz;
}


uint64_t Job::Info::pack( char* ptr ) const {

    using redux::util::pack;
    
    uint64_t count = pack( ptr, typeString );    // NB: the type-string has to be first in the packed block,
    count += pack( ptr+count, name );            //   it is used to identify which job-class to instantiate on the receiving side.
    count += pack( ptr+count, user );
    count += pack( ptr+count, host );
    count += pack( ptr+count, logFile );
    count += pack( ptr+count, id );
    count += pack( ptr+count, priority );
    count += pack( ptr+count, verbosity );
    count += pack( ptr+count, nThreads );
    count += pack( ptr+count, maxPartRetries );
    count += pack( ptr+count, step.load() );
    count += pack( ptr+count, state.load() );
    count += pack( ptr+count, to_time_t( submitTime ) );

    return count;

}


uint64_t Job::Info::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;
    
    uint64_t count = unpack( ptr, typeString );
    count += unpack( ptr+count, name );
    count += unpack( ptr+count, user );
    count += unpack( ptr+count, host );
    count += unpack( ptr+count, logFile );
    count += unpack( ptr+count, id, swap_endian );
    count += unpack( ptr+count, priority );
    count += unpack( ptr+count, verbosity );
    count += unpack( ptr+count, nThreads );
    count += unpack( ptr+count, maxPartRetries );
    uint8_t tmp;
    count += unpack( ptr+count, tmp );
    step.store( tmp );
    count += unpack( ptr+count, tmp );
    state.store( tmp );
    time_t timestamp;
    count += unpack( ptr+count, timestamp, swap_endian );
    submitTime = boost::posix_time::from_time_t( timestamp );

    return count;
}


std::string Job::Info::printHeader( void ) {
    string hdr = alignRight( "ID", 5 ) + alignCenter( "type", 10 ) + alignCenter( "submitted", 20 );
    hdr += alignCenter( "name", 15 ) + alignLeft( "user", 15 ) + alignCenter( "priority", 8 ) + alignCenter( "state", 8 );
    return hdr;
}


std::string Job::Info::print( void ) {
    string info = alignRight( std::to_string( id ), 5 ) + alignCenter( typeString, 10 );
    info += alignCenter( to_iso_extended_string( submitTime ), 20 );
    info += alignCenter( name, 15 ) + alignLeft( user + "@" + host, 15 ) + alignCenter( std::to_string( priority ), 8 ) + alignCenter( stateTag( state ), 3 ) + alignLeft( stepString( step ), 15 );
    return info;
}


void Job::parseProperties( po::variables_map&, bpt::ptree& tree ) {
    info.priority = tree.get<uint8_t>( "PRIORITY", 10 );
    info.logFile = tree.get<string>( "LOGFILE", "" );
    info.verbosity = tree.get<uint8_t>( "VERBOSITY", 0 );
    info.nThreads = tree.get<uint8_t>( "MAX_THREADS", 255 );
    info.maxPartRetries = tree.get<uint8_t>( "MAX_PART_RETRIES", 0 );
}


bpt::ptree Job::getPropertyTree( bpt::ptree* root ) {

    bpt::ptree tree;

    tree.put( "PRIORITY", info.priority );
    tree.put( "LOGFILE", info.logFile );
    tree.put( "VERBOSITY", info.verbosity );
    tree.put( "MAX_THREADS", info.nThreads );
    tree.put( "MAX_PART_RETRIES", info.maxPartRetries );

    if( root ) {
        root->push_back( bpt::ptree::value_type( "job", tree ) );
    }

    return tree;

}


Job::Job( void ) {
    info.user = getUname();
    info.host = boost::asio::ip::host_name();
}


Job::~Job( void ) {

}


size_t Job::size( void ) const {
    size_t sz = info.size();
    return sz;
}


uint64_t Job::pack( char* ptr ) const {
    return info.pack( ptr );
}


uint64_t Job::unpack( const char* ptr, bool swap_endian ) {
    return info.unpack( ptr, swap_endian );
}


bool Job::operator<( const Job& rhs ) {
    return ( info.id < rhs.info.id );
}


bool Job::operator!=( const Job& rhs ) {
    return ( info.id != rhs.info.id );
}
