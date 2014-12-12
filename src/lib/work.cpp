#include "redux/work.hpp"

#include "redux/util/datautil.hpp"
#include "redux/job.hpp"

using namespace redux;
using namespace std;

Part::Part() : id( 0 ), step( Job::JSTEP_QUEUED ), nRetries( 0 ) {

}

size_t Part::size( void ) const {
    return sizeof( id ) + sizeof( step ) + sizeof( nRetries ) ;
}


uint64_t Part::pack( char* ptr ) const {

    using redux::util::pack;
    
    uint64_t count = pack( ptr, id );
    count += pack( ptr+count, step );
    count += pack( ptr+count, nRetries );

    return count;

}


uint64_t Part::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;
    
    uint64_t count = unpack( ptr, id, swap_endian );
    count += unpack( ptr+count, step, swap_endian );
    count += unpack( ptr+count, nRetries, swap_endian );
    
    return count;
}


WorkInProgress::WorkInProgress( network::TcpConnection::Ptr c ) : connection( c ), nCompleted(0) {

}


size_t WorkInProgress::size( bool includeJob ) const {
    
    size_t sz = sizeof( size_t ) + sizeof(nCompleted) + 1 ; // nParts + includeJob
    if( includeJob ) {
        sz += job->size();
    }

    for( auto & it : parts ) {
        if( it ) {
            sz += it->size();
        }
    }
    
    return sz;
    
}


uint64_t WorkInProgress::pack( char* ptr, bool includeJob ) const {

    using redux::util::pack;
    
    uint64_t count(0);
    if( includeJob ) {
        count += pack( ptr, uint8_t(1) );
        count += job->pack( ptr+count );
    } else count += pack( ptr, uint8_t(0) );
    count += pack( ptr+count, parts.size() );
    for( auto & it : parts ) {
        if( it ) {
            count += it->pack( ptr+count );
        }
    }
    count += pack( ptr+count, nCompleted );
    
    return count;
    
}


uint64_t WorkInProgress::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;
    
    uint8_t hasJob;
    uint64_t count = unpack( ptr, hasJob );
    if( hasJob ) {
        string tmpS = string( ptr+count );
        job = Job::newJob( tmpS );
        if( job ) {
            count += job->unpack( ptr+count, swap_endian );
        }
        else throw invalid_argument( "Unrecognized Job tag: \"" + tmpS + "\"" );
    }
    if( job ) {
        count += job->unpackParts( ptr+count, parts, swap_endian );
    }
    else throw invalid_argument( "Can't unpack parts without a job instance..." );
    
    count += unpack( ptr+count, nCompleted, swap_endian );
    
    return count;
    
}


std::string WorkInProgress::print( void ) {
    
    string ret = "\"" + ( job ? job->info.name : string( "null" ) ) + "\"";
    if( parts.size() ) {
        ret += " (" + to_string( nCompleted ) + "/" + to_string( parts.size() ) + " part(s):";
        for( auto & it : parts ) {
            if( it ) {
                ret += " " + to_string( it->id );
            }
        }
        ret += ")";
    }
    
    return ret;
    
}
