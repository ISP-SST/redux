#include "redux/work.hpp"

#include "redux/util/datautil.hpp"
#include "redux/job.hpp"

using namespace redux;
using namespace std;

Part::Part() : id( 0 ), step( Job::JSTEP_QUEUED ), nRetries( 0 ) {

}

size_t Part::size( void ) const {
    size_t sz = sizeof( id ) + 2;
    return sz;
}


char* Part::pack( char* ptr ) const {

    using redux::util::pack;

    ptr = pack( ptr, id );
    ptr = pack( ptr, step );
    ptr = pack( ptr, nRetries );

    return ptr;

}


const char* Part::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;

    ptr = unpack( ptr, id, swap_endian );
    ptr = unpack( ptr, step, swap_endian );
    ptr = unpack( ptr, nRetries, swap_endian );

    return ptr;
}


WorkInProgress::WorkInProgress( network::Peer::Ptr p ) : peer( p ) {

}

size_t WorkInProgress::size( bool includeJob ) const {
    size_t sz = sizeof( size_t ); // nParts
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


char* WorkInProgress::pack( char* ptr, bool includeJob ) const {
    using redux::util::pack;
    if( includeJob ) {
        ptr = job->pack( ptr );
    }
    ptr = pack( ptr, parts.size() );
    for( auto & it : parts ) {
        if( it ) {
            ptr = it->pack( ptr );
        }
    }
    return ptr;
}


const char* WorkInProgress::unpack( const char* cptr, bool includeJob, bool swap_endian ) {
    using redux::util::unpack;
    if( includeJob ) {
        string tmpS = string( cptr );
        job = Job::newJob( tmpS );
        if( job ) {
            cptr = job->unpack( cptr, swap_endian );
        }
        else throw invalid_argument( "Unrecognized Job tag: \"" + tmpS + "\"" );
    }
    if( job ) {
        cptr = job->unpackParts( cptr, parts, swap_endian );
    }
    else throw invalid_argument( "Can't unpack parts without a job instance..." );
    return cptr;
}


std::string WorkInProgress::print( void ) {
    string ret = "\"" + ( job ? job->info.name : string( "null" ) ) + "\"";
    if( parts.size() ) {
        ret += " (" + to_string( parts.size() ) + " part(s):";
        for( auto & it : parts ) {
            if( it ) {
                ret += " " + to_string( it->id );
            }
        }
        ret += ")";
    }
    return ret;
}
