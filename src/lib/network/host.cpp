#include "redux/network/host.hpp"

#include "redux/util/arrayutil.hpp"
#include "redux/util/convert.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/logging/logger.hpp"
#include "redux/version.hpp"

#include <thread>
#include <unistd.h>
#include <sys/utsname.h>

#include <boost/date_time/posix_time/time_formatters.hpp>
#include <boost/date_time/posix_time/conversion.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
using boost::algorithm::iequals;

using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;

const std::string Host::StateNames[] = { "offline", "idle", "active", "error" };
const std::string Host::TypeNames[] = { "", "worker", "master", "m/w", "util", "u/w", "u/m", "u/m/w" };


namespace {
    
    set<uint64_t> host_ids;
    mutex gmtx;

    uint64_t getID( void ) {
        lock_guard<mutex> lock( gmtx );
        uint64_t id(0);
        while( host_ids.count( id ) ) id++;     // find first unused ID.
        host_ids.insert( id );
        return id;
    }
    
    void freeID( uint64_t id ) {
        lock_guard<mutex> lock( gmtx );
        host_ids.erase( id );
    }
    
    /*size_t idCount( void )  {
        lock_guard<mutex> lock( gmtx );
        return host_ids.size();
    }*/
    
}

Host::HostInfo::HostInfo( string username ) : peerType(0), connectPort(0), user(username) {

    int one = 1;
    littleEndian = *(char*)&one;
    reduxVersion = getVersionNumber();
    pid = getpid();
    nCores = std::thread::hardware_concurrency();
    startedAt = boost::posix_time::second_clock::universal_time();
    name = boost::asio::ip::host_name();
    struct utsname nm;
    if( uname( &nm ) < 0 ) {
        os = "-";
        arch = "-";
    } else {
        os = nm.sysname + string( " " ) + nm.release;
        arch = nm.machine;
    }

}


Host::HostStatus::HostStatus( void ) : currentJob( 0 ), maxThreads( std::thread::hardware_concurrency() ), listenPort(0),
    state( ST_IDLE ), progress( 0 ), statusString("idle") {
        
    lastSeen = boost::posix_time::second_clock::universal_time(); 
    lastActive = boost::posix_time::second_clock::universal_time();
    nThreads = maxThreads;
    load[0] = load[1] = 0;
}


Host::Host() : id(getID()), nConnections(0) {

}


Host::Host(const HostInfo& hi, uint64_t i) : info(hi), id(getID()), nConnections(0) {

}


Host::~Host() {
    freeID( id );
}


uint64_t Host::size(void) const {
    uint64_t sz = sizeof(uint64_t) + info.size() + status.size();
    return sz;
}


uint64_t Host::pack( char* ptr ) const {
    
    using redux::util::pack;
    uint64_t count = pack(ptr,id);
    count += info.pack(ptr+count);
    count += status.pack(ptr+count);
    
    return count;
}


uint64_t Host::unpack( const char* ptr, bool swap_endian ) {
    
    using redux::util::unpack;
    uint64_t count = unpack(ptr,id,swap_endian);
    count += info.unpack(ptr+count,swap_endian);
    count += status.unpack(ptr+count,swap_endian);
    
    return count;
}


void Host::touch(void) {

    status.lastSeen = boost::posix_time::second_clock::universal_time();
    
}


void Host::active(void) {

    status.lastActive = boost::posix_time::second_clock::universal_time();
    status.state = ST_ACTIVE;
    
}


void Host::limbo(void) {

    status.lastActive = boost::posix_time::second_clock::universal_time();
    status.state = ST_LIMBO;
    
}


void Host::idle(void) {

    status.state = ST_IDLE;
    status.statusString = "idle";

}


Host& Host::myInfo(void) {

    static Host singleton( getUname() );
    return singleton;
    
}


std::string Host::printHeader( int verbosity ) {
    string hdr = alignRight("ID",5) + alignCenter("NAME",25) + alignCenter("PID",7) + alignCenter("THREADS",10);
    hdr += alignLeft("VERSION",10) + alignLeft("Usage/Load",12) + alignCenter("Uptime",12) + alignCenter("Runtime",12);
    hdr += alignLeft("STATUS",18);
    if( verbosity ) hdr +=  alignCenter("port",6);
    return hdr;
}


std::string Host::print( int verbosity ) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::universal_time();
    boost::posix_time::time_duration elapsed = (now - status.lastActive);
    boost::posix_time::time_duration uptime = (now - info.startedAt);
    string elapsedString = "";
    if( (status.state != ST_IDLE) && !elapsed.is_not_a_date_time() ) elapsedString = to_simple_string(elapsed);
    string ret = alignRight(std::to_string(id),5) + alignCenter(info.name,25) + alignCenter(to_string(info.pid),7);
    ret += alignCenter(to_string(status.nThreads) + string("/") + to_string(info.nCores),10);
    ret += alignLeft(getVersionString(info.reduxVersion),10);
    ret += alignLeft(boost::str(boost::format("%.2f/%.2f") % status.load[0] % status.load[1]),12);
    ret += alignCenter(to_simple_string(uptime),12) + alignCenter(elapsedString,12);
    ret += alignLeft(status.statusString,18); 
    if( verbosity ) {
        if( id && status.listenPort ) {
            ret += alignRight( to_string(status.listenPort), 6 );
        }
    }
    return ret;
}


bool Host::operator>( const Host& rhs ) const {
//    if( info.peerType == rhs.info.peerType ) {
        if( iequals( info.name, rhs.info.name ) ) {
            return (info.pid > rhs.info.pid);
        } else return ( info.name.compare( rhs.info.name ) < 0 );
//    } else return (info.peerType > rhs.info.peerType);
}


bool Host::operator<( const Host& rhs ) const {
//    if( info.peerType == rhs.info.peerType ) {
        if( iequals( info.name, rhs.info.name ) ) return (info.pid < rhs.info.pid);
        return ( info.name < rhs.info.name );
//    } else return (info.peerType < rhs.info.peerType);
}


bool Host::operator==( const Host& rhs ) const {
    return (info == rhs.info);
}

/*
Host& Host::operator=( const Host& rhs ) {

   if( &rhs == this ) return *this;

    info.name = rhs.info.name;
    status.state = rhs.status.state;

    return *this;

}
*/

uint64_t Host::HostInfo::size(void) const {
    uint64_t sz = sizeof(littleEndian) + sizeof(reduxVersion) + sizeof(pid) + sizeof(peerType) + sizeof(nCores);
    sz += sizeof(time_t);   // startedAt is converted and transferred as time_t
    sz += name.length() + os.length() + arch.length() + user.size() + 4;
    return sz;
}


uint64_t Host::HostInfo::pack( char* ptr ) const {
    
    using redux::util::pack;
    
    *ptr = littleEndian;
    uint64_t count = 1;     // endianflag
    count += pack(ptr+count,reduxVersion);
    count += pack(ptr+count,pid);
    count += pack(ptr+count,peerType);
    count += pack(ptr+count,nCores);
    count += pack(ptr+count,redux::util::to_time_t( startedAt ));
    count += pack(ptr+count,name);
    count += pack(ptr+count,os);
    count += pack(ptr+count,arch);
    count += pack(ptr+count,user);
    
    return count;
    
}


uint64_t Host::HostInfo::unpack( const char* ptr, bool swap_endian ) {
    
    using redux::util::unpack;
    
    littleEndian = *ptr;
    uint64_t count = 1;     // endianflag
    count += unpack(ptr+count,reduxVersion, swap_endian);
    count += unpack(ptr+count,pid,swap_endian);
    count += unpack(ptr+count,peerType, swap_endian);
    count += unpack(ptr+count,nCores, swap_endian);
    time_t timestamp;
    count += unpack(ptr+count,timestamp,swap_endian);
    startedAt = boost::posix_time::from_time_t( timestamp );
    count += unpack(ptr+count,name);
    count += unpack(ptr+count,os);
    count += unpack(ptr+count,arch);
    count += unpack(ptr+count,user);
    
    return count;
    
}


bool Host::HostInfo::operator==(const HostInfo& rhs) const {
    return (name == rhs.name && pid == rhs.pid);
}


uint64_t Host::HostStatus::size(void) const {
    uint64_t sz = sizeof(currentJob) + sizeof(nThreads) + sizeof(maxThreads) + sizeof(listenPort);
    sz += sizeof(state) + sizeof(load) + sizeof(progress) + 2*sizeof(time_t);
    sz += statusString.length() + 1;
    return sz;
}


uint64_t Host::HostStatus::pack( char* ptr ) const {
    
    using redux::util::pack;
    
    uint64_t count = pack(ptr,nThreads);
    count += pack(ptr+count,maxThreads);
    count += pack(ptr+count,listenPort);
    count += pack(ptr+count,state);
    count += pack(ptr+count,currentJob);
    count += pack(ptr+count,load,2);
    count += pack(ptr+count,progress);
    count += pack(ptr+count,statusString);
    time_t tmpT(0);
    if( !lastSeen.is_not_a_date_time() ) tmpT = redux::util::to_time_t( lastSeen );
    count += pack(ptr+count, tmpT);
    tmpT = 0;
    if( !lastActive.is_not_a_date_time() ) tmpT = redux::util::to_time_t( lastActive );
    count += pack(ptr+count, tmpT);
    
    return count;
    
}


uint64_t Host::HostStatus::unpack( const char* ptr, bool swap_endian ) {
    
    using redux::util::unpack;
    
    lastSeen = boost::posix_time::ptime( boost::posix_time::not_a_date_time );
    lastActive = boost::posix_time::ptime( boost::posix_time::not_a_date_time );
    
    uint64_t count = unpack(ptr,nThreads, swap_endian);
    count += unpack(ptr+count,maxThreads, swap_endian);
    count += unpack(ptr+count,listenPort, swap_endian);
    count += unpack(ptr+count,state, swap_endian);
    count += unpack(ptr+count,currentJob, swap_endian);
    count += unpack(ptr+count,load, 2, swap_endian);
    count += unpack(ptr+count,progress, swap_endian);
    count += unpack(ptr+count,statusString, swap_endian);
    time_t timestamp;
    count += unpack(ptr+count,timestamp,swap_endian);
    if( timestamp ) lastSeen = boost::posix_time::from_time_t( timestamp );
    count += unpack(ptr+count,timestamp,swap_endian);
    if( timestamp ) lastActive = boost::posix_time::from_time_t( timestamp );
    
    return count;
    
}


TcpConnection& redux::network::operator<<(TcpConnection& conn, const Host::HostInfo& out) {
    
    uint64_t sz = out.size() + sizeof(uint64_t) + 1;
    shared_ptr<char> buf = rdx_get_shared<char>(sz);
    char* ptr = buf.get();
    memset( ptr, 0, sz );
    
    uint64_t count = pack( ptr, out.littleEndian );
    count += pack( ptr+count, out.size() );
    
    out.pack( ptr+count );
    
    size_t written = boost::asio::write( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( written != sz ) {
        //LOG_ERR << "TcpConnection << HostInfo: Failed to write buffer.";
    }
    
    return conn;
    
}


TcpConnection& redux::network::operator>>( TcpConnection& conn, Host::HostInfo& in ) {
    
    constexpr size_t hdrSize = sizeof(uint64_t)+1;
    std::array<char,hdrSize> hdr;

    size_t readBytes = boost::asio::read( conn.socket(), boost::asio::buffer(hdr.data(), hdrSize) );
    if( readBytes != hdrSize ) {
        //LOG_ERR << "TcpConnection >> HostInfo: Failed to read buffer.";
    }
    
    int one = 1;
    char* ptr = hdr.data();
    bool swap_endian = *((char*)&one) != *ptr;
    
    uint64_t sz(0);
    uint64_t count = unpack( ptr+1, sz, swap_endian );
    
    shared_ptr<char> buf = rdx_get_shared<char>(sz);

    count = boost::asio::read( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( count != sz ) {
        //LOG_ERR << "TcpConnection >> HostInfo: Failed to read buffer.";
    }
    
    in.unpack( buf.get(), swap_endian );

    return conn;
    
}
