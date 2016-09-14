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


Host::HostInfo::HostInfo( void ) : peerType(0), connectPort(0) {

    int one = 1;
    littleEndian = *(char*)&one;
    reduxVersion = getVersionNumber();
    pid = getpid();
    nCores = std::thread::hardware_concurrency();
    startedAt = boost::posix_time::second_clock::local_time();
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


Host::HostStatus::HostStatus( void ) : currentJob( 0 ), maxThreads( std::thread::hardware_concurrency() ), state( ST_IDLE ),
    loadAvg( 0 ), progress( 0 ), statusString("idle") {
        
    lastSeen = boost::posix_time::second_clock::local_time(); 
    lastActive = boost::posix_time::second_clock::local_time();
    nThreads = maxThreads;
    
}

Host::Host() : id(0), nConnections(0) {

}

Host::Host(const HostInfo& hi, uint64_t i) : info(hi), id(i), nConnections(0) {

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

    status.lastSeen = boost::posix_time::second_clock::local_time();
    
}


void Host::active(void) {

    status.lastActive = boost::posix_time::second_clock::local_time();
    status.state = ST_ACTIVE;
    
}


Host& Host::myInfo(void) {

    static Host singleton;
    return singleton;
    
}


std::string Host::printHeader(void) {
    string hdr = alignRight("ID",5) + alignCenter("NAME",25) + alignCenter("PID",7) + alignCenter("THREADS",9);
    hdr += alignLeft("VERSION",9) + alignCenter("loadavg",9) + alignCenter("Uptime",12) + alignCenter("Runtime",12) + "STATUS";
    return hdr;
}

std::string Host::print(void) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    boost::posix_time::time_duration elapsed = (now - status.lastActive);
    boost::posix_time::time_duration uptime = (now - info.startedAt);
    string ret = alignRight(std::to_string(id),5) + alignCenter(info.name,25) + alignCenter(to_string(info.pid),7);
    ret += alignCenter(to_string(status.nThreads) + string("/") + to_string(info.nCores),9);
    ret += alignCenter(getVersionString(info.reduxVersion),9);
    ret += alignCenter(boost::str(boost::format("%.2f") % status.loadAvg),9);
    ret += alignCenter(to_simple_string(uptime),12) + alignCenter(to_simple_string(elapsed),12);
    ret += status.statusString; 
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
    sz += name.length() + os.length() + arch.length() + 3;
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
    
    return count;
    
}


bool Host::HostInfo::operator==(const HostInfo& rhs) const {
    return (name == rhs.name && pid == rhs.pid);
}


uint64_t Host::HostStatus::size(void) const {
    uint64_t sz = sizeof(currentJob) + sizeof(nThreads) + sizeof(maxThreads);
    sz += sizeof(state) + sizeof(loadAvg) + sizeof(progress) + 2*sizeof(time_t);
    sz += statusString.length() + 1;
    return sz;
}


uint64_t Host::HostStatus::pack( char* ptr ) const {
    
    using redux::util::pack;
    
    uint64_t count = pack(ptr,nThreads);
    count += pack(ptr+count,maxThreads);
    count += pack(ptr+count,state);
    count += pack(ptr+count,currentJob);
    count += pack(ptr+count,loadAvg);
    count += pack(ptr+count,progress);
    count += pack(ptr+count,statusString);
    count += pack(ptr+count,redux::util::to_time_t( lastSeen ));
    count += pack(ptr+count,redux::util::to_time_t( lastActive ));
    
    return count;
    
}


uint64_t Host::HostStatus::unpack( const char* ptr, bool swap_endian ) {
    
    using redux::util::unpack;
    
    uint64_t count = unpack(ptr,nThreads, swap_endian);
    count += unpack(ptr+count,maxThreads, swap_endian);
    count += unpack(ptr+count,state, swap_endian);
    count += unpack(ptr+count,currentJob, swap_endian);
    count += unpack(ptr+count,loadAvg, swap_endian);
    count += unpack(ptr+count,progress, swap_endian);
    count += unpack(ptr+count,statusString, swap_endian);
    time_t timestamp;
    count += unpack(ptr+count,timestamp,swap_endian);
    lastSeen = boost::posix_time::from_time_t( timestamp );
    count += unpack(ptr+count,timestamp,swap_endian);
    lastActive = boost::posix_time::from_time_t( timestamp );
    
    return count;
    
}


TcpConnection& redux::network::operator<<(TcpConnection& conn, const Host::HostInfo& out) {
    
    uint64_t sz = out.size() + sizeof(uint64_t) + 1;
    shared_ptr<char> buf( new char[sz], []( char* p ){ delete[] p; } );
    char* ptr = buf.get();
    memset(ptr,0,sz);
    
    uint64_t count = pack(ptr,out.littleEndian);
    count += pack(ptr+count,out.size());
    
    out.pack(ptr+count);
    
    size_t written = boost::asio::write( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( written != sz ) {
        //LOG_ERR << "TcpConnection << HostInfo: Failed to write buffer.";
    }
    
    return conn;
    
}


TcpConnection& redux::network::operator>>(TcpConnection& conn, Host::HostInfo& in) {
    
    uint64_t sz = sizeof(uint64_t)+1;
    shared_ptr<char> buf( new char[sz], []( char* p ){ delete[] p; } );

    size_t readBytes = boost::asio::read( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( readBytes != sz ) {
        //LOG_ERR << "TcpConnection >> HostInfo: Failed to read buffer.";
    }
    
    int one = 1;
    char* ptr = buf.get();
    bool swap_endian = *((char*)&one) != *ptr;
    
    uint64_t count = unpack(ptr+1,sz,swap_endian);
    
    buf.reset( new char[sz], []( char* p ){ delete[] p; } );

    count = boost::asio::read( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( count != sz ) {
        //LOG_ERR << "TcpConnection >> HostInfo: Failed to read buffer.";
    }
    
    in.unpack(buf.get(), swap_endian);

    return conn;
    
}
