#include "redux/network/host.hpp"

#include "redux/util/arrayutil.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/logger.hpp"
#include "redux/version.hpp"

#include <thread>
#include <unistd.h>
#include <sys/utsname.h>

#include <boost/date_time/posix_time/time_formatters.hpp>
#include <boost/date_time/posix_time/conversion.hpp>
#include <boost/algorithm/string.hpp>
using boost::algorithm::iequals;

using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;

const std::string Host::StateNames[] = { "offline", "idle", "active", "error" };
const std::string Host::TypeNames[] = { "", "worker", "master", "m/w", "util", "u/w", "u/m", "u/m/w" };

#define lg Logger::lg
namespace {
    const std::string thisChannel = "net";
}

Host::HostInfo::HostInfo( void ) : peerType(0) {

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


Host::HostStatus::HostStatus( void ) : currentJob( 0 ), nThreads( std::thread::hardware_concurrency() ), state( ST_IDLE ),
    loadAvg( 0 ), progress( 0 ) {
        
    lastSeen = boost::posix_time::second_clock::local_time(); 

}

Host::Host() : id(0), nConnections(0) {

}

Host::Host(const HostInfo& hi, uint64_t i) : info(hi), id(i), nConnections(0) {

}


size_t Host::size(void) const {
    size_t sz = sizeof(uint64_t) + info.size() + status.size();
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

std::string Host::printHeader(void) {
    string hdr = alignRight("ID",5) + alignCenter("HOST",15) + alignCenter("PID",7) + alignCenter("STATE",9);
    hdr += alignLeft("VERSION",9) + alignCenter("ARCH",8) + alignCenter("OS",25) + alignCenter("GFLOPS",8);
    return hdr;
}

std::string Host::print(void) {
    string ret = alignRight(std::to_string(id),5) + alignCenter(info.name,15) + alignCenter(to_string(info.pid),7);
    ret += alignCenter(StateNames[status.state],9) + alignCenter(getVersionString(info.reduxVersion),9);
    ret += alignCenter(info.arch,8) + alignCenter(info.os,25) + alignCenter("-",8);
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
        if( iequals( info.name, rhs.info.name ) ) {
            return (info.pid < rhs.info.pid);
        } else return ( info.name.compare( rhs.info.name ) > 0 );
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

size_t Host::HostInfo::size(void) const {
    size_t sz = sizeof(littleEndian) + sizeof(reduxVersion) + sizeof(pid) + sizeof(peerType) + sizeof(nCores);
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
    count += pack(ptr+count,to_time_t( startedAt ));
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


size_t Host::HostStatus::size(void) const {
    size_t sz = sizeof(currentJob) + sizeof(nThreads) + sizeof(maxThreads);
    sz += sizeof(state) + sizeof(loadAvg) + sizeof(progress) + sizeof(time_t);
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
    count += pack(ptr+count,to_time_t( lastSeen ));
    
    return count;
    
}


uint64_t Host::HostStatus::unpack( const char* ptr, bool swap_endian ) {
    
    using redux::util::unpack;
    
    uint64_t count = unpack(ptr,nThreads);
    count += unpack(ptr+count,maxThreads);
    count += unpack(ptr+count,state, swap_endian);
    count += unpack(ptr+count,currentJob, swap_endian);
    count += unpack(ptr+count,loadAvg, swap_endian);
    count += unpack(ptr+count,progress, swap_endian);
    time_t timestamp;
    count += unpack(ptr+count,timestamp,swap_endian);
    lastSeen = boost::posix_time::from_time_t( timestamp );
    
    return count;
    
}


TcpConnection& redux::network::operator<<(TcpConnection& conn, const Host::HostInfo& out) {
    
    size_t sz = out.size() + sizeof(size_t) + 1;
    auto buf = sharedArray<char>(sz);
    char* ptr = buf.get();
    memset(ptr,0,sz);
    
    uint64_t count = pack(ptr,out.littleEndian);
    count += pack(ptr+count,out.size());
    
    out.pack(ptr+count);
    
    size_t written = boost::asio::write( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( written != sz ) {
        LOG_ERR << "TcpConnection << HostInfo: Failed to write buffer.";
    }
    
    return conn;
    
}


TcpConnection& redux::network::operator>>(TcpConnection& conn, Host::HostInfo& in) {
    
    size_t sz = sizeof(size_t)+1;
    auto buf = sharedArray<char>(sz);

    size_t readBytes = boost::asio::read( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( readBytes != sz ) {
        LOG_ERR << "TcpConnection >> HostInfo: Failed to read buffer.";
    }
    
    int one = 1;
    char* ptr = buf.get();
    bool swap_endian = *((char*)&one) != *ptr;
    
    uint64_t count = unpack(ptr+1,sz,swap_endian);
    
    buf = sharedArray<char>(sz);

    count = boost::asio::read( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( count != sz ) {
        LOG_ERR << "TcpConnection >> HostInfo: Failed to read buffer.";
    }
    
    in.unpack(buf.get(), swap_endian);

    return conn;
    
}
