#include "redux/network/peer.hpp"

#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/logger.hpp"
#include "redux/version.hpp"

#include <unistd.h>
#include <sys/utsname.h>

#include <boost/date_time/posix_time/time_formatters.hpp>
#include <boost/date_time/posix_time/conversion.hpp>

using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::lg
namespace {
    const std::string thisChannel = "net";
}


Peer::HostInfo::HostInfo( void ) {

    int one = 1;
    littleEndian = *(char*)&one;
    reduxVersion = getVersionNumber();
    pid = getpid();
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


Peer::PeerStatus::PeerStatus( void ) : currentJob( 0 ), nThreads( 1 ), state( PEER_IDLE ),
    loadAvg( 0 ), progress( 0 ) {

}

Peer::Peer(const HostInfo& hi, TcpConnection::ptr c, size_t i) : host(hi), id(i), conn(c) {

}


size_t Peer::size(void) const {
    size_t sz = sizeof(size_t) + sizeof(time_t) + host.size() + stat.size();
    return sz;
}


char* Peer::pack( char* ptr ) const {
    
    using redux::util::pack;
    
    ptr = pack(ptr,id);
    ptr = pack(ptr,to_time_t( lastSeen ));
    ptr = host.pack(ptr);
    ptr = stat.pack(ptr);
    
    return ptr;
}

const char* Peer::unpack( const char* ptr, bool needsSwap ) {
    
    using redux::util::unpack;
    
    ptr = unpack(ptr,id);
    time_t timestamp;
    ptr = unpack(ptr,timestamp);
    ptr = host.unpack(ptr,needsSwap);
    ptr = stat.unpack(ptr,needsSwap);
    
    if( needsSwap ) {
        swapEndian(id);
        swapEndian(timestamp);
    }
    lastSeen = boost::posix_time::from_time_t( timestamp );
    
    return ptr;
}


std::string Peer::printHeader(void) {
    string hdr = alignRight("ID",5) + alignCenter("HOST",15) + alignCenter("PID",7) + alignCenter("STATE",9);
    hdr += alignLeft("VERSION",9) + alignCenter("ARCH",8) + alignCenter("OS",25) + alignCenter("GFLOPS",8);
    return hdr;
}

std::string Peer::print(void) {
    string info = alignRight(std::to_string(id),5) + alignCenter(host.name,15) + alignCenter(to_string(host.pid),7);
    info += alignCenter(StateNames[stat.state],9) + alignCenter(getVersionString(host.reduxVersion),9);
    info += alignCenter(host.arch,8) + alignCenter(host.os,25) + alignCenter("-",8);
    return info;
}


bool Peer::operator>( const Peer& rhs ) const {
    //if( hostType != rhs.hostType ) return ( hostType > rhs.hostType );
    //if (IP != rhs.IP) return (IP > rhs.IP);
    //if (listeningPort != rhs.listeningPort) return (listeningPort > rhs.listeningPort);
    //if (PID != rhs.PID) return (PID > rhs.PID);
    return ( host.name.compare( rhs.host.name ) < 0 );
}


bool Peer::operator<( const Peer& rhs ) const {
    // if( hostType != rhs.hostType ) return ( hostType < rhs.hostType );
    //if (IP != rhs.IP) return (IP < rhs.IP);
    //if ( listeningPort != rhs.listeningPort ) return (listeningPort < rhs.listeningPort);
    //if ( hostName.compare( rhs.hostName ) == 0 ) return (PID < rhs.PID);
    return ( host.name.compare( rhs.host.name ) > 0 );
}


// use this one to compare the "configuration" of 2 hosts.
bool Peer::operator==( const Peer& rhs ) const {
    //if (IP && rhs.IP) return (IP == rhs.IP);
    //else
    return ( !host.name.compare( rhs.host.name ) );
}


// use this one to differentiate between different computers
bool Peer::operator!=( const Peer& rhs ) const {
    return ( host.name != rhs.host.name );
    uint tmp = 1; //IP;
    uint tmpRHS = 1; //rhs.IP;
    //if ( !tmp ) hostLookup( hostName, &tmp );
    //if ( !tmpRHS ) hostLookup( rhs.hostName, &tmpRHS );
    return ( tmp != tmpRHS );
    //if ( (hostName.size() > 0) && (rhs.hostName.size() > 0) ) return (bool)(hostName.compare( rhs.hostName ));
}


Peer& Peer::operator=( const Peer& rhs ) {
    if( &rhs == this ) return *this;

    host.name = rhs.host.name;
    stat.state = rhs.stat.state;


    return *this;

}


void Peer::HostInfo::print( void ) {
    
    cout << "Peer::HostInfo  littleEndian = " << ( int )littleEndian << endl;
    cout << "Peer::HostInfo  reduxVersion = " << reduxVersion << endl;
    cout << "Peer::HostInfo           pid = " << pid << endl;
    cout << "Peer::HostInfo     startedAt = " << to_iso_extended_string( startedAt ) << endl;
    cout << "Peer::HostInfo          name = " << name << endl;
    cout << "Peer::HostInfo            os = " << os << endl;
    cout << "Peer::HostInfo       machine = " << arch << endl;
    
}

size_t Peer::HostInfo::size(void) const {
    size_t sz = sizeof(littleEndian) + sizeof(reduxVersion) + sizeof(pid) + sizeof(time_t);
    sz += name.length() + os.length() + arch.length() + 3;
    return sz;
}


char* Peer::HostInfo::pack( char* ptr ) const {
    
    using redux::util::pack;
    
    *ptr++ = littleEndian;
    ptr = pack(ptr,reduxVersion);
    ptr = pack(ptr,pid);
    ptr = pack(ptr,to_time_t( startedAt ));
    ptr = pack(ptr,name);
    ptr = pack(ptr,os);
    ptr = pack(ptr,arch);
    
    return ptr;
}

const char* Peer::HostInfo::unpack( const char* ptr, bool needsSwap ) {
    
    using redux::util::unpack;
    
    littleEndian = *ptr++;
    ptr = unpack(ptr,reduxVersion);
    ptr = unpack(ptr,pid);
    time_t timestamp;
    ptr = unpack(ptr,timestamp);
    ptr = unpack(ptr,name);
    ptr = unpack(ptr,os);
    ptr = unpack(ptr,arch);
    
    if( needsSwap ) {
        swapEndian(reduxVersion);
        swapEndian(pid);
        swapEndian(timestamp);
    }
    startedAt = boost::posix_time::from_time_t( timestamp );
    
    return ptr;
}


bool Peer::HostInfo::operator==(const HostInfo& rhs) {
    return (name == rhs.name && pid == rhs.pid);
}


size_t Peer::PeerStatus::size(void) const {
    size_t sz = sizeof(size_t) + sizeof(State) + 2*sizeof(float) + 1;
    return sz;
}


char* Peer::PeerStatus::pack( char* ptr ) const {
    
    using redux::util::pack;
    
    *ptr++ = nThreads;
    ptr = pack(ptr,state);
    ptr = pack(ptr,currentJob);
    ptr = pack(ptr,loadAvg);
    ptr = pack(ptr,progress);
    
    return ptr;
}


const char* Peer::PeerStatus::unpack( const char* ptr, bool needsSwap ) {
    
    using redux::util::unpack;
    
    nThreads = *ptr++;
    ptr = unpack(ptr,state);
    ptr = unpack(ptr,currentJob);
    ptr = unpack(ptr,loadAvg);
    ptr = unpack(ptr,progress);
    
    if( needsSwap ) {
        swapEndian(currentJob);
        swapEndian(loadAvg);
        swapEndian(progress);
    }
    
    return ptr;
}


TcpConnection& redux::network::operator<<(TcpConnection& conn, const Peer::HostInfo& out) {
    
    size_t sz = out.size() + sizeof(size_t) + 1;
    std::unique_ptr<char[]> buf(new char[sz]);
    char* ptr = buf.get();
    memset(ptr,0,sz);
    
    *ptr++ = out.littleEndian;
    ptr = pack(ptr,out.size());
    
    out.pack(ptr);
    
    conn.writeAndCheck(buf.get(),sz);
    
    return conn;
    
}
TcpConnection& redux::network::operator>>(TcpConnection& conn, Peer::HostInfo& in) {
    
    size_t hdrSz = sizeof(size_t)+1;
    std::unique_ptr<char[]> buf(new char[hdrSz]);

    size_t count = boost::asio::read( conn.socket(), boost::asio::buffer(buf.get(), hdrSz) );
    if( count != hdrSz ) {
        LOG_ERR << "TcpConnection << HostInfo: Failed to read buffer.";
    }
    
    int one = 1;
    char* ptr = buf.get();
    bool needsSwap = *((char*)&one) != *reinterpret_cast<uint8_t*>(ptr++);
    
    size_t sz;
    unpack(ptr,sz);
    if(needsSwap) swapEndian(sz);
    
    buf.reset(new char[sz]);

    count = boost::asio::read( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( count != sz ) {
        LOG_ERR << "TcpConnection << HostInfo: Failed to read buffer.";
    }
    
    in.unpack(buf.get(), needsSwap);

    return conn;
    
}
