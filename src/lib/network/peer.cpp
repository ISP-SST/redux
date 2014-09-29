#include "redux/network/peer.hpp"

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

const std::string Peer::StateNames[] = { "offline", "idle", "active", "error" };
const std::string Peer::TypeNames[] = { "", "worker", "master", "m/w", "util", "u/w", "u/m", "u/m/w" };

#define lg Logger::lg
namespace {
    const std::string thisChannel = "net";
}

Peer::HostInfo::HostInfo( void ) : peerType(0) {

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


Peer::PeerStatus::PeerStatus( void ) : currentJob( 0 ), nThreads( std::thread::hardware_concurrency() ), state( PEER_IDLE ),
    loadAvg( 0 ), progress( 0 ) {

}

Peer::Peer() : id(0) {

}

Peer::Peer(const HostInfo& hi, TcpConnection::Ptr c, size_t i) : host(hi), id(i), conn(c) {

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

const char* Peer::unpack( const char* ptr, bool swap_endian ) {
    
    using redux::util::unpack;
    
    ptr = unpack(ptr,id,swap_endian);
    time_t timestamp;
    ptr = unpack(ptr,timestamp,swap_endian);
    ptr = host.unpack(ptr,swap_endian);
    ptr = stat.unpack(ptr,swap_endian);
    
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


shared_ptr<char> Peer::receiveBlock( size_t& blockSize, bool& swap_endian ) {

    using redux::util::unpack;
    shared_ptr<char> buf;
    
    try {
        char sz[sizeof( size_t )];
        size_t count = boost::asio::read( conn->socket(), boost::asio::buffer( sz, sizeof( size_t ) ) );

        if( count != sizeof( size_t ) ) {
            throw ios_base::failure( "blockSize: Only received " + to_string( count ) + "/" + to_string( sizeof( size_t ) ) + " bytes." );
        }
        
        int one = 1;
        char littleEndian = *(char*)&one;
        swap_endian = ( host.littleEndian != littleEndian );
        unpack( sz, blockSize, swap_endian );

        if(blockSize > 0) {
            buf.reset( new char[ blockSize ], []( char * p ) { delete[] p; } );
            count = boost::asio::read( conn->socket(), boost::asio::buffer( buf.get(), blockSize ) );

            if( count != blockSize ) {
                throw ios_base::failure( "Only received " + to_string( count ) + "/" + to_string( blockSize ) + " bytes." );
            }
        }
    }
    catch( const exception& e ) {
        LOG_ERR << "Failed to receive datablock: " << e.what();
        blockSize = 0;
        buf.reset();
    }
    return buf;

}


bool Peer::operator>( const Peer& rhs ) const {
    if( host.peerType == rhs.host.peerType ) {
        if( iequals( host.name, rhs.host.name ) ) {
            return (host.pid > rhs.host.pid);
        } else return ( host.name.compare( rhs.host.name ) < 0 );
    } else return (host.peerType > rhs.host.peerType);
}


bool Peer::operator<( const Peer& rhs ) const {
    if( host.peerType == rhs.host.peerType ) {
        if( iequals( host.name, rhs.host.name ) ) {
            return (host.pid < rhs.host.pid);
        } else return ( host.name.compare( rhs.host.name ) > 0 );
    } else return (host.peerType < rhs.host.peerType);
}


bool Peer::operator==( const Peer& rhs ) const {
    return (host == rhs.host);
}


Peer& Peer::operator=( const Peer& rhs ) {
    if( &rhs == this ) return *this;

    host.name = rhs.host.name;
    stat.state = rhs.stat.state;


    return *this;

}


size_t Peer::HostInfo::size(void) const {
    size_t sz = sizeof(littleEndian) + sizeof(reduxVersion) + sizeof(pid) + sizeof(uint16_t) + sizeof(time_t);
    sz += name.length() + os.length() + arch.length() + 4;
    return sz;
}


char* Peer::HostInfo::pack( char* ptr ) const {
    
    using redux::util::pack;
    
    *ptr++ = littleEndian;
    ptr = pack(ptr,reduxVersion);
    ptr = pack(ptr,pid);
    ptr = pack(ptr,peerType);
    ptr = pack(ptr,nCores);
    ptr = pack(ptr,to_time_t( startedAt ));
    ptr = pack(ptr,name);
    ptr = pack(ptr,os);
    ptr = pack(ptr,arch);
    
    return ptr;
}

const char* Peer::HostInfo::unpack( const char* ptr, bool swap_endian ) {
    
    using redux::util::unpack;
    
    littleEndian = *ptr++;
    ptr = unpack(ptr,reduxVersion, swap_endian);
    ptr = unpack(ptr,pid,swap_endian);
    ptr = unpack(ptr,peerType, swap_endian);
    ptr = unpack(ptr,nCores, swap_endian);
    time_t timestamp;
    ptr = unpack(ptr,timestamp,swap_endian);
    startedAt = boost::posix_time::from_time_t( timestamp );
    ptr = unpack(ptr,name);
    ptr = unpack(ptr,os);
    ptr = unpack(ptr,arch);
    
    return ptr;
}


bool Peer::HostInfo::operator==(const HostInfo& rhs) const {
    return (name == rhs.name && pid == rhs.pid);
}


size_t Peer::PeerStatus::size(void) const {
    size_t sz = sizeof(size_t) + sizeof(State) + 2*sizeof(float) + 2;
    return sz;
}


char* Peer::PeerStatus::pack( char* ptr ) const {
    
    using redux::util::pack;
    
    ptr = pack(ptr,nThreads);
    ptr = pack(ptr,maxThreads);
    ptr = pack(ptr,state);
    ptr = pack(ptr,currentJob);
    ptr = pack(ptr,loadAvg);
    ptr = pack(ptr,progress);
    
    return ptr;
}


const char* Peer::PeerStatus::unpack( const char* ptr, bool swap_endian ) {
    
    using redux::util::unpack;
    
    ptr = unpack(ptr,nThreads);
    ptr = unpack(ptr,maxThreads);
    ptr = unpack(ptr,state, swap_endian);
    ptr = unpack(ptr,currentJob, swap_endian);
    ptr = unpack(ptr,loadAvg, swap_endian);
    ptr = unpack(ptr,progress, swap_endian);
    
    return ptr;
}


TcpConnection& redux::network::operator<<(TcpConnection& conn, const Peer::HostInfo& out) {
    
    size_t sz = out.size() + sizeof(size_t) + 1;
    auto buf = sharedArray<char>(sz);
    char* ptr = buf.get();
    memset(ptr,0,sz);
    
    ptr = pack(ptr,out.littleEndian);
    ptr = pack(ptr,out.size());
    
    out.pack(ptr);
    
    size_t count = boost::asio::write( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( count != sz ) {
        LOG_ERR << "TcpConnection << HostInfo: Failed to write buffer.";
    }
    
    return conn;
    
}
TcpConnection& redux::network::operator>>(TcpConnection& conn, Peer::HostInfo& in) {
    
    size_t sz = sizeof(size_t)+1;
    auto buf = sharedArray<char>(sz);

    size_t count = boost::asio::read( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( count != sz ) {
        LOG_ERR << "TcpConnection >> HostInfo: Failed to read buffer.";
    }
    
    int one = 1;
    char* ptr = buf.get();
    bool swap_endian = *((char*)&one) != *ptr++;
    
    unpack(ptr,sz,swap_endian);
    
    buf = sharedArray<char>(sz);

    count = boost::asio::read( conn.socket(), boost::asio::buffer(buf.get(), sz) );
    if( count != sz ) {
        LOG_ERR << "TcpConnection >> HostInfo: Failed to read buffer.";
    }
    
    in.unpack(buf.get(), swap_endian);

    return conn;
    
}
