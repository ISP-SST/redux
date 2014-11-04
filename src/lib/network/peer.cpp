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


uint64_t Peer::pack( char* ptr ) const {
    
    using redux::util::pack;
    uint64_t count = pack(ptr,id);
    count += pack(ptr+count,to_time_t( lastSeen ));
    count += host.pack(ptr+count);
    count += stat.pack(ptr+count);
    
    return count;
}

uint64_t Peer::unpack( const char* ptr, bool swap_endian ) {
    
    using redux::util::unpack;
    uint64_t count = unpack(ptr,id,swap_endian);
    time_t timestamp;
    count += unpack(ptr+count,timestamp,swap_endian);
    count += host.unpack(ptr+count,swap_endian);
    count += stat.unpack(ptr+count,swap_endian);
    
    lastSeen = boost::posix_time::from_time_t( timestamp );
    
    return count;
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


uint64_t Peer::HostInfo::pack( char* ptr ) const {
    
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

uint64_t Peer::HostInfo::unpack( const char* ptr, bool swap_endian ) {
    
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


bool Peer::HostInfo::operator==(const HostInfo& rhs) const {
    return (name == rhs.name && pid == rhs.pid);
}


size_t Peer::PeerStatus::size(void) const {
    size_t sz = sizeof(size_t) + sizeof(State) + 2*sizeof(float) + 2;
    return sz;
}


uint64_t Peer::PeerStatus::pack( char* ptr ) const {
    
    using redux::util::pack;
    
    uint64_t count = pack(ptr,nThreads);
    count += pack(ptr+count,maxThreads);
    count += pack(ptr+count,state);
    count += pack(ptr+count,currentJob);
    count += pack(ptr+count,loadAvg);
    count += pack(ptr+count,progress);
    
    return count;
    
}


uint64_t Peer::PeerStatus::unpack( const char* ptr, bool swap_endian ) {
    
    using redux::util::unpack;
    uint64_t count = unpack(ptr,nThreads);
    count += unpack(ptr+count,maxThreads);
    count += unpack(ptr+count,state, swap_endian);
    count += unpack(ptr+count,currentJob, swap_endian);
    count += unpack(ptr+count,loadAvg, swap_endian);
    count += unpack(ptr+count,progress, swap_endian);
    
    return count;
    
}


TcpConnection& redux::network::operator<<(TcpConnection& conn, const Peer::HostInfo& out) {
    
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


TcpConnection& redux::network::operator>>(TcpConnection& conn, Peer::HostInfo& in) {
    
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
