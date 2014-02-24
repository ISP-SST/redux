#ifndef REDUX_NETWORK_PEER_HPP
#define REDUX_NETWORK_PEER_HPP

#include "redux/network/tcpconnection.hpp"

#include <memory>
#include <set>
#include <string>

#include <boost/date_time/posix_time/ptime.hpp>

namespace redux {

    namespace network {

        struct Peer {

            typedef std::shared_ptr<Peer> Ptr;
        
            enum State : uint8_t { PEER_OFFLINE=0, PEER_IDLE, PEER_ACTIVE, PEER_ERROR };
            enum Type : uint8_t { PEER_WORKER=1, PEER_MASTER, PEER_UTIL=4 };
            static const std::string StateNames[4];
            static const std::string TypeNames[8];
            
            struct HostInfo {       // "static" information, only exchanged on connect.
                char littleEndian;
                int reduxVersion;
                pid_t pid;
                uint8_t peerType;
                uint16_t nCores;
                boost::posix_time::ptime startedAt;
                std::string name, os, arch;
                HostInfo( void );
                size_t size(void) const;
                char* pack( char* ) const;
                const char* unpack( const char*, bool );
                bool operator==(const HostInfo&) const;
            } host;

            struct PeerStatus {     // volatile info, refreshed every now and then.
                size_t currentJob;
                uint8_t nThreads, maxThreads;
                State state;
                float loadAvg;
                float progress;
                PeerStatus( void );
                size_t size(void) const;
                char* pack( char* ) const;
                const char* unpack( const char*, bool );
            } stat;
            
            Peer(void);
            Peer(const HostInfo&, TcpConnection::Ptr c=nullptr, size_t i=0);

            size_t size(void) const;
            char* pack( char* ) const;
            const char* unpack( const char*, bool );
            static std::string printHeader(void);
            std::string print(void);
            
            std::shared_ptr<char> receiveBlock( size_t& blockSize, bool& swap_endian );
            
            Peer& operator=( const Peer& rhs );

            bool operator>( const Peer& rhs ) const;
            bool operator<( const Peer& rhs ) const;
            bool operator==( const Peer& rhs ) const;
            bool operator!=( const Peer& rhs ) const { return !(*this == rhs); };

            size_t id;
            boost::posix_time::ptime lastSeen;
            TcpConnection::Ptr conn;
            
        };

        TcpConnection& operator<<(TcpConnection&, const Peer::HostInfo&);
        TcpConnection& operator>>(TcpConnection&, Peer::HostInfo&);


    }   // network

}   // redux

#endif  // REDUX_NETWORK_PEER_HPP
