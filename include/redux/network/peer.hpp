#ifndef REDUX_NETWORK_PEER_HPP
#define REDUX_NETWORK_PEER_HPP

#include "redux/network/tcpconnection.hpp"

#include <string>

#include <boost/date_time/posix_time/ptime.hpp>

namespace redux {

    namespace network {

        struct Peer {

            typedef std::shared_ptr<Peer> ptr;
            enum State : uint8_t { PEER_OFFLINE = 0, PEER_IDLE, PEER_ACTIVE, PEER_ERROR };
            const std::string StateNames[4] = { "offline", "idle", "active", "fail" };
            
            bool operator>( const Peer& rhs ) const;
            bool operator<( const Peer& rhs ) const;
            bool operator==( const Peer& rhs ) const;
            bool operator!=( const Peer& rhs ) const;

            Peer& operator=( const Peer& rhs );

            template<class Archive> void serialize( Archive& ar, const unsigned int version );

            struct HostInfo {       // "static" information, only exchanged on connect.
                char littleEndian;
                int reduxVersion;
                pid_t pid;
                boost::posix_time::ptime startedAt;
                std::string name, os, arch;
                HostInfo( void );
                void print(void);
                size_t size(void) const;
                char* pack( char* ) const;
                const char* unpack( const char*, bool );
                bool operator==(const HostInfo&);
            } host;

            struct PeerStatus {     // volatile info, refreshed every now and then.
                size_t currentJob;
                char nThreads;
                State state;
                float loadAvg;
                float progress;
                PeerStatus( void );
                size_t size(void) const;
                char* pack( char* ) const;
                const char* unpack( const char*, bool );
                template<class Archive> void serialize( Archive& ar, const unsigned int version );
            } stat;
            
            Peer(void) {};
            Peer(const HostInfo&, TcpConnection::ptr c=nullptr, size_t i=0);
            
            size_t size(void) const;
            char* pack( char* ) const;
            const char* unpack( const char*, bool );
            static std::string printHeader(void);
            std::string print(void);
            
            size_t id;
            boost::posix_time::ptime lastSeen;
            TcpConnection::ptr conn;
            
        };

        TcpConnection& operator<<(TcpConnection&, const Peer::HostInfo&);
        TcpConnection& operator>>(TcpConnection&, Peer::HostInfo&);


    }   // network

}   // redux

#endif  // REDUX_NETWORK_PEER_HPP
