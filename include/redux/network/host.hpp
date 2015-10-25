#ifndef REDUX_NETWORK_HOST_HPP
#define REDUX_NETWORK_HOST_HPP

#include "redux/network/tcpconnection.hpp"

#include <memory>
#include <set>
#include <string>

#include <boost/date_time/posix_time/ptime.hpp>

namespace redux {

    namespace network {

        struct Host {

            typedef std::shared_ptr<Host> Ptr;
            struct Compare {
                bool operator()( const Ptr &a, const Ptr &b ) const { return ( (*a) < (*b) ); }
            };
            typedef std::set<Ptr, Compare> Set;
             
            enum State : uint8_t { ST_OFFLINE=0, ST_IDLE, ST_ACTIVE, ST_ERROR };
            enum Type : uint8_t { TP_WORKER=1, TP_MASTER, TP_UI=4 };
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
                uint64_t pack( char* ) const;
                uint64_t unpack( const char*, bool );
                bool operator==(const HostInfo&) const;
            } info;

            struct HostStatus {     // volatile info, refreshed every now and then.
                uint64_t currentJob;
                uint8_t nThreads, maxThreads;
                State state;
                float loadAvg;
                float progress;
                boost::posix_time::ptime lastSeen;
                HostStatus( void );
                size_t size(void) const;
                uint64_t pack( char* ) const;
                uint64_t unpack( const char*, bool );
            } status;
            
            Host(void);
            Host(const HostInfo&, uint64_t i=0);

            size_t size(void) const;
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );
            
            void touch(void);
            
            static std::string printHeader(void);
            std::string print(void);
            
            //Host& operator=( const Host& rhs );

            bool operator>( const Host& rhs ) const;
            bool operator<( const Host& rhs ) const;
            bool operator==( const Host& rhs ) const;
            bool operator!=( const Host& rhs ) const { return !(*this == rhs); };

            uint64_t id;
            int nConnections;
            
        };

        TcpConnection& operator<<(TcpConnection&, const Host::HostInfo&);
        TcpConnection& operator>>(TcpConnection&, Host::HostInfo&);


    }   // network

}   // redux

#endif  // REDUX_NETWORK_HOST_HPP
