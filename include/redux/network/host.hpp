#ifndef REDUX_NETWORK_HOST_HPP
#define REDUX_NETWORK_HOST_HPP

#include "redux/network/tcpconnection.hpp"
#include "redux/util/stringutil.hpp"

#include <atomic>
#include <memory>
#include <set>
#include <string>

#ifndef Q_MOC_RUN
#include <boost/date_time/posix_time/ptime.hpp>
#endif

namespace redux {

    namespace network {

        struct Host {

            typedef std::shared_ptr<Host> Ptr;
            struct Compare {
                bool operator()( const Ptr &a, const Ptr &b ) const {
                    if( !a || !b ) return a < b;
                    return ( (*a) < (*b) ); }
            };
            typedef std::set<Ptr, Compare> Set;
             
            enum State : uint8_t { ST_OFFLINE=0, ST_LIMBO, ST_IDLE, ST_ACTIVE, ST_ERROR };
            enum Type : uint8_t { TP_WORKER=1, TP_MASTER, TP_UI=4 };
            static const std::string StateNames[4];
            static const std::string TypeNames[8];
            
            struct HostInfo {       // "static" information, only exchanged on connect.
                char littleEndian;
                uint64_t reduxVersion;
                pid_t pid;
                uint8_t peerType;
                uint16_t nCores;
                uint16_t connectPort;
                boost::posix_time::ptime startedAt;
                std::string name, os, arch, connectName, user;
                HostInfo( std::string username="" );
                uint64_t size(void) const;
                uint64_t pack( char* ) const;
                uint64_t unpack( const char*, bool );
                bool operator==(const HostInfo&) const;
            } info;

            struct HostStatus {     // volatile info, refreshed every now and then.
                uint64_t currentJob;
                uint16_t nThreads, maxThreads;
                uint16_t listenPort;
                State state;
                float load[2];
                float progress;
                std::string statusString;
                boost::posix_time::ptime lastSeen;
                boost::posix_time::ptime lastActive;
                HostStatus( void );
                uint64_t size(void) const;
                uint64_t pack( char* ) const;
                uint64_t unpack( const char*, bool );
            } status;
            
            Host(void);
            Host(const HostInfo&, uint64_t i=0);
            ~Host();

            uint64_t size(void) const;
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );
            
            void touch(void);
            void idle(void);
            void active(void);
            void limbo(void);
            
            static Host& myInfo(void);
            static std::string printHeader( int verbosity=0 );
            std::string print( int verbosity=0 );
            
            //Host& operator=( const Host& rhs );

            bool operator>( const Host& rhs ) const;
            bool operator<( const Host& rhs ) const;
            bool operator==( const Host& rhs ) const;
            bool operator!=( const Host& rhs ) const { return !(*this == rhs); };

            uint64_t id;
            std::atomic<int> nConnections;
            
        };

        TcpConnection& operator<<(TcpConnection&, const Host::HostInfo&);
        TcpConnection& operator>>(TcpConnection&, Host::HostInfo&);


    }   // network

}   // redux

#endif  // REDUX_NETWORK_HOST_HPP
