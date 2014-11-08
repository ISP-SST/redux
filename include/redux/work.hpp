#ifndef REDUX_WORK_HPP
#define REDUX_WORK_HPP

#include "redux/network/tcpconnection.hpp"

#include <cstdint>
#include <memory>
#include <vector>

namespace redux {

    /*! @ingroup redux
     *  @{
     */

    class Job;

    struct Part {
        typedef std::shared_ptr<Part> Ptr;
        Part();
        virtual size_t size( void ) const;
        virtual uint64_t pack( char* ) const;
        virtual uint64_t unpack( const char*, bool swap_endian=false  );
        uint64_t id;
        uint8_t step, nRetries;
    };


    struct WorkInProgress {
        typedef std::shared_ptr<WorkInProgress> Ptr;
        WorkInProgress(network::TcpConnection::Ptr c=0);
        size_t size( bool includeJob ) const;
        uint64_t pack( char*, bool includeJob) const;
        uint64_t unpack( const char*, bool swap_endian=false );
        std::string print(void);
        std::shared_ptr<Job> job;
        std::vector<Part::Ptr> parts;
        network::TcpConnection::Ptr connection;            // Basically only used to separate remote/local jobs.
        uint16_t nCompleted;
    };



    /*! @} */

}  // redux


#endif // REDUX_WORK_HPP
