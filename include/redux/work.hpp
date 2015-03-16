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
        //Part(const Part&) = delete;
        virtual ~Part();
        virtual uint64_t size( void ) const;
        virtual uint64_t pack( char* ) const;
        virtual uint64_t unpack( const char*, bool swap_endian=false );
        bool operator==(const Part& rhs) const { return (id == rhs.id); }
        uint64_t id;
        uint8_t step, nRetries, partType;
    };


    struct WorkInProgress {
        typedef std::shared_ptr<WorkInProgress> Ptr;
        WorkInProgress(network::TcpConnection::Ptr c=0);
        WorkInProgress(const WorkInProgress&);
        ~WorkInProgress();
        uint64_t size(void) const;
        uint64_t pack(char*) const;
        uint64_t unpack(const char*, bool swap_endian=false);
        uint64_t workSize(void);
        uint64_t packWork(char*) const;
        uint64_t unpackWork(const char*, bool swap_endian=false);
        std::string print(void);
        std::shared_ptr<Job> job, previousJob;
        std::vector<Part::Ptr> parts;
        network::TcpConnection::Ptr connection;            // Basically only used to separate remote/local jobs.
        uint16_t nParts,nCompleted;
    };



    /*! @} */

}  // redux


#endif // REDUX_WORK_HPP
