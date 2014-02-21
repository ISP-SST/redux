#ifndef REDUX_WORK_HPP
#define REDUX_WORK_HPP

#include "redux/network/peer.hpp"

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
        virtual char* pack( char* ) const;
        virtual const char* unpack( const char*, bool swap_endian=false  );
        size_t id;
        uint8_t step, nRetries;
    };


    struct WorkInProgress {
        typedef std::shared_ptr<WorkInProgress> Ptr;
        WorkInProgress(network::Peer::Ptr p=0);
        size_t size( bool includeJob ) const;
        char* pack( char*, bool includeJob) const;
        const char* unpack( const char*, bool includeJob, bool swap_endian=false );
        std::string print(void);
        std::shared_ptr<Job> job;
        std::vector<Part::Ptr> parts;
        network::Peer::Ptr peer;            // Basically only used to separate remote/local jobs.
        size_t nTotal, nCompleted;
    };



    /*! @} */

}  // redux


#endif // REDUX_WORK_HPP
