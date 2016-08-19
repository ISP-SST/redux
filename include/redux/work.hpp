#ifndef REDUX_WORK_HPP
#define REDUX_WORK_HPP

#include "redux/network/tcpconnection.hpp"
#include "redux/util/cache.hpp"

#include <cstdint>
#include <memory>
#include <vector>

namespace redux {

    /*! @ingroup redux
     *  @{
     */

    class Job;

    struct Part : public redux::util::CacheItem {
        typedef std::shared_ptr<Part> Ptr;
        Part();
        //Part(const Part&) = delete;
        virtual ~Part();
        virtual inline uint64_t size( void ) const {
            static uint64_t fixed_sz = sizeof(id)+sizeof(step)+sizeof(nRetries)+sizeof(partType);
            return fixed_sz;
        }
        virtual uint64_t pack( char* ) const;
        virtual uint64_t unpack( const char*, bool swap_endian=false );
        size_t csize(void) const { return size(); };
        uint64_t cpack(char* p) const { return pack(p); };
        uint64_t cunpack(const char* p, bool e) { return unpack(p,e); };
        bool operator==(const Part& rhs) const { return (id == rhs.id); }
        uint64_t id;
        boost::posix_time::ptime partStarted;
        uint16_t step;
        uint8_t nRetries, partType;
    };


    struct WorkInProgress {
        typedef std::shared_ptr<WorkInProgress> Ptr;
        WorkInProgress(void);
        WorkInProgress(const WorkInProgress&);
        ~WorkInProgress();
        uint64_t size(void) const;
        uint64_t pack(char*) const;
        uint64_t unpack(const char*, bool swap_endian=false);
        void resetParts(void);
        uint64_t workSize(void);
        uint64_t packWork(char*) const;
        uint64_t unpackWork(const char*, bool swap_endian=false);
        void returnResults(void);
        std::string print(void);
        std::shared_ptr<Job> job, previousJob;
        std::vector<Part::Ptr> parts;
        boost::posix_time::ptime workStarted;
        bool isRemote;
        bool hasResults;
        uint16_t nParts,nCompleted;
    };



    /*! @} */

}  // redux


#endif // REDUX_WORK_HPP
