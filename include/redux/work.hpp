#ifndef REDUX_WORK_HPP
#define REDUX_WORK_HPP

#include "redux/network/tcpconnection.hpp"
#include "redux/util/cacheitem.hpp"

#ifdef RDX_TRACE_PARTS
#   include "redux/util/trace.hpp"
#endif

#include <cstdint>
#include <memory>
#include <vector>

namespace redux {

    /*! @ingroup redux
     *  @{
     */

    class Job;
    
    struct PackedData
#ifdef RDX_TRACE_PARTS
            : public redux::util::TraceObject<PackedData>
#endif
    {
        PackedData() : size(0), packedSize(0) {};
        void clear(void) { data.reset(); size = packedSize = 0; };
        std::shared_ptr<char> data;
        uint64_t size;
        uint64_t packedSize;
    };


    struct Part : public redux::util::CacheItem
#ifdef RDX_TRACE_PARTS
#ifndef RDX_TRACE_CACHE
            ,public redux::util::TraceObject<Part>
#endif
#endif
    {
        typedef std::shared_ptr<Part> Ptr;
        Part();
        //Part(const Part&) = delete;
        virtual ~Part();
        virtual inline uint64_t size( void ) const {
            static uint64_t fixed_sz = sizeof(id)+sizeof(step)+sizeof(nRetries)+sizeof(partType)
                                    +sizeof(nThreads)+sizeof(runtime_wall)+sizeof(runtime_cpu);
            return fixed_sz;
        }
        virtual uint64_t pack( char* ) const;
        virtual uint64_t unpack( const char*, bool swap_endian=false );
        virtual void load(void);
        virtual void unload(void);
        virtual void prePack( bool force=false ) {};
        size_t csize(void) const { return size(); };
        uint64_t cpack(char* p) const { return pack(p); };
        uint64_t cunpack(const char* p, bool e) { return unpack(p,e); };
        bool operator==(const Part& rhs) const { return (id == rhs.id); }
        bool operator<(const Part& rhs) const { return (id < rhs.id); }
        uint64_t id;
        boost::posix_time::ptime partStarted;
        uint16_t step;
        uint8_t nRetries, partType;
        uint16_t nThreads;
        float runtime_wall;
        float runtime_cpu;
        
        PackedData packed;
        
    };


    struct WorkInProgress : std::enable_shared_from_this<WorkInProgress>
#ifdef RDX_TRACE_PARTS
            ,public redux::util::TraceObject<WorkInProgress>
#endif
    {
        typedef std::shared_ptr<WorkInProgress> Ptr;
        WorkInProgress(void);
        WorkInProgress(const WorkInProgress&);
        ~WorkInProgress();
        uint64_t size(void) const;
        uint64_t pack(char*) const;
        uint64_t unpack(const char*, bool swap_endian=false);
        void reset(void);
        void resetParts(void);
        uint64_t workSize(void);
        uint64_t packWork(char*);
        uint64_t unpackWork(const char*, std::shared_ptr<Job>& tmpJob, bool swap_endian=false );
        void returnResults(void);
        bool operator<(const WorkInProgress& rhs) const;
        std::string print(void);
        std::weak_ptr<Job> job;
        uint32_t jobID;
        std::vector<Part::Ptr> parts;
        boost::posix_time::ptime workStarted;
        bool isRemote;
        bool hasResults;
        uint16_t nParts,nCompleted;
    };



    /*! @} */

}  // redux


#endif // REDUX_WORK_HPP
