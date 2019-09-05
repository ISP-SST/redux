# Global configuration for the gcc compiler
#

message(STATUS "Loading gcc-specific configuration. (${CMAKE_CURRENT_LIST_FILE})")

execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
if (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
    message(STATUS "C++11 activated.")
    add_definitions("-std=gnu++11")
elseif(GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
    message(STATUS "C++0x activated. If you get any errors, update to a compiler which fully supports C++11")
    add_definitions("-std=gnu++0x")
else ()
    message(FATAL_ERROR "C++11 needed. Therefore a gcc compiler with a version higher than 4.3 is needed.")   
endif()

if( RDX_TRACE )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -rdynamic -DRDX_TRACE")
    message(STATUS "Enabling trace (with symbol demangling)")
endif()

if( RDX_TRACE_ARRAY )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRDX_TRACE_ARRAY")
    message(STATUS "  Tracing Array objects")
endif()

if( RDX_TRACE_CACHE )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRDX_TRACE_CACHE")
    message(STATUS "  Tracing CacheItem objects")
endif()

if( RDX_TRACE_FILE )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRDX_TRACE_FILE")
    message(STATUS "  Tracing FileMeta/CachedFile objects")
endif()

if( RDX_TRACE_MEM )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRDX_TRACE_MEM")
    message(STATUS "  Tracing shared_ptr (de-)allocations")
endif()

if( RDX_TRACE_NET )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRDX_TRACE_NET")
    message(STATUS "  Tracing TcpConnection objects")
endif()

if( RDX_TRACE_JOB )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRDX_TRACE_JOB")
    message(STATUS "  Tracing Job objects")
endif()

if( RDX_TRACE_PARTS )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRDX_TRACE_PARTS")
    message(STATUS "  Tracing Part objects")
endif()
if( RDX_TRACE_PROC )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRDX_TRACE_PROC")
    message(STATUS "  Tracing Processing objects")
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -fopenmp -pedantic -Wall -fstrict-aliasing")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -ggdb3 -fvar-tracking -fvar-tracking-assignments")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -pg -fopenmp -DDEBUG_ -fno-omit-frame-pointer -ggdb3 -fvar-tracking -fvar-tracking-assignments")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -fopenmp -std=c99 -w -O3")
