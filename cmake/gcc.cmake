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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -pedantic -Wall -fstrict-aliasing")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -ggdb3 -fvar-tracking -fvar-tracking-assignments")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -pg -DDEBUG_ -fno-omit-frame-pointer -ggdb3 -fvar-tracking -fvar-tracking-assignments")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=c99 -w -O3")
