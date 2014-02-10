# Global configuration for the gcc compiler
#

message(STATUS "Loading gcc-specific configuration. (${CMAKE_CURRENT_LIST_FILE})")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -fstrict-aliasing -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -std=c++0x -pedantic -Wall -fstrict-aliasing -O3 -pg -D_DEBUG")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -w -O3")
