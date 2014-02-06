# Windows specific configuration
#
# Author: Tomas Hillberg

message(STATUS "Loading configuration for Windows. (${CMAKE_CURRENT_LIST_FILE})")
# define WIN64 for x86-64 builds (_WIN64 is defined internaly by the compiler)
if(CMAKE_CL_64)
    add_definitions( -DWIN64 )
endif(CMAKE_CL_64)

# Reduce the amount of defines and symbols brought in by including the windows.h header
add_definitions(-DWIN32_LEAN_AND_MEAN)

# Tell windows.h not to define the min/max macros
add_definitions(-DNOMINMAX)

