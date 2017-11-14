# Redux project config stuff
#  
# Author: Tomas Hillberg

# simple include guard
if(DEFINED RDX_CONFIGURATION_LOADED)
    return()
endif()
set(RDX_CONFIGURATION_LOADED True)

# Enable empty ELSE and ENDIF's
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# Add a '-d' to all debug targets (lib's)
set(CMAKE_DEBUG_POSTFIX "-d")

# Include debug libraries in install phase.
set(CMAKE_INSTALL_DEBUG_LIBRARIES "True")

# Allow for optional "lib" prefix when searching for libraries
# This is just to combine the two "defaults" for windows and apple/linux to be
# able to streamline the "find" script a bit.
set(CMAKE_FIND_LIBRARY_PREFIXES "" lib)

include(TestBigEndian)
test_big_endian(RDX_BIG_ENDIAN)
if (RDX_BIG_ENDIAN)
    add_definitions(-DRDX_BYTE_ORDER=RDX_BIG_ENDIAN)
else()
    add_definitions(-DRDX_BYTE_ORDER=RDX_LITTLE_ENDIAN)
endif()

option( RDX_SKIP_GUI "Should the gui/QT parts be skipped?" ON )
option( RDX_AUTO_REVISION "Automatically run the revision script on each build ?" ON )

option( RDX_TRACE_ALL "Enable tracing all supported objects" OFF )
option( RDX_TRACE_ARRAY "Enable tracing of Array/ArrayStats objects" OFF )
option( RDX_TRACE_CACHE "Enable tracing of CacheItem objects" OFF )
option( RDX_TRACE_FILE "Enable tracing of FileMeta/CachedFile objects" OFF )
option( RDX_TRACE_NET "Enable tracing of TcpConnection objects" OFF )
option( RDX_TRACE_JOB "Enable tracing of Job objects" OFF )
option( RDX_TRACE_PARTS "Enable tracing of Part/WorkInProgress objects" OFF )

set(RDX_TRACE_ALL ${RDX_TRACE_ALL} CACHE BOOL "Enable tracing all supported objects" FORCE)

if( RDX_TRACE_ALL )
    set( RDX_TRACE_ARRAY ON  )
    set( RDX_TRACE_CACHE ON )
    set( RDX_TRACE_FILE  ON )
    set( RDX_TRACE_NET   ON )
    set( RDX_TRACE_JOB   ON )
    set( RDX_TRACE_PARTS ON )
endif()

if( RDX_TRACE_ARRAY OR RDX_TRACE_CACHE OR RDX_TRACE_FILE OR RDX_TRACE_NET OR RDX_TRACE_JOB OR RDX_TRACE_PARTS )
    set( RDX_TRACE True )
endif()

# Load platform specific configuration
if(WIN32)
    include(${CMAKE_CURRENT_LIST_DIR}/windows.cmake)
elseif(APPLE)
    include(${CMAKE_CURRENT_LIST_DIR}/mac.cmake)
elseif(UNIX)
    include(${CMAKE_CURRENT_LIST_DIR}/unix.cmake)
else()
    message(FATAL_ERROR "Unsupported platform")
endif()

# Load compiler specific configuration
if(MSVC)
    include(${CMAKE_CURRENT_LIST_DIR}/msvc.cmake)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    include(${CMAKE_CURRENT_LIST_DIR}/clang.cmake)
elseif(CMAKE_COMPILER_IS_GNUCC)
    include(${CMAKE_CURRENT_LIST_DIR}/gcc.cmake)
else()
    message(FATAL_ERROR "Unsupported compiler")
endif()


include(${CMAKE_CURRENT_LIST_DIR}/dependencies.cmake)
file(TO_CMAKE_PATH "${PROJECT_SOURCE_DIR}/cmake/dependencies.cmake" PROJ_DEPS)
if(EXISTS "${PROJ_DEPS}")
    include( "${PROJ_DEPS}" )
endif()

include(${CMAKE_CURRENT_LIST_DIR}/macros.cmake)
    
# The thirdparty directory contains all thirdparty libraries used by REDUX.
# Try to get a default from the environment variable with the same name.
if( NOT IS_DIRECTORY "${THIRDPARTY_DIR}" AND IS_DIRECTORY "$ENV{THIRDPARTY_DIR}" )
    set( THIRDPARTY_DIR "$ENV{THIRDPARTY_DIR}/" )
endif()

if( THIRDPARTY_DIR )
    # Convert the path to cmake internal format (avoids '\'-related problems)
    file( TO_CMAKE_PATH "${THIRDPARTY_DIR}/" THIRDPARTY_DIR )
endif()

if( EXISTS "${THIRDPARTY_DIR}" )
    message(STATUS "Using thirdparty directory: \"${THIRDPARTY_DIR}\"")
endif()
set(THIRDPARTY_DIR "${THIRDPARTY_DIR}" CACHE PATH "Path to thirdparty software.")

get_filename_component(RDX_DIR ${RDX_DIR} REALPATH)

# Try to find cppcheck
# Rationale: running cppcheck outside the project means a lot of includes (e.g. thirdparty)
# will not be found, this way the tests will use the "real" include enviroment of the compiler.
find_program(RDX_CPPCHECK_EXECUTABLE NAMES "cppcheck" "cppcheck.exe" PATHS "C:\\Program Files\\cppcheck" "C:\\Program Files (x86)\\cppcheck")
if(NOT EXISTS "${RDX_CPPCHECK_EXECUTABLE}")
    message("Install cppcheck to enable target-specific source-code checking.")
else()
    option(RDX_CPPCHECK_TARGETS "Should the targets for CppCheck be included?" OFF)
endif()

include_directories(${RDX_DIR}/include)

message(STATUS "Redux will be installed into: \"${CMAKE_INSTALL_PREFIX}\"")

# Setup a variable to indicate endian:ness of current system.
include(TestBigEndian)
TEST_BIG_ENDIAN(RDX_BIG_ENDIAN)

# Create a variable containing a list of the REDUX libraries. This is usefull for specifying
# REDUX linkage with the target_link_libraries in projects using REDUX libraries.
# The order of the libraries should be optimized for linking (i.e. dependent libraries appear
# before the libraries they depend on).
# Note: reduxgui is not included in this list (to avoid that compilers adds it as a dependency).
#   If reduxgui is needed, it should be added to the link_libraries() call, and it should go
#   before RDX_LIBLIST because reduxgui depends on pretty much everything else.
set(RDX_LIBLIST 
    redux                         # depends on: none
    reduxgui                      # depends on: none
    CACHE INTERNAL "")
   

# To build an application using REDUX as a binary distribution from outside the source tree, 
# has to be run to make the project aware of the redux libraries. The reason is that
# CMake won't be able to handle the -d debug suffix correctly otherwise.
#file(GLOB BUILDING_REDUX "${RDX_DIR}/CMakeLists.txt")
#if (NOT BUILDING_REDUX AND NOT ${PROJECT_NAME} STREQUAL "REDUX")
#    foreach(LIB ${RDX_LIBLIST})
#    IMPORT_RDX_LIB(${LIB})
#    endforeach()
#endif()
