# Configuration specific for the Visual Studio compiler
#
# Author: Tomas Hillberg

message(STATUS "Loading configuration for MS Visual Studio. (${CMAKE_CURRENT_LIST_FILE})")


# Make Visual Studio stop complaining about using the unsecure (but non-standard)
# versions of e.g. printf, etc.
add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)
add_definitions(-D_SCL_SECURE_NO_WARNINGS)
add_definitions(-D_CRT_NONSTDC_NO_WARNINGS)

add_definitions(-DUSE_MSVC_THREAD)

set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /NODEFAULTLIB:msvcrt.lib")

# Turn of loading of the program only at its preferred base address. This is required for building
# dynamically loaded libraries and to enable profiling.
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FIXED:NO")


# Enable debugging of release builds
if (CMAKE_CL_64)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi /fp:fast /O2 /Oi /Ot")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /debug")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi /fp:fast /arch:SSE2 /O2 /Oi /Ot")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /debug")
endif()

