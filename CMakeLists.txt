cmake_minimum_required(VERSION  2.8)

project(REDUX)

if(NOT RDX_PROJECT)
    # Not included as a submodule to another project
    set(RDX_PROJECT "REDUX")
endif()

set(RDX_DIR ${CMAKE_CURRENT_LIST_DIR} CACHE PATH "Path to source tree.")
mark_as_advanced(RDX_DIR)

include(${RDX_DIR}/cmake/config.cmake)
include(${RDX_DIR}/cmake/doxygen.cmake)

include_directories(${RDX_DIR}/include)

add_subdirectory(src)
add_subdirectory(testsuite EXCLUDE_FROM_ALL)

