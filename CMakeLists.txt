project(REDUX)
cmake_minimum_required(VERSION  2.8)

if(NOT REDUX_PROJECT)
    # Not included as a submodule to another project
    set(REDUX_PROJECT "REDUX")
endif()

set(REDUX_DIR ${CMAKE_CURRENT_LIST_DIR} CACHE PATH "Path to source tree.")
mark_as_advanced(REDUX_DIR)

include(${REDUX_DIR}/cmake/config.cmake)
include(${REDUX_DIR}/cmake/doxygen.cmake)

include_directories(${REDUX_DIR}/include)

add_subdirectory(src)
add_subdirectory(testsuite)

# Installation (distribution) settings
install(DIRECTORY include/redux
        DESTINATION ${REDUX_INSTALL_DIR}/include/
        FILES_MATCHING PATTERN "*.hpp"
)

