cmake_minimum_required(VERSION  2.8)

IF(POLICY CMP0048)  # The ``project()`` command manages VERSION variables. cf. cmake --help-policy CMP0048
    CMAKE_POLICY(SET CMP0048 OLD)   # We manage versions manually
ENDIF()
IF(POLICY CMP0026)  # Disallow use of the LOCATION target property. cf. cmake --help-policy CMP0026
    CMAKE_POLICY(SET CMP0026 OLD)   # Some macros use this, should be updated to use the "$<TARGET_FILE:tgt>" syntax instead
ENDIF()

project(REDUX)

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
add_subdirectory(testsuite EXCLUDE_FROM_ALL)

# Installation (distribution) settings
install(DIRECTORY include/redux
        DESTINATION ${REDUX_INSTALL_DIR}/include/
        FILES_MATCHING PATTERN "*.hpp"
)

