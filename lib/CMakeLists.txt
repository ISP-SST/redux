#
# Collect code for the main library and create build target
#

USE_REDUX(redux)    # includes listed thirdparty dependencies and sets up paths etc.
#USE_EXTERNAL(fits ncdf)

set(REDUX_HDR_BASE "${REDUX_DIR}/include/redux")
add_definitions(-DBOOST_LOG_DYN_LINK)

##### revison info ######
set(REDUX_REVISION_SRC_FILE "${CMAKE_CURRENT_BINARY_DIR}/revision.cpp")
SET_SOURCE_FILES_PROPERTIES(${REDUX_REVISION_SRC_FILE} GENERATED)
CHECK_REVISION(REDUX "${CMAKE_CURRENT_LIST_DIR}/revision.cpp.template" "${REDUX_REVISION_SRC_FILE}")

file(GLOB REDUX_CPP "*.cpp")
file(GLOB REDUX_HPP "${REDUX_HDR_BASE}/*.hpp")

file(GLOB REDUX_MOMFBD_CPP "momfbd/*.cpp")
file(GLOB REDUX_MOMFBD_HPP "${REDUX_HDR_BASE}/momfbd/*.hpp")

file(GLOB REDUX_NETWORK_CPP "network/*.cpp")
file(GLOB REDUX_NETWORK_HPP "${REDUX_HDR_BASE}/network/*.hpp")

file(GLOB REDUX_FILE_CPP "file/*.cpp")
file(GLOB REDUX_FILE_HPP "${REDUX_HDR_BASE}/file/*.hpp")

file(GLOB REDUX_UTIL_CPP "util/*.cpp")
file(GLOB REDUX_UTIL_HPP "${REDUX_HDR_BASE}/util/*.hpp")

if (Fits_FOUND)
    option(REDUX_WITH_FITS "Build with FITS support" ON)
    if (REDUX_WITH_FITS)
        message(STATUS "Building with FITS support")
        add_definitions(-DREDUX_WITH_FITS)
    endif()
endif()

if (NetCDF_FOUND)
    option(REDUX_WITH_NCDF "Build with NetCDF support" ON)
    if (REDUX_WITH_NCDF)
        message(STATUS "Building with NetCDF support")
        add_definitions(-DREDUX_WITH_NCDF)
    endif()
endif()

REDUX_ADD_LIBRARY(redux ${REDUX_REVISION_SRC_FILE}
                    ${REDUX_CPP} ${REDUX_HPP} 
                    ${REDUX_MOMFBD_CPP} ${REDUX_MOMFBD_HPP}
                    ${REDUX_NETWORK_CPP} ${REDUX_NETWORK_HPP}
                    ${REDUX_FILE_CPP} ${REDUX_FILE_HPP}
                    ${REDUX_UTIL_CPP} ${REDUX_UTIL_HPP}
)

if(REDUX_AUTO_REVISION AND TARGET GIT_CHECK_REDUX)
    add_dependencies(redux GIT_CHECK_REDUX)
endif()
