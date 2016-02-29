#
# Collect code for the main library and create build target
#

USE_REDUX(redux)    # includes listed thirdparty dependencies and sets up paths etc.
#USE_EXTERNAL(fits ncdf)

set(REDUX_HDR_BASE "${REDUX_DIR}/include/redux")
add_definitions(-DBOOST_LOG_DYN_LINK)

file(GLOB REDUX_CPP "*.cpp")
file(GLOB REDUX_HPP "${REDUX_HDR_BASE}/*.hpp")

file(GLOB REDUX_ALGO_CPP "algo/*.cpp")
file(GLOB REDUX_ALGO_HPP "${REDUX_HDR_BASE}/algo/*.hpp")

file(GLOB REDUX_FILE_CPP "file/*.cpp")
file(GLOB REDUX_FILE_HPP "${REDUX_HDR_BASE}/file/*.hpp")

file(GLOB REDUX_IMAGE_CPP "image/*.cpp")
file(GLOB REDUX_IMAGE_HPP "${REDUX_HDR_BASE}/image/*.hpp")

file(GLOB REDUX_MATH_CPP "math/*.cpp")
file(GLOB REDUX_MATH_HPP "${REDUX_HDR_BASE}/math/*.hpp")

file(GLOB REDUX_MOMFBD_CPP "momfbd/*.cpp")
file(GLOB REDUX_MOMFBD_HPP "${REDUX_HDR_BASE}/momfbd/*.hpp")

file(GLOB REDUX_NETWORK_CPP "network/*.cpp")
file(GLOB REDUX_NETWORK_HPP "${REDUX_HDR_BASE}/network/*.hpp")

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

if (OpenCV_FOUND)
    option(REDUX_WITH_OPENCV "Build with OpenCV support" ON)
    if (REDUX_WITH_OPENCV)
        message(STATUS "Building with OpenCV support")
        add_definitions(-DREDUX_WITH_OPENCV)
        file(GLOB REDUX_UTIL_OPENCV_CPP "util/opencv/*.cpp")
    endif()
endif()

if( NOT REDUX_WITH_OPENCV )
    message(WARNING "Building without OpenCV support, some functionality might not be compiled.")
endif()


##### revison info ######
set(REDUX_REVISION_SRC_FILE "${CMAKE_CURRENT_BINARY_DIR}/revision.cpp")
SET_SOURCE_FILES_PROPERTIES(${REDUX_REVISION_SRC_FILE} GENERATED)
CHECK_REVISION(LIBREDUX "${CMAKE_CURRENT_LIST_DIR}/revision.cpp.template" "${REDUX_REVISION_SRC_FILE}")

REDUX_ADD_LIBRARY(redux ${REDUX_REVISION_SRC_FILE}
                    ${REDUX_CPP} ${REDUX_HPP} 
                    ${REDUX_ALGO_CPP} ${REDUX_ALGO_HPP}
                    ${REDUX_FILE_CPP} ${REDUX_FILE_HPP}
                    ${REDUX_IMAGE_CPP} ${REDUX_IMAGE_HPP}
                    ${REDUX_MATH_CPP} ${REDUX_MATH_HPP}
                    ${REDUX_MOMFBD_CPP} ${REDUX_MOMFBD_HPP}
                    ${REDUX_NETWORK_CPP} ${REDUX_NETWORK_HPP}
                    ${REDUX_SST_CPP} ${REDUX_SST_HPP}
                    ${REDUX_UTIL_CPP} ${REDUX_UTIL_OPENCV_CPP} ${REDUX_UTIL_HPP}
)

if(REDUX_AUTO_REVISION AND TARGET GIT_CHECK_LIBREDUX)
    add_dependencies(redux GIT_CHECK_LIBREDUX)
endif()


if (REDUX_SKIP_GUI)
    message(STATUS "REDUX_SKIP_GUI is set: not building reduxgui.")
    return()
endif()

USE_REDUX(reduxgui)

# Check that we found Qt4
if (NOT QT4_FOUND)
    message(WARNING "QT4 not found: can't build reduxgui!")
    return()
endif()

file(GLOB REDUX_GUI_CPP "gui/*.cpp")
file(GLOB REDUX_GUI_HPP "${REDUX_HDR_BASE}/gui/*.hpp")

set(CMAKE_AUTOMOC ON)

REDUX_ADD_LIBRARY(reduxgui #${REDUX_REVISION_SRC_FILE}
                    ${REDUX_GUI_CPP} ${REDUX_GUI_HPP} 
                    ${REDUX_GUI_MOC_SRC}
)


