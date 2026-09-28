#
# Collect code for the main library and create build target
#

USE_REDUX(redux)    # includes listed thirdparty dependencies and sets up paths etc.

set(RDX_HDR_BASE "${RDX_DIR}/include/redux")
add_definitions(-DBOOST_LOG_DYN_LINK)

file(GLOB RDX_CPP "*.cpp")
file(GLOB RDX_HPP "${RDX_HDR_BASE}/*.hpp")

#file(GLOB RDX_ALGO_CPP "algo/*.cpp")
#file(GLOB RDX_ALGO_HPP "${RDX_HDR_BASE}/algo/*.hpp")

file(GLOB RDX_FILE_CPP "file/*.cpp")
file(GLOB RDX_FILE_HPP "${RDX_HDR_BASE}/file/*.hpp")

file(GLOB RDX_IMAGE_CPP "image/*.cpp")
file(GLOB RDX_IMAGE_HPP "${RDX_HDR_BASE}/image/*.hpp")

file(GLOB RDX_LOGGING_CPP "logging/*.cpp")
file(GLOB RDX_LOGGING_HPP "${RDX_HDR_BASE}/logging/*.hpp")

file(GLOB RDX_MATH_CPP "math/*.cpp")
file(GLOB RDX_MATH_HPP "${RDX_HDR_BASE}/math/*.hpp")

file(GLOB RDX_MOMFBD_CPP "momfbd/*.cpp")
file(GLOB RDX_MOMFBD_HPP "${RDX_HDR_BASE}/momfbd/*.hpp")

file(GLOB RDX_NETWORK_CPP "network/*.cpp")
file(GLOB RDX_NETWORK_HPP "${RDX_HDR_BASE}/network/*.hpp")

file(GLOB RDX_UTIL_CPP "util/*.cpp")
file(GLOB RDX_UTIL_GSL_CPP "util/gsl/*.cpp")
file(GLOB RDX_UTIL_HPP "${RDX_HDR_BASE}/util/*.hpp")

set( RDX_MISSING_REQUIRED_DEPS "" )

if( NOT DEFINED Boost_FOUND )
    message(STATUS "Boost not found. Try your systems equivalent of \"apt-get install libboost-all-dev\"" )
    list(APPEND RDX_MISSING_REQUIRED_DEPS "Boost (>= 1.66)")
elseif( Boost_VERSION VERSION_LESS "1.41" )
    message(STATUS "Boost ${Boost_VERSION} was found, but redux is only tested with boost >= 1.41.")
    list(APPEND RDX_MISSING_REQUIRED_DEPS "Boost (found ${Boost_VERSION}, need >= 1.41)")
endif()

if( DEFINED Fits_FOUND )
    option(RDX_WITH_FITS "Build with FITS support" ON)
    if( RDX_WITH_FITS )
        message(STATUS "Building with FITS support")
        add_definitions(-DRDX_WITH_FITS)
    endif()
else()
    unset(RDX_WITH_FITS CACHE)
    message(STATUS "CFITSIO not found. FITS support is optional; try your systems equivalent of \"apt-get install libcfitsio-dev\" to enable it." )
endif()

if( DEFINED FFTW3_FOUND )
    message(STATUS "Building with FFTW3 support")
    set(RDX_WITH_FFTW3 ON CACHE INTERNAL "")
    add_definitions(-DRDX_WITH_FFTW3)
else()
    unset(RDX_WITH_FFTW3 CACHE)
    message(STATUS "FFTW3 not found. Try your systems equivalent of \"apt-get install libfftw3-dev\"" )
    list(APPEND RDX_MISSING_REQUIRED_DEPS "FFTW3")
endif()

if( DEFINED GSL_FOUND )
    message(STATUS "Building with GSL support")
    set(RDX_WITH_GSL ON CACHE INTERNAL "")
    add_definitions(-DRDX_WITH_GSL)
else()
    unset(RDX_WITH_GSL CACHE)
    message(STATUS "GSL not found (redux requires GSL >= 1.16). Try your systems equivalent of \"apt-get install libgsl-dev\"" )
    list(APPEND RDX_MISSING_REQUIRED_DEPS "GSL (>= 1.16)")
endif()

if( DEFINED OpenCV_FOUND )
    option(RDX_WITH_OPENCV "Build with OpenCV support" ON)
    if( RDX_WITH_OPENCV )
        message(STATUS "Building with OpenCV support")
        add_definitions(-DRDX_WITH_OPENCV)
        file(GLOB RDX_UTIL_OPENCV_CPP "util/opencv/*.cpp")
    else()
        message(FATAL_ERROR " redux has to be compiled with OpenCV support." )
    endif()
else()
    unset(RDX_WITH_OPENCV CACHE)
    message(FATAL_ERROR "OpenCV not found. Try your systems equivalent of \"apt-get install libopencv-*-dev\"" )
endif()

if( DEFINED ZLIB_FOUND )
    message(STATUS "Building with ZLIB support")
    set(RDX_WITH_ZLIB ON CACHE INTERNAL "")
    add_definitions(-DRDX_WITH_ZLIB)
else()
    unset(RDX_WITH_ZLIB CACHE)
    message(STATUS "ZLIB not found. Try your systems equivalent of \"apt-get install zlib1g-dev\"" )
    list(APPEND RDX_MISSING_REQUIRED_DEPS "ZLIB")
endif()

if( NOT Threads_FOUND )
    list(APPEND RDX_MISSING_REQUIRED_DEPS "Threads")
endif()

if( RDX_MISSING_REQUIRED_DEPS )
    string(REPLACE ";" ", " RDX_MISSING_REQUIRED_DEPS_STR "${RDX_MISSING_REQUIRED_DEPS}")
    message(FATAL_ERROR "redux cannot be built: missing or unsuitable required dependencies: ${RDX_MISSING_REQUIRED_DEPS_STR}. See the messages above for install hints.")
endif()

##### revison info ######
set(RDX_REVISION_SRC_FILE "${CMAKE_CURRENT_BINARY_DIR}/revision.cpp")
set_source_files_properties( ${RDX_REVISION_SRC_FILE} GENERATED )
CHECK_REVISION(LIBRDX "${CMAKE_CURRENT_LIST_DIR}/revision.cpp.template" "${RDX_REVISION_SRC_FILE}")

RDX_ADD_LIBRARY(redux ${RDX_REVISION_SRC_FILE}
                    ${RDX_CPP} ${RDX_HPP} 
                    ${RDX_ALGO_CPP} ${RDX_ALGO_HPP}
                    ${RDX_FILE_CPP} ${RDX_FILE_HPP}
                    ${RDX_IMAGE_CPP} ${RDX_IMAGE_HPP}
                    ${RDX_LOGGING_CPP} ${RDX_LOGGING_HPP}
                    ${RDX_MATH_CPP} ${RDX_MATH_HPP}
                    ${RDX_MOMFBD_CPP} ${RDX_MOMFBD_HPP}
                    ${RDX_NETWORK_CPP} ${RDX_NETWORK_HPP}
                    ${RDX_SST_CPP} ${RDX_SST_HPP}
                    ${RDX_UTIL_CPP} ${RDX_UTIL_OPENCV_CPP} ${RDX_UTIL_HPP}
                    ${RDX_UTIL_CPP} ${RDX_UTIL_GSL_CPP} ${RDX_UTIL_OPENCV_CPP} ${RDX_UTIL_HPP}
)

if(RDX_AUTO_REVISION AND TARGET GIT_CHECK_LIBRDX)
    add_dependencies(redux GIT_CHECK_LIBRDX)
endif()

# Installation (distribution) settings
install(DIRECTORY ${RDX_DIR}/include/redux
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/
        FILES_MATCHING PATTERN "*.hpp"
)


if( RDX_SKIP_GUI )
    message(STATUS "RDX_SKIP_GUI is set: not building reduxgui.")
    return()
endif()

USE_REDUX(reduxgui)

# Check that we found Qt4
if (NOT QT4_FOUND)
    message(WARNING "QT4 not found: can't build reduxgui!")
    return()
endif()

file(GLOB RDX_GUI_CPP "gui/*.cpp")
file(GLOB RDX_GUI_HPP "${RDX_HDR_BASE}/gui/*.hpp")

set(CMAKE_AUTOMOC ON)

RDX_ADD_LIBRARY(reduxgui #${RDX_REVISION_SRC_FILE}
                    ${RDX_GUI_CPP} ${RDX_GUI_HPP} 
                    ${RDX_GUI_MOC_SRC}
)


