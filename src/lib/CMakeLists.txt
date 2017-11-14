#
# Collect code for the main library and create build target
#

USE_REDUX(redux)    # includes listed thirdparty dependencies and sets up paths etc.
#USE_EXTERNAL(fits ncdf)

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

if( NOT DEFINED Boost_FOUND )
    message(STATUS "Boost not found. Try your systems equivalent of \"apt-get install libboost-all-dev\"" )
endif()

if( DEFINED Fits_FOUND )
    option(RDX_WITH_FITS "Build with FITS support" ON)
    if( RDX_WITH_FITS )
        message(STATUS "Building with FITS support")
        add_definitions(-DRDX_WITH_FITS)
    endif()
else()
    unset(RDX_WITH_FITS CACHE)
    message(STATUS "CFITSIO not found. Try your systems equivalent of \"apt-get install libcfitsio-dev\"" )
endif()

if( DEFINED NetCDF_FOUND )
    option(RDX_WITH_NCDF "Build with NetCDF support" ON)
    if( RDX_WITH_NCDF )
        message(STATUS "Building with NetCDF support")
        add_definitions(-DRDX_WITH_NCDF)
    endif()
endif()

if( DEFINED FFTW3_FOUND )
    option(RDX_WITH_FFTW3 "Build with FFTW3 support" ON)
    if( RDX_WITH_FFTW3 )
        message(STATUS "Building with FFTW3 support")
        add_definitions(-DRDX_WITH_FFTW3)
    endif()
else()
    unset(RDX_WITH_FFTW3 CACHE)
    message(STATUS "FFTW3 not found. Try your systems equivalent of \"apt-get install libfftw3-dev\"" )
endif()

if( DEFINED GSL_FOUND )
    option(RDX_WITH_GSL "Build with GSL support" ON)
    if( RDX_WITH_GSL )
        message(STATUS "Building with GSL support")
        add_definitions(-DRDX_WITH_GSL)
    endif()
else()
    unset(RDX_WITH_GSL CACHE)
    message(STATUS "GSL not found. Try your systems equivalent of \"apt-get install libgsl0-dev\"" )
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
    option(RDX_WITH_ZLIB "Build with ZLIB support" ON)
    if( RDX_WITH_ZLIB )
        message(STATUS "Building with ZLIB support")
        add_definitions(-DRDX_WITH_ZLIB)
    endif()
else()
    unset(RDX_WITH_GSL CACHE)
    message(STATUS "ZLIB not found. Try your systems equivalent of \"apt-get install zlib1g-dev\"" )
endif()

if( NOT RDX_WITH_OPENCV )
    message(STATUS "Building without OpenCV support, some functionality might not be compiled.")
endif()

if( NOT DEFINED Boost_FOUND OR NOT RDX_WITH_FFTW3 OR NOT RDX_WITH_GSL OR NOT RDX_WITH_ZLIB OR NOT Threads_FOUND )
    message(STATUS "The redux library will not be built (requires Boost, FFTW3, GSL, zlib and Threads).")
    return()
endif()

if( Boost_VERSION VERSION_LESS "1.41" )
    message(STATUS "The redux library will not be built (only tested with boost >= 1.41).")
    return()
endif()

##### revison info ######
set(RDX_REVISION_SRC_FILE "${CMAKE_CURRENT_BINARY_DIR}/revision.cpp")
set_source_files_properties( ${RDX_REVISION_SRC_FILE} GENERATED )
CHECK_REVISION(LIBREDUX "${CMAKE_CURRENT_LIST_DIR}/revision.cpp.template" "${RDX_REVISION_SRC_FILE}")

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

if(RDX_AUTO_REVISION AND TARGET GIT_CHECK_LIBREDUX)
    add_dependencies(redux GIT_CHECK_LIBREDUX)
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


