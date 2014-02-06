#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "Boost" )
#set( EXT_DEBUG "1" )

if (MSVC11)
    set (EXT_HINT "${THIRDPARTY_DIR}/boost/1.55.0")
    #set( LIB_SUBDIR "/vc110" )
    set( EXT_LIB_SUFFIXES "" "-vc110-mt-1_55" )
    set( EXT_LIB_DEBUG_SUFFIX "-vc110-mt-gd-1_55" )
else()
    #set (EXT_HINT "${THIRDPARTY_DIR}/vendor/boost/1.46.1")
    #set( EXT_LIB_SUFFIXES "" "-vc80-mt-1_46_1" )
    #set( EXT_LIB_DEBUG_SUFFIX "-vc80-mt-gd-1_46_1" )
    #set (EXT_HINT "${THIRDPARTY_DIR}/vendor/boost/1.51.0")
    set (EXT_HINT "/usr/local/")
    set( LIB_SUBDIR "/vc80" )
    set( EXT_LIB_SUFFIXES "" "-mt" )
    set( EXT_LIB_DEBUG_SUFFIX "-mt-d" )
endif()
# if(CMAKE_CL_64)
    # set(LIB_ARCH "64" )
# endif()

set( EXT_LIB_PREFIX boost_ )

set( EXT_LIBPATH_SUFFIXES "lib${LIB_ARCH}${LIB_SUBDIR}" )
set( EXT_COMPONENTS date_time filesystem locale program_options serialization log
                         system thread regex timer unit_test_framework )

set( EXT_HEADER_FILE boost/version.hpp )
set( EXT_VERSION_FILE boost/version.hpp )
set( EXT_MAJOR_REGEXP "define[ \t]+BOOST_VERSION[ \t]+" )
set( EXT_MINOR_REGEXP "Disabled-NotAvailable" )
set( EXT_PATCH_REGEXP "Disabled-NotAvailable" )

# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")

# From boost/version.hpp:
# Each component uses two digits, so:
#   major = BOOST_VERSION / 100000
#   minor = BOOST_VERSION / 100 % 1000
#   patch = BOOST_VERSION % 100

if(Boost_VERSION VERSION_GREATER 100)
    MATH(EXPR lb_major "${Boost_VERSION} / 100000")
    MATH(EXPR lb_minor "${Boost_VERSION} / 100 % 1000")
    MATH(EXPR lb_patch "${Boost_VERSION} % 100")
    set(Boost_VERSION "${lb_major}.${lb_minor}.${lb_patch}" CACHE STRING "Version" FORCE)
endif()


appendPaths()
