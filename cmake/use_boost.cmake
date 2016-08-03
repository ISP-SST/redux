#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "Boost" )
#set( EXT_DEBUG "1" )

set( EXT_LIB_PREFIX boost_ )
set( EXT_LIB_SUFFIXES "-mt" "" )
set( EXT_LIB_DEBUG_SUFFIX "-mt-d" )

set( EXT_LIBPATH_SUFFIXES "lib${LIB_ARCH}${LIB_SUBDIR}" )
set( EXT_COMPONENTS date_time filesystem program_options serialization system thread regex )

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
