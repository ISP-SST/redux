#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "GSL" )

set( USE_VERSION "1.16" )
set( EXT_COMPONENTS gsl gslcblas)

set( EXT_HEADER_FILE "gsl/gsl_version.h" )
set( EXT_VERSION_FILE "gsl/gsl_version.h" )
set( EXT_MAJOR_REGEXP "GSL_MAJOR_VERSION" )
set( EXT_MINOR_REGEXP "GSL_MINOR_VERSION" )
#set( EXT_PATCH_REGEXP "SUBMINOR_VERSION" )

set( EXT_HELPTEXT "Try your systems equivalent of \"apt-get install libgsl0-dev\"" )

# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")


appendPaths()
