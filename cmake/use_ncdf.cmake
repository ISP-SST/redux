#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "NetCDF" )

set( USE_VERSION "3.1" )
set( EXT_COMPONENTS netcdf )

set( EXT_HEADER_FILE "netcdf.h" )
#set( EXT_VERSION_FILE "netcdf.h" )
#set( EXT_MAJOR_REGEXP "MAJOR" )
#set( EXT_MINOR_REGEXP "MINOR" )
#set( EXT_PATCH_REGEXP "SUBMINOR_VERSION" )


# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")


appendPaths()
