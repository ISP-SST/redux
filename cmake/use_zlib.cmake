#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "ZLIB" )

set( EXT_COMPONENTS z )

set( EXT_HEADER_FILE "zlib.h" )
#set( EXT_VERSION_FILE "zlib.h" )
set( EXT_MAJOR_REGEXP "ZLIB_VER_MAJOR" )
set( EXT_MINOR_REGEXP "ZLIB_VER_MINOR" )
set( EXT_PATCH_REGEXP "ZLIB_VER_REVISION" )


# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")


appendPaths()
