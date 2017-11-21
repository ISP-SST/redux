#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "IDL" )
#set( EXT_DEBUG "1" )
#set( EXT_COMPONENTS idl )

if( IS_DIRECTORY "$ENV{IDL_DIR}" )   # always use the environment IDL directory if it exists.
    set( IDL_DIR "$ENV{IDL_DIR}" )
endif()

set( EXT_HINT ${IDL_DIR} "${IDL_DIR}/external/" "${IDL_DIR}/external/include/"
              "/usr/local/harris/idl/external/include/"
              "/opt/harris/idl/external/include/"
              "/usr/local/rsi/idl/external/include/"
              "/opt/rsi/external/include/"
              "/usr/local/exelis/idl/external/include/"
              "/opt/exelis/idl/external/include/"
              "/usr/local/itt/idl/external/include/"
              "/opt/itt/idl/external/include/"
)

set( EXT_HEADER_FILE "idl_export.h" )
set( EXT_VERSION_FILE "idl_export.h" )
set( EXT_MAJOR_REGEXP "VERSION_MAJOR" )
set( EXT_MINOR_REGEXP "VERSION_MINOR" )
set( EXT_PATCH_REGEXP "NotUsed" )


# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")

set( IDL_RDX 0 CACHE INTERNAL "" FORCE )
mark_as_advanced( IDL_RDX )
if( NOT DEFINED IDL_FOUND )     # no system installation found, use the header included with rdx
    set( EXT_HINT "${RDX_DIR}/include/IDL" )
    set( EXT_HEADER_FILE "idl_export.h" )
    set( EXT_VERSION_FILE "idl_export.h" )
    set( EXT_MAJOR_REGEXP "VERSION_MAJOR" )
    set( EXT_MINOR_REGEXP "VERSION_MINOR" )
    include( "${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake" )
    if( DEFINED IDL_FOUND )
        set( IDL_RDX 1 CACHE INTERNAL "" FORCE )
    endif()
endif()

if( DEFINED IDL_FOUND )
    message("Using IDL header file found in: ${IDL_DIR} (ver: ${IDL_VERSION})")
    message(" To use a different IDL version, set the IDL_DIR variable in your environment before running cmake.")
endif()


appendPaths()
