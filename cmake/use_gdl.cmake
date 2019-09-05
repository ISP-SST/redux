#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "GDL" )
#set( EXT_DEBUG "1" )
#set( EXT_COMPONENTS idl )

if( IS_DIRECTORY "$ENV{GDL_DIR}" )   # always use the environment GDL directory if it exists.
    set( GDL_DIR "$ENV{GDL_DIR}" )
elseif( IS_DIRECTORY "$ENV{GDL_PROJECT_DIR}" )
    set( GDL_DIR "$ENV{GDL_PROJECT_DIR}" )
endif()
#message(STATUS "GDL_DIR = ${GDL_DIR}  ENV{GDL_DIR} = $ENV{GDL_DIR}  ENV{GDL_PROJECT_DIR} = $ENV{GDL_PROJECT_DIR}")
set( GDL_DIR "/home.local/tomas/Projects/gdl/" )
set( EXT_HINT "${GDL_DIR}"
              "${GDL_DIR}/src/"
              "${GDL_DIR}/external/"
              "${GDL_DIR}/include/"
)

set( EXT_HEADER_FILE "envt.hpp" )


# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")

#set( GDL_RDX 0 CACHE INTERNAL "" FORCE )
#mark_as_advanced( GDL_RDX )

appendPaths()
