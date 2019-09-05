#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "EIGEN" )
#set( EXT_DEBUG "1" )

set( EXT_HEADER_FILE "Eigen/Core" )
#set( EXT_VERSION_FILE "Eigen/Core" )
#set( EXT_MAJOR_REGEXP "MAJOR_VERSION" )
#set( EXT_MINOR_REGEXP "MINOR_VERSION" )
#set( EXT_PATCH_REGEXP "SUBMINOR_VERSION" )
set( EXT_HINT "/usr/include/eigen"
              "/usr/include/eigen3"
)

# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")


appendPaths()
