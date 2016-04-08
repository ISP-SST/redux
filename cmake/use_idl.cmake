#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "IDL" )
set( EXT_REQUIRED_VERSION "7.1" )
#set( EXT_COMPONENTS idl )

set( EXT_HINT "/usr/local/itt/idl/external/include/"
              "/opt/itt/idl/external/include/"
              "${THIRDPARTY_DIR}/vendor/idl/${USE_VERSION}"
              "${REDUX_DIR}/include/IDL"
)

set( EXT_HEADER_FILE "idl_export.h" )
set( EXT_VERSION_FILE "idl_export.h" )
set( EXT_MAJOR_REGEXP "VERSION_MAJOR" )
set( EXT_MINOR_REGEXP "VERSION_MINOR" )
set( EXT_PATCH_REGEXP "VERSION_SUB" )


# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")


appendPaths()
