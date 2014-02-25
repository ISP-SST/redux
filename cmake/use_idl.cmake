#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "IDL" )
set( USE_VERSION "1.0" )
#set( EXT_COMPONENTS idl )

set( EXT_HINT "/usr/local/itt/idl/external/include/"
              "/opt/itt/idl/external/include/"
              "${THIRDPARTY_DIR}/vendor/idl/${USE_VERSION}"
              "${REDUX_DIR}/include/"
)

set( EXT_HEADER_FILE "idl_export.h" )
#set( EXT_VERSION_FILE "idl/idlv.h" )
#set( EXT_MAJOR_REGEXP "MAJOR_VERSION" )
#set( EXT_MINOR_REGEXP "MINOR_VERSION" )
#set( EXT_PATCH_REGEXP "SUBMINOR_VERSION" )


# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")


appendPaths()
