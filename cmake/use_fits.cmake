#
# Set input data for FindExternal.cmake
#

set(EXT_NAME "Fits")


if(WIN32)
    set(EXT_HINT "${THIRDPARTY_DIR}/vendor/cfitsio/3.360" )
    set( EXT_LIB_SUFFIXES "" "-3" )
endif()

set(EXT_HEADER_FILE fitsio.h)
set(EXT_VERSION_FILE fitsio.h)
set(EXT_LIB_DEBUG_SUFFIX d)
set(EXT_COMPONENTS cfitsio)

set(EXT_MAJOR_REGEXP "CFITSIO_MAJOR")
set(EXT_MINOR_REGEXP "CFITSIO_MINOR")
#set(EXT_PATCH_REGEXP "CV_SUBMINOR_VERSION")


# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")


appendPaths()
