#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "FFTW3" )

set( EXT_COMPONENTS fftw3 fftw3_threads )

if(WIN32)
    set (EXT_HINT "${THIRDPARTY_DIR}/vendor/fftw3/3.3.3")
    set( EXT_LIB_SUFFIXES "" "-3" )
endif()


set( EXT_HEADER_FILE "fftw3.h" )
#set( EXT_VERSION_FILE "fftw3.h" )
#set( EXT_MAJOR_REGEXP "MAJOR_VERSION" )
#set( EXT_MINOR_REGEXP "MINOR_VERSION" )
#set( EXT_PATCH_REGEXP "SUBMINOR_VERSION" )


# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")


appendPaths()
