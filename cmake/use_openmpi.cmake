#
# Set input data for FindExternal.cmake
#

set( EXT_NAME "OpenMPI" )

set( EXT_COMPONENTS mpi_cxx mpi dl hwloc)
set( EXT_INCLUDE_SUFFIXES mpi openmpi )

set( EXT_HEADER_FILE "mpi.h" )
set( EXT_VERSION_FILE "mpi.h" )

set( EXT_MAJOR_REGEXP "OMPI_MAJOR_VERSION" )
set( EXT_MINOR_REGEXP "OMPI_MINOR_VERSION" )
set( EXT_PATCH_REGEXP "OMPI_RELEASE_VERSION" )


# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")


appendPaths()
