USE_REDUX("")
USE_EXTERNAL( openmpi gsl threads )

if( DEFINED OpenMPI_FOUND )
    option(REDUX_WITH_OPENMPI "Build with OpenMPI_FOUND support" OFF)
    if( REDUX_WITH_OPENMPI )
        message(STATUS "Building with OpenMPI support")
        add_definitions(-DREDUX_WITH_OPENMPI)
    endif()
else()
    message(STATUS "OpenMPI not found. Try your systems equivalent of \"apt-get install libopenmpi-dev\"" )
endif()

if( NOT REDUX_WITH_OPENMPI OR NOT REDUX_WITH_GSL OR NOT Threads_FOUND )
    message(STATUS "The speckle calibration tool will not be built (requires OpenMPI, GSL and Threads).")
    return()
endif()

SET(MPI_CFLAGS "-DUSE_STDARG -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 -DUSE_STDARG=1 -DMALLOC_RET_VOID=1")

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CFLAGS})

REDUX_ADD_EXECUTABLE( stf_new stf.cpp zernike.cpp wam.cpp )
REDUX_ADD_EXECUTABLE( stf_newest stf_newest.cpp mpiwrapper.cpp zernike.cpp wam.cpp )


