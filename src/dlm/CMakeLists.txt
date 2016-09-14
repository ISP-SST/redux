#USE_REDUX(redux)
USE_EXTERNAL( boost fftw3 fits gsl idl opencv threads )


if (DEFINED IDL_FOUND)

    if( NOT IS_DIRECTORY "${IDL_DLM_DIR}" AND IS_DIRECTORY "$ENV{IDL_DLM_DIR}" )
        set( IDL_DLM_DIR "$ENV{IDL_DLM_DIR}/" )
    endif()
    
    if( NOT IDL_DLM_DIR )
        # if it is still blank, use ~/lib/dlm/ as default
        set(IDL_DLM_DIR "$ENV{HOME}/lib/dlm")
    endif()
    
    # Convert the path to cmake internal format (avoids '\'-related problems)
    file( TO_CMAKE_PATH "${IDL_DLM_DIR}/" IDL_DLM_DIR )
    
    # Store it in the cache (so it shows up in the cmake-gui)
    set(IDL_DLM_DIR "${IDL_DLM_DIR}" CACHE PATH "")
    message(STATUS "IDL DLM will be installed into: \"${IDL_DLM_DIR}\"")
    
    if( REDUX_WITH_FFTW3 )
        add_definitions(-DREDUX_WITH_FFTW3)
    endif()

    if( REDUX_WITH_FITS )
        add_definitions(-DREDUX_WITH_FITS)
    endif()

    if( REDUX_WITH_GSL )
        add_definitions(-DREDUX_WITH_GSL)
    endif()

    if( REDUX_WITH_OPENCV )
        add_definitions(-DREDUX_WITH_OPENCV)
    endif()
    
    add_definitions(-DBOOST_ALL_DYN_LINK)

    set(CMAKE_DEBUG_POSTFIX "")   # no debug-suffix for dlm libraries
    
    REDUX_SETDIRS()     # this will populate REDUX_CURRENT_LIBRARIES with the necessary dependencies
    
    set( REDUX_LIB_BASE "${REDUX_DIR}/src/lib/" )
    
    set( REDUX_ANA_DEPS
         "${REDUX_LIB_BASE}/file/fileana.cpp"
         "${REDUX_LIB_BASE}/file/anacompress.cpp"
         "${REDUX_LIB_BASE}/file/anadecompress.cpp"
    )
    add_library(ana SHARED ana.cpp ${REDUX_ANA_DEPS} )
    set(ANA_DLM_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libana.dlm")
    configure_file("libana.dlm.in" "${ANA_DLM_FILE_LOCATION}")
    target_link_libraries (ana ${REDUX_CURRENT_LIBRARIES})
    
    set( REDUX_MOMFBD_DEPS
         "${REDUX_LIB_BASE}/file/filemomfbd.cpp"
         "${REDUX_LIB_BASE}/math/functions.cpp"
    )
    add_library(momfbd SHARED momfbd.cpp ${REDUX_MOMFBD_DEPS} )
    set(MOMFBD_DLM_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libmomfbd.dlm")
    configure_file("libmomfbd.dlm.in" "${MOMFBD_DLM_FILE_LOCATION}")
    target_link_libraries (momfbd ${REDUX_CURRENT_LIBRARIES})
    
    add_library(gridmatch SHARED gridmatch.cpp )
    set(GRIDMATCH_DLM_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libgridmatch.dlm")
    configure_file("libgridmatch.dlm.in" "${GRIDMATCH_DLM_FILE_LOCATION}")
    
    install(TARGETS ana momfbd gridmatch DESTINATION ${IDL_DLM_DIR} )
    install(FILES ${ANA_DLM_FILE_LOCATION} ${MOMFBD_DLM_FILE_LOCATION} ${GRIDMATCH_DLM_FILE_LOCATION} DESTINATION ${IDL_DLM_DIR} )
        
    if( Boost_VERSION VERSION_LESS "1.41" OR NOT REDUX_WITH_FFTW3 OR NOT REDUX_WITH_GSL OR NOT Threads_FOUND )
        message(STATUS "The rdx DLM will not be built (requires Boost>=1.41, FFTW3, GSL and Threads).")
        return()
    endif()

        
    set( REDUX_RDX_DEPS idlutil.cpp imgtools.cpp modes.cpp string.cpp
         ${REDUX_ANA_DEPS} ${REDUX_MOMFBD_DEPS}
         "${REDUX_LIB_BASE}/file/fileio.cpp"
         "${REDUX_LIB_BASE}/image/descatter.cpp"
         "${REDUX_LIB_BASE}/image/grid.cpp"
         "${REDUX_LIB_BASE}/image/pupil.cpp"
         "${REDUX_LIB_BASE}/image/utils.cpp"
         "${REDUX_LIB_BASE}/image/zernike.cpp"
         "${REDUX_LIB_BASE}/math/linalg.cpp"
         "${REDUX_LIB_BASE}/momfbd/modes.cpp"
         "${REDUX_LIB_BASE}/util/arraystats.cpp"
         "${REDUX_LIB_BASE}/util/cache.cpp"
         "${REDUX_LIB_BASE}/util/progresswatch.cpp"
         "${REDUX_LIB_BASE}/util/stopwatch.cpp"
         "${REDUX_LIB_BASE}/util/stringutil.cpp"
         "${REDUX_LIB_BASE}/version.cpp"
         "${CMAKE_CURRENT_BINARY_DIR}/../lib/revision.cpp"
    )
    set_source_files_properties( "${CMAKE_CURRENT_BINARY_DIR}/../lib/revision.cpp" GENERATED )
    if( REDUX_WITH_FFTW3 )
        set( REDUX_RDX_DEPS ${REDUX_RDX_DEPS} fft.cpp "${REDUX_LIB_BASE}/image/fouriertransform.cpp")
    endif()
    if( REDUX_WITH_FITS )
        set( REDUX_RDX_DEPS ${REDUX_RDX_DEPS} "${REDUX_LIB_BASE}/file/filefits.cpp")
    endif()
    if( REDUX_WITH_OPENCV )
        file( GLOB REDUX_UTIL_OPENCV_CPP "${REDUX_LIB_BASE}/util/opencv/*.cpp")
        set( REDUX_RDX_DEPS ${REDUX_RDX_DEPS} cvutil.cpp ${REDUX_UTIL_OPENCV_CPP} )
    endif()
    
    add_library(rdx SHARED rdx.cpp ${REDUX_RDX_DEPS} )
    set(RDX_DLM_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/librdx.dlm")
    configure_file("librdx.dlm.in" "${RDX_DLM_FILE_LOCATION}")
    target_link_libraries (rdx ${REDUX_CURRENT_LIBRARIES})
    
    if(REDUX_AUTO_REVISION AND TARGET GIT_CHECK_LIBREDUX)
        add_dependencies(rdx GIT_CHECK_LIBREDUX)
    endif()
    
    install(TARGETS rdx DESTINATION ${IDL_DLM_DIR} )
    install(FILES ${RDX_DLM_FILE_LOCATION} DESTINATION ${IDL_DLM_DIR} )
   
else()
     message(WARNING "IDL headers not found, can not build the DLMs")
endif()
