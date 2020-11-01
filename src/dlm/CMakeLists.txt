#USE_REDUX(redux)
USE_EXTERNAL( boost fftw3 fits gsl idl opencv threads )


if (DEFINED IDL_FOUND)

    if( IDL_RDX )
        message( STATUS "Unable to find an IDL installation, using header from IDL 7.1." )
    endif()

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
    
    if( RDX_WITH_FFTW3 )
        add_definitions(-DRDX_WITH_FFTW3)
    endif()

    if( RDX_WITH_FITS )
        add_definitions(-DRDX_WITH_FITS)
    endif()

    if( RDX_WITH_GSL )
        add_definitions(-DRDX_WITH_GSL)
    endif()

    if( RDX_WITH_OPENCV )
        add_definitions(-DRDX_WITH_OPENCV)
    endif()
    
    add_definitions(-DBOOST_ALL_DYN_LINK)
    add_definitions(-DRDX_IDL_HDR_VER=${IDL_VERSION})

    set(CMAKE_DEBUG_POSTFIX "")   # no debug-suffix for dlm libraries
    
    RDX_SETDIRS()     # this will populate RDX_CURRENT_LIBRARIES with the necessary dependencies
    
    set( RDX_LIB_BASE "${RDX_DIR}/src/lib/" )
    
    set( RDX_ANA_DEPS
         "${RDX_LIB_BASE}/file/fileana.cpp"
         "${RDX_LIB_BASE}/file/anacompress.cpp"
         "${RDX_LIB_BASE}/file/anadecompress.cpp"
    )
    set( RDX_MOMFBD_DEPS
         "${RDX_LIB_BASE}/file/filefits.cpp"
         "${RDX_LIB_BASE}/util/stringutil.cpp"
         "${RDX_LIB_BASE}/file/filemomfbd.cpp"
         "${RDX_LIB_BASE}/math/functions.cpp"
    )
    set( RDX_RDX_DEPS filetools.cpp idlutil.cpp imgtools.cpp modes.cpp string.cpp
         ${RDX_ANA_DEPS} ${RDX_MOMFBD_DEPS}
         "${RDX_LIB_BASE}/file/fileio.cpp"
         "${RDX_LIB_BASE}/image/descatter.cpp"
         "${RDX_LIB_BASE}/image/grid.cpp"
         "${RDX_LIB_BASE}/image/pupil.cpp"
         "${RDX_LIB_BASE}/image/utils.cpp"
         "${RDX_LIB_BASE}/image/zernike.cpp"
         "${RDX_LIB_BASE}/math/linalg.cpp"
         "${RDX_LIB_BASE}/momfbd/modes.cpp"
         "${RDX_LIB_BASE}/util/arraystats.cpp"
         "${RDX_LIB_BASE}/util/cache.cpp"
         "${RDX_LIB_BASE}/util/progresswatch.cpp"
         "${RDX_LIB_BASE}/util/ricecompress.cpp"
         "${RDX_LIB_BASE}/util/stopwatch.cpp"
         "${RDX_LIB_BASE}/util/stringutil.cpp"
         "${RDX_LIB_BASE}/util/astro_comp.cpp"
         "${RDX_LIB_BASE}/util/solar_system_geo.cpp"
         "${RDX_LIB_BASE}/version.cpp"
         "${CMAKE_CURRENT_BINARY_DIR}/../lib/revision.cpp"
    )
    if( RDX_WITH_FFTW3 )
        set( RDX_RDX_DEPS ${RDX_RDX_DEPS} fft.cpp "${RDX_LIB_BASE}/image/fouriertransform.cpp")
    endif()
    if( RDX_WITH_FITS )
        set( RDX_RDX_DEPS ${RDX_RDX_DEPS} "${RDX_LIB_BASE}/file/filefits.cpp")
    endif()
    if( RDX_WITH_OPENCV )
        file( GLOB RDX_UTIL_OPENCV_CPP "${RDX_LIB_BASE}/util/opencv/*.cpp")
        set( RDX_RDX_DEPS ${RDX_RDX_DEPS} cvutil.cpp ${RDX_UTIL_OPENCV_CPP} )
    endif()
    if( RDX_TRACE )
        set(RDX_ANA_DEPS "${RDX_ANA_DEPS}" "${RDX_LIB_BASE}/util/trace.cpp")
        set(RDX_MOMFBD_DEPS "${RDX_MOMFBD_DEPS}" "${RDX_LIB_BASE}/util/trace.cpp")
        set(RDX_RDX_DEPS "${RDX_RDX_DEPS}" "${RDX_LIB_BASE}/util/trace.cpp")
    endif()

    add_library(ana SHARED ana.cpp ${RDX_ANA_DEPS} )
    set(ANA_DLM_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libana.dlm")
    configure_file("libana.dlm.in" "${ANA_DLM_FILE_LOCATION}")
    target_link_libraries (ana ${RDX_CURRENT_LIBRARIES})
    
    add_library(momfbd SHARED momfbd.cpp ${RDX_MOMFBD_DEPS} )
    set(MOMFBD_DLM_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libmomfbd.dlm")
    configure_file("libmomfbd.dlm.in" "${MOMFBD_DLM_FILE_LOCATION}")
    target_link_libraries (momfbd ${RDX_CURRENT_LIBRARIES})
    
    add_library(gridmatch SHARED gridmatch.cpp )
    set(GRIDMATCH_DLM_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libgridmatch.dlm")
    configure_file("libgridmatch.dlm.in" "${GRIDMATCH_DLM_FILE_LOCATION}")
    
    install(TARGETS ana momfbd gridmatch DESTINATION ${IDL_DLM_DIR} )
    install(FILES ${ANA_DLM_FILE_LOCATION} ${MOMFBD_DLM_FILE_LOCATION} ${GRIDMATCH_DLM_FILE_LOCATION} DESTINATION ${IDL_DLM_DIR} )

    
    if( Boost_VERSION VERSION_LESS "1.41" OR NOT RDX_WITH_FFTW3 OR NOT RDX_WITH_GSL OR NOT Threads_FOUND )
        message(STATUS "The rdx DLM will not be built (requires Boost>=1.41, FFTW3, GSL and Threads).")
        return()
    endif()
    
    set_source_files_properties( "${CMAKE_CURRENT_BINARY_DIR}/../lib/revision.cpp" GENERATED )
    
    add_library(rdx SHARED rdx.cpp convertlog.cpp ${RDX_RDX_DEPS} )
    set(RDX_DLM_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/librdx.dlm")
    configure_file("librdx.dlm.in" "${RDX_DLM_FILE_LOCATION}")
    target_link_libraries (rdx ${RDX_CURRENT_LIBRARIES})
    
    if(RDX_AUTO_REVISION AND TARGET GIT_CHECK_LIBRDX)
        add_dependencies(rdx GIT_CHECK_LIBRDX)
    endif()
    
    install(TARGETS rdx DESTINATION ${IDL_DLM_DIR} )
    install(FILES ${RDX_DLM_FILE_LOCATION} DESTINATION ${IDL_DLM_DIR} )
   
else()
     message(WARNING "IDL headers not found, can not build the DLMs")
endif()
