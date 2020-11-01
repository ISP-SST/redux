#USE_REDUX(redux)
#USE_EXTERNAL( boost fftw3 fits gsl gdl opencv threads )
USE_EXTERNAL( boost eigen fftw3 fits gdl )

message(STATUS "Building in GDL ...")

if ( DEFINED GDL_FOUND )

message(STATUS "GDL is defined.")
    if( NOT IS_DIRECTORY "${GDL_DLM_DIR}" AND IS_DIRECTORY "$ENV{GDL_DLM_DIR}" )
        set( GDL_DLM_DIR "$ENV{GDL_DLM_DIR}/" )
    endif()
    
    if( NOT GDL_DLM_DIR )
        # if it is still blank, use ~/lib/gdl/ as default
        set(GDL_DLM_DIR "$ENV{HOME}/lib/gdl")
    endif()
    
    # Convert the path to cmake internal format (avoids '\'-related problems)
    file( TO_CMAKE_PATH "${GDL_DLM_DIR}/" GDL_DLM_DIR )
    
    # Store it in the cache (so it shows up in the cmake-gui)
    set(GDL_DLM_DIR "${GDL_DLM_DIR}" CACHE PATH "")
    message(STATUS "GDL DLM will be installed into: \"${GDL_DLM_DIR}\"")
    
    add_definitions(-DBOOST_ALL_DYN_LINK)

    set(CMAKE_DEBUG_POSTFIX "")   # no debug-suffix for dlm libraries
    
    RDX_SETDIRS()     # this will populate RDX_CURRENT_LIBRARIES with the necessary dependencies
    
    if (DEFINED EIGEN_FOUND)
        add_definitions(-DUSE_EIGEN)    # needed so that the DLM alignment matches GDL (if built with eigen)
    endif()
    
    set( RDX_LIB_BASE "${RDX_DIR}/src/lib/" )
    
    set( RDX_RDX_DEPS 
#         "${RDX_LIB_BASE}/file/fileio.cpp"
         "${RDX_LIB_BASE}/version.cpp"
         "${CMAKE_CURRENT_BINARY_DIR}/../lib/revision.cpp"
    )
    add_library( rdx_gdl SHARED rdx.cpp two.cpp ${RDX_RDX_DEPS} )
    install(TARGETS rdx_gdl DESTINATION ${GDL_DLM_DIR} )
    get_property(fancy_lib_location TARGET rdx_gdl PROPERTY LOCATION)
    get_property(fancy_lib_locations TARGET rdx_gdl PROPERTY NAME)
    #get_target_property(Shared PROPERTIES LIBRARY_OUTPUT_NAME demo.so.1.2.3)

    message (STATUS "fancy_lib_location == ${fancy_lib_location}")
    message (STATUS "fancy_lib_location == ${fancy_lib_locations}")
    set( RDX_GDL_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/load_rdx_gdl.pro" )
    configure_file( "load_rdx_gdl_tpl.pro" "${RDX_GDL_FILE_LOCATION}" )
    target_link_libraries( rdx_gdl ${RDX_CURRENT_LIBRARIES} )

    
#    if( RDX_AUTO_REVISION AND TARGET GIT_CHECK_LIBREDUX )
#        add_dependencies( rdx_gdl GIT_CHECK_LIBREDUX )
#    endif()
    
    install(FILES ${RDX_GDL_FILE_LOCATION} DESTINATION ${GDL_DLM_DIR} )
   
else()
     message(WARNING "GDL headers not found, can not build the DLMs")
endif()
