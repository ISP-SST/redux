USE_REDUX(redux)
USE_EXTERNAL(idl)


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
    
    set(RDX_SRC "${REDUX_DIR}/src/lib")
    set(CMAKE_DEBUG_POSTFIX "")   # no debug-suffix for dlm libraries
    
    REDUX_SETDIRS()
    
    add_library(ana SHARED ana.cpp
                              ${RDX_SRC}/file/anacompress.cpp
                              ${RDX_SRC}/file/anadecompress.cpp
                              ${RDX_SRC}/file/fileana.cpp
                             # ${RDX_SRC}/util/stringutil.cpp
               )
    SET(ANA_DLM_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libana.dlm")
    configure_file("libana.dlm.in" "${ANA_DLM_FILE_LOCATION}")
    
    add_library(momfbd SHARED momfbd.cpp
                           ${RDX_SRC}/file/filemomfbd.cpp
                           ${RDX_SRC}/math/functions.cpp
                           ${RDX_SRC}/util/stringutil.cpp
               )
    SET(MOMFBD_DLM_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libmomfbd.dlm")
    configure_file("libmomfbd.dlm.in" "${MOMFBD_DLM_FILE_LOCATION}")
    
    
    add_library(rdx SHARED rdx.cpp
                           #${RDX_SRC}/util/stringutil.cpp
               )
    SET(RDX_DLM_FILE_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/librdx.dlm")
    configure_file("librdx.dlm.in" "${RDX_DLM_FILE_LOCATION}")
    
    
    install(TARGETS ana momfbd rdx DESTINATION ${IDL_DLM_DIR} )
    install(FILES ${ANA_DLM_FILE_LOCATION} ${MOMFBD_DLM_FILE_LOCATION} ${RDX_DLM_FILE_LOCATION} DESTINATION ${IDL_DLM_DIR} )
   
else()
     message(WARNING "IDL headers not found, can not build the DLMs")
endif()