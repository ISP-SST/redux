
find_package(Doxygen)
if(DOXYGEN_FOUND)
 
    set( DOX_SRCS "${REDUX_DIR}/doc ${REDUX_DIR}/include/redux ${REDUX_DIR}/lib")

    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY )
    add_custom_target( doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMENT "Generate documentation with Doxygen" VERBATIM )

    if(NOT REDUX_NO_AUTO_REVISION AND TARGET GIT_CHECK_REDUX)
        add_dependencies(doc GIT_CHECK_REDUX)
    endif()

endif()
