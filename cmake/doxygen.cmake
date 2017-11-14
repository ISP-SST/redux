
find_package(Doxygen)
if(DOXYGEN_FOUND)
 
    set( DOX_SRCS "${RDX_DIR}/doc ${RDX_DIR}/include/redux")

    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY )
    add_custom_target( doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMENT "Generate documentation with Doxygen" VERBATIM )

    if(NOT RDX_NO_AUTO_REVISION AND TARGET GIT_CHECK_LIBREDUX)
        add_dependencies(doc GIT_CHECK_LIBREDUX)
    endif()

endif()
