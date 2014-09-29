
##### revison info ######
set(REDUX_REVISION_SRC_FILE "${CMAKE_CURRENT_BINARY_DIR}/revision.cpp")
SET_SOURCE_FILES_PROPERTIES(${REDUX_REVISION_SRC_FILE} GENERATED)
CHECK_REVISION(LIBREDUX "${CMAKE_CURRENT_LIST_DIR}/revision.cpp.template" "${REDUX_REVISION_SRC_FILE}")

add_subdirectory(lib)
#add_subdirectory(bin)
add_subdirectory(dlm)

