# Get sha1 hash for current branch
execute_process( COMMAND "${GIT_EXECUTABLE}" describe "--always" "--abbrev=12" "--long" "--dirty=\ (Modified)"
                WORKING_DIRECTORY "${GIT_DIR}"
                OUTPUT_VARIABLE REVISION_ID
                OUTPUT_STRIP_TRAILING_WHITESPACE
               )

# Get timestamp of last commit on current branch
execute_process(COMMAND "${GIT_EXECUTABLE}" show "-s" "--format=%ci"
                WORKING_DIRECTORY "${GIT_DIR}"
                OUTPUT_VARIABLE COMMIT_TIME_tmp
                OUTPUT_STRIP_TRAILING_WHITESPACE
               )
if(COMMIT_TIME_tmp)
    string(REGEX REPLACE "([0-9]+)-([0-9]+)-([0-9]+) ([0-9:]+).*" "\\1/\\2/\\3 \\4" REVISION_COMMIT_TIME "${COMMIT_TIME_tmp}")
else()
    set(REVISION_COMMIT_TIME "unknown")
endif()

# Get commit message
execute_process(COMMAND "${GIT_EXECUTABLE}" show "-s" "--format=%s"
                WORKING_DIRECTORY "${GIT_DIR}"
                OUTPUT_VARIABLE REVISION_COMMIT_COMMENT
                OUTPUT_STRIP_TRAILING_WHITESPACE
               )

set(REVISION_DESCRIPTION "${REVISION_ID} - ${REVISION_COMMIT_COMMENT}")

# Ugly way to get current time, hopefully cmake will provide such functionality soon...
if (WIN32)
    execute_process(COMMAND "wmic" os get LocalDateTime OUTPUT_VARIABLE BUILD_TIME_tmp OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX MATCH "([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])([0-9][0-9])([0-9][0-9])([0-9][0-9])" BUILD_TIME_tmp2 ${BUILD_TIME_tmp})
    string(REGEX REPLACE "([0-9][0-9][0-9][0-9])([0-9][0-9])([0-9][0-9])([0-9][0-9])([0-9][0-9])([0-9][0-9])" "\\1-\\2-\\3 \\4:\\5:\\6" REVISION_BUILD_TIME ${BUILD_TIME_tmp2})
elseif (UNIX)
        execute_process(COMMAND "date" "+%Y-%m-%d %k:%M:%S" OUTPUT_VARIABLE REVISION_BUILD_TIME OUTPUT_STRIP_TRAILING_WHITESPACE)
else (WIN32)
    message(WARNING "date not implemented")
    set(REVISION_BUILD_TIME "2000-01-01 00:00:00")
endif (WIN32)


string( REGEX MATCHALL "[0-9.-]+g" _VERSION_STR "${REVISION_ID}" )
string( REGEX MATCHALL "[0-9]+" _VERSION_PARTS "${_VERSION_STR}" )

list(LENGTH _VERSION_PARTS _NPARTS)
set(_CMAKE_ARGS --no-warn-unused-cli)
if(_NPARTS GREATER 0)
    list(GET _VERSION_PARTS 0 ${REVNAME}_VERSION_MAJOR)
    if(_NPARTS GREATER 1)
        list(GET _VERSION_PARTS 1 ${REVNAME}_VERSION_MINOR)
    endif()
    if(_NPARTS GREATER 2)
        list(GET _VERSION_PARTS 2 ${REVNAME}_VERSION_PATCH)
    endif()
    if(_NPARTS GREATER 3)
        list(GET _VERSION_PARTS 3 ${REVNAME}_VERSION_COMMIT)
    endif()
endif()

# Set to 0 if undefined
if(NOT ${REVNAME}_VERSION_MAJOR)
    set(${REVNAME}_VERSION_MAJOR 0)
endif()
if(NOT ${REVNAME}_VERSION_MINOR)
    set(${REVNAME}_VERSION_MINOR 0)
endif()
if(NOT ${REVNAME}_VERSION_PATCH)
    set(${REVNAME}_VERSION_PATCH 0)
endif()
if(NOT ${REVNAME}_VERSION_COMMIT)
    set(${REVNAME}_VERSION_COMMIT 0)
endif()

list(APPEND _CMAKE_ARGS -D${REVNAME}_VERSION_MAJOR:STRING=${${REVNAME}_VERSION_MAJOR})
list(APPEND _CMAKE_ARGS -D${REVNAME}_VERSION_MINOR:STRING=${${REVNAME}_VERSION_MINOR})
list(APPEND _CMAKE_ARGS -D${REVNAME}_VERSION_PATCH:STRING=${${REVNAME}_VERSION_PATCH})
list(APPEND _CMAKE_ARGS -D${REVNAME}_VERSION_COMMIT:STRING=${${REVNAME}_VERSION_COMMIT})

# Update values in the cmake cache
execute_process(COMMAND ${CMAKE_COMMAND} ${_CMAKE_ARGS} ${BUILD_DIR})

# Generate revision-file from specified template
configure_file("${REVISION_TEMPLATE_LOCATION}" "${REVISION_FILE_LOCATION}")
