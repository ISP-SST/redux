###############################################################################
#
# Input:
#   These input variables can be defined before invoking FindExternal.cmake
#       EXT_NAME                           Name of the package
#       EXT_REQUIRED_VERSION               Minimum supported version number
#       EXT_PATHS                          Additional paths to search for package
#                                          By default the system paths are used.
#       EXT_HINT                           Prefered path,
#                                          Will be searched before system paths.
#       EXT_HEADER_FILE                    Name of file located in the main
#                                          include directory
#       EXT_INCLUDE_SUFFIXES               Search for EXT_HEADER_FILE in these subdirs
#       EXT_VERSION_FILE                   Name of file containing the version
#       EXT_VERSION_REGEXP                 Regexp for filtering out the version string
#                                          Disabled by default
#       EXT_MAJOR_REGEXP                   Regexp for filtering out the major version number
#                                          Default: "_MAJOR"
#       EXT_MINOR_REGEXP                   Regexp for filtering out the minor version number (default: 
#                                          Default: "_MINOR"
#       EXT_PATCH_REGEXP                   Regexp for filtering out the patch number
#                                          Default: "(_PATCH|_BUILD|_SUBMINOR|_RELEASE|_BETA|_REVISION)"
#       EXT_LIBPATH_SUFFIXES               Search for libraries in these subdirs.
#       EXT_LIB_PREFIX                     String to prepend to all component names
#       EXT_LIB_SUFFIXES                   Also search for alternate filenames
#                                          e.g. appending version number to the library names
#       EXT_LIB_DEBUG_SUFFIX               Find debug libraries with this suffix to the names.
#       EXT_COMPONENTS                     Required libraries
#       EXT_DEBUG                          Enable debug printouts
#                                          Helps to see why a component is not located.
#
#
# Output:
#   The following are set when a package is located
#       ${EXT_NAME}_FOUND
#       ${EXT_NAME}_VERSION                (only if EXT_VERSION_FILE was specified)
#       ${EXT_NAME}_INCLUDE_DIR
#       ${EXT_NAME}_LIBRARIES              (only if EXT_COMPONENTS were specified)
#       ${EXT_NAME}_LIB_DIRS               (only if EXT_COMPONENTS were specified)
#



######## Helper Functions #########

macro( appendPaths )
    # Append to current lists if found
    if( ${EXT_NAME}_FOUND )
        append_libs_unique( REDUX_CURRENT_LIBRARIES "${${EXT_NAME}_LIBRARIES}" )
        append_paths_unique( REDUX_CURRENT_INCLUDES ${${EXT_NAME}_INCLUDE_DIR} )
        append_paths_unique( REDUX_CURRENT_LIBDIRS ${${EXT_NAME}_LIB_DIRS} )
    endif()
endmacro()

string(TOUPPER ${EXT_NAME} EXT_UPPERNAME)

macro( storeVars )
    if(EXT_DEBUG)
        message(STATUS "FindExternal:\"${EXT_NAME}\"  storeVars()")
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_INCLUDE_DIR: \"${${EXT_NAME}_INCLUDE_DIR}\"")
    endif()
    set( ${EXT_NAME}_FOUND ${${EXT_NAME}_FOUND} CACHE INTERNAL "" )
    set( ${EXT_NAME}_DIR ${${EXT_NAME}_DIR} CACHE PATH "${EXT_NAME} location" FORCE )
    if(EXT_VERSION_FILE)
        set( ${EXT_NAME}_VERSION ${${EXT_NAME}_VERSION} CACHE STRING "Version" FORCE )
    endif()
    set( ${EXT_NAME}_INCLUDE_DIR ${${EXT_NAME}_INCLUDE_DIR} CACHE PATH "Location of header files" FORCE )
    mark_as_advanced( ${EXT_NAME}_DIR ${EXT_NAME}_VERSION ${EXT_NAME}_INCLUDE_DIR )
    if(EXT_COMPONENTS)
        set( ${EXT_NAME}_LIBRARIES ${${EXT_NAME}_LIBRARIES}  CACHE STRING "List of found libraries" FORCE )
        set( ${EXT_NAME}_LIB_DIRS ${${EXT_NAME}_LIB_DIRS} CACHE STRING "Locations of found libraries" FORCE )
    mark_as_advanced( ${EXT_NAME}_LIBRARIES ${EXT_NAME}_LIB_DIRS )
    endif()
endmacro()


macro( clearVars )
    if(EXT_DEBUG)
        message("FindExternal:\"${EXT_NAME}\"  clearVars()")
    endif()
    set( ${EXT_NAME}_DIR "" CACHE PATH "" )
    mark_as_advanced( ${EXT_NAME}_DIR )
    unset( ${EXT_NAME}_FOUND )
    unset( ${EXT_NAME}_VERSION CACHE )
    unset( ${EXT_NAME}_INCLUDE_DIR CACHE )
    unset( ${EXT_NAME}_LIBRARIES CACHE )
    unset( ${EXT_NAME}_LIB_DIRS CACHE )
endmacro()


macro( clearLocalVars )
    #if(EXT_DEBUG)
    #    message(STATUS "FindExternal:\"${EXT_NAME}\"  clearLocalVars()")
    #endif()
    set( EXT_DEBUG PARENT_SCOPE)
    unset( EXT_DEBUG )
    unset( EXT_DEBUG CACHE )
    unset( EXT_FOUND )
    unset( EXT_REQUIRED_VERSION )
    unset( EXT_LIB_PREFIX )
    unset( EXT_LIB_SUFFIXES )
    unset( EXT_LIB_DEBUG_SUFFIX )
    unset( EXT_INCLUDE_SUFFIXES )
    unset( EXT_LIBPATH_SUFFIXES )
    unset( EXT_HINT )
    unset( EXT_HEADER_FILE )
    unset( EXT_VERSION_FILE )
    unset( EXT_COMPONENTS )
    unset( EXT_PATHS )
    unset( EXT_VERSION_STR )
    unset( EXT_VERSION_MAJOR )
    unset( EXT_VERSION_MINOR )
    unset( EXT_VERSION_PATCH )
    unset( EXT_VERSION_REGEXP )
    unset( EXT_MAJOR_REGEXP )
    unset( EXT_MINOR_REGEXP )
    unset( EXT_PATCH_REGEXP )
endmacro()


# ------------------------------------------------------------------------
# Helper function to extract version tags of the form
# *VERSION_(MAJOR|MINOR|(BUILD|PATCH|SUBMINOR)) from a given file.
# Values are exported to parent scope as variables:
#    "VERSION_STR"
#    "VERSION_MAJOR"
#    "VERSION_MINOR"
#    "VERSION_PATCH"
# ------------------------------------------------------------------------
macro( EXTRACT_FILE_VERSION FILE_NAME )

    if(EXT_DEBUG)
        message(STATUS "FindExternal:\"${EXT_NAME}\"  parsing version file: \"${FILE_NAME}\"")
    endif()

    set( _VERSION_STR )
    set( _MAJOR )
    set( _MINOR )
    set( _PATCH )

    if(NOT EXT_VERSION_REGEXP)
        set(EXT_VERSION_REGEXP "DISABLED__TOO_GENERIC")
    endif()
    if(NOT EXT_MAJOR_REGEXP)
        set( EXT_MAJOR_REGEXP "_MAJOR" )
    endif()
    if(NOT EXT_MINOR_REGEXP)
        set( EXT_MINOR_REGEXP "_MINOR" )
    endif()
    if(NOT EXT_PATCH_REGEXP)
        set( EXT_PATCH_REGEXP "(_PATCH|_BUILD|_SUBMINOR|_RELEASE|_BETA|_REVISION)")
    endif()


    if( EXISTS ${FILE_NAME} )
        # parse file
        file( STRINGS ${FILE_NAME} STR1 REGEX "${EXT_MAJOR_REGEXP}" )
        file( STRINGS ${FILE_NAME} STR2 REGEX "${EXT_MINOR_REGEXP}" )
        file( STRINGS ${FILE_NAME} STR3 REGEX "${EXT_PATCH_REGEXP}" )
        file( STRINGS ${FILE_NAME} STR4 REGEX "${EXT_VERSION_REGEXP}" )

        if(EXT_DEBUG)
            message(STATUS "FindExternal:\"${EXT_NAME}\"  STR1: ${STR1}")
            message(STATUS "FindExternal:\"${EXT_NAME}\"  STR2: ${STR2}")
            message(STATUS "FindExternal:\"${EXT_NAME}\"  STR3: ${STR3}")
            message(STATUS "FindExternal:\"${EXT_NAME}\"  STR4: ${STR4}")
        endif()
        
        if(STR4)
            string( REGEX MATCHALL "[0-9]+[.0-9]*" _VERSION_STR ${STR4} )
            string( SUBSTRING "${_VERSION_STR}" 0 1 _MAJOR )
            string( SUBSTRING "${_VERSION_STR}" 2 1 _MINOR )
            string( SUBSTRING "${_VERSION_STR}" 4 1 _PATCH )
        else()
            if(STR1)
                string( REGEX MATCHALL "[0-9]+" _MAJOR ${STR1} )
                set( _VERSION_STR "${_MAJOR}" )
            endif()
            if(STR2)
                string( REGEX MATCHALL "[0-9]+" _MINOR ${STR2} )
                set( _VERSION_STR "${_VERSION_STR}.${_MINOR}" )
            endif()
            if(STR3)
                string( REGEX MATCHALL "[0-9]+" _PATCH ${STR3} )
                set( _VERSION_STR "${_VERSION_STR}.${_PATCH}" )
            endif()
        endif()
    endif()

    # export variables to parent scope
    set( EXT_VERSION_STR   ${_VERSION_STR} )
    set( EXT_VERSION_MAJOR ${_MAJOR} )
    set( EXT_VERSION_MINOR ${_MINOR} )
    set( EXT_VERSION_PATCH ${_PATCH} )

endmacro()


macro( notFound )
    if( ${EXT_NAME}_HELPTEXT )
        message( ${${EXT_NAME}_HELPTEXT} )
    endif()
    clearVars()
    clearLocalVars()
    return()
endmacro()

##### End of helper functions #####


# return immediately if we already located this package.
# if( ${EXT_NAME}_FOUND AND IS_DIRECTORY "${${EXT_NAME}_INCLUDE_DIR}" )
#     #if(EXT_DEBUG)
#     #    message(STATUS "FindExternal:\"${EXT_NAME}\"  already located, skipping.")
#     #endif()
#     clearLocalVars()
#     return()
# endif()

 # Clear variables
clearVars()


# Set some default subdirs to search for libraries in (usually lib64 and/or lib)
if( MSVC AND NOT CMAKE_CL_64 )
    set( DEFAULT_LIBPATH_SUFFIXES "lib" )
else()
    set( DEFAULT_LIBPATH_SUFFIXES "lib64" "lib" )
endif()
if( EXT_LIBPATH_SUFFIXES )
    set( DEFAULT_LIBPATH_SUFFIXES "${EXT_LIBPATH_SUFFIXES};${DEFAULT_LIBPATH_SUFFIXES}" )
endif()

set( DEFAULT_INCLUDE_SUFFIXES "" "include" "${EXT_NAME}" "${EXT_NAME}/include" "include/${EXT_NAME}" ".." )
if( EXT_INCLUDE_SUFFIXES )
    set( DEFAULT_INCLUDE_SUFFIXES "${DEFAULT_INCLUDE_SUFFIXES};${EXT_INCLUDE_SUFFIXES}" )
endif()

# Check if ${EXT_UPPERNAME}_DIR is defined in the environment
# skip if ${EXT_NAME}_DIR already is a valid directory
if( NOT IS_DIRECTORY "${${EXT_NAME}_DIR}" )
    if( IS_DIRECTORY "${${EXT_UPPERNAME}_DIR}" )
        set( ${EXT_NAME}_DIR "${${EXT_UPPERNAME}_DIR}" )
    elseif( DEFINED ENV{${EXT_UPPERNAME}_DIR} )
        set( ${EXT_NAME}_DIR "$ENV{${EXT_UPPERNAME}_DIR}" )
    endif()
endif()


if(EXT_DEBUG)
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_REQUIRED_VERSION: \"${EXT_REQUIRED_VERSION}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_LIB_PREFIX: \"${EXT_LIB_PREFIX}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_LIB_SUFFIXES: \"${EXT_LIB_SUFFIXES}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  DEFAULT_LIBPATH_SUFFIXES: \"${DEFAULT_LIBPATH_SUFFIXES}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  DEFAULT_INCLUDE_SUFFIXES: \"${DEFAULT_INCLUDE_SUFFIXES}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_HINT: \"${EXT_HINT}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_PATHS: \"${EXT_PATHS}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_HEADER_FILE: \"${EXT_HEADER_FILE}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_VERSION_FILE: \"${EXT_VERSION_FILE}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_COMPONENTS: \"${EXT_COMPONENTS}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_VERSION_REGEXP: \"${EXT_VERSION_REGEXP}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_MAJOR_REGEXP: \"${EXT_MAJOR_REGEXP}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_MINOR_REGEXP: \"${EXT_MINOR_REGEXP}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  EXT_PATCH_REGEXP: \"${EXT_PATCH_REGEXP}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  CMAKE_SYSTEM_PREFIX_PATH   \"${CMAKE_SYSTEM_PREFIX_PATH}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  CMAKE_SYSTEM_LIBRARY_PATH   \"${CMAKE_SYSTEM_LIBRARY_PATH}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  CMAKE_LIBRARY_PATH   \"${CMAKE_LIBRARY_PATH}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  CMAKE_INCLUDE_PATH   \"${CMAKE_INCLUDE_PATH}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  CMAKE_LIBRARY_PATH   \"$ENV{CMAKE_LIBRARY_PATH}\"")
    message(STATUS "FindExternal:\"${EXT_NAME}\"  CMAKE_INCLUDE_PATH   \"$ENV{CMAKE_INCLUDE_PATH}\"")
endif()



# search for EXT_HEADER_FILE to determine the include path.
if ( NOT EXISTS ${${EXT_NAME}_INCLUDE_DIR} )
    find_path(
        ${EXT_NAME}_INCLUDE_DIR "${EXT_HEADER_FILE}"
        HINTS "${${EXT_NAME}_DIR}" ${EXT_HINT} ${EXT_PATHS} ENV CMAKE_INCLUDE_PATH
        PATHS ${EXT_PATHS} "${EXT_HINT}" ENV CMAKE_INCLUDE_PATH
        PATH_SUFFIXES ${DEFAULT_INCLUDE_SUFFIXES}
    )
    if( EXISTS ${${EXT_NAME}_INCLUDE_DIR} )
        get_filename_component(${EXT_NAME}_INCLUDE_DIR ${${EXT_NAME}_INCLUDE_DIR} REALPATH)
        string( REGEX REPLACE "/include(.*)" "" ${EXT_NAME}_DIR "${${EXT_NAME}_INCLUDE_DIR}" )
    endif()
endif()

if(EXT_DEBUG)
    message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_INCLUDE_DIR: \"${${EXT_NAME}_INCLUDE_DIR}\"")
endif()

# set ${EXT_NAME}_DIR to the parent of the include path
if( EXISTS ${${EXT_NAME}_INCLUDE_DIR} )
    if( NOT EXISTS ${${EXT_NAME}_DIR} )
        string( REGEX REPLACE "/include(.*)" "" ${EXT_NAME}_DIR "${${EXT_NAME}_INCLUDE_DIR}" )
    endif()
else()
    unset( ${EXT_NAME}_INCLUDE_DIR CACHE )
endif()

# extract version
if( EXT_VERSION_FILE )

    find_file( ${EXT_NAME}_VERSION_FILE
        NAMES ${EXT_VERSION_FILE}
        HINTS "${${EXT_NAME}_DIR}" "${EXT_HINT}" ${EXT_PATHS} ENV CMAKE_INCLUDE_PATH
        PATHS "${${EXT_NAME}_DIR}" "${EXT_HINT}" "${${EXT_NAME}_INCLUDE_DIR}" "${${EXT_NAME}_INCLUDE_DIR}/.." ENV CMAKE_INCLUDE_PATH
        PATH_SUFFIXES ${DEFAULT_INCLUDE_SUFFIXES}
    )
    if(EXT_DEBUG)
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_DIR: \"${${EXT_NAME}_DIR}\"")
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_VERSION_FILE: \"${${EXT_NAME}_VERSION_FILE}\"")
    endif()

    if ( EXISTS ${${EXT_NAME}_VERSION_FILE} )
        extract_file_version( ${${EXT_NAME}_VERSION_FILE} )
        set( ${EXT_NAME}_VERSION "${EXT_VERSION_STR}" CACHE INTERNAL "" )
        set( ${EXT_NAME}_VERSION_MAJOR "${EXT_VERSION_MAJOR}" CACHE INTERNAL "" )
        set( ${EXT_NAME}_VERSION_MINOR "${EXT_VERSION_MINOR}" CACHE INTERNAL "" )
        set( ${EXT_NAME}_VERSION_PATCH "${EXT_VERSION_PATCH}" CACHE INTERNAL "" )
    endif()
    if(EXT_DEBUG)
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_VERSION: \"${${EXT_NAME}_VERSION}\"")
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_VERSION_MAJOR: \"${${EXT_NAME}_VERSION_MAJOR}\"")
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_VERSION_MINOR: \"${${EXT_NAME}_VERSION_MINOR}\"")
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_VERSION_PATCH: \"${${EXT_NAME}_VERSION_PATCH}\"")
    endif()
    unset( ${EXT_NAME}_VERSION_FILE CACHE )

    if( ${EXT_NAME}_VERSION AND EXT_REQUIRED_VERSION AND ${EXT_NAME}_VERSION VERSION_LESS EXT_REQUIRED_VERSION )
        message("FindExternal:\"${EXT_NAME}\"  Version mismatch.  Found: ${${EXT_NAME}_VERSION} < ${EXT_REQUIRED_VERSION}")
        notFound()
    endif()
    
endif()

if( EXT_COMPONENTS )
    set( ${EXT_NAME}_COMPONENTS CACHE INTERNAL "" FORCE )
    set( ${EXT_NAME}_REQ CACHE INTERNAL "" FORCE )
    set( ${EXT_NAME}_LIB_DIRS CACHE INTERNAL "" FORCE )

#     if(EXT_VERSION_MAJOR)
#         list( APPEND EXT_LIB_SUFFIXES "${EXT_VERSION_MAJOR}" )
#         if(EXT_VERSION_MINOR)
#             list( APPEND EXT_LIB_SUFFIXES "${EXT_VERSION_MAJOR}.${EXT_VERSION_MINOR}" )
#         endif()
#     endif()
#     
    foreach( __COMP ${EXT_COMPONENTS} )
        list( APPEND ${EXT_NAME}_COMPONENTS "${__COMP}" )
        set( ${EXT_NAME}_REQ "${${EXT_NAME}_REQ} ${EXT_NAME}_${__COMP}_LIBRARY" )
    endforeach()
    set( ${EXT_NAME}_FIND_COMPONENTS ${${EXT_NAME}_COMPONENTS} )

    # construct template for library names
    if(NOT DEFINED EXT_LIB_SUFFIXES)
        set(__LIBNAME_TEMPLATE ${EXT_LIB_PREFIX}__COMP__ )
    else()
        set( __LIBNAME_TEMPLATE )
        foreach( __SUFFIX ${EXT_LIB_SUFFIXES} )
            list( APPEND __LIBNAME_TEMPLATE ${EXT_LIB_PREFIX}__COMP__${__SUFFIX})
        endforeach()
        list( APPEND __LIBNAME_TEMPLATE ${EXT_LIB_PREFIX}__COMP__) # no prefix as last resort.
    endif()

    if(EXT_DEBUG)
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_REQ : \"${${EXT_NAME}_REQ}\"")
        message("FindExternal:\"${EXT_NAME}\"  __LIBNAME_TEMPLATE : \"${__LIBNAME_TEMPLATE}\"")
    endif()
    # if "EXT_LIB_DEBUG_SUFFIX" is set, also look for debug libraries
    set( __LIBNAME_DEBUG_TEMPLATE )
    if( EXT_LIB_DEBUG_SUFFIX )
        foreach( __NAME ${__LIBNAME_TEMPLATE} )
            list( APPEND __LIBNAME_DEBUG_TEMPLATE ${__NAME}${EXT_LIB_DEBUG_SUFFIX})
        endforeach()
        if(EXT_DEBUG)
            message("FindExternal:\"${EXT_NAME}\"  EXT_LIB_DEBUG_SUFFIX : \"${EXT_LIB_DEBUG_SUFFIX}\"")
        endif()
    endif()

    # find libraries for all components
    foreach( __COMP IN LISTS ${EXT_NAME}_COMPONENTS )

        set( ${EXT_NAME}_${__COMP}_LIBRARY_RELEASE )
        set( ${EXT_NAME}_${__COMP}_LIBRARY_DEBUG )

        # find release library
        string( REGEX REPLACE "__COMP__" "${__COMP}" __LIBNAMES "${__LIBNAME_TEMPLATE}" )
        find_library(
            ${EXT_NAME}_${__COMP}_LIBRARY_RELEASE
            NAMES ${__LIBNAMES}
            HINTS "${${EXT_NAME}_DIR}" "${${EXT_NAME}_DIR}/.." "${EXT_HINT}" ${EXT_PATHS} ENV CMAKE_LIBRARY_PATH
            PATH_SUFFIXES ${DEFAULT_LIBPATH_SUFFIXES}
            PATHS "${${EXT_NAME}_DIR}" "${${EXT_NAME}_DIR}/.." ${EXT_PATHS} "${EXT_HINT}" ENV CMAKE_LIBRARY_PATH
        )
        if(EXT_DEBUG)
            message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_${__COMP}_LIBRARY_RELEASE : \"${${EXT_NAME}_${__COMP}_LIBRARY_RELEASE}\"")
        endif()

        # find debug library
        if( __LIBNAME_DEBUG_TEMPLATE )
            string( REGEX REPLACE "__COMP__" "${__COMP}" __LIBNAMES "${__LIBNAME_DEBUG_TEMPLATE}" )
            find_library(
                ${EXT_NAME}_${__COMP}_LIBRARY_DEBUG
                NAMES ${__LIBNAMES}
                HINTS "${${EXT_NAME}_DIR}" "${${EXT_NAME}_DIR}/.." "${EXT_HINT}" ${EXT_PATHS} ENV CMAKE_LIBRARY_PATH
                PATH_SUFFIXES ${DEFAULT_LIBPATH_SUFFIXES}
                PATHS "${${EXT_NAME}_DIR}" "${${EXT_NAME}_DIR}/.." ${EXT_PATHS} "${EXT_HINT}" ENV CMAKE_LIBRARY_PATH
            )
            if(EXT_DEBUG)
                message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_${__COMP}_LIBRARY_DEBUG : \"${${EXT_NAME}_${__COMP}_LIBRARY_DEBUG}\"")
            endif()
        endif()

        unset( ${EXT_NAME}_${__COMP}_LIBRARY )
        unset( ${EXT_NAME}_${__COMP}_LIBRARY CACHE )

        if( ${EXT_NAME}_${__COMP}_LIBRARY_DEBUG )
            get_filename_component( LIBPATH "${${EXT_NAME}_${__COMP}_LIBRARY_DEBUG}" PATH )
            get_filename_component(${EXT_NAME}_${__COMP}_LIBRARY_DEBUG ${${EXT_NAME}_${__COMP}_LIBRARY_DEBUG} NAME_WE)
            string(REGEX REPLACE "^lib" "" ${EXT_NAME}_${__COMP}_LIBRARY_DEBUG ${${EXT_NAME}_${__COMP}_LIBRARY_DEBUG})
            set( ${EXT_NAME}_${__COMP}_LIBRARY ${${EXT_NAME}_${__COMP}_LIBRARY_DEBUG} )
            append_paths_unique( ${EXT_NAME}_LIB_DIRS "${LIBPATH}/" )
        endif()

        if( ${EXT_NAME}_${__COMP}_LIBRARY_RELEASE )
            get_filename_component( LIBPATH "${${EXT_NAME}_${__COMP}_LIBRARY_RELEASE}" PATH )
            get_filename_component(${EXT_NAME}_${__COMP}_LIBRARY_RELEASE ${${EXT_NAME}_${__COMP}_LIBRARY_RELEASE} NAME_WE)
            string(REGEX REPLACE "^lib" "" ${EXT_NAME}_${__COMP}_LIBRARY_RELEASE ${${EXT_NAME}_${__COMP}_LIBRARY_RELEASE})
            set( ${EXT_NAME}_${__COMP}_LIBRARY ${${EXT_NAME}_${__COMP}_LIBRARY_RELEASE} )
            append_paths_unique( ${EXT_NAME}_LIB_DIRS "${LIBPATH}/" )
        endif()

        # if both are found...
        if( ${EXT_NAME}_${__COMP}_LIBRARY_DEBUG AND ${EXT_NAME}_${__COMP}_LIBRARY_RELEASE )
            set( ${EXT_NAME}_${__COMP}_LIBRARY debug ${${EXT_NAME}_${__COMP}_LIBRARY_DEBUG} optimized ${${EXT_NAME}_${__COMP}_LIBRARY_RELEASE} )
        endif()

        if(EXT_DEBUG)
            message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_${__COMP}_LIBRARY : \"${${EXT_NAME}_${__COMP}_LIBRARY}\"")
        endif()
        mark_as_advanced( ${EXT_NAME}_${__COMP}_LIBRARY )
        unset( ${EXT_NAME}_${__COMP}_LIBRARY_RELEASE CACHE )
        unset( ${EXT_NAME}_${__COMP}_LIBRARY_DEBUG CACHE )

        # add to list of found libraries, and remove from component list
        if( ${EXT_NAME}_${__COMP}_LIBRARY )
            list( APPEND ${EXT_NAME}_LIBRARIES "${${EXT_NAME}_${__COMP}_LIBRARY}" )
            append_libs_unique( ${EXT_NAME}_LIBRARIES "${${EXT_NAME}_${__COMP}_LIBRARY}" )
            list( REMOVE_ITEM ${EXT_NAME}_COMPONENTS ${__COMP} )
        endif()
    endforeach()
    if(EXT_DEBUG)
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_LIBRARIES : \"${${EXT_NAME}_LIBRARIES}\"")
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_COMPONENTS : \"${${EXT_NAME}_COMPONENTS}\"")
    endif()

    
    # verify that all components were found
    set( ${EXT_NAME}_REQUIRED_COMPONENTS_FOUND TRUE CACHE INTERNAL "" FORCE )
    list( LENGTH ${EXT_NAME}_COMPONENTS REMAINING_COMPONENTS )
    if(EXT_DEBUG)
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_LIBRARIES : \"${${EXT_NAME}_LIBRARIES}\"")
        message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_COMPONENTS : \"${${EXT_NAME}_COMPONENTS}\"")
        message("FindExternal:\"${EXT_NAME}\"  REMAINING_COMPONENTS : \"${REMAINING_COMPONENTS}\"")
    endif()

    if( ${REMAINING_COMPONENTS} GREATER 0 )
        if(EXT_DEBUG)

            message("FindExternal:\"${EXT_NAME}\"  REMAINING_COMPONENTS > 0 : \"${REMAINING_COMPONENTS}\"")
        endif()
        set( ${EXT_NAME}_REQUIRED_COMPONENTS_FOUND FALSE  )
    endif()

    if( NOT ${EXT_NAME}_REQUIRED_COMPONENTS_FOUND )
        message( "The following required ${EXT_NAME} components"
                             " were not found: \"${${EXT_NAME}_COMPONENTS}\"" )
        notFound()
    endif()
else()
    unset( ${EXT_NAME}_LIBRARIES CACHE )
    unset( ${EXT_NAME}_LIB_DIRS CACHE )
endif()
 
if(EXT_DEBUG)
    message("FindExternal:\"${EXT_NAME}\"  ${EXT_NAME}_REQ: \"${${EXT_NAME}_REQ}\"")
else()
    set(${EXT_NAME}_FIND_QUIETLY TRUE)
endif()

include( FindPackageHandleStandardArgs )
# default handler for QUIET/REQUIRED arguments and sets ${EXT_NAME}_FOUND to TRUE if all listed variables are TRUE
if(EXT_COMPONENTS)
    find_package_handle_standard_args( ${EXT_NAME} REQUIRED_VARS ${EXT_NAME}_INCLUDE_DIR ${EXT_NAME}_REQ )
else()
    find_package_handle_standard_args( ${EXT_NAME} REQUIRED_VARS ${EXT_NAME}_INCLUDE_DIR )
endif()

if ( ${EXT_UPPERNAME}_FOUND )
    set(${EXT_UPPERNAME}_FOUND 1)
    if(EXT_DEBUG)
        message(STATUS "FindExternal:\"${EXT_NAME}\" Found.")
    endif()
    storeVars()
else()
    if(EXT_DEBUG)
        message(STATUS "FindExternal:\"${EXT_NAME}\" not found." )
    endif()
    clearVars()
endif()


clearLocalVars()

