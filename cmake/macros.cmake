# CMake convenience macros
#  
# Author: Tomas Hillberg


macro(APPEND_PATHS_UNIQUE LIST)

    list(LENGTH ${LIST} SIZE)
    foreach( elm ${ARGN} )
        #message(STATUS "APPEND_PATH_UNIQUE:  " ${elm})
        if (EXISTS ${elm})
            get_filename_component(elm2 ${elm} REALPATH)
            if (SIZE EQUAL 0 )               
                set(${LIST} ${elm2} CACHE INTERNAL "")
                set(SIZE 1)
            else()
                list(APPEND ${LIST} ${elm2})
            endif()
        endif()
    endforeach()

    list(LENGTH ${LIST} SIZE)
    if (SIZE GREATER 0 )
        list(REMOVE_DUPLICATES ${LIST})
    endif()
endmacro()


macro(APPEND_LIBS_UNIQUE LIST ELEMENT)
    set(prefix "")
    foreach( elm ${ELEMENT} ${ARGN} )
        if( "${elm}" STREQUAL "" )
        elseif( "${elm}" STREQUAL "optimized" OR "${elm}" STREQUAL "debug" )
            set(prefix "${elm};")
        else()
            list(FIND "${LIST}" "${elm}" FOUND)
            if (FOUND EQUAL -1)
                list(APPEND ${LIST} "${prefix}${elm}")
            endif()
            set(prefix "")
        endif()
    endforeach()
endmacro()


macro(ADD_TRANSLATION TARGETNAME ROOTDIR TSFILES)
    set(QM_FILE)

    add_custom_target (${TARGETNAME} DEPENDS ${QM_FILES})
        
    foreach(TSFILE ${${TSFILES}})
        # We have to do this in two steps due to bug in CMake regexp
        string(REGEX REPLACE "\"(.+)\"" "\\1" TSFILEPART_NOCITE ${TSFILE}) 
        string(REGEX REPLACE ".*/(.+)\\.ts" "\\1" TSFILEPART ${TSFILEPART_NOCITE})
        set(QM_FILE ${CMAKE_CURRENT_BINARY_DIR}/${TSFILEPART}.qm)

        add_custom_command(TARGET ${TARGETNAME}
            POST_BUILD COMMAND ${QT_LRELEASE_EXECUTABLE}
            ARGS ${TSFILE} -qm ${QM_FILE}
            DEPENDS ${TSFILE}
        )
        install(FILES ${QM_FILE} DESTINATION ${REDUX_INSTALL_DIR}/translations OPTIONAL)
    endforeach()
        
    string(REPLACE ";" " " TSFILESARG "${${TSFILES}}")
    add_custom_command(TARGET ${TARGETNAME}
        PRE_BUILD COMMAND ${QT_LUPDATE_EXECUTABLE}
        ARGS \"${ROOTDIR}/src\" \"${ROOTDIR}/include\" \"${ROOTDIR}/resources\" -ts ${TSFILESARG}
    )
endmacro()


# -----------------------------------------------------------------------------
# Macro for creating targets which will run CppCheck on the source files.
# Check the option "REDUX_CPPCHECK_TARGETS" (OFF by default) if you want/need these
# targets.
# -----------------------------------------------------------------------------
macro(MAKE_CPPCHECK tgt ...)
    if(EXISTS ${REDUX_CPPCHECK_EXECUTABLE} AND REDUX_CPPCHECK_TARGETS)
        get_property(IDS DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
        set(CPPCHECK_ARGS --force --enable=all -q -rp --std=c99 --std=posix --std=c++11)
        #set(CPPCHECK_ARGS --force --enable=all -q -rp --std=c99 --std=posix --std=c++11 --check-config)
        #list(APPEND CPPCHECK_ARGS "--template='{id}:{message}({file}:{line})'" --xml)
        list(APPEND CPPCHECK_ARGS "--template=gcc")
        #list(APPEND CPPCHECK_ARGS "--template=vs")
        #list(APPEND CPPCHECK_ARGS "--template='{file}:{line}: {severity}: {message}'")
        foreach(id ${IDS} ${CMAKE_SYSTEM_INCLUDE_PATH})
            string(REGEX MATCH ${PROJECT_SOURCE_DIR} FOUND ${id})
            if(EXISTS "${id}" AND FOUND)
                list(APPEND CPPCHECK_ARGS "-I${id}")
            else()
                list(APPEND CPPCHECK_ARGS "-i${id}")
            endif()
        endforeach()
        foreach(item ${ARGV})
            if(EXISTS "${item}")
               list(APPEND CPPCHECK_ARGS "${item}")
            elseif(EXISTS "${CMAKE_CURRENT_LIST_DIR}/${item}")
               list(APPEND CPPCHECK_ARGS "${CMAKE_CURRENT_LIST_DIR}/${item}")
            endif()
        endforeach()
        #list(APPEND CPPCHECK_ARGS "2>${CMAKE_BINARY_DIR}/cppcheck_${tgt}.log")
        add_custom_target( CPPCHECK_${tgt} ${REDUX_CPPCHECK_EXECUTABLE} ${CPPCHECK_ARGS}
            #COMMENT "Running CppCheck on target ${tgt}\noutput: \"${CMAKE_BINARY_DIR}/cppcheck_${tgt}.log\""
            COMMENT "Running CppCheck on target ${tgt}"
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        )
    endif()

endmacro()



# -----------------------------------------------------------------------------
# Macro to set up a custom installation script, which will scan the binary for
# dependencies and (try to) copy all those as well. This is not very well tested
# yet, but we'll see how it behaves...
# Note: The dynamic linking/dispatching of IPP can not be resolved like this
#  (at least I still haven't figure out how). Either copy all libs from IPP, or
#  use static linking (i.e. include specific header and compile for a specific cpu)
# -----------------------------------------------------------------------------
macro(INSTALL_DEPS TGT)
    get_target_property(TARGET_LOCATION ${TGT} LOCATION)
    get_target_property(TARGET_DEBUG_LOCATION ${TGT} DEBUG_LOCATION)
    get_property(TARGET_LINK_DIRECTORIES DIRECTORY PROPERTY LINK_DIRECTORIES)
    # Note: Seems like CMAKE_SYSTEM_PREFIX_PATH is not in cmake_path format by default
    # which causes some "unknown escape character" errors later. So clean it now.
    file(TO_CMAKE_PATH "${CMAKE_SYSTEM_PREFIX_PATH}" CMAKE_SYSTEM_PREFIX_PATH)
    # Replace VS variable $(OutDir) with CMake counterpart
    #get_filename_component(TMP_TARGET_PATH ${TARGET_LOCATION} REALPATH)
    string(REPLACE "$(OutDir)" "$${}{CMAKE_INSTALL_CONFIG_NAME}" TARGET_LOCATION ${TARGET_LOCATION})
    string(REPLACE "$(Configuration)" "$${}{CMAKE_INSTALL_CONFIG_NAME}" TARGET_LOCATION ${TARGET_LOCATION})
    # Ordinary install of the binary
    install(TARGETS ${TGT}
        RUNTIME DESTINATION ${REDUX_INSTALL_DIR}/bin
        LIBRARY DESTINATION ${REDUX_INSTALL_DIR}/lib
        ARCHIVE DESTINATION ${REDUX_INSTALL_DIR}/lib
        OPTIONAL
    )
    # Generate an install script which will check dependencies at install-time and copy those as well.
    configure_file(${REDUX_DIR}/cmake/install.cmake ${CMAKE_CURRENT_BINARY_DIR}/${TGT}_install.cmake)
    install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/${TGT}_install.cmake)
endmacro()




macro(USE_EXTERNAL ...)
    foreach(ext ${ARGV})
        file(TO_CMAKE_PATH "${REDUX_DIR}/cmake/use_${ext}.cmake" EXT_CONFIG)
        if(EXISTS "${EXT_CONFIG}")
            include( "${EXT_CONFIG}" )
        else()
            file(TO_CMAKE_PATH "${PROJECT_SOURCE_DIR}/cmake/use_${ext}.cmake" EXT_CONFIG2)
            if(EXISTS "${EXT_CONFIG2}")
                include( "${EXT_CONFIG2}" )
            else()
                message(WARNING "Can't locate configuration file for ${ext} (\"${EXT_CONFIG}\")")
            endif()
        endif()
    endforeach()
endmacro()


# -----------------------------------------------------------------------------
# This macro will accept an array and call the USE_EXTERNAL macro for each
# module-dependency listed in redux???_DEPS (see the dependencies.cmake file for all
# the module dependencies)
# It will first clear the REDUX_CURRENT_LIBRARIES variable, and also append the
# arguments to REDUX_CURRENT_LIBRARIES
# The common way to use these macros for an individual target:
#   USE_REDUX( reduxgui ${REDUX_LIBLIST} )    # will load all the dependencies for
#           all redux-modules, and also append the include/link paths needed.
#   USE_EXTERNAL( sapera photonfocus )  # Add additional dependencies
#   add_executable( "${MY_APPLICATION}" ${MY_SRC} ${MY_HPP} )
#   target_link_libraries( ${MY_APPLICATION} ${REDUX_CURRENT_LIBRARIES} )
# -----------------------------------------------------------------------------
macro(USE_REDUX ...)
    set(external_deps "")
    set(REDUX_CURRENT_LIBRARIES "")
    set(REDUX_CURRENT_INCLUDES "")
    set(REDUX_CURRENT_LIBDIRS "")
    foreach(lib ${ARGV})
        if(DEFINED ${lib}_DEPS)
            list(APPEND external_deps ${${lib}_DEPS})
        else()
            message(STATUS "${lib}_DEPS not defined.")
        endif()
        list(APPEND REDUX_CURRENT_LIBRARIES ${lib})
    endforeach()
    if(NOT external_deps STREQUAL "")
        USE_EXTERNAL(${external_deps})
    endif()
endmacro()



# -----------------------------------------------------------------------------
# Macro for delayloading DLLs on windows. Accepts multiple arguments.
# -----------------------------------------------------------------------------
macro(DELAYLOAD ...)
    if(WIN32)
        foreach( lib ${ARGV})
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /delayload:\"${lib}\"")
        endforeach()
    endif()
endmacro()



# -----------------------------------------------------------------------------
# A couple of convenience functions 
# -----------------------------------------------------------------------------
macro(REDUX_SETDIRS)
    if(DEFINED REDUX_CURRENT_INCLUDES)
        list(REMOVE_DUPLICATES REDUX_CURRENT_INCLUDES)
    endif()
    if(DEFINED REDUX_CURRENT_LIBDIRS)
        list(REMOVE_DUPLICATES REDUX_CURRENT_LIBDIRS)
    endif()
    include_directories( ${REDUX_CURRENT_INCLUDES} )
    link_directories( ${REDUX_CURRENT_LIBDIRS} )
endmacro()


macro(REDUX_ADD_EXECUTABLE target)
    REDUX_SETDIRS()
    add_executable(${target} ${ARGN})
    get_target_property(src ${target} SOURCES)
    MAKE_CPPCHECK(${target} ${src})
    target_link_libraries(${target} ${REDUX_CURRENT_LIBRARIES} ${REDUX_CURRENT_LIBRARIES}) # REDUX_LIBLIST is needed a 2nd time because of circular dependencies
    if(WIN32)   # allow delayloading DLLs - only for windows
        target_link_libraries(${target} delayimp)
    endif()
    set_target_properties(${target} PROPERTIES DEBUG_POSTFIX "-d")
    # install_deps needs some more tweaking...
    install_deps(${target})
endmacro()


macro(REDUX_ADD_LIBRARY target)
    REDUX_SETDIRS()
    add_library(${target} ${ARGN})
    get_target_property(src ${target} SOURCES)
    MAKE_CPPCHECK(${target} ${src})
    # TODO: Eventually the install-dir should really be user configurable.
    set_target_properties(${target} PROPERTIES DEBUG_POSTFIX "-d")
    install(TARGETS ${target} DESTINATION ${REDUX_INSTALL_DIR}/lib)
endmacro()


macro(FIX_CYGWIN_REPO)

    set(CYGWIN_GIT)
    
    if(IS_DIRECTORY "${GIT_WORKTREE}/.git")
        file(STRINGS "${GIT_WORKTREE}/.git/config" GITFILE_CONTENTS REGEX "submodule")
        string(REGEX REPLACE "\\[submodule \"([A-z0-9]+)\"\\]" "\\1" GIT_MODULE_LIST "${GITFILE_CONTENTS}")
        foreach(mod ${GIT_MODULE_LIST})
            file(STRINGS "${GIT_WORKTREE}/.git/modules/${mod}/config" GITFILE_CONTENTS REGEX "worktree(.)+ \\/cygdrive\\/")
            if(GITFILE_CONTENTS)
                set(CYGWIN_GIT "1")
            endif()
        endforeach()
    elseif(EXISTS "${GIT_WORKTREE}/.git") # submodule
        file(READ "${GIT_WORKTREE}/.git" GITFILE_CONTENTS)
        string(REGEX REPLACE "gitdir: (.*)\n$" "\\1" GIT_DIR "${GITFILE_CONTENTS}")
        if(${GIT_DIR} MATCHES "^/cygdrive/.*") # cygwin path
            set(CYGWIN_GIT "1")
        endif()
    endif()

    if(CYGWIN_GIT)
        message("This repo seems to have been created/cloned using cygwin-git, which uses cygwin specific paths.")
        message(STATUS "  -> attempting to use cygwin-git to retreive revision information")
        set(GIT_EXECUTABLE GIT_EXECUTABLE-NOTFOUND)     # cache has to be cleared to search again
        find_program(GIT_EXECUTABLE PATHS "c:/cygwin/bin/" NAMES "git" "git.exe" NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_PATH)
    else()
        set(GIT_EXECUTABLE GIT_EXECUTABLE-NOTFOUND)     # cache has to be cleared to search again
        find_program(GIT_EXECUTABLE PATHS "c:/Program Files (x86)/Git/bin/" NAMES "git" "git.exe" "tgit.exe" NO_CMAKE_PATH)
    endif()

endmacro()


# -----------------------------------------------------------------------------
# Macro for generic (git) versioning
# The template file will be configured with the following variables defined:
#     REVISION_COMMIT_TIME      (timestamp as string)
#     REVISION_COMMIT_COMMENT   (string)
#     REVISION_BUILD_TIME       (timestamp as string)
#     REVISION_ID               (string, e.g hash-tag)
#     REVISION_DESCRIPTION      (string, id + comment)
# -----------------------------------------------------------------------------
macro(CHECK_REVISION REVNAME REVISION_TEMPLATE_LOCATION REVISION_FILE_LOCATION) # optional 4th argument = path

    if(IS_DIRECTORY "${ARGN}") # default path is current project
        set(GIT_WORKTREE "${ARGN}")
    else()
        set(GIT_WORKTREE "${PROJECT_SOURCE_DIR}")
    endif()

    FIX_CYGWIN_REPO()
    if (EXISTS "${GIT_EXECUTABLE}" AND EXISTS "${GIT_WORKTREE}/.git")
        add_custom_target( GIT_CHECK_${REVNAME}
                        COMMAND ${CMAKE_COMMAND}
                            -DREVNAME:STRING="${REVNAME}"
                            -DGIT_DIR:STRING="${GIT_WORKTREE}"
                            -DREVISION_TEMPLATE_LOCATION="${REVISION_TEMPLATE_LOCATION}"
                            -DREVISION_FILE_LOCATION="${REVISION_FILE_LOCATION}"
                            -DBUILD_DIR="${CMAKE_BINARY_DIR}"
                            -P "${REDUX_DIR}/cmake/revision.cmake"
                            ${CMAKE_BINARY_DIR}
        )
        set_source_files_properties(${REVISION_FILE_LOCATION} PROPERTIES GENERATED 1)
    else()
        message(STATUS "Couldn't find git executable and/or .git directory")
        # Set default values for generating file without git data.
        set(REVISION_COMMIT_TIME "unknown")
        set(REVISION_BUILD_TIME "unknown")
        set(REVISION_ID "unknown")
        set(REVISION_COMMIT_COMMENT "No revision control")
        set(REVISION_DESCRIPTION "No revision control")
        set(${REVNAME}_VERSION_MAJOR 0)
        set(${REVNAME}_VERSION_PATCH 0)
        set(${REVNAME}_VERSION_COMMIT 0)
        configure_file("${REVISION_TEMPLATE_LOCATION}" "${REVISION_FILE_LOCATION}")
    endif()

endmacro()

