#
# Name of test
#
set(TEST_NAME "testsuite" CACHE INTERNAL "Name of test application")

USE_REDUX(redux)

enable_testing()

if( RDX_WITH_FFTW3 )
    add_definitions(-DRDX_WITH_FFTW3)
endif()

# Code for test modules
file(GLOB RDX_TEST_FILE_SRC  RELATIVE "${CMAKE_CURRENT_LIST_DIR}/" "file/*.cpp")
file(GLOB RDX_TEST_IMAGE_SRC RELATIVE "${CMAKE_CURRENT_LIST_DIR}/" "image/*.cpp")
#file(GLOB RDX_TEST_MATH_SRC RELATIVE "${CMAKE_CURRENT_LIST_DIR}/" "math/*.cpp")
#file(GLOB RDX_TEST_MOMFBD_SRC RELATIVE "${CMAKE_CURRENT_LIST_DIR}/" "momfbd/*.cpp")
#file(GLOB RDX_TEST_UTIL_SRC  RELATIVE "${CMAKE_CURRENT_LIST_DIR}/" "util/*.cpp")
set( RDX_TEST_SRC
        "testsuite.cpp"
        ${RDX_TEST_FILE_SRC}
        ${RDX_TEST_IMAGE_SRC}
        ${RDX_TEST_MATH_SRC}
        ${RDX_TEST_MOMFBD_SRC}
        ${RDX_TEST_UTIL_SRC}
)

file(GLOB RDX_TEST_HPP RELATIVE ${CMAKE_CURRENT_LIST_DIR} "include/*.hpp")

set(RDX_TESTDATA_DIR "${CMAKE_CURRENT_LIST_DIR}/testdata/")
add_definitions(-DRDX_TESTDATA_DIR="${RDX_TESTDATA_DIR}")
add_definitions(-DBOOST_TEST_DYN_LINK)


# Turn on support for SEH exceptions to let test catch a exe crash
if(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHa")
endif()

RDX_ADD_EXECUTABLE(${TEST_NAME} ${RDX_TEST_SRC} ${RDX_TEST_HPP})

set(TEST_CMD $<TARGET_FILE:${TEST_NAME}> --log_sink=test_log.xml
                              --log_format=XML
#                              --log_level=all 
                              --log_level=warning 
                              --report_sink=test_report.xml 
                              --report_format=XML 
                              --report_level=detailed )


# Fix linking problems by adding all the link-directories to the path
# before running the testsuite.
if(WIN32)
	get_property(TEST_LINK_DIRECTORIES DIRECTORY PROPERTY LINK_DIRECTORIES)
	file(TO_NATIVE_PATH "${TEST_LINK_DIRECTORIES}" TEST_LINK_DIRECTORIES)
	string(REPLACE ";" "$<SEMICOLON>" TEST_LINK_DIRECTORIES "${TEST_LINK_DIRECTORIES}")
	set(TEST_CMD set PATH=${TEST_LINK_DIRECTORIES}$<SEMICOLON>%PATH% COMMAND ${TEST_CMD})
endif()

add_custom_target( RUN_TESTSUITE
  ${TEST_CMD}
  DEPENDS ${TEST_NAME}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
