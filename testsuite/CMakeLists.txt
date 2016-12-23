#
# Name of test
#
set(TEST_NAME "testsuite" CACHE INTERNAL "Name of test application")

USE_REDUX(redux)

if( REDUX_WITH_FFTW3 )
    add_definitions(-DREDUX_WITH_FFTW3)
endif()

# Code for test module
file(GLOB REDUX_TEST_SRC RELATIVE ${CMAKE_CURRENT_LIST_DIR} "*.cpp")
file(GLOB REDUX_TEST_HPP RELATIVE ${CMAKE_CURRENT_LIST_DIR} "include/*.hpp")

set(REDUX_TESTDATA_DIR "${CMAKE_CURRENT_LIST_DIR}/testdata/")
add_definitions(-DREDUX_TESTDATA_DIR="${REDUX_TESTDATA_DIR}")

# Turn on support for SEH exceptions to let test catch a exe crash
if(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHa")
endif()

REDUX_ADD_EXECUTABLE(${TEST_NAME} ${REDUX_TEST_SRC} ${REDUX_TEST_HPP})

if("${CMAKE_BUILD_TYPE}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
	get_target_property(TEST_LOCATION ${TEST_NAME} DEBUG_LOCATION)
else()
	get_target_property(TEST_LOCATION ${TEST_NAME} LOCATION)
endif()

set(TEST_CMD ${TEST_LOCATION} --log_sink=test_log.xml
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
