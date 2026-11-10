#
# Name of test
#
set(TEST_NAME "testsuite" CACHE INTERNAL "Name of test application")

USE_REDUX(redux)

enable_testing()

if( RDX_WITH_FFTW3 )
    add_definitions(-DRDX_WITH_FFTW3)
endif()

if( RDX_WITH_FITS )
    add_definitions(-DRDX_WITH_FITS)
endif()

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

set( RDX_TEST_SRC
        "testsuite.cpp"
        "include/stringdiff.cpp"
        "file/ana.cpp"
        "file/file.cpp"
        "file/fits.cpp"
        "image/fourier.cpp"
        "image/grid.cpp"
        "image/util.cpp"
        "image/zernike.cpp"
        "util/array.cpp"
        "util/data.cpp"
        "util/string.cpp"
        "util/util.cpp"
        "math/differentiate.cpp"
        "math/helpers.cpp"
        "math/interval.cpp"
        "math/linalg.cpp"
        "momfbd/config.cpp"
        "momfbd/data.cpp"
        "momfbd/momfbd.cpp"
)


set(RDX_TESTDATA_DIR "${CMAKE_CURRENT_LIST_DIR}/testdata/")
add_definitions(-DRDX_TESTDATA_DIR="${RDX_TESTDATA_DIR}")


# Turn on support for SEH exceptions to let test catch a exe crash
if(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHa")
endif()

RDX_ADD_EXECUTABLE( ${TEST_NAME} ${RDX_TEST_SRC} )

set(TEST_CMD $<TARGET_FILE:${TEST_NAME}> -r junit -o test_log.xml )
# Detailed per-assertion XML report (can be GB size!)
#set(TEST_CMD ${TEST_CMD} COMMAND $<TARGET_FILE:${TEST_NAME}> -r xml -s -o test_report.xml )


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
