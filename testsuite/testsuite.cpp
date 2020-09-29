/*
#ifdef _MSC_VER                 // Microsoft Visual C++ Compiler
# pragma warning(push)	        // Disable for boost test only
#  pragma warning(disable:4244) // Disable warning about time_t being truncated 
                                // when time_t is 64 bits (The time_t is just 
                                // used to get a random seed)
#endif
#include <boost/test/included/unit_test_framework.hpp>
#ifdef _MSC_VER                 // Microsoft Visual C++ Compiler
# pragma warning(pop)  	        // Restore original warning level
#endif
*/

#include <boost/test/included/unit_test.hpp>
using namespace boost::unit_test;


namespace testsuite  {
    namespace file   { void add_tests( test_suite* ts ); }
    namespace image  { void add_tests( test_suite* ts ); }
    namespace math   { void add_tests( test_suite* ts ); }
    namespace momfbd { void add_tests( test_suite* ts ); }
    namespace util   { void add_tests( test_suite* ts ); }
}


test_suite* init_unit_test_suite( int , char* [] ) {

    framework::master_test_suite().p_name.value = "Redux Testsuite";

    test_suite* file_tests = BOOST_TEST_SUITE( "File" );
    testsuite::file::add_tests( file_tests );
    framework::master_test_suite().add( file_tests );

    test_suite* image_tests = BOOST_TEST_SUITE( "Image" );
    testsuite::image::add_tests( image_tests );
    framework::master_test_suite().add( image_tests );

    test_suite* math_tests = BOOST_TEST_SUITE( "Math" );
    testsuite::math::add_tests( math_tests );
    framework::master_test_suite().add( math_tests );

/*    test_suite* momfbd_tests = BOOST_TEST_SUITE( "MOMFBD" );
    testsuite::momfbd::add_tests( momfbd_tests );
    framework::master_test_suite().add( momfbd_tests );
*/
    test_suite* util_tests = BOOST_TEST_SUITE( "Util" );
    testsuite::util::add_tests( util_tests );
    framework::master_test_suite().add( util_tests );


    return 0;
}
