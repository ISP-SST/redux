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

using boost::unit_test_framework::test_suite;

#include <iostream>
using namespace std;

namespace testsuite {
    namespace file { void fileTest(void); }
    namespace image { void imageTest(void); }
    namespace math { void mathTest(void); }
    namespace momfbd { void momfbdTest(void); }
    namespace util { void utilTest(void); }
}

test_suite* init_unit_test_suite(int, char* []) {

	//test_suite* reduxTestSuite = BOOST_TEST_SUITE("REDUX Tests");
    cout << "muu " << __LINE__ << endl;
    //testsuite::file::fileTest();
    cout << "muu " << __LINE__ << endl;
    testsuite::image::imageTest();
    cout << "muu " << __LINE__ << endl;
    //testsuite::math::mathTest();
    cout << "muu " << __LINE__ << endl;
    //testsuite::momfbd::momfbdTest();
    cout << "muu " << __LINE__ << endl;
    //testsuite::util::utilTest();
    cout << "muu " << __LINE__ << endl;

    return nullptr;
}


*/

/*
#define BOOST_TEST_MODULE rdx_tests
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_SUITE( rdx_test_main )


BOOST_AUTO_TEST_CASE( rdx_ts_case1 )
{
  int a = 13, b = 12;
  BOOST_TEST(a == b);
  BOOST_TEST(a < b);
  BOOST_TEST(a - 1 < b);
  BOOST_TEST(b > a - 1);
}

BOOST_AUTO_TEST_SUITE_END()
*/



#include <boost/test/included/unit_test.hpp>
using namespace boost::unit_test;



namespace testsuite {
    namespace file { void add_tests( test_suite* ts ); }
    namespace image { void add_tests( test_suite* ts ); }
    namespace math { void add_tests( test_suite* ts ); }
    namespace momfbd { void add_tests( test_suite* ts ); }
    namespace util { void add_tests( test_suite* ts ); }
}


test_suite* init_unit_test_suite( int , char* [] ) {

    framework::master_test_suite().p_name.value = "Redux Testsuite";

    test_suite* file_tests = BOOST_TEST_SUITE( "redux::file" );
    testsuite::file::add_tests(file_tests);
    framework::master_test_suite().add( file_tests );

    test_suite* image_tests = BOOST_TEST_SUITE( "redux::image" );
    testsuite::image::add_tests(image_tests);
    framework::master_test_suite().add( image_tests );


    return 0;
}
