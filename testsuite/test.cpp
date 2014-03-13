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

namespace testsuite {
    namespace file { void fileTest(void); }
    namespace image { void imageTest(void); }
    namespace util { void utilTest(void); }
}

test_suite* init_unit_test_suite(int /*argc*/, char* /*argv*/[]) {

	//test_suite* sstTestSuite = BOOST_TEST_SUITE("REDUX Tests");
    
    testsuite::file::fileTest();
    testsuite::image::imageTest();
    testsuite::util::utilTest();

    return nullptr;
}
