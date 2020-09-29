#include <boost/test/unit_test.hpp>

using namespace boost::unit_test_framework;

namespace testsuite {

    namespace util {
        
        void add_array_tests( test_suite* ts );     // defined in array.cpp
        void add_data_tests( test_suite* ts );      // defined in data.cpp
        void add_string_tests( test_suite* ts );    // defined in string.cpp

        void add_tests( test_suite* ts ) {

            add_array_tests( ts );
            add_data_tests( ts );
            add_string_tests( ts );

        }

    }

}

