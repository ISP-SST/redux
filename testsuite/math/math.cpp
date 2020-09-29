#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;

namespace testsuite {

    namespace math {
        
        void add_differentiate_tests( test_suite* ts );     // defined in differentiate.cpp
        void add_helpers_tests( test_suite* ts );           // defined in helpers.cpp
        void add_interval_tests( test_suite* ts );          // defined in interval.cpp
        void add_linalg_tests( test_suite* ts );            // defined in linalg.cpp

        void add_tests( test_suite* ts ) {

            add_differentiate_tests( ts );
            add_helpers_tests( ts );
            add_interval_tests( ts );
            add_linalg_tests( ts );

        }

    }

}


