
#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;

namespace testsuite {
    namespace file {
        
        void add_ana_tests( test_suite* ts );       // defined in ana.cpp

        void add_tests( test_suite* ts ) {

            add_ana_tests( ts );

        }

    }
}
