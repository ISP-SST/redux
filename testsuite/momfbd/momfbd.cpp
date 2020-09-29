
#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;

namespace testsuite {

    namespace momfbd {
        
        //srand(time(NULL));
        
        void add_config_tests( test_suite* ts );    // defined in config.cpp
        void add_data_tests( test_suite* ts );      // defined in data.cpp

        void add_tests( test_suite* ts ) {

            add_config_tests( ts );
            add_data_tests( ts );

        }

    }

}
