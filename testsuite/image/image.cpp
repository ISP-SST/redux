
#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;

namespace testsuite {

    namespace image {
        
        //srand(time(NULL));
        
        void add_fourier_tests( test_suite* ts );   // defined in fourier.cpp
        void add_grid_tests( test_suite* ts );      // defined in grid.cpp
        void add_zernike_tests( test_suite* ts );   // defined in grid.cpp

        void add_tests( test_suite* ts ) {

            add_fourier_tests( ts );
            add_grid_tests( ts );
            add_zernike_tests( ts );

        }

    }

}

