
#include "redux/image/utils.hpp"
#include "redux/util/array.hpp"


#include <boost/test/unit_test.hpp>

using namespace redux::image;
using namespace redux::util;

using namespace std;


namespace testsuite {

    namespace image {

        void test_plane( void ) {
            
            size_t pSize(100);
            double coeffs[3] = { 0, 0, 0 };
            double fit_coeffs[3] = { 0, 0, 0 };
            double chisq(0);
            Array<float> plane( pSize, pSize );
            
            {
                // make a plane with slope in y-direction and constant in x.
                coeffs[0] = 1;
                makePlane( plane, coeffs, 1 );
                // fit and compare
                fitPlane( plane.ptr(), pSize, pSize, fit_coeffs, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs[0], 0.0000001 );
                BOOST_CHECK_SMALL( fit_coeffs[1], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs[2], pSize/2, 0.0000001 );

                // make a plane with slope in x-direction and constant in y.
                coeffs[0] = 0;
                coeffs[1] = 1;
                makePlane( plane, coeffs, 2 );
                // fit and compare
                fitPlane( plane.ptr(), pSize, pSize, fit_coeffs, &chisq );
                BOOST_CHECK_SMALL( fit_coeffs[0], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs[1], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs[2], pSize/2, 0.0000001 );
                
                // make a plane with slope in both x/y.
                coeffs[0] = 1;
                makePlane( plane, coeffs, 2 );
                // fit and compare
                fitPlane( plane.ptr(), pSize, pSize, fit_coeffs, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs[0], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs[1], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs[2], pSize, 0.0000001 );
                
                // make a plane with slope in both x/y and a specified mean value.
                coeffs[0] = 3.333;
                coeffs[1] = 12.34;
                coeffs[2] = -123;
                makePlane( plane, coeffs, 3 );
                // fit and compare
                fitPlane( plane.ptr(), pSize, pSize, fit_coeffs, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs[0], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs[1], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs[2], coeffs[2], 0.0000001 );
                
                // Check that the version returning a plane also gives right results.
                double fit_coeffs2[3] = { 0, 0, 0 };
                Array<float> plane2 = fitPlane( plane, false, fit_coeffs2, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs2[0], coeffs[0], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs2[1], coeffs[1], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs2[2], coeffs[2], 0.0000001 );
 
               // Verification that the returned plane is accurate
                fitPlane( plane2.ptr(), pSize, pSize, fit_coeffs, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs[0], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs[1], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs[2], coeffs[2], 0.0000001 );
                
                // same, but subtract average.
                plane2 = fitPlane( plane, true, fit_coeffs2, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs2[0], coeffs[0], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs2[1], coeffs[1], 0.0000001 );
                BOOST_CHECK_SMALL( fit_coeffs2[2], 0.0000001 );
                
                // Verification that the returned plane is accurate
                fitPlane( plane2.ptr(), pSize, pSize, fit_coeffs, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs[0], 0.0000001 );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs[1], 0.0000001 );
                BOOST_CHECK_SMALL( fit_coeffs[2], 0.0000001 );
                
            }

        }

        void util_tests( void ) {
            
            test_plane();

        }


        using namespace boost::unit_test;
        void add_util_tests( test_suite* ts ) {

            ts->add( BOOST_TEST_CASE_NAME( &util_tests, "Image/Util" ) );

        }

    }

}
