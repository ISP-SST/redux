
#include "redux/image/utils.hpp"
#include "redux/util/array.hpp"


#include <boost/test/unit_test.hpp>

using namespace redux::image;
using namespace redux::util;

using namespace std;


namespace testsuite {

    namespace image {

        void test_plane( void ) {
            
            double EPS = 1E-5;
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
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs[0], EPS );
                BOOST_CHECK_SMALL( fit_coeffs[1], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[2], pSize/2, EPS );

                // make a plane with slope in x-direction and constant in y.
                coeffs[0] = 0;
                coeffs[1] = 1;
                makePlane( plane, coeffs, 2 );
                // fit and compare
                fitPlane( plane.ptr(), pSize, pSize, fit_coeffs, &chisq );
                BOOST_CHECK_SMALL( fit_coeffs[0], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs[1], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[2], pSize/2, EPS );
                
                // make a plane with slope in both x/y.
                coeffs[0] = 1;
                makePlane( plane, coeffs, 2 );
                // fit and compare
                fitPlane( plane.ptr(), pSize, pSize, fit_coeffs, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs[0], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs[1], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[2], pSize, EPS );
                
                // make a plane with slope in both x/y and a specified mean value.
                coeffs[0] = 3.333;
                coeffs[1] = 12.34;
                coeffs[2] = -123;
                makePlane( plane, coeffs, 3 );
                // fit and compare
                fitPlane( plane.ptr(), pSize, pSize, fit_coeffs, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs[0], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs[1], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[2], coeffs[2], EPS );
                
                // Check that the version returning a plane also gives right results.
                double fit_coeffs2[3] = { 0, 0, 0 };
                Array<float> plane2 = fitPlane( plane, false, fit_coeffs2, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs2[0], coeffs[0], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs2[1], coeffs[1], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs2[2], coeffs[2], EPS );
 
               // Verification that the returned plane is accurate
                fitPlane( plane2.ptr(), pSize, pSize, fit_coeffs, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs[0], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs[1], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[2], coeffs[2], EPS );
                
                // same, but subtract average.
                plane2 = fitPlane( plane, true, fit_coeffs2, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs2[0], coeffs[0], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs2[1], coeffs[1], EPS );
                BOOST_CHECK_SMALL( fit_coeffs2[2], EPS );
                
                // Verification that the returned plane is accurate
                fitPlane( plane2.ptr(), pSize, pSize, fit_coeffs, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs[0], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs[1], EPS );
                BOOST_CHECK_SMALL( fit_coeffs[2], EPS );
                
                // 2 superposed planes, separated by masking
                double coeffs2[3] = { -1.23, 9.45, -23.456 };
                makePlane( plane2, coeffs2, 3 );
                size_t nPixels = pSize*pSize;
                float* pPtr = plane.ptr();
                float* pPtr2 = plane2.ptr();
                Array<bool> mask1(pSize,pSize);
                Array<float> mask2(pSize,pSize);
                bool* mPtr1 = mask1.ptr();
                float* mPtr2 = mask2.ptr();
                for( size_t i(0); i<nPixels; ++i ) {
                    if( i%4 ) {     // select some subset to be replaced with plane2
                        pPtr[i] = pPtr2[i];
                        mPtr2[i] = 1.0;
                        mPtr1[i] = false;
                    } else {
                        mPtr2[i] = 0.0;
                        mPtr1[i] = true;
                    }
                }
                // fit using mask1, to get plane
                fitPlane( pPtr, pSize, pSize, mPtr1, fit_coeffs, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs[0], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs[1], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[2], coeffs[2], EPS );
                // fit using mask2, to get plane
                fitPlane( pPtr, pSize, pSize, mPtr2, fit_coeffs, &chisq );
                BOOST_CHECK_CLOSE( fit_coeffs[0], coeffs2[0], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[1], coeffs2[1], EPS );
                BOOST_CHECK_CLOSE( fit_coeffs[2], coeffs2[2], EPS );
                
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
