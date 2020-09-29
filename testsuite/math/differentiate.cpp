#include "redux/math/differentiate.hpp"
#include "redux/util/bitoperations.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/test/unit_test.hpp>

using namespace redux::math;
using namespace redux::util;


namespace testsuite {

    namespace math {
        
        typedef boost::numeric::ublas::vector<double> Vector;
        typedef boost::numeric::ublas::matrix<double> Matrix;

        Vector complicatedFunction(const Vector &v) {
            double x = v(0), y = v(1), z = v(2);
            Vector out(6);

            out(0) = x;
            out(1) = y;
            out(2) = z;
            out(3) = x*x + y*y + z*z;
            out(4) = cos(x) + sin(y) + atan(z);
            out(5) = cos(x)*sin(y)*atan(z);
            
            return out;
        }

        Matrix complicatedJacobian(const Vector &v) {
            
            double x = v(0), y = v(1), z = v(2);

            Matrix j(6, 3);
            j(0,0) = 1; j(0,1) = 0; j(0,2) = 0;
            j(1,0) = 0; j(1,1) = 1; j(1,2) = 0;
            j(2,0) = 0; j(2,1) = 0; j(2,2) = 1;

            j(3,0) = 2*x; j(3,1) = 2*y; j(3,2) = 2*z;

            j(4,0) = -sin(x); j(4,1) = cos(y); j(4,2) = 1/(1+z*z);

            j(5,0) = -sin(x)*sin(y)*atan(z);
            j(5,1) = cos(x)*cos(y)*atan(z);
            j(5,2) = cos(x)*sin(y)*1/(1+z*z);

            return j;
        }


        void numericDifferentiatorTest(void) {
            NumericDifferentiator nd(complicatedFunction);

            srand( redux::util::mix(clock(), time(NULL), getpid()) );
            for (size_t i = 0; i < 100; ++i) {
                Vector xyz(3);
                xyz(0) = ((rand() % 1001) - 500) / 10.;
                xyz(1) = ((rand() % 1001) - 500) / 10.;
                xyz(2) = ((rand() % 1001) - 500) / 10.;

                Matrix m1 = nd(xyz);
                Matrix m2 = complicatedJacobian(xyz);

                BOOST_REQUIRE_EQUAL(m1.size1(), m2.size1());
                BOOST_REQUIRE_EQUAL(m1.size2(), m2.size2());
                for (size_t r = 0; r < m1.size1(); ++r) {
                    for (size_t c = 0; c < m1.size2(); ++c) {
                        BOOST_CHECK_CLOSE(m1(r,c), m2(r,c), 0.001);
                    }
                }
            }
        }

        
        using namespace boost::unit_test;
        void add_differentiate_tests( test_suite* ts ) {
            
            // Numerical differentiation
            ts->add( BOOST_TEST_CASE_NAME( &numericDifferentiatorTest, "Numerical differentiation" ) );


        }

    }

}


