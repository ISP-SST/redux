#include <boost/test/unit_test.hpp>
//#include <boost/test/floating_point_comparison.hpp>
//#include <sstream>

using namespace boost::unit_test_framework;

#include "redux/constants.hpp"
#include "redux/math/differentiate.hpp"
#include "redux/math/interval.hpp"
#include "redux/util/bitoperations.hpp"
#include "redux/util/stringutil.hpp"

#include <iostream>

using namespace redux::math;
using namespace redux::util;
using namespace std;

namespace {
    
    template <typename T>
    void typedIntervalTest(void) {

        Interval<T> interval(3, 9);
        BOOST_CHECK_EQUAL(interval.start(), 3);
        BOOST_CHECK_EQUAL(interval.end(), 9);

        BOOST_CHECK(interval.contains(5));
        BOOST_CHECK(!interval.contains(1));
        BOOST_CHECK(!interval.contains(20));

        BOOST_CHECK(interval.surrounds(Interval<T>(4, 8)));
        BOOST_CHECK(!interval.inside(Interval<T>(4, 8)));
        BOOST_CHECK(interval.inside(Interval<T>(1, 10)));
        BOOST_CHECK(!interval.surrounds(Interval<T>(1, 10)));

        BOOST_CHECK(interval.intersects(Interval<T>(0, 3)));
        BOOST_CHECK(interval.intersects(Interval<T>(9, 15)));
        BOOST_CHECK(interval.intersects(Interval<T>(1, 15)));
        BOOST_CHECK(!interval.intersects(Interval<T>(10, 15)));
        BOOST_CHECK(!interval.intersects(Interval<T>(0, 2)));

        BOOST_CHECK((interval == Interval<T>(3, 9)));
        BOOST_CHECK(!(interval == Interval<T>(2, 9)));

        BOOST_CHECK(!(interval != Interval<T>(3, 9)));
        BOOST_CHECK((interval != Interval<T>(2, 9)));

        BOOST_CHECK((interval > Interval<T>(0, 2)));
        BOOST_CHECK(!(interval > Interval<T>(10, 20)));
        BOOST_CHECK(!(interval < Interval<T>(0, 2)));
        BOOST_CHECK((interval < Interval<T>(10, 20)));

        BOOST_CHECK_THROW(interval.join(Interval<T>(10, 20)), std::logic_error);
        BOOST_CHECK_THROW(interval.joined(Interval<T>(10, 20)), std::logic_error);
        BOOST_CHECK_THROW(interval.intersection(Interval<T>(10, 20)), std::logic_error);

        BOOST_CHECK((interval.from(5) == Interval<T>(5, 9)));
        BOOST_CHECK((interval.until(5) == Interval<T>(3, 5)));

        BOOST_CHECK_THROW(interval.from(10), std::logic_error);
        BOOST_CHECK_THROW(interval.until(10), std::logic_error);

    }
    
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
    
    template <typename T>
    void minMaxMeanTest(size_t n) {
        BOOST_CHECK_EQUAL(0,0);
        std::unique_ptr<T[]> data( new T[n] );
        BOOST_CHECK_EQUAL(0,0);
        T* ptr = data.get();
        BOOST_CHECK( std::isnan( mean(ptr,0) ) );
        memset(ptr,0,n*sizeof(T));
        BOOST_CHECK_EQUAL(min(ptr,n),0);
        BOOST_CHECK_EQUAL(max(ptr,n),0);
        BOOST_CHECK_CLOSE(mean(ptr,n),0.0,1E-20);
        ptr[0] = std::numeric_limits<T>::lowest();
        ptr[n-1] = std::numeric_limits<T>::max();
        BOOST_CHECK_EQUAL(min(ptr,n),ptr[0]);
        BOOST_CHECK_EQUAL(max(ptr,n),ptr[n-1]);
        T mn = 13;      // just something that is representable in all integer/float types
        for(int i=0; i<n; ++i) ptr[i] = mn + (i%2)*2 - 1;   // mn ± 1
        if( n%2 ) ptr[0] = mn;                              // if n is odd, set one value to mn
        BOOST_CHECK_CLOSE( mean(ptr,n), mn, 1E-20);
        if( n > 1 ) {       // compare rms/stddev against naĩve implementation
            double rms(0), stddev(0);
            for(int i=0; i<n; ++i) {
                rms += ptr[i]*ptr[i];
                stddev += (ptr[i]-mn)*(ptr[i]-mn);
            }
            rms = sqrt(rms/n);
            stddev = sqrt(stddev/(n-1));
            double rms2, stddev2;
            BOOST_CHECK_CLOSE( rmsStddev(ptr,n,rms2,stddev2), mn, 1E-20);
            BOOST_CHECK_CLOSE( rms2, rms, 1E-20);
            BOOST_CHECK_CLOSE( stddev2, stddev, 1E-20);
            rmsStddev(ptr,n,mn,rms2,stddev2);
            BOOST_CHECK_CLOSE( rms2, rms, 1E-20);
            BOOST_CHECK_CLOSE( stddev2, stddev, 1E-20);

        }
        T mn2,mx;
        double avg;
        minMaxMean(ptr, n, mn2, mx, avg);
        BOOST_CHECK_EQUAL(mn2,mn-1);
        BOOST_CHECK_EQUAL(mx,mn+1);
        BOOST_CHECK_CLOSE( avg, mn, 1E-20);
    }
    
}


void helperTest(void) {
    
    // Basic tests
    BOOST_CHECK( almostEqual(0, 0) );
    BOOST_CHECK( almostEqual(0, 0, 0) );
    BOOST_CHECK( almostEqual(0, 0, 5) );
    BOOST_CHECK( !almostEqual(1.0, 0) );
    BOOST_CHECK( !almostEqual(0, 1.0) );

    // Some more interesting tests
    BOOST_CHECK( almostEqual(10000.0, 10000.0) );
    BOOST_CHECK( almostEqual(10000.0, 10000.0, 0) );
    BOOST_CHECK( !almostEqual(10000.0f, 10000.000977f, 0) );      // One float tic away
    BOOST_CHECK( almostEqual(10000.0f, 10000.000977f, 1) );       // One float tic away
    BOOST_CHECK( !almostEqual(-10000.0f, -10000.000977f, 0) );    // One float tic away
    BOOST_CHECK( almostEqual(-10000.0f, -10000.000977f, 1) );     // One float tic away
    
    // min, max, mean of some arbitrary types/sizes
    minMaxMeanTest<int>(42);
    minMaxMeanTest<int8_t>(150);
    minMaxMeanTest<uint8_t>(111);
    minMaxMeanTest<float>(77);
    minMaxMeanTest<double>(1008);
    minMaxMeanTest<size_t>(1009);
    
    
}


void rangeTest(void) {
    
    // try for signed/unsigned/floats
    typedIntervalTest<uint8_t>();
    typedIntervalTest<int>();
    typedIntervalTest<size_t>();
    typedIntervalTest<double>();

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

namespace testsuite {

    namespace math {

        void mathTest(void) {
            test_suite* ts = BOOST_TEST_SUITE("MATH");

            ts->add(BOOST_TEST_CASE(&helperTest));
            ts->add(BOOST_TEST_CASE(&rangeTest));
            ts->add(BOOST_TEST_CASE(&numericDifferentiatorTest));

            framework::master_test_suite().add(ts);
        }

    }

}
