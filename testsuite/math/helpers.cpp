#include "redux/math/helpers.hpp"

#include <boost/test/unit_test.hpp>

using namespace redux::math;


namespace testsuite {

    namespace math {
        
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
            for(unsigned int i=0; i<n; ++i) ptr[i] = mn + (i%2)*2 - 1;   // mn ± 1
            if( n%2 ) ptr[0] = mn;                              // if n is odd, set one value to mn
            BOOST_CHECK_CLOSE( mean(ptr,n), mn, 1E-20);
            if( n > 1 ) {       // compare rms/stddev against naĩve implementation
                double rms(0), stddev(0);
                for(unsigned int i=0; i<n; ++i) {
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


        using namespace boost::unit_test;
        void add_helpers_tests( test_suite* ts ) {
            
            ts->add( BOOST_TEST_CASE_NAME( &helperTest, "Math utilities" ) );


        }
        
    }

}

