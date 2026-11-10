#include "redux/math/helpers.hpp"

#include "catch.hpp"

#include <cstring>
#include <limits>
#include <memory>

using namespace redux::math;

namespace testsuite {

    namespace math {

        template <typename T>
        void minMaxMeanTest( size_t n ) {
            std::unique_ptr<T[]> data( new T[n] );
            T* ptr = data.get();
            CHECK( std::isnan( mean( ptr, 0 ) ) );
            memset( ptr, 0, n*sizeof(T) );
            CHECK( min(ptr,n) == 0 );
            CHECK( max(ptr,n) == 0 );
            CHECK( mean(ptr,n) == Approx(0.0).margin(1E-20) );
            ptr[0] = std::numeric_limits<T>::lowest();
            ptr[n-1] = std::numeric_limits<T>::max();
            CHECK( min(ptr,n) == ptr[0] );
            CHECK( max(ptr,n) == ptr[n-1] );
            T mn = 13;
            for( unsigned int i=0; i<n; ++i ) ptr[i] = mn + (i%2)*2 - 1;
            if( n%2 ) ptr[0] = mn;
            CHECK( mean(ptr,n) == Approx( (double)mn ).epsilon(1E-20) );
            if( n > 1 ) {
                double rms(0), stddev(0);
                for( unsigned int i=0; i<n; ++i ) {
                    rms += ptr[i]*ptr[i];
                    stddev += (ptr[i]-mn)*(ptr[i]-mn);
                }
                rms = sqrt(rms/n);
                stddev = sqrt(stddev/(n-1));
                double rms2, stddev2;
                CHECK( rmsStddev(ptr,n,rms2,stddev2) == Approx( (double)mn ).epsilon(1E-20) );
                CHECK( rms2 == Approx(rms).epsilon(1E-20) );
                CHECK( stddev2 == Approx(stddev).epsilon(1E-20) );
                rmsStddev(ptr,n,mn,rms2,stddev2);
                CHECK( rms2 == Approx(rms).epsilon(1E-20) );
                CHECK( stddev2 == Approx(stddev).epsilon(1E-20) );
            }
            T mn2, mx;
            double avg;
            minMaxMean( ptr, n, mn2, mx, avg );
            CHECK( mn2 == mn-1 );
            CHECK( mx == mn+1 );
            CHECK( avg == Approx( (double)mn ).epsilon(1E-20) );
        }

        TEST_CASE( "min/max/mean across integer and float types", "[math]" ) {
            SECTION( "int" )      { minMaxMeanTest<int>(42); }
            SECTION( "int8_t" )   { minMaxMeanTest<int8_t>(150); }
            SECTION( "uint8_t" )  { minMaxMeanTest<uint8_t>(111); }
            SECTION( "float" )    { minMaxMeanTest<float>(77); }
            SECTION( "double" )   { minMaxMeanTest<double>(1008); }
            SECTION( "size_t" )   { minMaxMeanTest<size_t>(1009); }
        }

        TEST_CASE( "almostEqual: exact and trivial cases", "[math]" ) {
            CHECK( almostEqual( 0, 0 ) );
            CHECK( almostEqual( 0, 0, 0 ) );
            CHECK( almostEqual( 0, 0, 5 ) );
            CHECK_FALSE( almostEqual( 1.0, 0 ) );
            CHECK_FALSE( almostEqual( 0, 1.0 ) );
        }

        TEST_CASE( "almostEqual: ULP tolerance at a realistic magnitude", "[math]" ) {
            CHECK( almostEqual( 10000.0, 10000.0 ) );
            CHECK( almostEqual( 10000.0, 10000.0, 0 ) );
            CHECK_FALSE( almostEqual( 10000.0f, 10000.000977f, 0 ) );
            CHECK( almostEqual( 10000.0f, 10000.000977f, 1 ) );
            CHECK_FALSE( almostEqual( -10000.0f, -10000.000977f, 0 ) );
            CHECK( almostEqual( -10000.0f, -10000.000977f, 1 ) );
        }

    }

}
