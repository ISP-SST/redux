#include "redux/util/endian.hpp"
#include "redux/util/datautil.hpp"

#include "catch.hpp"

#include <cstdlib>
#include <ctime>

using namespace redux::util;

using namespace std;

namespace testsuite {

    namespace util {

        namespace {

            template <typename T, typename U=T>
            void testRestrict() {

                CHECK( (restrict<T,U>( 0, 5 )) == T(0) );
                CHECK( (restrict<T,U>( 5, 5 )) == T(5) );
                CHECK( (restrict<T,U>( 3, 5 )) == T(3) );

                CHECK( (restrict<T,U>( 6, 5 )) == T(5) );
                CHECK( (restrict<T,U>( 1, 5, 2 )) == T(2) );

                CHECK( (restrict<T,U>( 6, 0, 5 )) == T(5) );
                CHECK( (restrict<T,U>( 1, 2, 5 )) == T(2) );

            }

            template <typename T, typename U=T>
            void testRestrictPeriodic() {

                CHECK( (restrict_periodic<T,U>( 0, 5 )) == T(0) );
                CHECK( (restrict_periodic<T,U>( 5, 5 )) == T(0) );
                CHECK( (restrict_periodic<T,U>( 3, 5 )) == T(3) );

                CHECK( (restrict_periodic<T,U>( 6, 5 )) == T(1) );
                CHECK( (restrict_periodic<T,U>( 1, 5, 2 )) == T(4) );
                CHECK( (restrict_periodic<T,U>( 13, 5 )) == T(3) );
                CHECK( (restrict_periodic<T,U>( 6, 0, 5 )) == T(1) );
                CHECK( (restrict_periodic<T,U>( 1, 2, 5 )) == T(4) );
                CHECK( (restrict_periodic<T,U>( 13, 0, 5 )) == T(3) );

                CHECK( (restrict_periodic<T,U>( 5, 10, 5 )) == T(5) );
                CHECK( (restrict_periodic<T,U>( 10, 10, 5 )) == T(5) );
                CHECK( (restrict_periodic<T,U>( 7, 10, 5 )) == T(7) );
                CHECK( (restrict_periodic<T,U>( 11, 10, 5 )) == T(6) );
                CHECK( (restrict_periodic<T,U>( 4, 10, 5 )) == T(9) );
                CHECK( (restrict_periodic<T,U>( 5, 5, 10 )) == T(5) );
                CHECK( (restrict_periodic<T,U>( 10, 5, 10 )) == T(5) );
                CHECK( (restrict_periodic<T,U>( 7, 5, 10 )) == T(7) );
                CHECK( (restrict_periodic<T,U>( 11, 5, 10 )) == T(6) );
                CHECK( (restrict_periodic<T,U>( 4, 5, 10 )) == T(9) );

            }

            template <typename T, typename U=T>
            void testRestrictReflected() {

                CHECK( (restrict_reflected<T,U>( 0, 5 )) == T(0) );
                CHECK( (restrict_reflected<T,U>( 5, 5 )) == T(5) );
                CHECK( (restrict_reflected<T,U>( 3, 5 )) == T(3) );
                CHECK( (restrict_reflected<T,U>( 13, 5 )) == T(3) );

                CHECK( (restrict_reflected<T,U>( 6, 5 )) == T(4) );
                CHECK( (restrict_reflected<T,U>( 1, 5, 2 )) == T(3) );

                CHECK( (restrict_reflected<T,U>( 5, 10, 5 )) == T(5) );
                CHECK( (restrict_reflected<T,U>( 10, 10, 5 )) == T(10) );
                CHECK( (restrict_reflected<T,U>( 7, 10, 5 )) == T(7) );
                CHECK( (restrict_reflected<T,U>( 11, 10, 5 )) == T(9) );
                CHECK( (restrict_reflected<T,U>( 4, 10, 5 )) == T(6) );

            }

        }

        TEST_CASE( "Abs across signed/unsigned/float types", "[util]" ) {
            srand(time(NULL));

            uint8_t tmpu8( rand() );
            CHECK( Abs(tmpu8) == tmpu8 );
            uint16_t tmpu16( rand() );
            CHECK( Abs(tmpu16) == tmpu16 );
            uint32_t tmpu32( rand() );
            CHECK( Abs(tmpu32) == tmpu32 );
            uint64_t tmpu64( rand() );
            CHECK( Abs(tmpu64) == tmpu64 );

            int8_t tmp8(122);
            CHECK( Abs(-tmp8) == tmp8 );
            int16_t tmp16(5432);
            CHECK( Abs(-tmp16) == tmp16 );
            int32_t tmp32(65422345);
            CHECK( Abs(-tmp32) == tmp32 );
            int64_t tmp64(23452524523);
            CHECK( Abs(-tmp64) == tmp64 );

            float tmpF(1.234567E-13);
            CHECK( Abs(-tmpF) == tmpF );
            double tmpD(1.234567E-54);
            CHECK( Abs(-tmpD) == tmpD );
        }

        TEST_CASE( "restrict/restrict_periodic/restrict_reflected across type combinations", "[util]" ) {

            SECTION( "restrict" ) {
                testRestrict<int8_t,int>();
                testRestrict<uint8_t,int>();
                testRestrict<int16_t,int>();
                testRestrict<uint16_t,int>();
                testRestrict<int,int8_t>();
                testRestrict<int,uint8_t>();
                testRestrict<int,int16_t>();
                testRestrict<int,uint16_t>();
                testRestrict<int,int>();
                testRestrict<int,float>();
                testRestrict<float,int>();
                testRestrict<float,float>();
            }

            SECTION( "restrict_periodic" ) {
                testRestrictPeriodic<int8_t,int>();
                testRestrictPeriodic<uint8_t,int>();
                testRestrictPeriodic<int16_t,int>();
                testRestrictPeriodic<uint16_t,int>();
                testRestrictPeriodic<int,int8_t>();
                testRestrictPeriodic<int,uint8_t>();
                testRestrictPeriodic<int,int16_t>();
                testRestrictPeriodic<int,uint16_t>();
                testRestrictPeriodic<int,int>();
                testRestrictPeriodic<int,float>();
                testRestrictPeriodic<float,int>();
                testRestrictPeriodic<float,float>();
            }

            SECTION( "restrict_reflected" ) {
                testRestrictReflected<int8_t,int>();
                testRestrictReflected<uint8_t,int>();
                testRestrictReflected<int16_t,int>();
                testRestrictReflected<uint16_t,int>();
                testRestrictReflected<int,int8_t>();
                testRestrictReflected<int,uint8_t>();
                testRestrictReflected<int,int16_t>();
                testRestrictReflected<int,uint16_t>();
                testRestrictReflected<int,int>();
                testRestrictReflected<int,float>();
                testRestrictReflected<float,int>();
                testRestrictReflected<float,float>();
            }

        }

        TEST_CASE( "array-wrapper swap()", "[util]" ) {
            srand(time(NULL));

            uint32_t a[10], b[10], aa[10], bb[10];
            for( uint16_t i( 0 ); i < 10; i++ ) {
                a[i] = aa[i] = rand();
                b[i] = bb[i] = rand();
            }

            swap( a, b, 10 );
            for( uint16_t i( 0 ); i < 10; i++ ) {
                CHECK( b[i] == aa[i] );
                CHECK( a[i] == bb[i] );
            }
        }

        TEST_CASE( "endian swap", "[util]" ) {

            union {
                uint32_t i;
                uint8_t c[4];
                bool b[4];
            } a, aa;

            a.i = 1;
            CHECK( a.b[0] == (RDX_BYTE_ORDER == RDX_LITTLE_ENDIAN) );

            for( int i(0); i < 4; i++ ) {
                a.c[i] = aa.c[i] = i;
            }
            swapEndian( a.i );
            for( int i(0); i < 4; i++ ) {
                CHECK( a.c[i] == 3 - i );
            }

            uint32_t b[] = { aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i };
            swapEndian( b, 10 );
            for( int i(0); i < 10; i++ ) {
                CHECK( b[i] == a.i );
            }

        }

    }

}
