

#include "redux/util/endian.hpp"
#include "redux/util/datautil.hpp"

#include <boost/test/unit_test.hpp>

using namespace redux::util;

using namespace std;
using namespace boost::unit_test_framework;

namespace testsuite {

    namespace util {

        template <typename T, typename U=T>
        void testRestrict( void ) {
            
            // within allowed range
            BOOST_CHECK_EQUAL( (restrict<T,U>( 0, 5 )), T(0) );
            BOOST_CHECK_EQUAL( (restrict<T,U>( 5, 5 )), T(5) );     // high endpoint maps back to low side
            BOOST_CHECK_EQUAL( (restrict<T,U>( 3, 5 )), T(3) );     // inside range -> no change
            
            // outside range
            BOOST_CHECK_EQUAL( (restrict<T,U>( 6, 5 )), T(5) );
            BOOST_CHECK_EQUAL( (restrict<T,U>( 1, 5, 2 )), T(2) );

            // outside reversed range
            BOOST_CHECK_EQUAL( (restrict<T,U>( 6, 0, 5 )), T(5) );
            BOOST_CHECK_EQUAL( (restrict<T,U>( 1, 2, 5 )), T(2) );

        }
        
        template <typename T, typename U=T>
        void testRestrictPeriodic( void ) {
            
            // within allowed range
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 0, 5 )), T(0) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 5, 5 )), T(0) );    // high endpoint maps back to low side
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 3, 5 )), T(3) );    // inside range -> no change
            
            // outside range
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 6, 5 )), T(1) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 1, 5, 2 )), T(4) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 13, 5 )), T(3) );   // same as 13 mod 5
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 6, 0, 5 )), T(1) ); // repeat the above for reversed range
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 1, 2, 5 )), T(4) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 13, 0, 5 )), T(3) );
            
            // same but for shifted ranges
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 5, 10, 5 )), T(5) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 10, 10, 5 )), T(5) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 7, 10, 5 )), T(7) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 11, 10, 5 )), T(6) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 4, 10, 5 )), T(9) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 5, 5, 10 )), T(5) ); // repeat the above for reversed range
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 10, 5, 10 )), T(5) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 7, 5, 10 )), T(7) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 11, 5, 10 )), T(6) );
            BOOST_CHECK_EQUAL( (restrict_periodic<T,U>( 4, 5, 10 )), T(9) );
            
        }
        
        template <typename T, typename U=T>
        void testRestrictReflected( void ) {
            
            // within allowed range
            BOOST_CHECK_EQUAL( (restrict_reflected<T,U>( 0, 5 )), T(0) );
            BOOST_CHECK_EQUAL( (restrict_reflected<T,U>( 5, 5 )), T(5) );     // high endpoint maps back to low side
            BOOST_CHECK_EQUAL( (restrict_reflected<T,U>( 3, 5 )), T(3) );     // inside range -> no change
            BOOST_CHECK_EQUAL( (restrict_reflected<T,U>( 13, 5 )), T(3) );    // same as 13 mod 5
            
            // outside range
            BOOST_CHECK_EQUAL( (restrict_reflected<T,U>( 6, 5 )), T(4) );
            BOOST_CHECK_EQUAL( (restrict_reflected<T,U>( 1, 5, 2 )), T(3) );
            
            // same but for shifted ranges
            BOOST_CHECK_EQUAL( (restrict_reflected<T,U>( 5, 10, 5 )), T(5) );
            BOOST_CHECK_EQUAL( (restrict_reflected<T,U>( 10, 10, 5 )), T(10) );
            BOOST_CHECK_EQUAL( (restrict_reflected<T,U>( 7, 10, 5 )), T(7) );
            BOOST_CHECK_EQUAL( (restrict_reflected<T,U>( 11, 10, 5 )), T(9) );
            BOOST_CHECK_EQUAL( (restrict_reflected<T,U>( 4, 10, 5 )), T(6) );
            
        }
        
        void dataTest( void ) {
            
            {   // test the Abs template for some types
                uint8_t tmpu8( rand() );
                BOOST_CHECK_EQUAL( Abs(tmpu8), tmpu8 );
                uint16_t tmpu16( rand() );
                BOOST_CHECK_EQUAL( Abs(tmpu16), tmpu16 );
                uint32_t tmpu32( rand() );
                BOOST_CHECK_EQUAL( Abs(tmpu32), tmpu32 );
                uint64_t tmpu64( rand() );
                BOOST_CHECK_EQUAL( Abs(tmpu64), tmpu64 );
                
                int8_t tmp8(122);
                BOOST_CHECK_EQUAL( Abs(-tmp8), tmp8 );
                int16_t tmp16(5432);
                BOOST_CHECK_EQUAL( Abs(-tmp16), tmp16 );
                int32_t tmp32(65422345);
                BOOST_CHECK_EQUAL( Abs(-tmp32), tmp32 );
                int64_t tmp64(23452524523);
                BOOST_CHECK_EQUAL( Abs(-tmp64), tmp64 );
                
                float tmpF(1.234567E-13);
                BOOST_CHECK_EQUAL( Abs(-tmpF), tmpF );
                double tmpD(1.234567E-54);
                BOOST_CHECK_EQUAL( Abs(-tmpD), tmpD );
                
            }
        
            // test the restrict functions
            
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
            

            
            uint32_t a[10], b[10], aa[10], bb[10];
            for( uint16_t i( 0 ); i < 10; i++ ) {
                a[i] = aa[i] = rand();
                b[i] = bb[i] = rand();
            }

            // check the array wrapper for std::swap
            swap( a, b, 10 );
            for( uint16_t i( 0 ); i < 10; i++ ) {
                BOOST_CHECK_EQUAL( b[i], aa[i] );
                BOOST_CHECK_EQUAL( a[i], bb[i] );
            }

        }


        void endianTest( void ) {

            union {
                uint32_t i;
                uint8_t c[4];
                bool b[4];
            } a, aa;

            a.i = 1;
            BOOST_CHECK( a.b[0] == (RDX_BYTE_ORDER == RDX_LITTLE_ENDIAN) );

            for( int i(0); i < 4; i++ ) {
                a.c[i] = aa.c[i] = i;
            }
            swapEndian( a.i );
            for( int i(0); i < 4; i++ ) {
                BOOST_CHECK_EQUAL( a.c[i], 3 - i ); // verify
            }

            uint32_t b[] = { aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i };
            swapEndian( b, 10 );
            for( int i(0); i < 10; i++ ) {
                BOOST_CHECK_EQUAL( b[i], a.i ); // verify
            }

        }


        void add_data_tests( test_suite* ts ) {
            
            srand (time(NULL));
            
            ts->add( BOOST_TEST_CASE_NAME( &dataTest, "Data tools"  ) );
            ts->add( BOOST_TEST_CASE_NAME( &endianTest, "Endian tests"  ) );

        }

    }

}

