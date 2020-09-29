

#include "redux/util/bitoperations.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/boundvalue.hpp"
#include "redux/util/datautil.hpp"

#include <boost/test/unit_test.hpp>

using namespace redux::util;

using namespace std;
using namespace boost::unit_test_framework;

namespace testsuite {

    namespace util {


        uint64_t naiveBitCount( uint64_t v ) {
            uint64_t cnt( 0 );
            do {
                if( v & 1 ) cnt++;
            }
            while( v >>= 1 );
            return cnt;
        }



        void bitTest( void ) {

            // deBruijn based log2
            for( uint8_t i( 1 ); i < 32; ++i ) {
                BOOST_CHECK_EQUAL( redux::util::log2( ( 1 << i ) + 1 ), i );
            }

            // deBruijn LSB detector
            for( uint8_t i( 0 ); i < 32; ++i ) {
                BOOST_CHECK_EQUAL( findLSB( 1 << i ), i );
            }

            // deBruijn LSB detector, 1-based index
            for( uint8_t i( 0 ); i < 32; ++i ) {
                BOOST_CHECK_EQUAL( findLSB1( 1 << i ), i + 1 );
            }

            // 64-bit deBruijn LSB detector, 1-based index
            uint64_t tmp64( 1 );
            for( uint8_t i( 0 ); i < 64; ++i ) {
                BOOST_CHECK_EQUAL( findLSB64( tmp64 << i ), i + 1 );
            }

            // next power of two
            for( uint8_t i( 0 ); i < 31; ++i ) {
                BOOST_CHECK_EQUAL( nextPowerOfTwo( ( 1 << i ) + 1 ), ( 1 << ( i + 1 ) ) );
            }

            uint32_t mask( ( 1 << 16 ) - 1 );
            // swap first 2 bytes of a.i[0] with the first 2 bytes of a.i[1]
            union {
                uint64_t l;
                uint32_t i[2];
                uint16_t s[4];
                uint8_t b[8];
            } a, b;
            for( uint16_t i( 0 ); i < 10; i++ ) {
                a.i[0] = b.i[0] = rand();
                a.i[1] = b.i[1] = rand();
                swapBits( a.i[0], a.i[1], mask );
                BOOST_CHECK_EQUAL( a.s[0], b.s[2] ); // verify
                BOOST_CHECK_EQUAL( a.s[2], b.s[0] );
                BOOST_CHECK_EQUAL( a.s[1], b.s[1] ); // check that the unmasked bits are unchanged.
                BOOST_CHECK_EQUAL( a.s[3], b.s[3] );
            }

            // swap first 2 bytes of each integer in array with matching bytes in another array
            uint32_t aa[10], bb[10], cc[10];
            memset( bb, 0, 10 * sizeof( uint32_t ) );
            for( uint16_t i( 0 ); i < 10; i++ ) {
                aa[i] = cc[i] = rand();
            }
            swapBits( aa, bb, mask, 10 );
            for( uint16_t i( 0 ); i < 10; i++ ) {
                BOOST_CHECK_EQUAL( cc[i] & ( mask << 16 ), aa[i] );
                BOOST_CHECK_EQUAL( cc[i]&mask, bb[i] );
            }

            memcpy( aa, cc, 10 * sizeof( uint32_t ) );
            // flip the first 16 bits in aa[10]
            flipBits( aa, mask, 10 );
            for( uint16_t i( 0 ); i < 10; i++ ) {
                BOOST_CHECK_EQUAL( cc[i] ^ mask, aa[i] );
            }

            for( int i( 0 ); i < 15; ++i ) {
                a.i[0] = rand();
                a.i[1] = rand();
                //BOOST_CHECK_EQUAL( countBits(a.b[0]), naiveBitCount(b.s[0]) );
                //BOOST_CHECK_EQUAL( countBits(a.s[0]), naiveBitCount(a.s[0]) );
                BOOST_CHECK_EQUAL( countBits( a.i[0] ), naiveBitCount( a.i[0] ) );
                BOOST_CHECK_EQUAL( countBits( a.l ), naiveBitCount( a.l ) );
            }


        }


        void boundValueTest( void ) {

            typedef BoundValue<int, redux::util::detail::UNDEFINEDTRIM> invalidBV;
            BOOST_CHECK_THROW( invalidBV(), invalid_argument );

            BoundValue<int> intVal( 3, -9, 9 );                 // ordinary construction / assignment
            BOOST_CHECK_EQUAL( intVal, 3 );
            BOOST_CHECK_EQUAL( intVal.min(), -9 );
            BOOST_CHECK_EQUAL( intVal.max(),  9 );
            intVal = 14;
            BOOST_CHECK_EQUAL( intVal, 9 );
            intVal = -14.0;                                     // test implicit arithmetic cast.
            BOOST_CHECK_EQUAL( intVal, -9 );

            BoundValue<double> dblVal( 3, 9, -9 );              // reverse min/max order
            BOOST_CHECK_EQUAL( dblVal, 3.0 );
            BOOST_CHECK_EQUAL( dblVal.min(), -9.0 );
            BOOST_CHECK_EQUAL( dblVal.max(),  9.0 );
            dblVal = 14;
            BOOST_CHECK_EQUAL( dblVal, 9.0 );
            dblVal = ( int ) - 14;                              // test implicit arithmetic cast.
            BOOST_CHECK_EQUAL( dblVal, -9.0 );

            BoundValue<double> dblVal2( dblVal );               // copy construction of the same type
            BOOST_CHECK_EQUAL( dblVal2, -9.0 );
            BOOST_CHECK_EQUAL( dblVal2.min(), -9.0 );
            BOOST_CHECK_EQUAL( dblVal2.max(),  9.0 );

            BoundValue<int> intVal2( dblVal );                  // copy construction from different type
            BOOST_CHECK_EQUAL( intVal2, -9 );
            BOOST_CHECK_EQUAL( intVal2.min(), -9 );
            BOOST_CHECK_EQUAL( intVal2.max(),  9 );


            BoundValue<int, redux::util::detail::WRAP> wrappedVal( 0, 0, 5 ); // range [0,5)  ->  equivalent to modulo 5
            wrappedVal = 6;
            BOOST_CHECK_EQUAL( wrappedVal, 1 );
            wrappedVal = -1;
            BOOST_CHECK_EQUAL( wrappedVal, 4 );
            wrappedVal = -126;
            BOOST_CHECK_EQUAL( wrappedVal, 4 );
            wrappedVal = 126;
            BOOST_CHECK_EQUAL( wrappedVal, 1 );


            double delta( 0.05 );
            double epsilon = 1E-9;

            BoundValue<double, redux::util::detail::WRAP> periodicVal( 0, 0, 2 * M_PI );
            periodicVal = 3 * M_PI + delta;
            BOOST_CHECK_CLOSE( ( double )periodicVal, M_PI + delta, epsilon );
            periodicVal = -3 * M_PI / 2 + delta;
            BOOST_CHECK_CLOSE( ( double )periodicVal, M_PI / 2 + delta, epsilon );


            BoundValue<double, redux::util::detail::REFLECT> inclVal( 0, 0, M_PI ); // e.g. inclination = PI/2+\delta -> PI/2 - \delta
            inclVal = M_PI + delta;
            BOOST_CHECK_CLOSE( ( double )inclVal, M_PI - delta, epsilon );
            inclVal = 2 * M_PI + delta;
            BOOST_CHECK_CLOSE( ( double )inclVal, delta, epsilon );
            inclVal = 3 * M_PI + delta;
            BOOST_CHECK_CLOSE( ( double )inclVal, M_PI - delta, epsilon );
            inclVal = -delta;
            BOOST_CHECK_CLOSE( ( double )inclVal, delta, epsilon );
            inclVal = -M_PI - delta;
            BOOST_CHECK_CLOSE( ( double )inclVal, M_PI - delta, epsilon );
            inclVal = -2 * M_PI - delta;
            BOOST_CHECK_CLOSE( ( double )inclVal, delta, epsilon );


            BoundValue<int, redux::util::detail::REFLECT> reflectedVal( 0, 0, 4 );
            reflectedVal = 0;
            BOOST_CHECK_EQUAL( reflectedVal, 0 );
            reflectedVal = 5;
            BOOST_CHECK_EQUAL( reflectedVal, 3 );
            reflectedVal = 9;
            BOOST_CHECK_EQUAL( reflectedVal, 1 );
            reflectedVal = 11;
            BOOST_CHECK_EQUAL( reflectedVal, 3 );
            reflectedVal.setMax( 5 );
            reflectedVal = -3;
            BOOST_CHECK_EQUAL( reflectedVal, 3 );
            reflectedVal = -5;
            BOOST_CHECK_EQUAL( reflectedVal, 5 );
            reflectedVal = -9;
            BOOST_CHECK_EQUAL( reflectedVal, 1 );
            reflectedVal = -11;
            BOOST_CHECK_EQUAL( reflectedVal, 1 );

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
            
            ts->add( BOOST_TEST_CASE_NAME( &bitTest, "Bit manipulations"  ) );
            ts->add( BOOST_TEST_CASE_NAME( &boundValueTest, "BoundValue"  ) );
            ts->add( BOOST_TEST_CASE_NAME( &dataTest, "Data tools"  ) );
            ts->add( BOOST_TEST_CASE_NAME( &endianTest, "Endian tests"  ) );

        }

    }

}

