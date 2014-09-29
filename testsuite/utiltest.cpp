#include <boost/test/unit_test.hpp>


#include "redux/constants.hpp"
#include "redux/util/array.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/bitoperations.hpp"
#include "redux/util/boundvalue.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <iostream>
#include <memory>

using namespace redux::util;

using namespace std;
using namespace boost::unit_test_framework;

namespace {

    uint64_t naiveCount( uint64_t v ) {
        uint64_t cnt( 0 );
        do {
            if( v & 1 ) cnt++;
        }
        while( v >>= 1 );
        return cnt;
    }

}

void arrayTest( void ) {
    
    Array<int> array4x5( 4, 5 );
    Array<int> array3x3;

    // check that set throws for nValues != nElements;
    BOOST_CHECK_THROW( array4x5.set( 1, 2, 3 ), logic_error );
    BOOST_CHECK_THROW( array4x5.set( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 ), logic_error );

    array4x5.set( 1,  2,  3,  4,  5.0,   // all values will be cast to int in set
                  6,  7,  8,  9, 10,
                  11, 12, 13, 14, 15,
                  16, 17, 18, 19, 20 );

    // check the iterator is accessing the right elements.
    Array<int>::iterator it = array4x5.begin();
    BOOST_CHECK_EQUAL( it.pos(), 0 );
    it = array4x5.end();
    BOOST_CHECK_EQUAL( it.pos(), 20 );
    BOOST_CHECK_EQUAL( *--it, 20 );
    it = array4x5.pos(1,1);
    BOOST_CHECK_EQUAL( it.pos(), 6 );
    BOOST_CHECK_EQUAL( it.step(1).pos(), 7 );
    BOOST_CHECK_EQUAL( it.step(-1).pos(), 6 );
    BOOST_CHECK_EQUAL( it.step(1,1).pos(), 12 );
    BOOST_CHECK_EQUAL( it.step(-1,-1).pos(), 6 );
    BOOST_CHECK_EQUAL( it.step(1,-1).pos(), 10 );
    BOOST_CHECK_EQUAL( it.step(-1,1).pos(), 6 );
   
    // check index out of bounds when using at()
    BOOST_CHECK_THROW( array4x5.at( 1, 5 ), out_of_range );
    BOOST_CHECK_THROW( array4x5.at( 5, 1 ), out_of_range );
    BOOST_CHECK_THROW( array4x5.at( -1, 1 ), out_of_range );
    BOOST_CHECK_THROW( array4x5.at( 1, -1 ), out_of_range );
    
    // test iterator returned by Array::pos() for some values
    BOOST_CHECK_EQUAL( *array4x5.pos( 0, 0 ), 1 );
    BOOST_CHECK_EQUAL( *array4x5.pos( 1, 1 ), 7 );
    BOOST_CHECK_EQUAL( *array4x5.pos( 3, 4 ), 20 );

   
    // check values
    for( int i = 0; i < 4; ++i ) {
        for( int j = 0; j < 5; ++j ) {
            BOOST_CHECK_EQUAL( array4x5( i, j ), 5 * i + j + 1 );
        }
    }

    // test iterators
    it = array4x5.begin();
    Array<int>::iterator endit = array4x5.end();
    Array<int>::const_iterator cit = array4x5.begin();
    Array<int>::const_iterator cendit = array4x5.end();

    // prefix ++
    for( int i = 1; i <= 20; ++i ) {
        BOOST_CHECK_EQUAL( *it, i );
        BOOST_CHECK_EQUAL( *cit, i );
        ++it;
        ++cit;
    }

    // at end()
    BOOST_CHECK( it == endit );
    BOOST_CHECK( cit == cendit );


    // prefix --
    for( int i = 20; i > 0; --i ) {
        BOOST_CHECK_EQUAL( *--it, i );
        BOOST_CHECK_EQUAL( *--cit, i );
    }

    // postfix ++
    for( int i = 1; i <= 20; ++i ) {
        BOOST_CHECK_EQUAL( *it++, i );
        BOOST_CHECK_EQUAL( *cit++, i );
    }

    // at end()
    BOOST_CHECK( it == endit );
    BOOST_CHECK( cit == cendit );

    // postfix --
    for( int i = 20; i > 0; --i ) {
        it--;
        cit--;
        BOOST_CHECK_EQUAL( *it, i );
        BOOST_CHECK_EQUAL( *cit, i );
    }

    // test with std auto iterator (by reference)
    int i = 1;
    for( auto & ait : array4x5 ) {
        BOOST_CHECK_EQUAL( ait, i++ );
    }

    // shallow copy (shared data)
    array3x3 = array4x5;
    BOOST_CHECK( array3x3 == array4x5 );
    // verify that data is shared
    BOOST_CHECK( array3x3.ptr() == array4x5.ptr() );

    // deep copy
    array3x3 = array4x5.copy();
    BOOST_CHECK( array3x3 == array4x5 );
    // verify that data is not shared
    BOOST_CHECK( array3x3.ptr() != array4x5.ptr() );

    // test sub-array constructor
    Array<int> subarray( array4x5, 1, 3, 1, 3 ); // 3x3 sub-array, with offset (1,1)

    // verify that data is shared
    BOOST_CHECK( subarray.ptr() == array4x5.ptr() );
    // define an array equal to the center 3x3 elements of the 5x5 array above for reference
    array3x3.resize( 3, 3 );
    array3x3.set( 7,  8,  9,
                  12, 13, 14,
                  17, 18, 19 );
    // verify that data matches
    BOOST_CHECK( subarray == array3x3 );


    // test subiterator::step() to step along any dimension
    BOOST_CHECK_EQUAL( *(subarray.pos( 1, 1 ).step()), 14 );
    BOOST_CHECK_EQUAL( *(subarray.pos( 1, 1 ).step(1)), 14 );
    BOOST_CHECK_EQUAL( *(subarray.pos( 1, 1 ).step(1,-1)), 17 );
    BOOST_CHECK_EQUAL( *(subarray.pos( 1, 1 ).step(-1,-1)), 7 );
    BOOST_CHECK_EQUAL( *(subarray.pos( 1, 1 ).step(-1,1)), 9 );

    // test assigning to sub-array
    array3x3.set( 70,  80,  90,
                  120, 130, 140,
                  170, 180, 190 );
    subarray.set( array3x3 );

    // array4x5 should now be:
    //     1,   2,   3,   4,  5,
    //     6,  70,  80,  90, 10,
    //    11, 120, 130, 140, 15,
    //    16, 170, 180, 190, 20 );
    //
    // check values
    for( int i = 0; i < 4; ++i ) {
        for( int j = 0; j < 5; ++j ) {
            int val = 5 * i + j + 1;
            if( i >= 1 && i <= 3 && j >= 1 && j <= 3 ) val *= 10;
            BOOST_CHECK_EQUAL( array4x5( i, j ), val );
        }
    }

    Array<int> array4x5x6( 4, 5, 6 );
    size_t cnt( 0 );
    for( auto & it : array4x5x6 ) {
        it = ++cnt;
    }
    BOOST_CHECK_EQUAL( cnt, 120 );
    
    {
        // arbitrary subarray.
        Array<int> subarray( array4x5x6, 1, 2, 2, 3, 3, 4 );
        
        // check that the subarray is accessing the right elements.
        for( size_t i = 0; i < 2; ++i ) {
            for( size_t j = 0; j < 2; ++j ) {
                for( size_t k = 0; k < 2; ++k ) {
                    BOOST_CHECK_EQUAL( array4x5x6( i + 1, j + 2, k + 3 ), subarray( i, j, k ) );
                }
            }
        }

        // assign to subarray
        subarray = 999;
        cnt = 0;
        for( auto cit=subarray.begin(); cit != subarray.end(); ++cit ) {
            BOOST_CHECK_EQUAL( *cit, 999 );
            cnt++;
        }
        BOOST_CHECK_EQUAL( cnt, 8 );
        
        cnt = 0;
        for( auto cit=subarray.end(); cit != subarray.begin();  ) {
            BOOST_CHECK_EQUAL( *--cit, 999 );
            cnt++;
        }
        BOOST_CHECK_EQUAL( cnt, 8 );

    }
    
    {
        // test access as raw multidimensional array
        cnt=0;
        for( auto & it : array4x5x6 ) {
            it = ++cnt;
        }
        shared_ptr<int**> raiiArray = array4x5x6.get(4,5,6);
        int*** rawArray = raiiArray.get();
        cnt = 0;
        for( size_t i=0; i<4; ++i) {
            for( size_t j=0; j<5; ++j) {
                for( size_t k=0; k<6; ++k) {
                    BOOST_CHECK_EQUAL( rawArray[i][j][k], ++cnt );
                }
            }
        }
    }
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
        //BOOST_CHECK_EQUAL( countBits(a.b[0]), naiveCount(b.s[0]) );
        //BOOST_CHECK_EQUAL( countBits(a.s[0]), naiveCount(a.s[0]) );
        BOOST_CHECK_EQUAL( countBits( a.i[0] ), naiveCount( a.i[0] ) );
        BOOST_CHECK_EQUAL( countBits( a.l ), naiveCount( a.l ) );
    }


}


void boundValueTest( void ) {

    typedef BoundValue<int, detail::UNDEFINEDTRIM> invalidBV;
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


    BoundValue<int, detail::WRAP> wrappedVal( 0, 0, 5 ); // range [0,5)  ->  equivalent to modulo 5
    wrappedVal = 6;
    BOOST_CHECK_EQUAL( wrappedVal, 1 );
    wrappedVal = -1;
    BOOST_CHECK_EQUAL( wrappedVal, 4 );
    wrappedVal = -126;
    BOOST_CHECK_EQUAL( wrappedVal, 4 );
    wrappedVal = 126;
    BOOST_CHECK_EQUAL( wrappedVal, 1 );


    using redux::PI;
    double delta( 0.05 );
    double epsilon = 1E-9;

    BoundValue<double, detail::WRAP> periodicVal( 0, 0, 2 * PI );
    periodicVal = 3 * PI + delta;
    BOOST_CHECK_CLOSE( ( double )periodicVal, PI + delta, epsilon );
    periodicVal = -3 * PI / 2 + delta;
    BOOST_CHECK_CLOSE( ( double )periodicVal, PI / 2 + delta, epsilon );


    BoundValue<double, detail::REFLECT> inclVal( 0, 0, PI ); // e.g. inclination = PI/2+\delta -> PI/2 - \delta
    inclVal = PI + delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, PI - delta, epsilon );
    inclVal = 2 * PI + delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, delta, epsilon );
    inclVal = 3 * PI + delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, PI - delta, epsilon );
    inclVal = -delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, delta, epsilon );
    inclVal = -PI - delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, PI - delta, epsilon );
    inclVal = -2 * PI - delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, delta, epsilon );


    BoundValue<int, detail::REFLECT> reflectedVal( 0, 0, 4 );
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
    BOOST_CHECK( a.b[0] == ( REDUX_BYTE_ORDER == REDUX_LITTLE_ENDIAN ) );

    for( int i( 0 ); i < 4; i++ ) {
        a.c[i] = aa.c[i] = i;
    }
    swapEndian( a.i );
    for( int i( 0 ); i < 4; i++ ) {
        BOOST_CHECK_EQUAL( a.c[i], 3 - i ); // verify
    }

    uint32_t b[] = { aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i };
    swapEndian( b, 10 );
    for( int i( 0 ); i < 4; i++ ) {
        BOOST_CHECK_EQUAL( b[i], a.i ); // verify
    }

}


void stringTest( void ) {
    
    // checks
    BOOST_CHECK( onlyDigits( "1234567890" ) );
    BOOST_CHECK( !onlyDigits( "1234d567890" ) );
    BOOST_CHECK( onlyAlpha( "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
    BOOST_CHECK( !onlyAlpha( "abcdefghijklmnopqrstuvwxyz7ABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
    BOOST_CHECK( onlyAlnum( "1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
    BOOST_CHECK( !onlyAlnum( "1234567890abcdefghijklmnopqrstuvw@xyzABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
    BOOST_CHECK( onlyHex( "1234567890abcdefABCDEF" ) );
    BOOST_CHECK( !onlyHex( "1234567890aQbcdefABCDEF" ) );
    BOOST_CHECK( isInteger( "1234567890" ) );
    BOOST_CHECK( isInteger( "0x123" ) );
    BOOST_CHECK( isInteger( "-19" ) );
    BOOST_CHECK( !isInteger( "-19.4" ) );
    BOOST_CHECK( !isInteger( "" ) );


    // align
    BOOST_CHECK_EQUAL( alignLeft( "string", 10, '-' ),   "string----" );
    BOOST_CHECK_EQUAL( alignRight( "string", 10, '-' ),  "----string" );
    BOOST_CHECK_EQUAL( alignCenter( "string", 10, '|' ), "||string||" );
    BOOST_CHECK_EQUAL( alignCenter( "string", 11, '|' ), "|||string||" ); // for odd padding, the extra char goes on the left side


    // hexString
    BOOST_CHECK_EQUAL( hexString( uint8_t( 0xc ) ), "0xc" );
    BOOST_CHECK_EQUAL( hexString( int16_t( 0x301 ) ), "0x301" );
    BOOST_CHECK_EQUAL( hexString( int32_t( 0xf5a3b1 ) ), "0xf5a3b1" );
    BOOST_CHECK_EQUAL( hexString( size_t( 0xff0af5a3b1 ) ), "0xff0af5a3b1" );
    BOOST_CHECK_EQUAL( hexString( uint8_t( 0xc ), false ), "c" );
    BOOST_CHECK_EQUAL( hexString( int16_t( 0x301 ), false ), "301" );
    BOOST_CHECK_EQUAL( hexString( int32_t( 0xf5a3b1 ), false ), "f5a3b1" );
    BOOST_CHECK_EQUAL( hexString( size_t( 0xff0af5a3b1 ), false ), "ff0af5a3b1" );


    // bitString
    BOOST_CHECK_EQUAL( bitString( uint8_t( 0xc ) ), "00001100" );
    BOOST_CHECK_EQUAL( bitString( int16_t( 0x301 ) ), "00000011 00000001" );
    uint32_t tmp = 0xf5a3b1;
    BOOST_CHECK_EQUAL( bitString( tmp ), "00000000 11110101 10100011 10110001" );
    BOOST_CHECK_EQUAL( bitString( ( uint8_t* )&tmp, 4 ), "00000000 11110101 10100011 10110001" );


}

namespace testsuite {

    namespace util {

        void utilTest( void ) {

            test_suite* ts = BOOST_TEST_SUITE( "UTIL" );

            ts->add( BOOST_TEST_CASE( &arrayTest ) );
            ts->add( BOOST_TEST_CASE( &bitTest ) );
            ts->add( BOOST_TEST_CASE( &boundValueTest ) );
            ts->add( BOOST_TEST_CASE( &dataTest ) );
            ts->add( BOOST_TEST_CASE( &endianTest ) );
            ts->add( BOOST_TEST_CASE( &stringTest ) );

            framework::master_test_suite().add( ts );

        }

    }

}

