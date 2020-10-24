#include <boost/test/unit_test.hpp>

#include "redux/util/bitoperations.hpp"
#include "redux/util/boundvalue.hpp"
#include "redux/util/point.hpp"
#include "redux/util/region.hpp"

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
        
        template<typename T>
        void typedBVTest( void ) {
            
            BoundValue<T> bv( T(3), T(2), T(7) );                 // default/truncated restriction
            
            BOOST_CHECK_EQUAL( bv, T(3) );
            BOOST_CHECK_EQUAL( bv.min(),  T(2) );
            BOOST_CHECK_EQUAL( bv.max(),  T(7) );
            BOOST_CHECK_EQUAL( bv.span(), T(5) );
            
            // assign
            bv = T(6); BOOST_CHECK_EQUAL( bv, T(6) );
            bv = T(1); BOOST_CHECK_EQUAL( bv, bv.min() );
            bv = T(9); BOOST_CHECK_EQUAL( bv, bv.max() );
            
            // +=, -=
            bv -= T(4); BOOST_CHECK_EQUAL( bv, T(3) );
            bv -= T(2); BOOST_CHECK_EQUAL( bv, bv.min() );
            bv += T(5); BOOST_CHECK_EQUAL( bv, T(7) );
            bv += T(7); BOOST_CHECK_EQUAL( bv, bv.max() );
            
            // *=, /=
            bv = T(2);
            bv *= T(2); BOOST_CHECK_EQUAL( bv, T(4) );
            bv /= T(2); BOOST_CHECK_EQUAL( bv, T(2) );
            
            // ++, --
            bv = T(5);
            BOOST_CHECK_EQUAL( bv++, T(5) );
            BOOST_CHECK_EQUAL( ++bv, T(7) );
            BOOST_CHECK_EQUAL( ++bv, T(7) );    // should not move outside boundary
            bv = T(4);
            BOOST_CHECK_EQUAL( bv--, T(4) );
            BOOST_CHECK_EQUAL( --bv, T(2) );
            BOOST_CHECK_EQUAL( --bv, T(2) );    // should not move outside boundary
            
            // <<, >>
            stringstream ss;
            ss << bv;
            bv++;       // destroy value and then restore it from ss.
            ss >> bv;
            BOOST_CHECK_EQUAL( bv, T(2) );
            
            bv.setRestrictType( detail::WRAP );         // periodic restriction
            bv = T(6); BOOST_CHECK_EQUAL( bv, T(6) );
            bv = T(2); BOOST_CHECK_EQUAL( bv, T(2) );
            bv = T(7); BOOST_CHECK_EQUAL( bv, T(2) );
            bv = T(1); BOOST_CHECK_EQUAL( bv, T(6) );
            bv = T(9); BOOST_CHECK_EQUAL( bv, T(4) );
            
            bv.setRestrictType( detail::REFLECT );      // reflection at boundaries
            bv = T(6); BOOST_CHECK_EQUAL( bv, T(6) );
            bv = T(2); BOOST_CHECK_EQUAL( bv, T(2) );
            bv = T(7); BOOST_CHECK_EQUAL( bv, T(7) );
            bv = T(1); BOOST_CHECK_EQUAL( bv, T(3) );
            bv = T(9); BOOST_CHECK_EQUAL( bv, T(5) );
            
        }
        
        void LimitBVTest( int val, int mn, int mx ) {
            BoundValue<int> bv( val, mn, mx );
            BOOST_CHECK_EQUAL( bv, val );
            BOOST_CHECK_EQUAL( bv.max(), std::max<int>(mn,mx) );
            BOOST_CHECK_EQUAL( bv.min(), std::min<int>(mn,mx) );
            BOOST_CHECK_EQUAL( bv.span(), Abs(mx-mn) );
        }

        void boundValueTest( void ) {

            typedef BoundValue<int, detail::UNDEFINEDTRIM> invalidBV;
            BOOST_CHECK_THROW( invalidBV(), invalid_argument );
            
            typedBVTest<int8_t>();
            typedBVTest<uint8_t>();
            typedBVTest<int16_t>();
            typedBVTest<uint16_t>();
            typedBVTest<int>();
            typedBVTest<uint64_t>();
            typedBVTest<float>();
            typedBVTest<double>();
            
            LimitBVTest( 3, -9, 9 );
            LimitBVTest( 3, 9, -9 );

            BoundValue<int> intBV( 3, -9, 9 );
            BOOST_CHECK_NO_THROW( (BoundValue<float>(intBV)) );
            BoundValue<float> floatBV( intBV );
            BOOST_CHECK_EQUAL( floatBV, 3 );
            floatBV = 2;
            BOOST_CHECK_EQUAL( floatBV, 2 );
            intBV = floatBV;
            BOOST_CHECK_EQUAL( intBV, 2 );
            

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


        template <typename T>
        void typedPointTest( void ) {
            
            PointType<T> p;                     // empty constructor
            BOOST_CHECK_EQUAL( p.y, T(0) );
            BOOST_CHECK_EQUAL( p.x, T(0) );
            
            p = T(1);                           // assign value
            BOOST_CHECK_EQUAL( p.y, T(1) );
            BOOST_CHECK_EQUAL( p.x, T(1) );
            
            p = PointType<T>(2,2);              // assign Point
            BOOST_CHECK_EQUAL( p.y, T(2) );
            BOOST_CHECK_EQUAL( p.x, T(2) );
            
            p = PointType<uint8_t>(3,3);        // assign different Point
            BOOST_CHECK_EQUAL( p.y, T(3) );
            BOOST_CHECK_EQUAL( p.x, T(3) );
            
            T a(10), b(20);
            p = PointType<T>(a);        // single-value constructor
            BOOST_CHECK_EQUAL( p.y, a );
            BOOST_CHECK_EQUAL( p.x, T(0) );
            
            p = PointType<T>(a,b);      // two-value constructor
            BOOST_CHECK_EQUAL( p.y, a );
            BOOST_CHECK_EQUAL( p.x, b );
            
            BOOST_CHECK_EQUAL( p.max(), std::max<T>(a,b) );
            BOOST_CHECK_EQUAL( p.min(), std::min<T>(a,b) );
            
            PointI r = p.round();
            BOOST_CHECK_EQUAL( r.y, lround(p.y) );
            BOOST_CHECK_EQUAL( r.x, lround(p.x) );
            
            // min/max
            BOOST_CHECK_EQUAL( p.min(), std::min<T>(a,b) );
            BOOST_CHECK_EQUAL( p.max(), std::max<T>(a,b) );
            BOOST_CHECK_EQUAL( p.min(PointI(a,b)), PointType<T>(a,b) );         // same
            BOOST_CHECK_EQUAL( p.min(PointI(a+1,b+1)), PointType<T>(a,b) );
            BOOST_CHECK_EQUAL( p.min(PointI(a-1,b)), PointType<T>(a-1,b) );
            BOOST_CHECK_EQUAL( p.min(PointI(a,b-1)), PointType<T>(a,b-1) );
            BOOST_CHECK_EQUAL( p.min(PointI(a-1,b-1)), PointType<T>(a-1,b-1) );
            BOOST_CHECK_EQUAL( p.max(PointI(a,b)), PointType<T>(a,b) );         // same
            BOOST_CHECK_EQUAL( p.max(PointI(a-1,b-1)), PointType<T>(a,b) );
            BOOST_CHECK_EQUAL( p.max(PointI(a+1,b)), PointType<T>(a+1,b) );
            BOOST_CHECK_EQUAL( p.max(PointI(a,b+1)), PointType<T>(a,b+1) );
            BOOST_CHECK_EQUAL( p.max(PointI(a+1,b+1)), PointType<T>(a+1,b+1) );

            
            // ==, !=, <, >
            p = PointType<T>(3,4);          // first value (y or "row") takes precedence in comparisons.
            BOOST_CHECK_EQUAL( p, PointType<T>(3,4) );
            BOOST_CHECK_NE( p, PointType<T>(2,4) );
            BOOST_CHECK_NE( p, PointType<T>(3,5) );
            BOOST_CHECK_LT( p, PointType<T>(3,5) );
            BOOST_CHECK_LT( p, PointType<T>(4,4) );
            BOOST_CHECK_GT( p, PointType<T>(3,3) );
            BOOST_CHECK_GT( p, PointType<T>(2,4) );
            p = 5;
            BOOST_CHECK_EQUAL( p, 5 );
            BOOST_CHECK_NE( p, 4 );
            
            // +=, -=
            p = PointType<T>(a,b);
            p += 2;
            BOOST_CHECK_EQUAL( p, PointType<T>(a+2,b+2) );
            p -= 2;
            BOOST_CHECK_EQUAL( p, PointType<T>(a,b) );
            p += PointI(2,2);
            BOOST_CHECK_EQUAL( p, PointType<T>(a+2,b+2) );
            p -= PointI(2,2);
            BOOST_CHECK_EQUAL( p, PointType<T>(a,b) );
            
            // +, -
            BOOST_CHECK_EQUAL( (p+2), PointType<T>(a+2,b+2) );
            BOOST_CHECK_EQUAL( (p-2), PointType<T>(a-2,b-2) );
            BOOST_CHECK_EQUAL( (p+PointI(2,3)), PointType<T>(a+2,b+3) );
            BOOST_CHECK_EQUAL( (p-PointI(3,2)), PointType<T>(a-3,b-2) );

            // *=, /=
            p = PointType<T>(a,b);
            BOOST_CHECK_EQUAL( (p*=2), PointType<T>(a*2,b*2) );
            BOOST_CHECK_EQUAL( (p/=2), PointType<T>(a,b) );
            BOOST_CHECK_EQUAL( (p*=PointI(2,3)), PointType<T>(a*2,b*3) );
            BOOST_CHECK_EQUAL( (p/=PointI(2,3)), PointType<T>(a,b) );
            
            // *, /
            BOOST_CHECK_EQUAL( (p*2), PointType<T>(a*2,b*2) );
            BOOST_CHECK_EQUAL( (p/2), PointType<T>(a/2,b/2) );
            BOOST_CHECK_EQUAL( (p*PointI(2,3)), PointType<T>(a*2,b*3) );
            BOOST_CHECK_EQUAL( (p/PointI(5,2)), PointType<T>(a/5,b/2) );
            
            // pack/unpack
            uint64_t sz = p.size();
            char buf[32];      // large enough for all basic types
            uint64_t count = p.pack( buf );
            BOOST_CHECK_EQUAL( count, sz );
            p = 0;
            BOOST_CHECK_EQUAL( p, PointType<T>(0,0) );
            count = p.unpack( buf, false );
            BOOST_CHECK_EQUAL( count, sz );
            BOOST_CHECK_EQUAL( p, PointType<T>(a,b) );
            
        }
        
        void pointTest( void ) {
            
            // test for various types
            typedPointTest<int8_t>();
            typedPointTest<uint8_t>();
            typedPointTest<int16_t>();
            typedPointTest<uint16_t>();
            typedPointTest<int>();
            typedPointTest<int64_t>();
            typedPointTest<float>();
            
            // cast to (string)
            PointI ip(2,3);
            BOOST_CHECK_EQUAL( (string)ip, ("("+to_string(ip.y)+","+to_string(ip.x)+")") );
            
            // Abs, absMax
            BOOST_CHECK_EQUAL( ip.Abs(), ip );
            BOOST_CHECK_EQUAL( PointI(-2,3).Abs(), ip );
            BOOST_CHECK_EQUAL( PointI(2,-3).Abs(), ip );
            BOOST_CHECK_EQUAL( PointI(-2,-3).Abs(), ip );
            BOOST_CHECK_EQUAL( ip.maxAbs(), 3 );
            PointD dp(2,3);
            BOOST_CHECK_EQUAL( dp.Abs(), dp );
            BOOST_CHECK_EQUAL( PointD(-2,3).Abs(), dp );
            BOOST_CHECK_EQUAL( PointD(2,-3).Abs(), dp );
            BOOST_CHECK_EQUAL( PointD(-2,-3).Abs(), dp );
            BOOST_CHECK_EQUAL( dp.maxAbs(), 3 );
           
            // round
            BOOST_CHECK_EQUAL( PointD(1.5,3.49).round(), ip );
            BOOST_CHECK_EQUAL( PointD(2.49,2.50).round(), ip );
          
        }
        
        template <typename T>
        void typedRegionTest( void ) {
            
            RegionType<T> r;                        // empty constructor
            BOOST_CHECK_EQUAL( r.first, T(0) );
            BOOST_CHECK_EQUAL( r.last, T(0) );
            
            r = RegionType<T>(2,3);                 // construct from upper limit only
            BOOST_CHECK_EQUAL( r.first, T(0) );
            BOOST_CHECK_EQUAL( r.last, PointType<T>(2,3) );
            r = RegionType<T>(PointType<T>(3,2));
            BOOST_CHECK_EQUAL( r.first, T(0) );
            BOOST_CHECK_EQUAL( r.last, PointType<T>(3,2) );
            
            r = RegionType<T>(1,2,3,4);             // construct from upper & lower limits
            BOOST_CHECK_EQUAL( r.first, PointType<T>(1,2) );
            BOOST_CHECK_EQUAL( r.last, PointType<T>(3,4) );
            r = RegionType<T>( PointType<T>(2,3), PointType<T>(4,5) );
            BOOST_CHECK_EQUAL( r.first, PointType<T>(2,3) );
            BOOST_CHECK_EQUAL( r.last, PointType<T>(4,5) );
            
            // shift
            r.shift( 2, 3 );
            BOOST_CHECK_EQUAL( r.first, PointType<T>(4,6) );
            BOOST_CHECK_EQUAL( r.last, PointType<T>(6,8) );
            r.shift( PointType<T>(2,3) );
            BOOST_CHECK_EQUAL( r.first, PointType<T>(6,9) );
            BOOST_CHECK_EQUAL( r.last, PointType<T>(8,11) );
            r.shiftY( -4 );
            BOOST_CHECK_EQUAL( r.first, PointType<T>(2,9) );
            BOOST_CHECK_EQUAL( r.last, PointType<T>(4,11) );
            r.shiftX( -6 );
            BOOST_CHECK_EQUAL( r.first, PointType<T>(2,3) );
            BOOST_CHECK_EQUAL( r.last, PointType<T>(4,5) );
            
            // normalize
            r = RegionType<T>(10,3,3,10);
            BOOST_CHECK_EQUAL( r.first, PointType<T>(10,3) );
            BOOST_CHECK_EQUAL( r.last, PointType<T>(3,10) );
            r.normalize();
            BOOST_CHECK_EQUAL( r.first, PointType<T>(3,3) );
            BOOST_CHECK_EQUAL( r.last, PointType<T>(10,10) );
            
            // isInside
            r = RegionType<T>(10,20,30,40);
            BOOST_CHECK( r.isInside(PointType<T>(10,20)) );
            BOOST_CHECK( !r.isInside(PointType<T>(9,20)) );
            BOOST_CHECK( !r.isInside(PointType<T>(10,19)) );
            BOOST_CHECK( r.isInside(PointType<T>(30,40)) );
            BOOST_CHECK( !r.isInside(PointType<T>(31,40)) );
            BOOST_CHECK( !r.isInside(PointType<T>(30,41)) );
            BOOST_CHECK( r.isInside(PointType<T>(20,30)) );
            
            // isNormal
            BOOST_CHECK( RegionType<T>(10,20,30,40).isNormal() );
            BOOST_CHECK( !RegionType<T>(30,20,10,40).isNormal() );
            BOOST_CHECK( !RegionType<T>(10,40,30,20).isNormal() );
            BOOST_CHECK( !RegionType<T>(30,40,10,20).isNormal() );
            
            // diff
            BOOST_CHECK_EQUAL( r.diff(PointI(10,20)), PointI() );
            BOOST_CHECK_EQUAL( r.diff(PointI(30,40)), PointI() );
            BOOST_CHECK_EQUAL( r.diff(PointI(9,30)), PointI(-1,0) );
            BOOST_CHECK_EQUAL( r.diff(PointI(20,19)), PointI(0,-1) );
            BOOST_CHECK_EQUAL( r.diff(PointI(31,22)), PointI(1,0) );
            BOOST_CHECK_EQUAL( r.diff(PointI(31,55)), PointI(1,15) );
            
            // intersection
            RegionType<T> rr = r.intersection(RegionType<T>(10,20,30,40));  // the same as r
            BOOST_CHECK_EQUAL( rr, r );
            rr = r.intersection(RegionType<T>(15,25,25,35));        // completely within r
            BOOST_CHECK_EQUAL( rr, RegionType<T>(15,25,25,35) );
            rr = r.intersection(RegionType<T>(2,2,9,9));            // completely outside r
            BOOST_CHECK_EQUAL( rr, RegionType<T>() );
            rr = r.intersection(RegionType<T>(15,25,85,95));        // first inside r, last outside
            BOOST_CHECK_EQUAL( rr, RegionType<T>(15,25,30,40) );
            rr = r.intersection(RegionType<T>(2,3,25,35));          // last inside r, first outside
            BOOST_CHECK_EQUAL( rr, RegionType<T>(10,20,25,35) );
            
            // enclosure
            rr = r.enclosure(RegionType<T>(10,20,30,40));           // the same as r
            BOOST_CHECK_EQUAL( rr, r );
            rr = r.enclosure(RegionType<T>(15,25,25,35));           // completely within r
            BOOST_CHECK_EQUAL( rr, r );
            rr = r.enclosure(RegionType<T>(2,2,9,9));               // completely outside r
            BOOST_CHECK_EQUAL( rr, RegionType<T>(2,2,30,40) );
            rr = r.enclosure(RegionType<T>(15,25,85,95));           // first inside r, last outside
            BOOST_CHECK_EQUAL( rr, RegionType<T>(10,20,85,95) );
            rr = r.enclosure(RegionType<T>(2,3,25,35));             // last inside r, first outside
            BOOST_CHECK_EQUAL( rr, RegionType<T>(2,3,30,40) );
            
            // expand/shrink
            rr = r;     // = (10,20,30,40)
            rr.expand(PointType<T>(4,3));
            BOOST_CHECK_EQUAL( rr, RegionType<T>(6,17,34,43) );
            rr.expand(2);
            BOOST_CHECK_EQUAL( rr, RegionType<T>(4,15,36,45) );
            rr.shrink(PointType<T>(4,3));
            BOOST_CHECK_EQUAL( rr, RegionType<T>(8,18,32,42) );
            rr.shrink(2);
            BOOST_CHECK_EQUAL( rr, RegionType<T>(10,20,30,40) );
            rr.shrinkSigned( PointI(-1,0) );
            BOOST_CHECK_EQUAL( rr, RegionType<T>(11,20,30,40) );
            rr.shrinkSigned( PointI(0,-1) );
            BOOST_CHECK_EQUAL( rr, RegionType<T>(11,21,30,40) );
            rr.shrinkSigned( PointI(-1,-1) );
            BOOST_CHECK_EQUAL( rr, RegionType<T>(12,22,30,40) );
            rr.shrinkSigned( PointI(1,0) );
            BOOST_CHECK_EQUAL( rr, RegionType<T>(12,22,29,40) );
            rr.shrinkSigned( PointI(0,1) );
            BOOST_CHECK_EQUAL( rr, RegionType<T>(12,22,29,39) );
            rr.shrinkSigned( PointI(1,1) );
            BOOST_CHECK_EQUAL( rr, RegionType<T>(12,22,28,38) );
            
            // include
            rr = r;     // = (10,20,30,40)
            rr.include(PointType<T>(40,30));                        // y > last.y
            BOOST_CHECK_EQUAL( rr, RegionType<T>(10,20,40,40) );
            rr.include(30,50);                                      // x > last.x
            BOOST_CHECK_EQUAL( rr, RegionType<T>(10,20,40,50) );
            rr.include(100,110);                                    // both > last
            BOOST_CHECK_EQUAL( rr, RegionType<T>(10,20,100,110) );
            rr.include(9,30);                                       // y < first.y
            BOOST_CHECK_EQUAL( rr, RegionType<T>(9,20,100,110) );
            rr.include(30,7);                                       // x < first.x
            BOOST_CHECK_EQUAL( rr, RegionType<T>(9,7,100,110) );
            rr.include(2,3);                                        // both < first
            BOOST_CHECK_EQUAL( rr, RegionType<T>(2,3,100,110) );
            rr.include(20,30);                                      // both inside -> no change
            BOOST_CHECK_EQUAL( rr, RegionType<T>(2,3,100,110) );
           
            // grow
            rr = r;     // = (10,20,30,40)
            rr.grow(RegionType<T>(10,20,30,40));                        // the same as r
            BOOST_CHECK_EQUAL( rr, r );
            rr = r; rr.grow(RegionType<T>(15,25,25,35));                // completely within r
            BOOST_CHECK_EQUAL( rr, r );
            rr = r; rr.grow(RegionType<T>(2,2,9,9));                    // completely outside r
            BOOST_CHECK_EQUAL( rr, RegionType<T>(2,2,30,40) );
            rr = r; rr.grow(RegionType<T>(15,25,85,95));                // first inside r, last outside
            BOOST_CHECK_EQUAL( rr, RegionType<T>(10,20,85,95) );
            rr = r; rr.grow(RegionType<T>(2,3,25,35));                  // last inside r, first outside
            BOOST_CHECK_EQUAL( rr, RegionType<T>(2,3,30,40) );
            rr = r; rr.grow(PointType<T>(25,27));                       // grow by point inside
            BOOST_CHECK_EQUAL( rr, RegionType<T>(10,20,30,40) );
            rr = r; rr.grow(PointType<T>(5,7));                         // grow by point outside
            BOOST_CHECK_EQUAL( rr, RegionType<T>(5,7,30,40) );
            rr = r; rr.grow(PointType<T>(55,77));                       // grow by point outside
            BOOST_CHECK_EQUAL( rr, RegionType<T>(10,20,55,77) );
            
            // getAsGlobal
            BOOST_CHECK_EQUAL( r.getAsGlobal(PointI(1,2)), PointI(11,22) );
            BOOST_CHECK_EQUAL( r.getAsGlobal(PointI(-1,-2)), PointI(9,18) );
            BOOST_CHECK_EQUAL( r.getAsGlobal(r.last-r.first), r.last );
            BOOST_CHECK_EQUAL( RegionType<T>(30,20,10,40).getAsGlobal(PointI(1,2)), PointI(29,22) );
            BOOST_CHECK_EQUAL( RegionType<T>(10,40,30,20).getAsGlobal(PointI(1,2)), PointI(11,38) );
            BOOST_CHECK_EQUAL( RegionType<T>(30,40,10,20).getAsGlobal(PointI(1,2)), PointI(29,38) );
            BOOST_CHECK_EQUAL( RegionType<T>(30,20,10,40).getAsGlobal(PointI(-1,-2)), PointI(31,18) );
            BOOST_CHECK_EQUAL( RegionType<T>(10,40,30,20).getAsGlobal(PointI(-1,-2)), PointI(9,42) );
            BOOST_CHECK_EQUAL( RegionType<T>(30,40,10,20).getAsGlobal(PointI(-1,-2)), PointI(31,42) );
            BOOST_CHECK_EQUAL( RegionType<T>(30,40,10,20).getAsGlobal(PointI(20,20)), PointI(10,20) );
            
            // pack/unpack
            r = RegionType<T>(2,3,4,5); 
            uint64_t sz = r.size();
            char buf[64];      // large enough for all basic types
            uint64_t count = r.pack( buf );
            BOOST_CHECK_EQUAL( count, sz );
            r = 0;
            BOOST_CHECK_EQUAL( r, RegionType<T>(0,0) );
            count = r.unpack( buf, false );
            BOOST_CHECK_EQUAL( count, sz );
            BOOST_CHECK_EQUAL( r, RegionType<T>(2,3,4,5) );
            
    
        }
        
        void regionTest( void ) {
            
            // test for various types
            typedRegionTest<int8_t>();
            typedRegionTest<uint8_t>();
            typedRegionTest<int16_t>();
            typedRegionTest<uint16_t>();
            typedRegionTest<int>();
            typedRegionTest<uint64_t>();
            typedRegionTest<float>();
            typedRegionTest<double>();
            
            // cast to (string)
            RegionType<int> ir(2,3,4,5);
            BOOST_CHECK_EQUAL( (string)ir, ("[2:4,3:5]") );
          
        }
        
        void add_array_tests( test_suite* ts );     // defined in array.cpp
        void add_data_tests( test_suite* ts );      // defined in data.cpp
        void add_string_tests( test_suite* ts );    // defined in string.cpp

        void add_tests( test_suite* ts ) {

            ts->add( BOOST_TEST_CASE_NAME( &bitTest, "Bit manipulations"  ) );
            ts->add( BOOST_TEST_CASE_NAME( &boundValueTest, "BoundValue"  ) );
            ts->add( BOOST_TEST_CASE_NAME( &pointTest, "Point struct" ) );
            ts->add( BOOST_TEST_CASE_NAME( &regionTest, "Region struct" ) );

            add_array_tests( ts );
            add_data_tests( ts );
            add_string_tests( ts );

        }

    }

}

