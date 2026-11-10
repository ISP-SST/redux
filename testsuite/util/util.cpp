#include "catch.hpp"

#include "redux/util/bitoperations.hpp"
#include "redux/util/boundvalue.hpp"
#include "redux/util/point.hpp"
#include "redux/util/region.hpp"

#include <cstdlib>
#include <cstring>
#include <sstream>

using namespace redux::util;

using namespace std;

namespace testsuite {

    namespace util {

        namespace {
            uint64_t naiveBitCount( uint64_t v ) {
                uint64_t cnt( 0 );
                do {
                    if( v & 1 ) cnt++;
                }
                while( v >>= 1 );
                return cnt;
            }
        }

        TEST_CASE( "Bit manipulations", "[util]" ) {

            // deBruijn based log2
            for( uint8_t i( 1 ); i < 32; ++i ) {
                CHECK( (redux::util::log2( ( 1 << i ) + 1 )) == (i) );
            }

            // deBruijn LSB detector
            for( uint8_t i( 0 ); i < 32; ++i ) {
                CHECK( (findLSB( 1 << i )) == (i) );
            }

            // deBruijn LSB detector, 1-based index
            for( uint8_t i( 0 ); i < 32; ++i ) {
                CHECK( (findLSB1( 1 << i )) == (i + 1) );
            }

            // 64-bit deBruijn LSB detector, 1-based index
            uint64_t tmp64( 1 );
            for( uint8_t i( 0 ); i < 64; ++i ) {
                CHECK( (findLSB64( tmp64 << i )) == (i + 1) );
            }

            // next power of two
            for( uint8_t i( 0 ); i < 31; ++i ) {
                CHECK( (nextPowerOfTwo( ( 1 << i ) + 1 )) == (( 1 << ( i + 1 ) )) );
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
                CHECK( (a.s[0]) == (b.s[2]) ); // verify
                CHECK( (a.s[2]) == (b.s[0]) );
                CHECK( (a.s[1]) == (b.s[1]) ); // check that the unmasked bits are unchanged.
                CHECK( (a.s[3]) == (b.s[3]) );
            }

            // swap first 2 bytes of each integer in array with matching bytes in another array
            uint32_t aa[10], bb[10], cc[10];
            memset( bb, 0, 10 * sizeof( uint32_t ) );
            for( uint16_t i( 0 ); i < 10; i++ ) {
                aa[i] = cc[i] = rand();
            }
            swapBits( aa, bb, mask, 10 );
            for( uint16_t i( 0 ); i < 10; i++ ) {
                CHECK( (cc[i] & ( mask << 16 )) == (aa[i]) );
                CHECK( (cc[i]&mask) == (bb[i]) );
            }

            memcpy( aa, cc, 10 * sizeof( uint32_t ) );
            // flip the first 16 bits in aa[10]
            flipBits( aa, mask, 10 );
            for( uint16_t i( 0 ); i < 10; i++ ) {
                CHECK( (cc[i] ^ mask) == (aa[i]) );
            }

            for( int i( 0 ); i < 15; ++i ) {
                a.i[0] = rand();
                a.i[1] = rand();
                //CHECK( (countBits(a.b[0])) == (naiveBitCount(b.s[0])) );
                //CHECK( (countBits(a.s[0])) == (naiveBitCount(a.s[0])) );
                CHECK( (countBits( a.i[0] )) == (naiveBitCount( a.i[0] )) );
                CHECK( (countBits( a.l )) == (naiveBitCount( a.l )) );
            }

        }

        TEMPLATE_TEST_CASE( "BoundValue across types", "[util]",
                            int8_t, uint8_t, int16_t, uint16_t, int, uint64_t, float, double ) {

            using T = TestType;
            BoundValue<T> bv( T(3), T(2), T(7) );                 // default/truncated restriction

            CHECK( (bv) == (T(3)) );
            CHECK( (bv.min()) == (T(2)) );
            CHECK( (bv.max()) == (T(7)) );
            CHECK( (bv.span()) == (T(5)) );

            // assign
            bv = T(6); CHECK( (bv) == (T(6)) );
            bv = T(1); CHECK( (bv) == (bv.min()) );
            bv = T(9); CHECK( (bv) == (bv.max()) );

            // +=, -=
            bv -= T(4); CHECK( (bv) == (T(3)) );
            bv -= T(2); CHECK( (bv) == (bv.min()) );
            bv += T(5); CHECK( (bv) == (T(7)) );
            bv += T(7); CHECK( (bv) == (bv.max()) );

            // *=, /=
            bv = T(2);
            bv *= T(2); CHECK( (bv) == (T(4)) );
            bv /= T(2); CHECK( (bv) == (T(2)) );

            // ++, --
            bv = T(5);
            CHECK( (bv++) == (T(5)) );
            CHECK( (++bv) == (T(7)) );
            CHECK( (++bv) == (T(7)) );    // should not move outside boundary
            bv = T(4);
            CHECK( (bv--) == (T(4)) );
            CHECK( (--bv) == (T(2)) );
            CHECK( (--bv) == (T(2)) );    // should not move outside boundary

            // <<, >>
            stringstream ss;
            ss << bv;
            bv++;       // destroy value and then restore it from ss.
            ss >> bv;
            CHECK( (bv) == (T(2)) );

            bv.setRestrictType( detail::WRAP );         // periodic restriction
            bv = T(6); CHECK( (bv) == (T(6)) );
            bv = T(2); CHECK( (bv) == (T(2)) );
            bv = T(7); CHECK( (bv) == (T(2)) );
            bv = T(1); CHECK( (bv) == (T(6)) );
            bv = T(9); CHECK( (bv) == (T(4)) );

            bv.setRestrictType( detail::REFLECT );      // reflection at boundaries
            bv = T(6); CHECK( (bv) == (T(6)) );
            bv = T(2); CHECK( (bv) == (T(2)) );
            bv = T(7); CHECK( (bv) == (T(7)) );
            bv = T(1); CHECK( (bv) == (T(3)) );
            bv = T(9); CHECK( (bv) == (T(5)) );

        }

        namespace {
            void LimitBVTest( int val, int mn, int mx ) {
                BoundValue<int> bv( val, mn, mx );
                CHECK( (bv) == (val) );
                CHECK( (bv.max()) == (std::max<int>(mn,mx)) );
                CHECK( (bv.min()) == (std::min<int>(mn,mx)) );
                CHECK( (bv.span()) == (Abs(mx-mn)) );
            }
        }

        TEST_CASE( "BoundValue", "[util]" ) {

            typedef BoundValue<int, detail::UNDEFINEDTRIM> invalidBV;
            CHECK_THROWS_AS( invalidBV(), invalid_argument );

            LimitBVTest( 3, -9, 9 );
            LimitBVTest( 3, 9, -9 );

            BoundValue<int> intBV( 3, -9, 9 );
            CHECK_NOTHROW( (BoundValue<float>(intBV)) );
            BoundValue<float> floatBV( intBV );
            CHECK( (floatBV) == (3.0f) );
            floatBV = 2;
            CHECK( (floatBV) == (2.0f) );
            intBV = floatBV;
            CHECK( (intBV) == (2) );


            double delta( 0.05 );
            double epsilon = 1E-9;

            BoundValue<double, redux::util::detail::WRAP> periodicVal( 0, 0, 2 * M_PI );
            periodicVal = 3 * M_PI + delta;
            CHECK( (( double )periodicVal) == Approx(M_PI + delta).epsilon((epsilon)/100.0) );
            periodicVal = -3 * M_PI / 2 + delta;
            CHECK( (( double )periodicVal) == Approx(M_PI / 2 + delta).epsilon((epsilon)/100.0) );


            BoundValue<double, redux::util::detail::REFLECT> inclVal( 0, 0, M_PI ); // e.g. inclination = PI/2+\delta -> PI/2 - \delta
            inclVal = M_PI + delta;
            CHECK( (( double )inclVal) == Approx(M_PI - delta).epsilon((epsilon)/100.0) );
            inclVal = 2 * M_PI + delta;
            CHECK( (( double )inclVal) == Approx(delta).epsilon((epsilon)/100.0) );
            inclVal = 3 * M_PI + delta;
            CHECK( (( double )inclVal) == Approx(M_PI - delta).epsilon((epsilon)/100.0) );
            inclVal = -delta;
            CHECK( (( double )inclVal) == Approx(delta).epsilon((epsilon)/100.0) );
            inclVal = -M_PI - delta;
            CHECK( (( double )inclVal) == Approx(M_PI - delta).epsilon((epsilon)/100.0) );
            inclVal = -2 * M_PI - delta;
            CHECK( (( double )inclVal) == Approx(delta).epsilon((epsilon)/100.0) );


            BoundValue<int, redux::util::detail::REFLECT> reflectedVal( 0, 0, 4 );
            reflectedVal = 0;
            CHECK( (reflectedVal) == (0) );
            reflectedVal = 5;
            CHECK( (reflectedVal) == (3) );
            reflectedVal = 9;
            CHECK( (reflectedVal) == (1) );
            reflectedVal = 11;
            CHECK( (reflectedVal) == (3) );
            reflectedVal.setMax( 5 );
            reflectedVal = -3;
            CHECK( (reflectedVal) == (3) );
            reflectedVal = -5;
            CHECK( (reflectedVal) == (5) );
            reflectedVal = -9;
            CHECK( (reflectedVal) == (1) );
            reflectedVal = -11;
            CHECK( (reflectedVal) == (1) );

        }


        TEMPLATE_TEST_CASE( "Point across types", "[util]",
                            int8_t, uint8_t, int16_t, uint16_t, int, int64_t, float ) {

            using T = TestType;
            PointType<T> p;                     // empty constructor
            CHECK( (p.y) == (T(0)) );
            CHECK( (p.x) == (T(0)) );

            p = T(1);                           // assign value
            CHECK( (p.y) == (T(1)) );
            CHECK( (p.x) == (T(1)) );

            p = PointType<T>(2,2);              // assign Point
            CHECK( (p.y) == (T(2)) );
            CHECK( (p.x) == (T(2)) );

            p = PointType<uint8_t>(3,3);        // assign different Point
            CHECK( (p.y) == (T(3)) );
            CHECK( (p.x) == (T(3)) );

            T a(10), b(20);
            p = PointType<T>(a);        // single-value constructor
            CHECK( (p.y) == (a) );
            CHECK( (p.x) == (T(0)) );

            p = PointType<T>(a,b);      // two-value constructor
            CHECK( (p.y) == (a) );
            CHECK( (p.x) == (b) );

            CHECK( (p.max()) == (std::max<T>(a,b)) );
            CHECK( (p.min()) == (std::min<T>(a,b)) );

            PointI r = p.round();
            CHECK( (r.y) == (lround(p.y)) );
            CHECK( (r.x) == (lround(p.x)) );

            // min/max
            CHECK( (p.min()) == (std::min<T>(a,b)) );
            CHECK( (p.max()) == (std::max<T>(a,b)) );
            CHECK( (p.min(PointI(a,b))) == (PointType<T>(a,b)) );         // same
            CHECK( (p.min(PointI(a+1,b+1))) == (PointType<T>(a,b)) );
            CHECK( (p.min(PointI(a-1,b))) == (PointType<T>(a-1,b)) );
            CHECK( (p.min(PointI(a,b-1))) == (PointType<T>(a,b-1)) );
            CHECK( (p.min(PointI(a-1,b-1))) == (PointType<T>(a-1,b-1)) );
            CHECK( (p.max(PointI(a,b))) == (PointType<T>(a,b)) );         // same
            CHECK( (p.max(PointI(a-1,b-1))) == (PointType<T>(a,b)) );
            CHECK( (p.max(PointI(a+1,b))) == (PointType<T>(a+1,b)) );
            CHECK( (p.max(PointI(a,b+1))) == (PointType<T>(a,b+1)) );
            CHECK( (p.max(PointI(a+1,b+1))) == (PointType<T>(a+1,b+1)) );


            // ==, !=, <, >
            p = PointType<T>(3,4);          // first value (y or "row") takes precedence in comparisons.
            CHECK( (p) == (PointType<T>(3,4)) );
            CHECK( (p) != (PointType<T>(2,4)) );
            CHECK( (p) != (PointType<T>(3,5)) );
            CHECK( (p) < (PointType<T>(3,5)) );
            CHECK( (p) < (PointType<T>(4,4)) );
            CHECK( (p) > (PointType<T>(3,3)) );
            CHECK( (p) > (PointType<T>(2,4)) );
            p = 5;
            CHECK( (p) == (5) );
            CHECK( (p) != (4) );

            // +=, -=
            p = PointType<T>(a,b);
            p += 2;
            CHECK( (p) == (PointType<T>(a+2,b+2)) );
            p -= 2;
            CHECK( (p) == (PointType<T>(a,b)) );
            p += PointI(2,2);
            CHECK( (p) == (PointType<T>(a+2,b+2)) );
            p -= PointI(2,2);
            CHECK( (p) == (PointType<T>(a,b)) );

            // +, -
            CHECK( ((p+2)) == (PointType<T>(a+2,b+2)) );
            CHECK( ((p-2)) == (PointType<T>(a-2,b-2)) );
            CHECK( ((p+PointI(2,3))) == (PointType<T>(a+2,b+3)) );
            CHECK( ((p-PointI(3,2))) == (PointType<T>(a-3,b-2)) );

            // *=, /=
            p = PointType<T>(a,b);
            CHECK( ((p*=2)) == (PointType<T>(a*2,b*2)) );
            CHECK( ((p/=2)) == (PointType<T>(a,b)) );
            CHECK( ((p*=PointI(2,3))) == (PointType<T>(a*2,b*3)) );
            CHECK( ((p/=PointI(2,3))) == (PointType<T>(a,b)) );

            // *, /
            CHECK( ((p*2)) == (PointType<T>(a*2,b*2)) );
            CHECK( ((p/2)) == (PointType<T>(a/2,b/2)) );
            CHECK( ((p*PointI(2,3))) == (PointType<T>(a*2,b*3)) );
            CHECK( ((p/PointI(5,2))) == (PointType<T>(a/5,b/2)) );

            // pack/unpack
            uint64_t sz = p.size();
            char buf[32];      // large enough for all basic types
            uint64_t count = p.pack( buf );
            CHECK( (count) == (sz) );
            p = 0;
            CHECK( (p) == (PointType<T>(0,0)) );
            count = p.unpack( buf, false );
            CHECK( (count) == (sz) );
            CHECK( (p) == (PointType<T>(a,b)) );

        }

        TEST_CASE( "Point struct", "[util]" ) {

            // cast to (string)
            PointI ip(2,3);
            CHECK( ((string)ip) == (("("+to_string(ip.y)+","+to_string(ip.x)+")")) );

            // Abs, absMax
            CHECK( (ip.Abs()) == (ip) );
            CHECK( (PointI(-2,3).Abs()) == (ip) );
            CHECK( (PointI(2,-3).Abs()) == (ip) );
            CHECK( (PointI(-2,-3).Abs()) == (ip) );
            CHECK( (ip.maxAbs()) == (3) );
            PointD dp(2,3);
            CHECK( (dp.Abs()) == (dp) );
            CHECK( (PointD(-2,3).Abs()) == (dp) );
            CHECK( (PointD(2,-3).Abs()) == (dp) );
            CHECK( (PointD(-2,-3).Abs()) == (dp) );
            CHECK( (dp.maxAbs()) == (3) );

            // round
            CHECK( (PointD(1.5,3.49).round()) == (ip) );
            CHECK( (PointD(2.49,2.50).round()) == (ip) );

        }

        TEMPLATE_TEST_CASE( "Region across types", "[util]",
                            int8_t, uint8_t, int16_t, uint16_t, int, uint64_t, float, double ) {

            using T = TestType;
            RegionType<T> r;                        // empty constructor
            CHECK( (r.first) == (T(0)) );
            CHECK( (r.last) == (T(0)) );

            r = RegionType<T>(2,3);                 // construct from upper limit only
            CHECK( (r.first) == (T(0)) );
            CHECK( (r.last) == (PointType<T>(2,3)) );
            r = RegionType<T>(PointType<T>(3,2));
            CHECK( (r.first) == (T(0)) );
            CHECK( (r.last) == (PointType<T>(3,2)) );

            r = RegionType<T>(1,2,3,4);             // construct from upper & lower limits
            CHECK( (r.first) == (PointType<T>(1,2)) );
            CHECK( (r.last) == (PointType<T>(3,4)) );
            r = RegionType<T>( PointType<T>(2,3), PointType<T>(4,5) );
            CHECK( (r.first) == (PointType<T>(2,3)) );
            CHECK( (r.last) == (PointType<T>(4,5)) );

            // shift
            r.shift( 2, 3 );
            CHECK( (r.first) == (PointType<T>(4,6)) );
            CHECK( (r.last) == (PointType<T>(6,8)) );
            r.shift( PointType<T>(2,3) );
            CHECK( (r.first) == (PointType<T>(6,9)) );
            CHECK( (r.last) == (PointType<T>(8,11)) );
            r.shiftY( -4 );
            CHECK( (r.first) == (PointType<T>(2,9)) );
            CHECK( (r.last) == (PointType<T>(4,11)) );
            r.shiftX( -6 );
            CHECK( (r.first) == (PointType<T>(2,3)) );
            CHECK( (r.last) == (PointType<T>(4,5)) );

            // normalize
            r = RegionType<T>(10,3,3,10);
            CHECK( (r.first) == (PointType<T>(10,3)) );
            CHECK( (r.last) == (PointType<T>(3,10)) );
            r.normalize();
            CHECK( (r.first) == (PointType<T>(3,3)) );
            CHECK( (r.last) == (PointType<T>(10,10)) );

            // isInside
            r = RegionType<T>(10,20,30,40);
            CHECK(  r.isInside(PointType<T>(10,20))  );
            CHECK(  !r.isInside(PointType<T>(9,20))  );
            CHECK(  !r.isInside(PointType<T>(10,19))  );
            CHECK(  r.isInside(PointType<T>(30,40))  );
            CHECK(  !r.isInside(PointType<T>(31,40))  );
            CHECK(  !r.isInside(PointType<T>(30,41))  );
            CHECK(  r.isInside(PointType<T>(20,30))  );

            // isNormal
            CHECK(  RegionType<T>(10,20,30,40).isNormal()  );
            CHECK(  !RegionType<T>(30,20,10,40).isNormal()  );
            CHECK(  !RegionType<T>(10,40,30,20).isNormal()  );
            CHECK(  !RegionType<T>(30,40,10,20).isNormal()  );

            // diff
            CHECK( (r.diff(PointI(10,20))) == (PointI()) );
            CHECK( (r.diff(PointI(30,40))) == (PointI()) );
            CHECK( (r.diff(PointI(9,30))) == (PointI(-1,0)) );
            CHECK( (r.diff(PointI(20,19))) == (PointI(0,-1)) );
            CHECK( (r.diff(PointI(31,22))) == (PointI(1,0)) );
            CHECK( (r.diff(PointI(31,55))) == (PointI(1,15)) );

            // intersection
            RegionType<T> rr = r.intersection(RegionType<T>(10,20,30,40));  // the same as r
            CHECK( (rr) == (r) );
            rr = r.intersection(RegionType<T>(15,25,25,35));        // completely within r
            CHECK( (rr) == (RegionType<T>(15,25,25,35)) );
            rr = r.intersection(RegionType<T>(2,2,9,9));            // completely outside r
            CHECK( (rr) == (RegionType<T>()) );
            rr = r.intersection(RegionType<T>(15,25,85,95));        // first inside r, last outside
            CHECK( (rr) == (RegionType<T>(15,25,30,40)) );
            rr = r.intersection(RegionType<T>(2,3,25,35));          // last inside r, first outside
            CHECK( (rr) == (RegionType<T>(10,20,25,35)) );

            // enclosure
            rr = r.enclosure(RegionType<T>(10,20,30,40));           // the same as r
            CHECK( (rr) == (r) );
            rr = r.enclosure(RegionType<T>(15,25,25,35));           // completely within r
            CHECK( (rr) == (r) );
            rr = r.enclosure(RegionType<T>(2,2,9,9));               // completely outside r
            CHECK( (rr) == (RegionType<T>(2,2,30,40)) );
            rr = r.enclosure(RegionType<T>(15,25,85,95));           // first inside r, last outside
            CHECK( (rr) == (RegionType<T>(10,20,85,95)) );
            rr = r.enclosure(RegionType<T>(2,3,25,35));             // last inside r, first outside
            CHECK( (rr) == (RegionType<T>(2,3,30,40)) );

            // expand/shrink
            rr = r;     // = (10,20,30,40)
            rr.expand(PointType<T>(4,3));
            CHECK( (rr) == (RegionType<T>(6,17,34,43)) );
            rr.expand(2);
            CHECK( (rr) == (RegionType<T>(4,15,36,45)) );
            rr.shrink(PointType<T>(4,3));
            CHECK( (rr) == (RegionType<T>(8,18,32,42)) );
            rr.shrink(2);
            CHECK( (rr) == (RegionType<T>(10,20,30,40)) );
            rr.shrinkSigned( PointI(-1,0) );
            CHECK( (rr) == (RegionType<T>(11,20,30,40)) );
            rr.shrinkSigned( PointI(0,-1) );
            CHECK( (rr) == (RegionType<T>(11,21,30,40)) );
            rr.shrinkSigned( PointI(-1,-1) );
            CHECK( (rr) == (RegionType<T>(12,22,30,40)) );
            rr.shrinkSigned( PointI(1,0) );
            CHECK( (rr) == (RegionType<T>(12,22,29,40)) );
            rr.shrinkSigned( PointI(0,1) );
            CHECK( (rr) == (RegionType<T>(12,22,29,39)) );
            rr.shrinkSigned( PointI(1,1) );
            CHECK( (rr) == (RegionType<T>(12,22,28,38)) );

            // include
            rr = r;     // = (10,20,30,40)
            rr.include(PointType<T>(40,30));                        // y > last.y
            CHECK( (rr) == (RegionType<T>(10,20,40,40)) );
            rr.include(30,50);                                      // x > last.x
            CHECK( (rr) == (RegionType<T>(10,20,40,50)) );
            rr.include(100,110);                                    // both > last
            CHECK( (rr) == (RegionType<T>(10,20,100,110)) );
            rr.include(9,30);                                       // y < first.y
            CHECK( (rr) == (RegionType<T>(9,20,100,110)) );
            rr.include(30,7);                                       // x < first.x
            CHECK( (rr) == (RegionType<T>(9,7,100,110)) );
            rr.include(2,3);                                        // both < first
            CHECK( (rr) == (RegionType<T>(2,3,100,110)) );
            rr.include(20,30);                                      // both inside -> no change
            CHECK( (rr) == (RegionType<T>(2,3,100,110)) );

            // grow
            rr = r;     // = (10,20,30,40)
            rr.grow(RegionType<T>(10,20,30,40));                        // the same as r
            CHECK( (rr) == (r) );
            rr = r; rr.grow(RegionType<T>(15,25,25,35));                // completely within r
            CHECK( (rr) == (r) );
            rr = r; rr.grow(RegionType<T>(2,2,9,9));                    // completely outside r
            CHECK( (rr) == (RegionType<T>(2,2,30,40)) );
            rr = r; rr.grow(RegionType<T>(15,25,85,95));                // first inside r, last outside
            CHECK( (rr) == (RegionType<T>(10,20,85,95)) );
            rr = r; rr.grow(RegionType<T>(2,3,25,35));                  // last inside r, first outside
            CHECK( (rr) == (RegionType<T>(2,3,30,40)) );
            rr = r; rr.grow(PointType<T>(25,27));                       // grow by point inside
            CHECK( (rr) == (RegionType<T>(10,20,30,40)) );
            rr = r; rr.grow(PointType<T>(5,7));                         // grow by point outside
            CHECK( (rr) == (RegionType<T>(5,7,30,40)) );
            rr = r; rr.grow(PointType<T>(55,77));                       // grow by point outside
            CHECK( (rr) == (RegionType<T>(10,20,55,77)) );

            // getAsGlobal
            CHECK( (r.getAsGlobal(PointI(1,2))) == (PointI(11,22)) );
            CHECK( (r.getAsGlobal(PointI(-1,-2))) == (PointI(9,18)) );
            CHECK( (r.getAsGlobal(r.last-r.first)) == (r.last) );
            CHECK( (RegionType<T>(30,20,10,40).getAsGlobal(PointI(1,2))) == (PointI(29,22)) );
            CHECK( (RegionType<T>(10,40,30,20).getAsGlobal(PointI(1,2))) == (PointI(11,38)) );
            CHECK( (RegionType<T>(30,40,10,20).getAsGlobal(PointI(1,2))) == (PointI(29,38)) );
            CHECK( (RegionType<T>(30,20,10,40).getAsGlobal(PointI(-1,-2))) == (PointI(31,18)) );
            CHECK( (RegionType<T>(10,40,30,20).getAsGlobal(PointI(-1,-2))) == (PointI(9,42)) );
            CHECK( (RegionType<T>(30,40,10,20).getAsGlobal(PointI(-1,-2))) == (PointI(31,42)) );
            CHECK( (RegionType<T>(30,40,10,20).getAsGlobal(PointI(20,20))) == (PointI(10,20)) );

            // pack/unpack
            r = RegionType<T>(2,3,4,5);
            uint64_t sz = r.size();
            char buf[64];      // large enough for all basic types
            uint64_t count = r.pack( buf );
            CHECK( (count) == (sz) );
            r = 0;
            CHECK( (r) == (RegionType<T>(0,0)) );
            count = r.unpack( buf, false );
            CHECK( (count) == (sz) );
            CHECK( (r) == (RegionType<T>(2,3,4,5)) );

        }

        TEST_CASE( "Region struct", "[util]" ) {

            // cast to (string)
            RegionType<int> ir(2,3,4,5);
            CHECK( ((string)ir) == (("[2:4,3:5]")) );

        }

    }

}
