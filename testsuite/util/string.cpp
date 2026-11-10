#include "redux/util/stringutil.hpp"

#include "catch.hpp"

#include <boost/numeric/conversion/cast.hpp>

#include <cstdlib>
#include <limits>
#include <string>
#include <type_traits>

using namespace redux::util;
using namespace std;

namespace testsuite {

    namespace util {

        TEST_CASE( "character-class checks", "[string]" ) {
            CHECK( onlyDigits( "1234567890" ) );
            CHECK_FALSE( onlyDigits( "1234d567890" ) );
            CHECK( onlyAlpha( "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
            CHECK_FALSE( onlyAlpha( "abcdefghijklmnopqrstuvwxyz7ABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
            CHECK( onlyAlnum( "1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
            CHECK_FALSE( onlyAlnum( "1234567890abcdefghijklmnopqrstuvw@xyzABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
            CHECK( onlyHex( "1234567890abcdefABCDEF" ) );
            CHECK_FALSE( onlyHex( "1234567890aQbcdefABCDEF" ) );
            CHECK( isInteger( "1234567890" ) );
            CHECK( isInteger( "0x123" ) );
            CHECK( isInteger( "-19" ) );
            CHECK_FALSE( isInteger( "-19.4" ) );
            CHECK_FALSE( isInteger( "" ) );
        }

        TEST_CASE( "string alignment/padding", "[string]" ) {
            SECTION( "left" )   { CHECK( alignLeft( "string", 10, '-' ) == "string----" ); }
            SECTION( "right" )  { CHECK( alignRight( "string", 10, '-' ) == "----string" ); }
            SECTION( "center, even padding" ) {
                CHECK( alignCenter( "string", 10, '|' ) == "||string||" );
            }
            SECTION( "center, odd padding goes left" ) {
                CHECK( alignCenter( "string", 11, '|' ) == "|||string||" );
            }
        }

        TEST_CASE( "hexString formatting", "[string]" ) {
            CHECK( hexString( uint8_t( 0xc ) ) == "0xc" );
            CHECK( hexString( int16_t( 0x301 ) ) == "0x301" );
            CHECK( hexString( int32_t( 0xf5a3b1 ) ) == "0xf5a3b1" );
            CHECK( hexString( size_t( 0xff0af5a3b1 ) ) == "0xff0af5a3b1" );
            CHECK( hexString( uint8_t( 0xc ), false ) == "c" );
            CHECK( hexString( int16_t( 0x301 ), false ) == "301" );
            CHECK( hexString( int32_t( 0xf5a3b1 ), false ) == "f5a3b1" );
            CHECK( hexString( size_t( 0xff0af5a3b1 ), false ) == "ff0af5a3b1" );
        }

        TEST_CASE( "bitString formatting", "[string]" ) {
            CHECK( bitString( uint8_t( 0xc ) ) == "00001100" );
            CHECK( bitString( int16_t( 0x301 ) ) == "00000011 00000001" );
            uint32_t tmp = 0xf5a3b1;
            CHECK( bitString( tmp ) == "00000000 11110101 10100011 10110001" );
            CHECK( bitString( (uint8_t*)&tmp, 4 ) == "00000000 11110101 10100011 10110001" );
        }

        TEST_CASE( "stringTo<bool> parsing", "[string]" ) {
            CHECK( stringTo<bool>( " T " ) );
            CHECK( stringTo<bool>( "true" ) );
            CHECK( stringTo<bool>( "yes" ) );
            CHECK_FALSE( stringTo<bool>( " f " ) );
            CHECK_FALSE( stringTo<bool>( "fAlse" ) );
            CHECK_FALSE( stringTo<bool>( "no" ) );
            CHECK( stringTo<bool>( " 153" ) );
            CHECK_FALSE( stringTo<bool>( " 0 " ) );
        }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"

        namespace {

            string decimalString( __int128 value ) {
                bool neg = value < 0;
                unsigned __int128 uval = neg ? -(unsigned __int128)value : (unsigned __int128)value;
                string s;
                do {
                    s.insert( s.begin(), '0' + int(uval % 10) );
                    uval /= 10;
                } while( uval > 0 );
                if( neg ) s.insert( s.begin(), '-' );
                return s;
            }

        }

#pragma GCC diagnostic pop

        TEMPLATE_TEST_CASE( "stringTo<T> integer cast + overflow", "[string]",
                            int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t ) {
            using T = TestType;
            T tMax = numeric_limits<T>::max();
            T tMin = numeric_limits<T>::min();
            T tMid = rand() % tMax;
            CHECK( stringTo<T>( " " + to_string(tMid) + " " ) == tMid );
            CHECK( stringTo<T>( " " + to_string(tMax) + " " ) == tMax );
            CHECK( stringTo<T>( " " + to_string(tMin) + " " ) == tMin );
            __int128 tmp = tMax;
            CHECK_THROWS_AS( stringTo<T>( " " + decimalString(tmp+1) + " " ), boost::numeric::positive_overflow );
            tmp = tMin;
            CHECK_THROWS_AS( stringTo<T>( " " + decimalString(tmp-1) + " " ), boost::numeric::negative_overflow );
            if( !is_signed<T>::value ) {
                CHECK_THROWS_AS( stringTo<T>( " -1 " ), boost::numeric::negative_overflow );
            }
        }

        TEST_CASE( "stringTo<float/double> parsing", "[string]" ) {
            CHECK( stringTo<float>( " 0 " ) == 0.0f );
            CHECK( stringTo<float>( " -1234.5 " ) == -1234.5f );
            CHECK( stringTo<double>( " 1.5E4 " ) == 1.5E4 );
            CHECK( stringTo<double>( " -1.5E-4 " ) == -1.5E-4 );

            float tmpF = stringTo<float>( " nan " );
            CHECK( std::isnan( tmpF ) );
            CHECK_FALSE( std::signbit( tmpF ) );
            tmpF = stringTo<float>( " -NaN " );
            CHECK( std::isnan( tmpF ) );
            CHECK( std::signbit( tmpF ) );
            tmpF = stringTo<float>( " inf " );
            CHECK( std::isinf( tmpF ) );
            CHECK_FALSE( std::signbit( tmpF ) );
            tmpF = stringTo<float>( " -Inf " );
            CHECK( std::isinf( tmpF ) );
            CHECK( std::signbit( tmpF ) );

            double tmpD = stringTo<double>( " nan " );
            CHECK( std::isnan( tmpD ) );
            CHECK_FALSE( std::signbit( tmpD ) );
            tmpD = stringTo<double>( " -NaN " );
            CHECK( std::isnan( tmpD ) );
            CHECK( std::signbit( tmpD ) );
            tmpD = stringTo<double>( " inf " );
            CHECK( std::isinf( tmpD ) );
            CHECK_FALSE( std::signbit( tmpD ) );
            tmpD = stringTo<double>( " -Inf " );
            CHECK( std::isinf( tmpD ) );
            CHECK( std::signbit( tmpD ) );
            CHECK( std::isnan( stringTo<double>( " NaN " ) ) );
        }

        TEST_CASE( "stringTo<T> throws on invalid input", "[string]" ) {
            CHECK_THROWS_AS( stringTo<int>( " 1234.5 " ), boost::bad_lexical_cast );
            CHECK_THROWS_AS( stringTo<bool>( " 1234.5 " ), boost::bad_lexical_cast );
        }

    }

}
