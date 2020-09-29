#include "redux/util/stringutil.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/multiprecision/cpp_int.hpp>

using namespace redux::util;
using namespace std;

using namespace boost::unit_test_framework;

namespace testsuite {

    namespace util {

    
        template <typename T>
        void testIntCast( void ) {
            T tMax = std::numeric_limits<T>::max();
            T tMin = std::numeric_limits<T>::min();
            T tMid = rand()%tMax;       // some value in the allowed range
            BOOST_CHECK( stringTo<T>(" "+to_string(tMid)+" ") == tMid );
            BOOST_CHECK( stringTo<T>(" "+to_string(tMax)+" ") == tMax );
            BOOST_CHECK( stringTo<T>(" "+to_string(tMin)+" ") == tMin );
            boost::multiprecision::int128_t tmp = tMax;
            BOOST_CHECK_THROW( stringTo<T>(" " + boost::lexical_cast<string>( tmp+1 ) + " "), boost::numeric::positive_overflow );
            tmp = tMin;
            BOOST_CHECK_THROW( stringTo<T>(" " + boost::lexical_cast<string>( tmp-1 ) + " "), boost::numeric::negative_overflow );
            if( !is_signed<T>::value ) {
                BOOST_CHECK_THROW( stringTo<T>(" -1 "), boost::numeric::negative_overflow );
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

            // stringTo (lexical casting)
            // some variants of bool cast
            BOOST_CHECK( stringTo<bool>(" T ") );
            BOOST_CHECK( stringTo<bool>("true") );
            BOOST_CHECK( stringTo<bool>("yes") );
            BOOST_CHECK( !stringTo<bool>(" f ") );
            BOOST_CHECK( !stringTo<bool>("fAlse") );
            BOOST_CHECK( !stringTo<bool>("no") );
            // integers should be cast as 0 = false, rest = true
            BOOST_CHECK( stringTo<bool>(" 153") );
            BOOST_CHECK( !stringTo<bool>(" 0 ") );

            
            // test cast & overflow for integer types
            testIntCast<int8_t>();
            testIntCast<uint8_t>();
            testIntCast<int16_t>();
            testIntCast<uint16_t>();
            testIntCast<int32_t>();
            testIntCast<uint32_t>();
            testIntCast<int64_t>();
            testIntCast<uint64_t>();

            
            // floating point types
            BOOST_CHECK( stringTo<float>(" 0 ") == 0.0f );
            BOOST_CHECK( stringTo<float>(" -1234.5 ") == -1234.5 );
            BOOST_CHECK( stringTo<double>(" 1.5E4 ") == 1.5E4 );
            BOOST_CHECK( stringTo<double>(" -1.5E-4 ") == -1.5E-4 );
            
            // special values    
            float tmpF = stringTo<float>(" nan ");
            BOOST_CHECK( std::isnan( tmpF ) );
            BOOST_CHECK( !std::signbit( tmpF ) );
            tmpF = stringTo<float>(" -NaN ");
            BOOST_CHECK( std::isnan( tmpF ) );
            BOOST_CHECK( std::signbit( tmpF ) );
            tmpF = stringTo<float>(" inf ");
            BOOST_CHECK( std::isinf( tmpF ) );
            BOOST_CHECK( !std::signbit( tmpF ) );
            tmpF = stringTo<float>(" -Inf ");
            BOOST_CHECK( std::isinf( tmpF ) );
            BOOST_CHECK( std::signbit( tmpF ) );
            
            float tmpD = stringTo<double>(" nan ");
            BOOST_CHECK( std::isnan( tmpD ) );
            BOOST_CHECK( !std::signbit( tmpD ) );
            tmpD = stringTo<double>(" -NaN ");
            BOOST_CHECK( std::isnan( tmpD ) );
            BOOST_CHECK( std::signbit( tmpD ) );
            tmpD = stringTo<double>(" inf ");
            BOOST_CHECK( std::isinf( tmpD ) );
            BOOST_CHECK( !std::signbit( tmpD ) );
            tmpD = stringTo<double>(" -Inf ");
            BOOST_CHECK( std::isinf( tmpD ) );
            BOOST_CHECK( std::signbit( tmpD ) );
            BOOST_CHECK( std::isnan( stringTo<double>(" NaN ") ) );

            
            // Test some invalid casts
            BOOST_CHECK_THROW( stringTo<int>(" 1234.5 "), boost::bad_lexical_cast );
            BOOST_CHECK_THROW( stringTo<bool>(" 1234.5 "), boost::bad_lexical_cast );

            

        }

        
        void add_string_tests( test_suite* ts ) {

            ts->add( BOOST_TEST_CASE_NAME( &stringTest, "String manipulations"  ) );

        }

    }

}

