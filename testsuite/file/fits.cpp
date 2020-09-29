
#include "redux/file/filefits.hpp"
#include "redux/image/image.hpp"

#include "testsuite.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/test/unit_test.hpp>

using namespace redux::file;
using namespace redux::image;
using namespace redux::util;
using namespace testsuite::file;

using namespace std;
using namespace boost::unit_test;


#ifndef RDX_TESTDATA_DIR
#error RDX_TESTDATA_DIR not set
#endif


namespace testsuite {

    namespace file {

                            
        const string fitsFiles[] = {"gradient_8u_4x5.fits",
                                    "gradient_16s_4x5.fits",
                                    "gradient_32s_4x5.fits",
                                    "gradient_32f_4x5.fits",
                                    "gradient_64f_4x5.fits",
                                    "gradient_32s_40x50.fits"
                                };
                                
        const char testFileFits[] = "testsuite.fits";
        const char testFileHeader[] = "header.fits";

        
        template <typename T>
        void writeAndVerify( const string& filename, const T& indata ) {
            writeFile( filename, indata );
            T data;
            readFile( filename, data );
            BOOST_TEST( data.nDimensions() == indata.nDimensions() );
            BOOST_TEST( data.dimSize( 0 ) == indata.dimSize( 0 ) );
            BOOST_TEST( data.dimSize( 1 ) == indata.dimSize( 1 ) );
            for( size_t j = 0; j < data.dimSize( 0 ); ++j ) {
                for( size_t k = 0; k < data.dimSize( 1 ); ++k ) {
                    BOOST_TEST( data( j, k ) == indata( j, k ) );
                }
            }
        }

        
        template <typename T>
        void readFitsAs( void ) {
            T data;
            for( size_t i = 0; i < 5; ++i ) {
                readFile( RDX_TESTDATA_DIR + fitsFiles[i], data );
                BOOST_TEST( data.nDimensions() == 2 );
                BOOST_TEST( data.dimSize( 0 ) == 4 );
                BOOST_TEST( data.dimSize( 1 ) == 5 );
                for( size_t j = 0; j < data.dimSize( 0 ); ++j ) {
                    for( size_t k = 0; k < data.dimSize( 1 ); ++k ) {
                        BOOST_TEST( data( j, k ) == j + k );
                    }
                }
            }
        }
        
        
        void testFitsCards( void ) {
            
            const string c("comment");
            const double d(-1234.567890123456789);
            const string  ss("strVal");                                                                                                         // short string
            const string  ls("Too long string to fit inside a single card, let's see if it will be truncated in the expected place etc.");      // long string
            const string  ls2("Too long string to fit inside a single card, let''s see if it will be truncated in the expected place etc.");    // escapedlong string
            
            {
                string key0 = "key";
                string val0 = "value";
                const string allowed_key_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ01923456789-_";
                for( int i(0); i<256; ++i ) {
                    string key = key0;
                    string value = val0;
                    key[1] = char(i);       // insert character in the middle of key/value (to avoid being trimmed away)
                    value[2] = char(i);
                    boost::trim( key );                                                 // reamove leading/trailing spaces
                    boost::to_upper( key );                                             // make uppercase
                    
                    size_t found = key.find_first_not_of( allowed_key_chars );          // look for illegal characters
                    if( found != string::npos ) {                                       // makeKey should trow if 
                        BOOST_CHECK_THROW( Fits::makeCard( key ), std::domain_error );
                        BOOST_CHECK_THROW( Fits::makeCard( key, "", true ), std::domain_error );
                    } else {
                        BOOST_CHECK_NO_THROW( Fits::makeCard( key ) );
                        BOOST_CHECK_NO_THROW( Fits::makeCard( key, "", true ) );
                    }

                    if( i<32 || i>126 ) {                                               // makeStringValue sohuld throw if ascii value is outside of [32, 126]
                        BOOST_CHECK_THROW( Fits::makeCard( key0, value ), std::domain_error );
                        BOOST_CHECK_THROW( Fits::makeCard( key0, key0, value ), std::domain_error );
                        BOOST_CHECK_THROW( Fits::makeCard( key0, value, "", true ), std::domain_error );
                        BOOST_CHECK_THROW( Fits::makeCard( key0, key0, value, true ), std::domain_error );
                    } else {
                        BOOST_CHECK_NO_THROW( Fits::makeCard( key0, value ) );
                        BOOST_CHECK_NO_THROW( Fits::makeCard( key0, key0, value ) );
                        BOOST_CHECK_NO_THROW( Fits::makeCard( key0, value, "", true ) );
                        BOOST_CHECK_NO_THROW( Fits::makeCard( key0, key0, value, true ) );
                    }
                    
                }
                
                BOOST_CHECK_NO_THROW( Fits::makeCard( "TooLongKeyN@me" ) );     // This should not throw, since key is truncated to 8 characters.
                
            }
            
            
            {
                using namespace boost::posix_time;
                ptime pt( boost::gregorian::date( 2020, 9, 24 ) );
                ptime ptt = pt + time_duration(23,03,37);
                
                const char* vStrings[] = { "                   T",
                                        "                   F",
                                        "                  -8",
                                        "                 248",
                                        "                 -16",
                                        "               65520",  // 5
                                        "                 -32",
                                        "          4294967264",
                                        "                 -64",
                                        "18446744073709551552",
                                        "'2020-09-24'",          // 10
                                        "'2020-09-24T23:03:37'"
                };

                BOOST_CHECK( compare_strings( Fits::makeValue ( true ),           vStrings[0] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue ( false ),          vStrings[1] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue ( int8_t(-8) ),     vStrings[2] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue ( uint8_t(-8) ),    vStrings[3] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue ( int16_t(-16) ),   vStrings[4] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue ( uint16_t(-16) ),  vStrings[5] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue ( int32_t(-32) ),   vStrings[6] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue ( uint32_t(-32) ),  vStrings[7] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue ( int64_t(-64) ),   vStrings[8] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue ( uint64_t(-64) ),  vStrings[9] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue ( pt ),             vStrings[10] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue ( ptt ),            vStrings[11] ) );
                BOOST_CHECK( compare_strings( Fits::makeValue( "a string  " ),   "'a string'" ) );    // trailing spaces should be stripped
                BOOST_CHECK( compare_strings( Fits::makeValue( "  a string" ),   "'  a string'" ) );  // ...but not leading.
                
                BOOST_CHECK(  Fits::getValue<bool>( vStrings[0] ) );
                BOOST_CHECK( !Fits::getValue<bool>( vStrings[1] ) );
                BOOST_CHECK( Fits::getValue<int8_t>( vStrings[2] ) == -8 );
                BOOST_CHECK( Fits::getValue<uint8_t>( vStrings[3] ) == 248 );
                BOOST_CHECK( Fits::getValue<int16_t>( vStrings[4] ) == -16 );
                BOOST_CHECK( Fits::getValue<uint16_t>( vStrings[5] ) == 65520 );
                BOOST_CHECK( Fits::getValue<int32_t>( vStrings[6] ) == -32 );
                BOOST_CHECK( Fits::getValue<uint32_t>( vStrings[7] ) == 4294967264 );
                BOOST_CHECK( Fits::getValue<int64_t>( vStrings[8] ) == -64 );
                BOOST_CHECK( Fits::getValue<uint64_t>( vStrings[9] ) == 18446744073709551552U );
                BOOST_CHECK( Fits::getValue<ptime>( vStrings[10] ) == pt );
                BOOST_CHECK( Fits::getValue<ptime>( vStrings[11] ) == ptt );
                BOOST_CHECK( Fits::getValue<string>( "'a string'" ) == "a string" );
                BOOST_CHECK( Fits::getValue<string>( "'  a string'" ) == "  a string" );

            }

        
            
            {
                const char* cards[] = {
                //    CardRuler: |0       1|0       2|0       3|0       4|0       5|0       6|0       7|0       8|0
                                "BOOLT   =                    T / comment                                        ",
                                "BOOLF   =                    F / comment                                        ",
                                "INT8S   =                   -8 / comment                                        ",
                                "INT8U   =                  248 / comment                                        ",
                                "INT16S  =                  -16 / comment                                        ",
                                "INT16U  =                65520 / comment                                        ", // 5
                                "INT32S  =                  -32 / comment                                        ",
                                "INT32U  =           4294967264 / comment                                        ",
                                "INT64S  =                  -64 / comment                                        ",
                                "INT64U  = 18446744073709551552 / comment                                        ",
                                "FLOAT   =         -1234.567871 / comment                                        ", // 10
                                "DOUBLE  =         -1234.567890 / comment                                        ",
                                "END                                                                             ",
                                "LONGCOMMToo long string to fit inside a single card, let's see if it will be tru",
                                "LONGCOMMToo long string to fit inside a single card, let's see if it will be truncated in the expected place etc.",
                                "TOOLONGK= 'strVal'             / comment                                        ", // 15
                                "LONGSTR = 'Too long string to fit inside a single card, let''s see if it will b'",
                                "LONGSTR = 'Too long string to fit inside a single card, let''s see if it will be truncated in the expected place etc.' / comment",
                                "LONGCOMM= 'strVal'             / Too long string to fit inside a single card, le",
                                "LONGCOMM= 'strVal'             / Too long string to fit inside a single card, let's see if it will be truncated in the expected place etc.",
                                "NOVALUE Too long string to fit inside a single card, let's see if it will be tru", // 20
                                "NOVALUE Too long string to fit inside a single card, let's see if it will be truncated in the expected place etc."
                };

                // Test makeCard
                BOOST_CHECK( compare_strings( Fits::makeCard("boolt",  true, c ),                cards[0] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("boolt",  true, c, true),           cards[0] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("boolf",  false, c ),               cards[1] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("boolf",  false, c , true),         cards[1] ) );
                
                BOOST_CHECK( compare_strings( Fits::makeCard("int8s",  int8_t(-8), c ),          cards[2] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int8s",  int8_t(-8), c, true ),    cards[2] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int8u",  uint8_t(-8), c ),         cards[3] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int8u",  uint8_t(-8), c, true ),   cards[3] ) );
                
                BOOST_CHECK( compare_strings( Fits::makeCard("int16s", int16_t(-16), c ),        cards[4] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int16s", int16_t(-16), c, true ),  cards[4] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int16u", uint16_t(-16), c ),       cards[5] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int16u", uint16_t(-16), c, true ), cards[5] ) );
                
                BOOST_CHECK( compare_strings( Fits::makeCard("int32s", int32_t(-32), c ),        cards[6] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int32s", int32_t(-32), c, true ),  cards[6] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int32u", uint32_t(-32), c ),       cards[7] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int32u", uint32_t(-32), c, true ), cards[7] ) );
                
                BOOST_CHECK( compare_strings( Fits::makeCard("int64s", int64_t(-64), c ),        cards[8] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int64s", int64_t(-64), c, true ),  cards[8] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int64u", uint64_t(-64), c ),       cards[9] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("int64u", uint64_t(-64), c, true ), cards[9] ) );
                
                BOOST_CHECK( compare_strings( Fits::makeCard("float",  float(d), c ),            cards[10] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("float",  float(d), c, true ),      cards[10] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("double", double(d), c ),           cards[11] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("double", double(d), c, true ),     cards[11] ) );
                
                BOOST_CHECK( compare_strings( Fits::makeCard("END"),                             cards[12] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("END", "", true),                   cards[12] ) );
                
                BOOST_CHECK( compare_strings( Fits::makeCard("LongComment", ls),                 cards[13] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("LongComment", ls, true),           cards[14] ) );
                
                BOOST_CHECK( compare_strings( Fits::makeCard("TooLongKey", ss, c),               cards[15] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("TooLongKey", ss, c, true),         cards[15] ) );
                
                BOOST_CHECK( compare_strings( Fits::makeCard("LongSTR", ls, c),                  cards[16] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("LongSTR", ls, c, true),            cards[17] ) );
                
                BOOST_CHECK( compare_strings( Fits::makeCard("LongComment", ss, ls),             cards[18] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("LongComment", ss, ls, true),       cards[19] ) );
                
                BOOST_CHECK( compare_strings( Fits::makeCard("NoValue", string(), ls),           cards[20] ) );
                BOOST_CHECK( compare_strings( Fits::makeCard("NoValue", string(), ls, true),     cards[21] ) );
            }

            //CardRuler: |0       1|0       2|0       3|0       4|0       5|0       6|0       7|0       8|0
            BOOST_TEST( "END                                                                             " == Fits::makeCard("END") );
            BOOST_TEST( "LONGCOMMToo long string to fit inside a single card, let's see if it will be tru" == Fits::makeCard("LongComment", ls) );
            // makeCard should truncate too long keys/values/comments. TODO: a makeCardS routine to split long strings in multiple lines
            BOOST_TEST( "TOOLONGK= 'strVal'             / comment                                        " == Fits::makeCard("TooLongKey", ss, c ) );
            BOOST_TEST( "LONGSTR = 'Too long string to fit inside a single card, let''s see if it will b'" == Fits::makeCard("LongSTR", ls, c ) );
            BOOST_TEST( "LONGCOMM= 'strVal'             / Too long string to fit inside a single card, le" == Fits::makeCard("LongComment", ss, ls ) );
            BOOST_TEST( "NOVALUE Too long string to fit inside a single card, let's see if it will be tru" == Fits::makeCard("NoValue", string(), ls ) );

            

            //CardRuler: |0       1|0       2|0       3|0       4|0       5|0       6|0       7|0       8|0


            // TODO lots of card-tests missing. Modification, change to/from long-cards etc.
            
        }
        
        
        template <typename T>
        void writeAndVerifyFitsHdr( const string& filename, const Array<T>& indata, const vector<string>& cards ) {

            shared_ptr<Fits> inHdr( new Fits() );
            inHdr->primaryHDU.cards = cards;

            Fits::write( filename, indata, inHdr );
            Array<T> data;
            shared_ptr<Fits> outHdr;
            Fits::read( filename, data, outHdr );

            BOOST_CHECK( data == indata );
            BOOST_TEST( data.nDimensions() == indata.nDimensions() );
            BOOST_TEST( data.dimSize( 0 ) == indata.dimSize( 0 ) );
            BOOST_TEST( data.dimSize( 1 ) == indata.dimSize( 1 ) );
            for( size_t j = 0; j < data.dimSize( 0 ); ++j ) {
                for( size_t k = 0; k < data.dimSize( 1 ); ++k ) {
                    BOOST_TEST( data( j, k ) == indata( j, k ) );
                }
            }
            
            const vector<string>& outCards = outHdr->primaryHDU.cards;
            // Verify no cards got lost in the write/read
            for( const auto& ic: inHdr->primaryHDU.cards ) {
                auto it = std::find( outCards.begin(), outCards.end(), ic );
                if( it == outCards.end() ) {
                    cout << "Output does not contain card: \"" << ic << "\"" << endl;   // TODO proper test, not just message
                } //else cout << "Card: \"" << ic << "\" -> OK!" << endl;
                //BOOST_ASSERT_MSG( it != outCards.end(), "Output does not contain card: \"" + ic + "\"" );
                //BOOST_TEST( std::find( outCards.begin(), outCards.end(), ic) != outCards.end() );
            }
        }



        void fits_test( void ) {

            // test reading/casting to various types.
            readFitsAs<Array<uint8_t>>();
            readFitsAs<Array<int16_t>>();
            readFitsAs<Array<int32_t>>();
            readFitsAs<Array<int64_t>>();
            readFitsAs<Array<float>>();
            readFitsAs<Array<double>>();

            readFitsAs<Image<uint8_t>>();
            readFitsAs<Image<int16_t>>();
            readFitsAs<Image<int32_t>>();
            readFitsAs<Image<int64_t>>();
            readFitsAs<Image<float>>();
            readFitsAs<Image<double>>();
            
            // read 32-bit signed integer file 
            Array<int32_t> array;
            readFile( RDX_TESTDATA_DIR + fitsFiles[5], array );
            BOOST_TEST( array.nDimensions() == 2 );
            BOOST_TEST( array.dimSize( 0 ) == 40 );
            BOOST_TEST( array.dimSize( 1 ) == 50 );
            for( size_t j = 0; j < array.dimSize( 0 ); ++j ) {
                for( size_t k = 0; k < array.dimSize( 1 ); ++k ) {
                    BOOST_TEST( array( j, k ) == j + k );
                }
            }

            writeAndVerify( testFileFits, array.copy<uint8_t>() );
            //writeAndVerify( string("dbg_")+testFileFits, array.copy<int8_t>() );  // FIXME: causes illegal bitpix value on close.
            writeAndVerify( testFileFits, array.copy<int16_t>() );
            writeAndVerify( testFileFits, array.copy<uint16_t>() );
            writeAndVerify( testFileFits, array );                       // int32_t
            writeAndVerify( testFileFits, array.copy<uint32_t>() );
            writeAndVerify( testFileFits, array.copy<int64_t>() );
            writeAndVerify( testFileFits, array.copy<uint64_t>() );
            writeAndVerify( testFileFits, array.copy<float>() );
            writeAndVerify( testFileFits, array.copy<double>() );
            
            testFitsCards();

            vector<string> cards;
            Fits::addCard( cards, Fits::makeCard("Blaha",1,"No comment") );
            for( size_t i(0); i<1; ++i ) {
                writeAndVerifyFitsHdr( testFileFits, array, cards );
            }
        }


        void add_fits_tests( test_suite* ts ) {

            ts->add( BOOST_TEST_CASE_NAME( &fits_test, "File FITS" ) );

        }

    }

}
