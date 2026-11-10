#include "redux/file/filefits.hpp"
#include "redux/image/image.hpp"

#include "stringdiff.hpp"

#include <boost/algorithm/string.hpp>

using namespace redux::file;
using namespace redux::image;
using namespace redux::util;
using namespace testsuite::file;

using namespace std;


#ifndef RDX_TESTDATA_DIR
#error RDX_TESTDATA_DIR not set
#endif


namespace testsuite {

    namespace file {

        namespace {

            const string fitsFiles[] = {"gradient_8u_4x5.fits",
                                        "gradient_16s_4x5.fits",
                                        "gradient_32s_4x5.fits",
                                        "gradient_32f_4x5.fits",
                                        "gradient_64f_4x5.fits",
                                        "gradient_32s_40x50.fits"
                                    };

            const char testFileFits[] = "testsuite.fits";

            template <typename T>
            void writeAndVerifyFits( const string& filename, const T& indata ) {
                writeFile( filename, indata );
                T data;
                readFile( filename, data );
                CHECK(  data.nDimensions() == indata.nDimensions()  );
                CHECK(  data.dimSize( 0 ) == indata.dimSize( 0 )  );
                CHECK(  data.dimSize( 1 ) == indata.dimSize( 1 )  );
                for( size_t j = 0; j < data.dimSize( 0 ); ++j ) {
                    for( size_t k = 0; k < data.dimSize( 1 ); ++k ) {
                        CHECK(  data( j, k ) == indata( j, k )  );
                    }
                }
            }

            template <typename T>
            void writeAndVerifyFitsHdr( const string& filename, const Array<T>& indata, const vector<string>& cards ) {

                shared_ptr<Fits> inHdr( new Fits() );
                inHdr->primaryHDU.cards = cards;

                Fits::write( filename, indata, inHdr );
                Array<T> data;
                shared_ptr<Fits> outHdr;
                Fits::read( filename, data, outHdr );

                CHECK(  data == indata  );
                CHECK(  data.nDimensions() == indata.nDimensions()  );
                CHECK(  data.dimSize( 0 ) == indata.dimSize( 0 )  );
                CHECK(  data.dimSize( 1 ) == indata.dimSize( 1 )  );
                for( size_t j = 0; j < data.dimSize( 0 ); ++j ) {
                    for( size_t k = 0; k < data.dimSize( 1 ); ++k ) {
                        CHECK(  data( j, k ) == indata( j, k )  );
                    }
                }

                const vector<string>& outCards = outHdr->primaryHDU.cards;
                // Verify no cards got lost in the write/read
                for( const auto& ic: inHdr->primaryHDU.cards ) {
                    CHECK( std::find( outCards.begin(), outCards.end(), ic ) != outCards.end() );
                }
            }

        }

        TEMPLATE_TEST_CASE( "FITS read/cast to various container types", "[file]",
                            Array<uint8_t>, Array<int16_t>, Array<int32_t>, Array<int64_t>, Array<float>, Array<double>,
                            Image<uint8_t>, Image<int16_t>, Image<int32_t>, Image<int64_t>, Image<float>, Image<double> ) {
            TestType data;
            for( size_t i = 0; i < 5; ++i ) {
                readFile( RDX_TESTDATA_DIR + fitsFiles[i], data );
                CHECK(  data.nDimensions() == 2  );
                CHECK(  data.dimSize( 0 ) == 4  );
                CHECK(  data.dimSize( 1 ) == 5  );
                for( size_t j = 0; j < data.dimSize( 0 ); ++j ) {
                    for( size_t k = 0; k < data.dimSize( 1 ); ++k ) {
                        CHECK(  data( j, k ) == j + k  );
                    }
                }
            }
        }

        TEST_CASE( "FITS read/write roundtrip across types", "[file]" ) {

            // read 32-bit signed integer file
            Array<int32_t> array;
            readFile( RDX_TESTDATA_DIR + fitsFiles[5], array );
            CHECK(  array.nDimensions() == 2  );
            CHECK(  array.dimSize( 0 ) == 40  );
            CHECK(  array.dimSize( 1 ) == 50  );
            for( size_t j = 0; j < array.dimSize( 0 ); ++j ) {
                for( size_t k = 0; k < array.dimSize( 1 ); ++k ) {
                    CHECK(  array( j, k ) == j + k  );
                }
            }

            writeAndVerifyFits( testFileFits, array.copy<uint8_t>() );
            writeAndVerifyFits( testFileFits, array.copy<int16_t>() );
            writeAndVerifyFits( testFileFits, array.copy<uint16_t>() );
            writeAndVerifyFits( testFileFits, array );                       // int32_t
            writeAndVerifyFits( testFileFits, array.copy<uint32_t>() );
            writeAndVerifyFits( testFileFits, array.copy<int64_t>() );
            writeAndVerifyFits( testFileFits, array.copy<uint64_t>() );
            writeAndVerifyFits( testFileFits, array.copy<float>() );
            writeAndVerifyFits( testFileFits, array.copy<double>() );

            vector<string> cards;
            Fits::addCard( cards, Fits::makeCard("Blaha",1,"No comment") );
            writeAndVerifyFitsHdr( testFileFits, array, cards );
        }

        TEST_CASE( "FITS cards", "[file]" ) {

            const string c("comment");
            const double d(-1234.567890123456789);
            const string  ss("strVal");                                                                                                         // short string
            const string  ls("Too long string to fit inside a single card, let's see if it will be truncated in the expected place etc.");      // long string

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
                        CHECK_THROWS_AS( Fits::makeCard( key ), std::domain_error );
                        CHECK_THROWS_AS( Fits::makeCard( key, "", true ), std::domain_error );
                    } else {
                        CHECK_NOTHROW( Fits::makeCard( key ) );
                        CHECK_NOTHROW( Fits::makeCard( key, "", true ) );
                    }

                    if( i<32 || i>126 ) {                                               // makeStringValue sohuld throw if ascii value is outside of [32, 126]
                        CHECK_THROWS_AS( Fits::makeCard( key0, value ), std::domain_error );
                        CHECK_THROWS_AS( Fits::makeCard( key0, key0, value ), std::domain_error );
                        CHECK_THROWS_AS( Fits::makeCard( key0, value, "", true ), std::domain_error );
                        CHECK_THROWS_AS( Fits::makeCard( key0, key0, value, true ), std::domain_error );
                    } else {
                        CHECK_NOTHROW( Fits::makeCard( key0, value ) );
                        CHECK_NOTHROW( Fits::makeCard( key0, key0, value ) );
                        CHECK_NOTHROW( Fits::makeCard( key0, value, "", true ) );
                        CHECK_NOTHROW( Fits::makeCard( key0, key0, value, true ) );
                    }

                }

                CHECK_NOTHROW( Fits::makeCard( "TooLongKeyN@me" ) );     // This should not throw, since key is truncated to 8 characters.

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

                CHECK_STRINGS_EQUAL( Fits::makeValue ( true ), vStrings[0] );
                CHECK_STRINGS_EQUAL( Fits::makeValue ( false ), vStrings[1] );
                CHECK_STRINGS_EQUAL( Fits::makeValue ( int8_t(-8) ), vStrings[2] );
                CHECK_STRINGS_EQUAL( Fits::makeValue ( uint8_t(-8) ), vStrings[3] );
                CHECK_STRINGS_EQUAL( Fits::makeValue ( int16_t(-16) ), vStrings[4] );
                CHECK_STRINGS_EQUAL( Fits::makeValue ( uint16_t(-16) ), vStrings[5] );
                CHECK_STRINGS_EQUAL( Fits::makeValue ( int32_t(-32) ), vStrings[6] );
                CHECK_STRINGS_EQUAL( Fits::makeValue ( uint32_t(-32) ), vStrings[7] );
                CHECK_STRINGS_EQUAL( Fits::makeValue ( int64_t(-64) ), vStrings[8] );
                CHECK_STRINGS_EQUAL( Fits::makeValue ( uint64_t(-64) ), vStrings[9] );
                CHECK_STRINGS_EQUAL( Fits::makeValue ( pt ), vStrings[10] );
                CHECK_STRINGS_EQUAL( Fits::makeValue ( ptt ), vStrings[11] );
                CHECK_STRINGS_EQUAL( Fits::makeValue( "a string  " ), "'a string'" );    // trailing spaces should be stripped
                CHECK_STRINGS_EQUAL( Fits::makeValue( "  a string" ), "'  a string'" );  // ...but not leading.

                CHECK(   Fits::getValue<bool>( vStrings[0] )  );
                CHECK(  !Fits::getValue<bool>( vStrings[1] )  );
                CHECK(  Fits::getValue<int8_t>( vStrings[2] ) == -8  );
                CHECK(  Fits::getValue<uint8_t>( vStrings[3] ) == 248  );
                CHECK(  Fits::getValue<int16_t>( vStrings[4] ) == -16  );
                CHECK(  Fits::getValue<uint16_t>( vStrings[5] ) == 65520  );
                CHECK(  Fits::getValue<int32_t>( vStrings[6] ) == -32  );
                CHECK(  Fits::getValue<uint32_t>( vStrings[7] ) == 4294967264  );
                CHECK(  Fits::getValue<int64_t>( vStrings[8] ) == -64  );
                CHECK(  Fits::getValue<uint64_t>( vStrings[9] ) == 18446744073709551552U  );
                CHECK(  Fits::getValue<ptime>( vStrings[10] ) == pt  );
                CHECK(  Fits::getValue<ptime>( vStrings[11] ) == ptt  );
                CHECK(  Fits::getValue<string>( "'a string'" ) == "a string"  );
                CHECK(  Fits::getValue<string>( "'  a string'" ) == "  a string"  );

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
                CHECK_STRINGS_EQUAL( Fits::makeCard("boolt",  true, c ), cards[0] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("boolt",  true, c, true), cards[0] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("boolf",  false, c ), cards[1] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("boolf",  false, c , true), cards[1] );

                CHECK_STRINGS_EQUAL( Fits::makeCard("int8s",  int8_t(-8), c ), cards[2] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int8s",  int8_t(-8), c, true ), cards[2] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int8u",  uint8_t(-8), c ), cards[3] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int8u",  uint8_t(-8), c, true ), cards[3] );

                CHECK_STRINGS_EQUAL( Fits::makeCard("int16s", int16_t(-16), c ), cards[4] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int16s", int16_t(-16), c, true ), cards[4] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int16u", uint16_t(-16), c ), cards[5] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int16u", uint16_t(-16), c, true ), cards[5] );

                CHECK_STRINGS_EQUAL( Fits::makeCard("int32s", int32_t(-32), c ), cards[6] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int32s", int32_t(-32), c, true ), cards[6] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int32u", uint32_t(-32), c ), cards[7] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int32u", uint32_t(-32), c, true ), cards[7] );

                CHECK_STRINGS_EQUAL( Fits::makeCard("int64s", int64_t(-64), c ), cards[8] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int64s", int64_t(-64), c, true ), cards[8] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int64u", uint64_t(-64), c ), cards[9] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("int64u", uint64_t(-64), c, true ), cards[9] );

                CHECK_STRINGS_EQUAL( Fits::makeCard("float",  float(d), c ), cards[10] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("float",  float(d), c, true ), cards[10] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("double", double(d), c ), cards[11] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("double", double(d), c, true ), cards[11] );

                CHECK_STRINGS_EQUAL( Fits::makeCard("END"), cards[12] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("END", "", true), cards[12] );

                CHECK_STRINGS_EQUAL( Fits::makeCard("LongComment", ls), cards[13] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("LongComment", ls, true), cards[14] );

                CHECK_STRINGS_EQUAL( Fits::makeCard("TooLongKey", ss, c), cards[15] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("TooLongKey", ss, c, true), cards[15] );

                CHECK_STRINGS_EQUAL( Fits::makeCard("LongSTR", ls, c), cards[16] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("LongSTR", ls, c, true), cards[17] );

                CHECK_STRINGS_EQUAL( Fits::makeCard("LongComment", ss, ls), cards[18] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("LongComment", ss, ls, true), cards[19] );

                CHECK_STRINGS_EQUAL( Fits::makeCard("NoValue", string(), ls), cards[20] );
                CHECK_STRINGS_EQUAL( Fits::makeCard("NoValue", string(), ls, true), cards[21] );
            }

            //CardRuler: |0       1|0       2|0       3|0       4|0       5|0       6|0       7|0       8|0
            CHECK(  "END                                                                             " == Fits::makeCard("END")  );
            CHECK(  "LONGCOMMToo long string to fit inside a single card, let's see if it will be tru" == Fits::makeCard("LongComment", ls)  );
            // makeCard should truncate too long keys/values/comments. TODO: a makeCardS routine to split long strings in multiple lines
            CHECK(  "TOOLONGK= 'strVal'             / comment                                        " == Fits::makeCard("TooLongKey", ss, c )  );
            CHECK(  "LONGSTR = 'Too long string to fit inside a single card, let''s see if it will b'" == Fits::makeCard("LongSTR", ls, c )  );
            CHECK(  "LONGCOMM= 'strVal'             / Too long string to fit inside a single card, le" == Fits::makeCard("LongComment", ss, ls )  );
            CHECK(  "NOVALUE Too long string to fit inside a single card, let's see if it will be tru" == Fits::makeCard("NoValue", string(), ls )  );

            // TODO lots of card-tests missing. Modification, change to/from long-cards etc.

        }

    }

}
