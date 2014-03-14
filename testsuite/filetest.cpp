#include <boost/test/unit_test.hpp>


#include "redux/file/fileio.hpp"
#include "redux/file/fileana.hpp"
#include "redux/image/image.hpp"
#include "redux/util/array.hpp"
#include "redux/util/stringutil.hpp"


using namespace redux::file;
using namespace redux::image;
using namespace redux::util;

using namespace std;
using namespace boost::unit_test_framework;

#ifndef REDUX_TESTDATA_DIR
#error REDUX_TESTDATA_DIR not set
#endif


namespace {

    const string anaFiles[] = { "gradient_8u_4x5.f0",
                                "gradient_16s_4x5_le.f0",
                                "gradient_32s_4x5_le.f0",
                                "gradient_32f_4x5_le.f0",
                                "gradient_64f_4x5_le.f0",
                                "gradient_32s_40x50_le.f0",
                                "gradient_32s_40x50_be.f0",
                                "gradient_32s_40x50_le.fz",
                                "gradient_32s_40x50_be.fz"
                              };
                              
    const string testFileAna = "testsuite_ana.f0";

    template <typename T>
    void readAnaAs( void ) {
        T data;
        for( size_t i = 0; i < 5; ++i ) {
            readFile( REDUX_TESTDATA_DIR + anaFiles[i], data );
            BOOST_CHECK_EQUAL( data.nDimensions(), 2 );
            BOOST_CHECK_EQUAL( data.dimSize( 0 ), 4 );
            BOOST_CHECK_EQUAL( data.dimSize( 1 ), 5 );
            for( size_t j = 0; j < data.dimSize( 0 ); ++j ) {
                for( size_t k = 0; k < data.dimSize( 1 ); ++k ) {
                    BOOST_CHECK_EQUAL( data( j, k ), j + k );
                }
            }
        }
    }
    
    template <typename T>
    void writeAndVerifyAna( const T& indata ) {
        redux::file::Ana::write(testFileAna,indata);
        T data;
        readFile( testFileAna, data );
        BOOST_CHECK_EQUAL( data.nDimensions(), indata.nDimensions() );
        BOOST_CHECK_EQUAL( data.dimSize( 0 ), indata.dimSize( 0 ) );
        BOOST_CHECK_EQUAL( data.dimSize( 1 ), indata.dimSize( 1 ) );
        for( size_t j = 0; j < data.dimSize( 0 ); ++j ) {
            for( size_t k = 0; k < data.dimSize( 1 ); ++k ) {
                BOOST_CHECK_EQUAL( data( j, k ), indata( j, k ) );
            }
        }
    }

    template <typename T>
    void writeAndVerifyCompressedAna( const Image<T>& indata ) {
        for(int sliceSize=1; sliceSize<8*sizeof(T); ++sliceSize) {       // test all values of sliceSize
            redux::file::Ana::write(testFileAna,indata,sliceSize);
            Image<T> data;
            readFile( testFileAna, data );
            auto hdr = static_pointer_cast<redux::file::Ana>(data.hdr);
            if( hdr->m_Header.subf&1 ) {
                BOOST_CHECK_EQUAL( hdr->m_CompressedHeader.slice_size, sliceSize );
            }
            BOOST_CHECK_EQUAL( data.nDimensions(), indata.nDimensions() );
            BOOST_CHECK_EQUAL( data.dimSize( 0 ), indata.dimSize( 0 ) );
            BOOST_CHECK_EQUAL( data.dimSize( 1 ), indata.dimSize( 1 ) );
            for( size_t j = 0; j < data.dimSize( 0 ); ++j ) {
                for( size_t k = 0; k < data.dimSize( 1 ); ++k ) {
                    BOOST_CHECK_EQUAL( data( j, k ), indata( j, k ) );
                }
            }
        }
    }


}


void anaTest( void ) {

    // test reading/casting to various types.
    readAnaAs<Array<uint8_t>>();
    readAnaAs<Array<int16_t>>();
    readAnaAs<Array<int32_t>>();
    readAnaAs<Array<int64_t>>();
    readAnaAs<Array<float>>();
    readAnaAs<Array<double>>();

    readAnaAs<Image<uint8_t>>();
    readAnaAs<Image<int16_t>>();
    readAnaAs<Image<int32_t>>();
    readAnaAs<Image<int64_t>>();
    readAnaAs<Image<float>>();
    readAnaAs<Image<double>>();
    
    // test reading file saved on little-endian machine:
    Array<int32_t> array;
    readFile( REDUX_TESTDATA_DIR + anaFiles[5], array );
    BOOST_CHECK_EQUAL( array.nDimensions(), 2 );
    BOOST_CHECK_EQUAL( array.dimSize( 0 ), 40 );
    BOOST_CHECK_EQUAL( array.dimSize( 1 ), 50 );
    for( size_t j = 0; j < array.dimSize( 0 ); ++j ) {
        for( size_t k = 0; k < array.dimSize( 1 ); ++k ) {
            BOOST_CHECK_EQUAL( array( j, k ), j + k );
        }
    }

    // test reading file saved on big-endian machine:
    readFile( REDUX_TESTDATA_DIR + anaFiles[6], array );
    BOOST_CHECK_EQUAL( array.nDimensions(), 2 );
    BOOST_CHECK_EQUAL( array.dimSize( 0 ), 40 );
    BOOST_CHECK_EQUAL( array.dimSize( 1 ), 50 );
    for( size_t j = 0; j < array.dimSize( 0 ); ++j ) {
        for( size_t k = 0; k < array.dimSize( 1 ); ++k ) {
            BOOST_CHECK_EQUAL( array( j, k ), j + k );
        }
    }

    // test reading compressed file saved on little-endian machine:
    Image<int32_t> image;
    readFile( REDUX_TESTDATA_DIR + anaFiles[7], image );
    BOOST_CHECK_EQUAL( image.nDimensions(), 2 );
    BOOST_CHECK_EQUAL( image.dimSize( 0 ), 40 );
    BOOST_CHECK_EQUAL( image.dimSize( 1 ), 50 );
    for( size_t j = 0; j < image.dimSize( 0 ); ++j ) {
        for( size_t k = 0; k < image.dimSize( 1 ); ++k ) {
            BOOST_CHECK_EQUAL( image( j, k ), j + k );
        }
    }

    // test reading compressed file saved on big-endian machine:
    // TODO: this test-file is broken, fix it.
//     readFile( REDUX_TESTDATA_DIR + anaFiles[8], array );
//     BOOST_CHECK_EQUAL( array.nDimensions(), 2 );
//     BOOST_CHECK_EQUAL( array.dimSize( 0 ), 40 );
//     BOOST_CHECK_EQUAL( array.dimSize( 1 ), 50 );
//     for( size_t j = 0; j < array.dimSize( 0 ); ++j ) {
//         for( size_t k = 0; k < array.dimSize( 1 ); ++k ) {
//             BOOST_CHECK_EQUAL( array( j, k ), j + k );
//         }
//     }

    writeAndVerifyAna(array.copy<uint8_t>());
    writeAndVerifyAna(array.copy<int16_t>());
    writeAndVerifyAna(array);                       // int32_t
    //writeAndVerifyAna(array.copy<int64_t>());     // TODO: int64 fails, fix it.
    writeAndVerifyAna(array.copy<float>());
    writeAndVerifyAna(array.copy<double>());

    writeAndVerifyAna(image.copy<uint8_t>());
    writeAndVerifyAna(image.copy<int16_t>());
    writeAndVerifyAna(image);                       // int32_t
    //writeAndVerifyAna(image.copy<int64_t>());     // TODO: int64 fails, fix it.
    writeAndVerifyAna(image.copy<float>());
    writeAndVerifyAna(image.copy<double>());


    writeAndVerifyCompressedAna(image.copy<uint8_t>());
    writeAndVerifyCompressedAna(image.copy<int16_t>());
    writeAndVerifyCompressedAna(image);                       // int32_t
    
    
    // test writing subimage
    image.resize(7,7);
    int cnt=0;
    for(auto& it: image) it = ++cnt;
    Image<int32_t> subimage(image,1,5,1,5);
    writeAndVerifyAna(subimage);            
    writeAndVerifyCompressedAna(subimage);


}


namespace testsuite {

    namespace file {

        void fileTest( void ) {

            test_suite* ts = BOOST_TEST_SUITE( "FILE" );

            ts->add( BOOST_TEST_CASE( &anaTest ) );

            framework::master_test_suite().add( ts );

        }

    }

}
