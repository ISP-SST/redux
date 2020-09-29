
#include "redux/file/fileana.hpp"
#include "redux/image/image.hpp"
#include "redux/util/array.hpp"


#include "testsuite.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/test/unit_test.hpp>

using namespace redux::file;
using namespace redux::image;
using namespace redux::util;

using namespace std;
using namespace boost::unit_test;


#ifndef RDX_TESTDATA_DIR
#error RDX_TESTDATA_DIR not set
#endif


namespace testsuite {

    namespace file {


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
                                
        const char testFileAna[]  = "testsuite.f0";


        template <typename T>
        void readAnaAs( void ) {
            T data;
            for( size_t i = 0; i < 5; ++i ) {
                readFile( RDX_TESTDATA_DIR + anaFiles[i], data );
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
        
        template <typename T>
        void writeAndVerify( const string& filename, const T& indata ) {
            writeFile( testFileAna, indata );
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
        void writeAndVerifyCompressedAna( const Image<T>& indata ) {
            for(uint8_t sliceSize=1; sliceSize<8*sizeof(T); ++sliceSize) {       // test all values of sliceSize
                redux::file::Ana::write( testFileAna, indata, sliceSize );
                Image<T> data;
                readFile( testFileAna, data );
                auto hdr = static_pointer_cast<redux::file::Ana>(data.meta);
                if( hdr->m_Header.subf&1 ) {
                    BOOST_TEST( hdr->m_CompressedHeader.slice_size == sliceSize );
                }
                BOOST_TEST( data.nDimensions() == indata.nDimensions() );
                BOOST_TEST( data.dimSize( 0 ) == indata.dimSize( 0 ) );
                BOOST_TEST( data.dimSize( 1 ) == indata.dimSize( 1 ) );
                for( size_t j = 0; j < data.dimSize( 0 ); ++j ) {
                    for( size_t k = 0; k < data.dimSize( 1 ); ++k ) {
                        BOOST_TEST( data( j, k ) == indata( j, k ) );
                    }
                }
            }
        }


        void ana_test( void ) {

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
            readFile( RDX_TESTDATA_DIR + anaFiles[5], array );
            BOOST_TEST( array.nDimensions() == 2 );
            BOOST_TEST( array.dimSize( 0 ) == 40 );
            BOOST_TEST( array.dimSize( 1 ) == 50 );
            for( size_t j = 0; j < array.dimSize( 0 ); ++j ) {
                for( size_t k = 0; k < array.dimSize( 1 ); ++k ) {
                    BOOST_TEST( array( j, k ) == j + k );
                }
            }

            // test reading file saved on big-endian machine:
            readFile( RDX_TESTDATA_DIR + anaFiles[6], array );
            BOOST_TEST( array.nDimensions() == 2 );
            BOOST_TEST( array.dimSize( 0 ) == 40 );
            BOOST_TEST( array.dimSize( 1 ) == 50 );
            for( size_t j = 0; j < array.dimSize( 0 ); ++j ) {
                for( size_t k = 0; k < array.dimSize( 1 ); ++k ) {
                    BOOST_TEST( array( j, k ) == j + k );
                }
            }

            // test reading compressed file saved on little-endian machine:
            Image<int32_t> image;
            readFile( RDX_TESTDATA_DIR + anaFiles[7], image );
            BOOST_TEST( image.nDimensions(), 2 );
            BOOST_TEST( image.dimSize( 0 ) == 40 );
            BOOST_TEST( image.dimSize( 1 ) == 50 );
            for( size_t j = 0; j < image.dimSize( 0 ); ++j ) {
                for( size_t k = 0; k < image.dimSize( 1 ); ++k ) {
                    BOOST_TEST( image( j, k ) == j + k );
                }
            }

            // test reading compressed file saved on big-endian machine:
            // TODO: this test-file is broken, fix it.
        //     readFile( RDX_TESTDATA_DIR + anaFiles[8], array );
        //     BOOST_TEST( array.nDimensions(), 2 );
        //     BOOST_TEST( array.dimSize( 0 ), 40 );
        //     BOOST_TEST( array.dimSize( 1 ), 50 );
        //     for( size_t j = 0; j < array.dimSize( 0 ); ++j ) {
        //         for( size_t k = 0; k < array.dimSize( 1 ); ++k ) {
        //             BOOST_TEST( array( j, k ), j + k );
        //         }
        //     }

            writeAndVerify( testFileAna, array.copy<uint8_t>() );
            writeAndVerify( testFileAna, array.copy<int16_t>() );
            writeAndVerify( testFileAna, array );                       // int32_t
            //writeAndVerify( testFileAna, array.copy<int64_t>() );     // TODO: int64 fails, fix it.
            writeAndVerify( testFileAna, array.copy<float>() );
            writeAndVerify( testFileAna, array.copy<double>() );

            writeAndVerify( testFileAna, image.copy<uint8_t>() );
            writeAndVerify( testFileAna, image.copy<int16_t>() );
            writeAndVerify( testFileAna, image);                        // int32_t
            //writeAndVerify( testFileAna, image.copy<int64_t>() );     // TODO: int64 fails, fix it.
            writeAndVerify( testFileAna, image.copy<float>() );
            writeAndVerify( testFileAna, image.copy<double>() );


            writeAndVerifyCompressedAna( image.copy<uint8_t>() );
            writeAndVerifyCompressedAna( image.copy<int16_t>() );
            writeAndVerifyCompressedAna( image );                       // int32_t
            
            
            // test writing subimage
            image.resize(7,7);
            int cnt=0;
            for(auto& value: image) value = ++cnt;
            Image<int32_t> subimage(image,1,5,1,5);
            writeAndVerify( testFileAna, subimage );            
            writeAndVerifyCompressedAna( subimage );

            // test reading into subimage
            Image<int32_t> imagecopy = image.copy();
            Image<int32_t> subimagecopy(imagecopy,1,5,1,5);
            subimage *= 10;
            BOOST_CHECK( imagecopy != image );
            redux::file::Ana::write( testFileAna, subimage );
            redux::file::Ana::read( testFileAna, subimagecopy );
            BOOST_CHECK( imagecopy == image );
            
        }


        void add_ana_tests( test_suite* ts ) {

            ts->add( BOOST_TEST_CASE_NAME( &ana_test, "File ANA" ) );

        }

    }

}

