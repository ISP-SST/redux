#include <boost/test/unit_test.hpp>


#include "redux/image/image.hpp"
#include "redux/util/stringutil.hpp"


using namespace redux::image;
using namespace redux::util;

using namespace std;
using namespace boost::unit_test_framework;

namespace {
}

void imageTests( void ) {

    size_t xSize = 7;
    size_t ySize = 7;
    size_t firstX = 1;
    size_t lastX = 5;
    size_t firstY = 1;
    size_t lastY = 4;
    size_t xs = lastX - firstX + 1;
    size_t ys = lastY - firstY + 1;

    // create/resize
    Image<int> image( 10, 15 );
    BOOST_CHECK_EQUAL( image.dimSize(0), 10 );
    BOOST_CHECK_EQUAL( image.dimSize(1), 15 );
    image.resize( ySize, xSize );
    BOOST_CHECK_EQUAL( image.dimSize(0), ySize );
    BOOST_CHECK_EQUAL( image.dimSize(1), xSize );

    {
        // test comparison and assignment
        Image<int16_t> image16;
        image16 = image = 3;      // set all values to 3
        BOOST_CHECK_EQUAL( image16.dimSize(0), ySize );
        BOOST_CHECK_EQUAL( image16.dimSize(1), xSize );
        BOOST_CHECK( image16 == image );
        image( 1, 1 ) = 123;      // 1 wrong value
        BOOST_CHECK( !( image16 == image ) );

        // test that images of different dimensions compare as different.
        image.resize( ySize, xSize + 1 );
        image = 3;      // all values equal
        BOOST_CHECK( !( image16 == image ) );
        image.resize( ySize - 1, xSize );
        image = 3;      // all values equal
        BOOST_CHECK( !( image16 == image ) );

    }

    // test operators with scalars
    image = 1;
    for( auto it : image ) { BOOST_CHECK_EQUAL( it, 1 ); }
    image += 10;
    for( auto it : image ) { BOOST_CHECK_EQUAL( it, 11 ); }
    image -= 1;
    for( auto it : image ) { BOOST_CHECK_EQUAL( it, 10 ); }
    image *= 2;
    for( auto it : image ) { BOOST_CHECK_EQUAL( it, 20 ); }
    image /= 7;
    for( auto it : image ) { BOOST_CHECK_EQUAL( it, 2 ); }
    image.zero();
    for( auto it : image ) { BOOST_CHECK_EQUAL( it, 0 ); }


    /***** test iterators *****/
    {

        image.resize( ySize, xSize );
        image = 1;
        Image<int> subimage( image, firstY, lastY, firstX, lastX );
        BOOST_CHECK_EQUAL( subimage.dimSize(0), lastY-firstY+1 );
        BOOST_CHECK_EQUAL( subimage.dimSize(1), lastX-firstX+1 );
        for( size_t y = 0; y < ySize; ++y ) {
            for( size_t x = 0; x < xSize; ++x ) {
                image( y, x ) = y * xSize + x + 1;        // set some values: 1,2,3...
            }
        }

        Image<int>::iterator it = image.begin();
        Image<int>::iterator endit = image.end();
        Image<int>::const_iterator cit = image.begin();
        Image<int>::const_iterator cendit = image.end();

        // postfix ++
        for( size_t i = 1; i <= ySize * xSize; ++i ) {
            BOOST_CHECK_EQUAL( *it++, i );
            BOOST_CHECK_EQUAL( *cit++, i );
        }

        // at end()
        BOOST_CHECK( it == endit );
        BOOST_CHECK( cit == cendit );

        // prefix --
        for( size_t i = ySize * xSize; i > 0; --i ) {
            BOOST_CHECK_EQUAL( *--it, i );
            BOOST_CHECK_EQUAL( *--cit, i );
        }

        // at begin()
        BOOST_CHECK( it == image.begin() );
        BOOST_CHECK( cit == image.begin() );

        // prefix ++
        for( size_t i = 1; i <= ySize * xSize; ++i ) {
            BOOST_CHECK_EQUAL( *it, i );
            BOOST_CHECK_EQUAL( *cit, i );
            ++it;
            ++cit;
        }


        // at end()
        BOOST_CHECK( it == endit );
        BOOST_CHECK( cit == cendit );


        // postfix --
        for( int i = ySize * xSize; i > 0; --i ) {
            it--;
            cit--;
            BOOST_CHECK_EQUAL( *it, i );
            BOOST_CHECK_EQUAL( *cit, i );
        }


        // set/check values via auto.
        int cnt = 0;
        for( auto & ait : image ) {
            ait = ++cnt;
        }

        cnt = 0;
        for( auto ait : image ) {
            BOOST_CHECK_EQUAL( ait, ++cnt );
        }
    }
    /**************************/

    image.resize( ySize, xSize );

    Image<int> subimage( image, firstY/*1*/, lastY/*4*/, firstX/*1*/, lastX/*5*/ );

    // compare pointers to verify that the data is shared
    BOOST_CHECK( subimage.ptr() == image.ptr() );

    image = 2;
    for( size_t y = 0; y < ys; ++y ) {
        for( size_t x = 0; x < xs; ++x ) {
            int val = image( firstY + y, firstX + x );
            BOOST_CHECK_EQUAL( subimage( y, x ), val );
            subimage( y, x ) *= 2;
            BOOST_CHECK_EQUAL( image( firstY + y, firstX + x ), 2 * val );  // verify that the shared data has been modified
        }
    }

    // assign to subimage
    image = 2;
    subimage = 13;
    for( size_t y = 0; y < ySize; ++y ) {
        for( size_t x = 0; x < xSize; ++x ) {
            if( y >= firstY && y <= lastY && x >= firstX && x <= lastX ) {
                BOOST_CHECK_EQUAL( image( y, x ), 13 );
            }
            else {
                BOOST_CHECK_EQUAL( image( y, x ), 2 );
            }
        }
    }

    image = 5;
    Image<int16_t> image2(ySize,xSize);
    image2 = 2;
    image += image2;
    for( auto ait : image ) {
        BOOST_CHECK_EQUAL( ait, 7 );
    }

    image *= image2;
    for( auto ait : image ) {
        BOOST_CHECK_EQUAL( ait, 14 );
    }

    image -= image2;
    for( auto ait : image ) {
        BOOST_CHECK_EQUAL( ait, 12 );
    }

    image /= image2;
    for( auto ait : image ) {
        BOOST_CHECK_EQUAL( ait, 6 );
    }

}


namespace testsuite {

    namespace image {

        void imageTest( void ) {

            test_suite* ts = BOOST_TEST_SUITE( "IMAGE" );

            ts->add( BOOST_TEST_CASE( &imageTests ) );

            framework::master_test_suite().add( ts );

        }

    }

}
