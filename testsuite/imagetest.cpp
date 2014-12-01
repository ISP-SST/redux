#include <boost/test/unit_test.hpp>


#include "redux/image/image.hpp"
#include "redux/image/utils.hpp"
#include "redux/util/stringutil.hpp"


using namespace redux::image;
using namespace redux::util;

using namespace std;
using namespace boost::unit_test_framework;

namespace {
    
}

void imageTests( void ) {

    // size of test image
    size_t xSize = 9;
    size_t ySize = 7;
    
    // location of subimage
    size_t subFirstX = 2;
    size_t subLastX = 5;
    size_t subFirstY = 3;
    size_t subLastY = 4;
    
    size_t subXSize = subLastX - subFirstX + 1;
    size_t subYSize = subLastY - subFirstY + 1;

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
        Image<int> subimage( image, subFirstY, subLastY, subFirstX, subLastX );
        BOOST_CHECK_EQUAL( subimage.dimSize(0), subLastY-subFirstY+1 );
        BOOST_CHECK_EQUAL( subimage.dimSize(1), subLastX-subFirstX+1 );
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
    Image<int> subimage( image, subFirstY, subLastY, subFirstX, subLastX );

    // compare pointers to verify that the data is shared
    BOOST_CHECK( subimage.ptr() == image.ptr() );
    
    // check that original is dense, and subimage not
    BOOST_CHECK( image.dense() );
    BOOST_CHECK( !subimage.dense() );

    // global assign to subimage
    image = 2;
    subimage = 13;
    for( size_t y = 0; y < ySize; ++y ) {
        for( size_t x = 0; x < xSize; ++x ) {
            if( y >= subFirstY && y <= subLastY && x >= subFirstX && x <= subLastX ) {
                BOOST_CHECK_EQUAL( image( y, x ), 13 );
            }
            else {
                BOOST_CHECK_EQUAL( image( y, x ), 2 );
            }
        }
    }

    // elementwise modification
    image = 2;
    for( size_t y = 0; y < subYSize; ++y ) {
        for( size_t x = 0; x < subXSize; ++x ) {
            int val = image( subFirstY + y, subFirstX + x );
            BOOST_CHECK_EQUAL( subimage( y, x ), val );
            subimage( y, x ) *= 2;
            BOOST_CHECK_EQUAL( image( subFirstY + y, subFirstX + x ), 2 * val );  // verify that the shared data has been modified
        }
    }
    
    
    // shift subimage
    {
        int cnt=1;
        for( size_t y = 0; y < ySize; ++y ) {
            for( size_t x = 0; x < xSize; ++x ) {
                image( y, x ) = cnt++;
            }
        }
        
        // check that the shift is limited by the borders of the underlying array.
        BOOST_CHECK_EQUAL( subimage.shift(0,-1000), -subFirstY);
        BOOST_CHECK_EQUAL( subimage.shift(1,-1000), -subFirstX);
        // now the subimage should start at (0,0) and end at (subYSize-1,subYSize-1)
        BOOST_CHECK_EQUAL( subimage( 0, 0 ), image( 0, 0 ) );
        BOOST_CHECK_EQUAL( subimage( subYSize-1, subXSize-1 ), image( subYSize-1, subXSize-1) );
        
        BOOST_CHECK_EQUAL( subimage.shift(0,1000), ySize-subYSize);
        BOOST_CHECK_EQUAL( subimage.shift(1,1000), xSize-subXSize);
        // now the subimage should start at (ySize-subYSize,xSize-subXSize) and end at (ySize-1,xSize-1)
        BOOST_CHECK_EQUAL( subimage( 0, 0 ), image( ySize-subYSize, xSize-subXSize ) );
        BOOST_CHECK_EQUAL( subimage( subYSize-1, subXSize-1 ), image( ySize-1, xSize-1) );
        
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
    
    
    // Grids
    // TODO: Better testcase, this is basically the same math as int the Grid-constructor.
    {
        size_t nPoints = 50;
        float originX = 11;
        float originY = 11.5;
        Grid grid(nPoints,originY,originX);
        float** distPtr = grid.distance.get();
        float** anglePtr = grid.angle.get();
        for( size_t i = 0; i < nPoints; ++i ) {
            double yDist = i - originY;
            for( size_t j = 0; j < nPoints; ++j ) {
                double xDist = j - originX;
                BOOST_CHECK_CLOSE( distPtr[i][j], (float)sqrt(yDist*yDist+xDist*xDist), 1E-10);
                if(yDist || xDist) {
                    BOOST_CHECK_CLOSE( anglePtr[i][j], (float)atan2(yDist, xDist), 1E-5 );
                } else  BOOST_CHECK_EQUAL( anglePtr[i][j], 0.0 );
            }
        }

        
    }

}

void transformTest( void ) {
    
    size_t sizeX = 50;
    size_t sizeY = 14;
    size_t sizeZ = 5;
    Array<int> array( sizeZ, sizeY, sizeX );
    
    int cnt = 0;
    for(auto& it: array) it = ++cnt;
  
    // check values
    for( size_t i = 0; i < sizeZ; ++i ) {
        for( size_t j = 0; j < sizeY; ++j ) {
            for( size_t k = 0; k < sizeX; ++k ) {
                BOOST_CHECK_EQUAL( array( i, j, k ), sizeX*sizeY*i + j*sizeX + k + 1 );
            }
        }
    }
    
    shared_ptr<int**> sharedPtr = array.get(sizeZ,sizeY,sizeX);
    int*** rawPtr = sharedPtr.get();
    // check reverseX
    for( size_t i = 0; i < sizeZ; ++i ) reverseX(rawPtr[i],sizeY,sizeX);
    for( size_t i = 0; i < sizeZ; ++i ) {
        for( size_t j = 0; j < sizeY; ++j ) {
            for( size_t k = 0; k < sizeX; ++k ) {
                BOOST_CHECK_EQUAL( array( i, j, k ), sizeX*sizeY*i + j*sizeX + (sizeX-k-1) + 1 );
            }
        }
    }
    for( size_t i = 0; i < sizeZ; ++i ) reverseX(rawPtr[i],sizeY,sizeX);        // flip back
    
    // check reverseY
    for( size_t i = 0; i < sizeZ; ++i ) reverseY(rawPtr[i],sizeY,sizeX);
    for( size_t i = 0; i < sizeZ; ++i ) {
        for( size_t j = 0; j < sizeY; ++j ) {
            for( size_t k = 0; k < sizeX; ++k ) {
                BOOST_CHECK_EQUAL( array( i, j, k ), sizeX*sizeY*i + (sizeY-j-1)*sizeX + k + 1 );
            }
        }
    }
    for( size_t i = 0; i < sizeZ; ++i ) reverseY(rawPtr[i],sizeY,sizeX);        // flip back
    
    Array<int> array2 = array.copy();           // make a copy for comparison

    // check transpose
    for( size_t i = 0; i < sizeZ; ++i ) transpose(*rawPtr[i],sizeY,sizeX);
    array.permuteDimensions(1,2);
    for( size_t i = 0; i < sizeZ; ++i ) {
        for( size_t j = 0; j < sizeY; ++j ) {
            for( size_t k = 0; k < sizeX; ++k ) {
                BOOST_CHECK_EQUAL( array(i,k,j), array2(i,j,k) );
            }
        }
    }
    
}



namespace testsuite {

    namespace image {

        void imageTest( void ) {

            test_suite* ts = BOOST_TEST_SUITE( "IMAGE" );

            ts->add( BOOST_TEST_CASE( &imageTests ) );
            ts->add( BOOST_TEST_CASE( &transformTest ) );

            framework::master_test_suite().add( ts );

        }

    }

}
