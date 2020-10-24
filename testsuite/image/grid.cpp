// #include "redux/file/fileana.hpp"
// #include "redux/image/fouriertransform.hpp"
// #include "redux/image/image.hpp"
#include "redux/image/grid.hpp"
//#include "redux/util/cache.hpp"
// #include "redux/image/utils.hpp"
// #include "redux/math/functions.hpp"
// #include "redux/util/stringutil.hpp"
// 
// #include <iostream>
// 
// //#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>
// 
// using namespace redux::file;
using namespace redux::image;
using namespace redux::util;
// using namespace redux;
using namespace std;


namespace testsuite {

    namespace image {

        
        void testGrid( int xSize, int ySize=0, bool randomOrigin=true ) {
            
            if( ySize == 0 ) ySize = xSize;
            
            double originX(0);      // if origin is integer, the generated dist/angle will be centered on the border between the specified pixels.
            double originY(0);      // if origin is integer+0.5, the generated dist/angle will be centered at the middle of the specified pixel.
                                    // The generated dist/angle values are calculated for the middle of the pixels.

            if( randomOrigin ) {
                originX = (rand()%(300*xSize) - 100*xSize)/300.0;     // pick a random origin within [-nPoints, 2*nPoints] area, with 0.01 resolution
                originY = (rand()%(300*ySize) - 100*ySize)/300.0;
            }
            shared_ptr<Grid> grid = Grid::get( ySize, xSize, originY, originX );
            auto dist2D = grid->dist2D();
            auto angle2D = grid->angle2D();
            double** distPtr = dist2D.get();
            double** anglePtr = angle2D.get();

            for( int y(0); y < ySize; ++y ) {
                double yDist = y + 0.5 - originY;           // where the origin "should" be, by the design of the Grid-class
                for( int x(0); x < xSize; ++x ) {
                    double xDist = x + 0.5 - originX;       // where the origin "should" be, by the design of the Grid-class
                    BOOST_CHECK_SMALL( distPtr[y][x] - sqrt(yDist*yDist+xDist*xDist), 0.0001 );
                    if( yDist || xDist ) {
                        BOOST_CHECK_SMALL( anglePtr[y][x] - atan2(yDist, xDist), 0.0001 );
                    } else {
                        BOOST_CHECK_SMALL( anglePtr[y][x], 0.0001 );
                    }
                }
            }
            
        }

        void grid_test( void ) {

            // Test generated grids for some different sizes (even/odd/large/small/...)
            vector<int> variousSizes = { 50, 49, 500, 497, 9, 8 };
            size_t nVS = variousSizes.size();
            
            for( auto& i: variousSizes ) testGrid( i );
            
            size_t nGrids = Cache::size<Grid::ID,shared_ptr<Grid>>();   // All the grids should be cached
            BOOST_TEST( nGrids == nVS );
            
            Cache::clear<Grid::ID,shared_ptr<Grid>>();                  // Remove the grids from the cache.
            nGrids = Cache::size<Grid::ID,shared_ptr<Grid>>();
            BOOST_TEST( nGrids == 0 );

            // Same as above, but origin fixed at (0,0)
            for( auto& i: variousSizes ) testGrid( i, 0, false );
            
            nGrids = Cache::size<Grid::ID,shared_ptr<Grid>>();          // All the grids should be cached
            BOOST_TEST( nGrids == nVS );
            
            Cache::clear<Grid::ID,shared_ptr<Grid>>();                  // Remove the grids from the cache.
            nGrids = Cache::size<Grid::ID,shared_ptr<Grid>>();
            BOOST_TEST( nGrids == 0 );


            // Same as above, but with x/y-sizes different
            vector<PointI> variousSizes2D = { PointI(50, 48), PointI(49,48), PointI(500, 497), PointI(497, 333), PointI(9,9), PointI(13,5) };
            nVS = variousSizes2D.size();
            
            for( auto& i: variousSizes2D ) testGrid( i.x, i.y );
            
            nGrids = Cache::size<Grid::ID,shared_ptr<Grid>>();          // All the grids should be cached
            BOOST_TEST( nGrids == nVS );
            
            Cache::clear<Grid::ID,shared_ptr<Grid>>();                  // Remove the grids from the cache.
            nGrids = Cache::size<Grid::ID,shared_ptr<Grid>>();
            BOOST_TEST( nGrids == 0 );

            // Check if the generator throws when called with 0-size.
            shared_ptr<Grid> tmpGrid;
            BOOST_CHECK_THROW( tmpGrid=Grid::get( 0, 1.0, 2.2 ), logic_error );         // Should throw a logic_error if called with size=0
            BOOST_CHECK_THROW( tmpGrid=Grid::get( 13, 0, 1.0, 2.2 ), logic_error );     // ...same but for 2D-sizes.
            
            nGrids = Cache::size<Grid::ID,shared_ptr<Grid>>();                          // No new Grids should have been created/added to Cache
            BOOST_TEST( nGrids == 0 );
            
        }


        using namespace boost::unit_test;
        void add_grid_tests( test_suite* ts ) {

            srand(time(NULL));
            
            // Grid, Pupil & Modes
            ts->add( BOOST_TEST_CASE_NAME( &grid_test, "Grid" ) );

        }

    }

}
