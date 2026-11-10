#include "redux/image/grid.hpp"

#include "catch.hpp"

#include <cstdlib>
#include <ctime>

using namespace redux::image;
using namespace redux::util;
using namespace std;


namespace testsuite {

    namespace image {

        namespace {

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
                        CHECK( std::fabs(distPtr[y][x] - sqrt(yDist*yDist+xDist*xDist)) < (0.001) );
                        if( yDist || xDist ) {
                            CHECK( std::fabs(anglePtr[y][x] - atan2(yDist, xDist)) < (0.001) );
                        } else {
                            CHECK( std::fabs(anglePtr[y][x]) < (0.001) );
                        }
                    }
                }

            }

            // after generating grids via testGrid(), checks the cache grew to the
            // expected count, then clears it and confirms it's actually empty afterwards.
            void checkAndClearGridCache( size_t expected ) {
                size_t nGrids = Cache::size<Grid::ID,shared_ptr<Grid>>();
                CHECK( nGrids == expected );
                Cache::clear<Grid::ID,shared_ptr<Grid>>();
                nGrids = Cache::size<Grid::ID,shared_ptr<Grid>>();
                CHECK( nGrids == 0 );
            }

        }

        TEST_CASE( "Grid", "[image]" ) {
            srand(time(NULL));

            // Test generated grids for some different sizes (even/odd/large/small/...)
            vector<int> variousSizes = { 50, 49, 500, 497, 9, 8 };

            for( auto& i: variousSizes ) testGrid( i );
            checkAndClearGridCache( variousSizes.size() );     // all grids should have been cached, then removed

            // Same as above, but origin fixed at (0,0)
            for( auto& i: variousSizes ) testGrid( i, 0, false );
            checkAndClearGridCache( variousSizes.size() );

            // Same as above, but with x/y-sizes different
            vector<PointI> variousSizes2D = { PointI(50, 48), PointI(49,48), PointI(500, 497), PointI(497, 333), PointI(9,9), PointI(13,5) };

            for( auto& i: variousSizes2D ) testGrid( i.x, i.y );
            checkAndClearGridCache( variousSizes2D.size() );

            // Check if the generator throws when called with 0-size.
            shared_ptr<Grid> tmpGrid;
            CHECK_THROWS_AS( tmpGrid=Grid::get( 0, 1.0, 2.2 ), logic_error );         // Should throw a logic_error if called with size=0
            CHECK_THROWS_AS( tmpGrid=Grid::get( 13, 0, 1.0, 2.2 ), logic_error );     // ...same but for 2D-sizes.

            checkAndClearGridCache( 0 );     // no new Grids should have been created/added to Cache
        }

    }

}
