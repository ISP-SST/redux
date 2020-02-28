#include "redux/image/grid.hpp"

#include "redux/util/arrayutil.hpp"
#include "redux/util/cache.hpp"

#include <mutex>

using namespace redux::image;
using namespace redux::util;

using namespace std;

namespace {
    
    mutex mtx;
    
}

Grid::ID::ID(uint32_t nPoints) : size(nPoints, nPoints), origin (nPoints / 2, nPoints / 2) {
    if (! (nPoints % 2)) {  // for nPoints even, place origin between mid-points
        origin.x += 0.5;
        origin.y += 0.5;
    }
}


Grid::ID::ID(uint32_t nPoints, float originY, float originX) : size (nPoints, nPoints), origin (originY, originX) {
    
}


Grid::ID::ID(uint32_t nPointsY, uint32_t nPointsX, float originY, float originX) : size (nPointsY, nPointsX), origin (originY, originX) {

}


bool Grid::ID::operator<( const Grid::ID& rhs ) const {
    if( size == rhs.size ) {
        return (origin < rhs.origin);
    }
    return ( size < rhs.size );
}


void Grid::init (void) {

    distance = rdx_get_shared<double>( id.size.y*id.size.x );
    angle = rdx_get_shared<double>( id.size.y*id.size.x );
    auto d2D = dist2D();
    auto a2D = angle2D();
    double** distPtr = d2D.get();
    double** anglePtr = a2D.get();
    PointF shiftedOrigin = id.origin - 0.5;                // (0,0) means cetered in the lower left corner, (0.5,0.5) means on the first pixel
    for (unsigned int y = 0; y < id.size.y; ++y) {
        long double yDist = y - shiftedOrigin.y;
        long double y2 = yDist*yDist;
        for (unsigned int x = 0; x < id.size.x; ++x) {
            long double xDist = x - shiftedOrigin.x;
            if (yDist || xDist) {
                long double x2 = xDist*xDist;
                distPtr[y][x] = sqrtl( y2 + x2 );
                anglePtr[y][x] = atan2l( yDist, xDist );      // note: slow index is Y, fast is X
            } else distPtr[y][x] = anglePtr[y][x] = 0;      // this pixel is at the origin -> set the angle to 0.
        }
    }

}


shared_ptr<Grid> Grid::get( const Grid::ID& id ) {
    
    if( (id.size.x == 0) || (id.size.y == 0) ) {
        throw logic_error("Grid::ID::size can not be 0!");
    }
    
    shared_ptr<Grid>& grid = Cache::get<ID,shared_ptr<Grid>>( id, nullptr );
    unique_lock<mutex> lock(mtx);
    if( !grid ) grid.reset( new Grid() );
    if( grid->id.size == 0 ) {
        grid->id = id;
        grid->init();
    }
    return grid;
}


void Grid::clear(void) {
    
    Cache::clear<ID,shared_ptr<Grid>>();
    
}
