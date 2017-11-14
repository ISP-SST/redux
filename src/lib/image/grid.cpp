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

    distance = sharedArray<float> (id.size.y, id.size.x);
    angle = sharedArray<float> (id.size.y, id.size.x);
    float** distPtr = distance.get();
    float** anglePtr = angle.get();
    for (unsigned int y = 0; y < id.size.y; ++y) {
        double yDist = y - id.origin.y;
        double y2 = yDist*yDist;
        for (unsigned int x = 0; x < id.size.x; ++x) {
            double xDist = x - id.origin.x;
            if (yDist || xDist) {
                double x2 = xDist*xDist;
                distPtr[y][x] = sqrt (y2 + x2);
                anglePtr[y][x] = atan2 (yDist, xDist);      // note: slow index is Y, fast is X
            } else distPtr[y][x] = anglePtr[y][x] = 0;      // this pixel is at the origin -> set the angle to 0.
        }
    }

}


shared_ptr<Grid> Grid::get(const Grid::ID& id) {
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

/*
double Zernike::covariance( int32_t i, int32_t j ) {
    
    if(i > j) std::swap(i, j);      // it's symmetric, so only store 1
    
    double& cov = Cache::get<PairID,double>( PairID(i,j), std::numeric_limits<double>::infinity() );
    unique_lock<mutex> lock(get().mtx);
    if( ! isfinite(cov) ) { // not calulated yet
        cov = calcZernikeCovariance(i, j);
    }
    return cov;

}
*/