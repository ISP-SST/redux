#ifndef REDUX_IMAGE_GRID_HPP
#define REDUX_IMAGE_GRID_HPP

#include "redux/util/point.hpp"

#include <memory>

namespace redux {

    namespace image {
        

        /*! Container for an equidistant grid. Distances to origin (in pixels) and angles (in radians) are stored in
         *  shared arrays. By default the origin is centered in the grid (between points for a grid with even points).
         */
        struct Grid {
            struct ID {
                explicit ID( uint32_t nPoints );
                ID( uint32_t nPoints, float originY, float originX );
                ID( uint32_t nPoints, redux::util::PointF origin ) : ID( nPoints, origin.y, origin.x ) {};
                ID( uint32_t nPointsY, uint32_t nPointsX, float originY, float originX );
                redux::util::Point size;
                redux::util::PointF origin;
                bool operator<( const ID& rhs ) const;
            } id;
            std::shared_ptr<float*> distance;
            std::shared_ptr<float*> angle;
            Grid(void) : id(0) { };
            void init(void);
            bool operator<( const Grid& rhs ) const { return ( id < rhs.id ); }
            static std::shared_ptr<Grid> get(const ID&);
            static std::shared_ptr<Grid> get(uint32_t n) { return get( ID(n) ); };
            static std::shared_ptr<Grid> get(uint32_t n, float y, float x) { return get( ID(n,y,x) ); };
            static std::shared_ptr<Grid> get(uint32_t ny, uint32_t nx, float y, float x) { return get( ID(ny,nx,y,x) ); };
            static void clear(void);
        };

    }   // image

}   // redux


#endif  // REDUX_IMAGE_GRID_HPP
