#ifndef REDUX_IMAGE_UTIL_HPP
#define REDUX_IMAGE_UTIL_HPP

#include "redux/image/image.hpp"
#include "redux/util/array.hpp"
#include "redux/types.hpp"

#include <algorithm>
#include <functional>
#include <iostream>

namespace redux {

    namespace image {
        

        /*! Container for an equidistant grid. Distances to origin (in pixels) and angles (in radians) are stored in
         *  shared arrays. By default the origin is centered in the grid (between points for a grid with even points).
         */
        struct Grid {
            Point size;
            PointF origin;
            std::shared_ptr<float*> distance;
            std::shared_ptr<float*> angle;
            Grid( uint32_t nPoints );
            Grid( uint32_t nPoints, float originY, float originX );
            Grid( uint32_t nPointsY, uint32_t nPointsX, float originY, float originX );
            void init(void);
            bool operator<( const Grid& rhs ) const { if(size == rhs.size ) return (origin < rhs.origin); return ( size < rhs.size ); }
        };
        
        float makePupil( util::Array<float>& aperture, uint32_t nPoints, float radius, bool normalize=false);
        float makePupil_mvn( float** data, uint32_t nPoints, float radius );

        template <typename T>
        void reverseX( T** data, size_t sizeY, size_t sizeX ) {
            for( size_t y = 0; y < sizeY; ++y ) {
                std::reverse( data[y], data[y] + sizeX );
            }
        }
        

        template <typename T>
        void reverseY( T** data, size_t sizeY, size_t sizeX ) {
            size_t halfY = sizeY / 2;
            size_t nBytes = sizeX * sizeof( T );
            T* tmp = new T[sizeX];
            for( size_t y = 0; y < halfY; ++y ) {
                memcpy( tmp, data[y], nBytes );
                memcpy( data[y], data[sizeY - y - 1], nBytes );
                memcpy( data[sizeY - y - 1], tmp, nBytes );
            }
            delete[] tmp;
        }


        template <typename T>
        double total( const redux::util::Array<T>& in ) {
            double sum( 0 );
            for( auto it : in ) {
                sum += it;
            }
            return sum;
        }


        template <typename T, typename U>
        double total( const redux::util::Array<T>& in, const redux::util::Array<U>& weight ) {
            if( !in.sameSizes( weight ) ) {
                throw std::logic_error( "Array dimensions does not match." );
            }
            double sum( 0 );
            typename redux::util::Array<U>::const_iterator wit = weight.begin();
            for( auto it : in ) {
                sum += it * *wit;
                ++wit;
            }
            return sum;
        }


        template <typename T>
        double mean( const redux::util::Array<T>& in ) {
            size_t nElements = in.nElements();
            if( nElements == 0 ) return 0;
            double sum( 0 );
            for( auto it : in ) {
                sum += it;
            }
            return sum / static_cast<double>( nElements );
        }


        template <typename T>
        T median( std::vector<T> in ) {
            std::nth_element( in.begin(), in.begin() + in.size() / 2, in.end() );
            return *( in.begin() + in.size() / 2 );
        }


        template <typename T>
        T median( const redux::util::Array<T>& in ) {
            redux::util::Array<T> tmp = in.copy();  // nth_element is destructive, so make a deep copy.
            std::nth_element( tmp.begin(), tmp.mid(), tmp.end() );
            return *tmp.mid();
        }


        template <typename T>
        double inverseDistanceWeight( T**, size_t sizeY, size_t sizeX, size_t posY, size_t posX );

        double inv_dist_wght( float **a, size_t sizeY, size_t sizeX, size_t posY, size_t posX );

        template <typename T>
        double horizontalInterpolation( T**, size_t sizeY, size_t sizeX, size_t posY, size_t posX );


        template <typename T, typename U>
        double chisq( const redux::util::Array<T>& a, const redux::util::Array<U>& b ) {
            if( !a.sameSizes( b ) ) {
                throw std::logic_error( "Array dimensions does not match." );
            }
            size_t nElements = a.nElements();
            if( nElements == 0 ) return 0.0;
            double chisq = 0;
            typename redux::util::Array<U>::const_iterator bit = b.begin();
            for( auto ait : a ) {
                double tmp = ait - *bit++;
                chisq += tmp * tmp;
            }
            return chisq / static_cast<double>( nElements );
        }


        template <typename T, typename U, typename V>
        double chisq( const redux::util::Array<T>& a, const redux::util::Array<U>& b, const redux::util::Array<V>& weight ) {
            if( !a.sameSizes( b ) || !a.sameSizes( weight ) ) {
                throw std::logic_error( "Array dimensions does not match." );
            }
            double tmp, chisq = 0;
            typename redux::util::Array<U>::const_iterator bit = b.begin();
            typename redux::util::Array<V>::const_iterator wit = weight.begin()--;
            size_t count(0);
            for( auto ait : a ) {
                if( *++wit ) {
                    tmp = ( ait - *bit++ ) * ( *wit );
                    chisq += tmp * tmp;
                    ++count;
                }
            }
            if( count == 0 ) return 0.0;
            return chisq / static_cast<double>( count );
        }


        template <typename T, typename Predicate>
        void fillPixels( redux::util::Array<T>& array, T fillValue, Predicate predicate = std::bind2nd( std::less_equal<T>(), 0 ) ) {
            for( auto & it : array ) {
                if( predicate( it ) ) it = fillValue;
            }
        }


        template <typename T, typename Predicate>
        void fillPixels( T** array, size_t sy, size_t sx, std::function<double( size_t, size_t )> filler, Predicate predicate = std::bind2nd( std::less_equal<T>(), 0 ) ) {
            std::map<size_t, T> tmp;
            T* ptr = *array;
            size_t offset = 0;
            for( size_t y = 0; y < sy; ++y ) {
                for( size_t x = 0; x < sx; ++x ) {
                    if( predicate( array[y][x] ) ) tmp.insert( std::pair<size_t, T>( offset, filler( y, x ) ) );
                    ++offset;
                }
            }
            size_t cnt = 0;
            for( auto it : tmp ) {
                ptr[it.first] = it.second;
                ++cnt;
            }
        }


        template <typename T, typename FillFunction, typename Predicate>
        void fillPixels( redux::util::Array<T>& array, FillFunction* filler, Predicate predicate = std::bind2nd( std::less_equal<T>(), 0 ) ) {
            for( auto it = array.begin(); it != array.end(); ++it ) {
                if( predicate( *it ) ) *it = filler( it );
            }
        }

        template <typename T>
        void apodize( redux::util::Array<T>& array, size_t blendRegion );
        
        template <typename T>
        void checkIfMultiFrames( redux::image::Image<T>& img );
        
        /*! Fit and subtract the backscatter patter for semi-transparent CCDs 
         *  @param data Image to be cleaned
         *  @param ccdGain Image containing the backscatter pattern
         *  @param psf PSF for the system
         *  @param maxIterations Exit condition
         *  @param minImprovement Exit if the improvement in one step is less than this factor (1 means no improvement)
         *  @param epsilon Exit if the metric reaches this value.
         */
        template <typename T, typename U>
        void descatter(redux::util::Array<T>& data, const redux::util::Array<U>& ccdgain, const redux::util::Array<U>& psf,
                       int maxIterations=50, double minImprovement=1, double epsilon=1E-10 );

    }   // image

}   // redux


#endif  // REDUX_IMAGE_UTIL_HPP
