#ifndef REDUX_IMAGE_UTIL_HPP
#define REDUX_IMAGE_UTIL_HPP

#include "redux/util/array.hpp"

#include <algorithm>
#include <functional>
#include <iostream>

namespace redux {

    namespace image {

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


        /*! @fn size_t nextLocation( int i, size_t M, size_t sizeX )
        *   @brief Find the new location of element "i" after transposing an MxN matrix.
        *   @param i Original location (as linear offset from the first element in the matrix)
        *   @param sizeY Size of the first dimension ( "row", slow )
        *   @param sizeX Size of the second dimension ( "column", fast )
        *   @returns The new location
        */
        inline size_t nextLocation( int i, size_t sizeY, size_t sizeX ) {
            return ( i % sizeX ) * sizeY + i / sizeX;
        }

        /*! @fn void transpose( T* data, size_t sizeY, size_t sizeX )
         *  @brief Transpose a 2D array (MxN matrix) of arbitrary type and size.
         *  @details Performs an in-place transpose of the matrix. The Method works by traversing through the permutaion-cycles, swapping the elements pairwise.
         *  @param data Input 2D Array (matrix)
         *  @param sizeY Size of the first dimension ("row",slow)
         *  @param sizeX Size of the second dimension ("column",fast)
         */
        template <class T>
        void transpose( T* data, size_t sizeY, size_t sizeX ) {

            if( sizeY && sizeX && data ) {
                size_t stillToMove = sizeY * sizeX;
                for( size_t i = 0; stillToMove; ++i ) {
                    size_t j, k;
                    for( j = nextLocation( i, sizeY, sizeX ); j > i; j = nextLocation( j, sizeY, sizeX ) ) ; //cycle.push_back(j+1);
                    if( j < i ) continue; // If true, we already traversed this cycle earlier.
                    // Note: j=nextLocation(i,sizeY,sizeX) => i=nextLocation(j,sizeX,sizeY) (cf. transposing MxN vs. NxM)
                    // We need to traverse the cycle backwards to get the elements in the right place by simple swapping, so interchange sizeY & sizeX
                    for( k = i, j = nextLocation( i, sizeX, sizeY ); j != i; k = j, j = nextLocation( j, sizeX, sizeY ) ) {
                        std::swap( data[k], data[j] );
                        --stillToMove;
                    }
                    --stillToMove;
                }
            }
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
            double tmp, chisq = 0;
            typename redux::util::Array<U>::const_iterator bit = b.begin();
            for( auto ait : a ) {
                tmp = ait - *bit++;
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
            size_t count;
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
                    //std::cout << "x = " << x << "    y = " << y << "    arr = " << array[y][x] << "    pred = " << (int)predicate(array[y][x]) << std::endl;
                    if( predicate( array[y][x] ) ) tmp.insert( std::pair<size_t, T>( offset, filler( y, x ) ) );
                    ++offset;
                }
            }
            size_t cnt = 0;
            for( auto it : tmp ) {
                //std::cout << "first = " << it.first << "    second = " << it.second << "    orig = " << ptr[it.first] << std::endl;
                ptr[it.first] = it.second;
                ++cnt;
            }
            std::cout << "cnt = " << cnt << std::endl;
            std::cout << "fillPixels   sy = " << sy << "   sx = " << sx << std::endl;
        }


        template <typename T, typename FillFunction, typename Predicate>
        void fillPixels( redux::util::Array<T>& array, FillFunction* filler, Predicate predicate = std::bind2nd( std::less_equal<T>(), 0 ) ) {
            for( auto it = array.begin(); it != array.end(); ++it ) {
                if( predicate( *it ) ) *it = filler( it );
            }
        }


    }   // image

}   // redux


#endif  // REDUX_IMAGE_UTIL_HPP
