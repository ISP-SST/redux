#ifndef REDUX_IMAGE_UTIL_HPP
#define REDUX_IMAGE_UTIL_HPP

#include "redux/image/image.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/types.hpp"

#include <atomic>
#include <functional>
#include <iostream>
#include <map>
#include <thread>

namespace redux {

    namespace image {
        
        double makePupil_thi( double** pupil, uint32_t nPoints, double radius);
        double makePupil_mvn( double** pupil, int nPoints, double radius );
        template <typename T>
        double makePupil_old( util::Array<T>& pupil, uint32_t nPoints, double radius) {
            pupil.resize(nPoints,nPoints);
            T **ptr = redux::util::makePointers(pupil.get(),nPoints,nPoints);
            double area = makePupil_mvn(ptr,nPoints,radius); // FIXME: temporarily using MvN pupilmaker for easier debugging.
            redux::util::delPointers(ptr);
            return area;
        }

        void makeZernike_thi( double** mode, int j, uint32_t nPoints, double radius, double angle=0 );
        void makeZernike_mvn( double** mode, int j, uint32_t nPoints, double radius, double angle=0 );
        template <typename T>
        void makeZernike( util::Array<T>& mode, int j, uint32_t nPoints, double radius, double angle=0) {
            mode.resize(nPoints,nPoints);
            T **ptr = redux::util::makePointers(mode.get(),nPoints,nPoints);
            makeZernike_mvn(ptr,nPoints,radius,angle); // FIXME: temporarily using MvN Z-maker for easier debugging.
            redux::util::delPointers(ptr);
        }

        
        template <typename T>
        redux::util::Array<T> fitPlane( const redux::util::Array<T>& in, bool subtract_mean=false );

        template <typename T>
        double total( const redux::util::Array<T>& in ) {
            double sum( 0 );
            for( auto &value : in ) {
                sum += value;
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
            for( auto &value : in ) {
                sum += value * *wit;
                ++wit;
            }
            return sum;
        }


        template <typename T>
        double mean( const redux::util::Array<T>& in ) {
            size_t nElements = in.nElements();
            if( nElements == 0 ) return 0;
            return total( in ) / static_cast<double>( nElements );
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
        void connectedRegion(T** data, uint8_t** mask, size_t sizeY, size_t sizeX, unsigned int y, unsigned int x, T threshold=T(0));
        template <typename T>
        void connectedRegion(T** data, size_t sizeY, size_t sizeX, unsigned int y, unsigned int x, T threshold=T(0)) {
            uint8_t** mask = redux::util::newArray<uint8_t>(sizeY,sizeX);
            connectedRegion(data, mask, sizeY, sizeX, y, x, threshold);
            std::transform(*data, *data+sizeY*sizeX, *mask, *data, redux::math::multiply<T,uint8_t>());
            redux::util::delArray(mask);
        }
        template <typename T>
        void connectedRegion(T* data, size_t sizeY, size_t sizeX, unsigned int y, unsigned int x, T threshold=T(0)) {
            T** ptr = redux::util::makePointers(data, sizeY, sizeX);
            connectedRegion(ptr, sizeY, sizeX, y, x, threshold);
            redux::util::delPointers(ptr);
        }
       
        template <typename T>
        void smooth(T** data, size_t sizeY, size_t sizeX, size_t nY, size_t nX);
        template <typename T>
        void smooth(T** data, size_t sizeY, size_t sizeX, size_t n) { smooth(data, sizeY, sizeX, n, n); }
        template <typename T>
        void smooth(T* data, size_t sizeY, size_t sizeX, size_t nY, size_t nX) {
            T** ptr = redux::util::makePointers(data, sizeY, sizeX);
            smooth(ptr, sizeY, sizeX, nY, nX);
            redux::util::delPointers(ptr);
        }
        template <typename T>
        void smooth(T* data, size_t sizeY, size_t sizeX, size_t n) { smooth(data, sizeY, sizeX, n, n); }
        
        template <typename T>
        void ScharmerFilter (T** data, double** q2_inv, size_t sizeY, size_t sizeX, double noise_power, double frequency_cutoff);
        template <typename T>
        void ScharmerFilter(T* data, double* q2_inv, size_t sizeY, size_t sizeX, double noise_power, double frequency_cutoff) {
            T** ptr = redux::util::makePointers(data, sizeY, sizeX);
            double** qPtr = redux::util::makePointers(q2_inv, sizeY, sizeX);
            ScharmerFilter(ptr, qPtr, sizeY, sizeX, noise_power, frequency_cutoff);
            redux::util::delPointers(ptr);
            redux::util::delPointers(qPtr);
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
            for( auto &avalue : a ) {
                double tmp = avalue - *bit++;
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
            for( auto &avalue : a ) {
                if( *++wit ) {
                    tmp = ( avalue - *bit++ ) * ( *wit );
                    chisq += tmp * tmp;
                    ++count;
                }
            }
            if( count == 0 ) return 0.0;
            return chisq / static_cast<double>( count );
        }

        template <typename T>
        void clipImage( redux::util::Array<T>& img, const std::vector<int16_t> clip, bool symmetricClip=false )  {
            
            std::vector<int16_t> thisClip = clip;
            if( thisClip.size() == 2 ) {        // apply same values to both dimensions.
                thisClip.insert( thisClip.end(), clip.begin(), clip.end() );
            }
            
            if( thisClip.size() == 4 ) {
                bool flipX = false, flipY = false;
                // we have the y (row/slow) dimension first, momfbd cfg-files (and thus alignClip) has x first.
                if ( thisClip[0] > thisClip[1] ) {
                    std::swap( thisClip[0], thisClip[1] );
                    flipX = true;
                }
                if ( thisClip[2] > thisClip[3] ) {
                    std::swap( thisClip[2], thisClip[3] );
                    flipY = true;
                }
                for( auto & index : thisClip )
                    --index;       // NOTE: momfbd cfg files uses 1-based indexes, internally we start with 0.
                size_t sy = thisClip[3] - thisClip[2] + 1;
                size_t sx = thisClip[1] - thisClip[0] + 1;
                if( symmetricClip ) {
                    const std::vector<size_t>& dims = img.dimensions();
                    int skewY = (dims[0] - sy) / 2  - thisClip[2];
                    int skewX = (dims[1] - sx) / 2  - thisClip[0];
                    thisClip[0] += skewX;
                    thisClip[1] += skewX;
                    thisClip[2] += skewY;
                    thisClip[3] += skewY;
                }
                img.setLimits( thisClip[2], thisClip[3], thisClip[0], thisClip[1] );
                img.trim();

                if( flipX || flipY ) {
                    std::shared_ptr<T*> arrayPtr = img.reshape(sy, sx);
                    T** imgPtr = arrayPtr.get();
                    if (flipX) redux::util::reverseX(imgPtr, sy, sx);
                    if (flipY) redux::util::reverseY(imgPtr, sy, sx);
                }
            }
        }
        
        template <typename T>
        void clipImage( redux::image::Image<T>& arr, const std::vector<int16_t> clip, bool symmetricClip=false )  {
            clipImage( reinterpret_cast<redux::util::Array<T>&>(arr), clip, symmetricClip );
        }

        template <typename T, typename Predicate>
        void fillPixels( redux::util::Array<T>& array, T fillValue, Predicate predicate = std::bind2nd( std::less_equal<T>(), 0 ) ) {
            for( auto & value : array ) {
                if( predicate( value ) ) value = fillValue;
            }
        }


        template <typename T, typename Predicate, typename MaskType=uint8_t>
        void fillPixels( T** array, size_t sy, size_t sx, std::function<double( size_t, size_t )> filler, Predicate predicate = std::bind2nd( std::less_equal<T>(), 0 ), MaskType** mask=nullptr  ) {
            std::map<size_t, T> tmp;
            T* ptr = *array;
            size_t offset = 0;
            for( size_t y = 0; y < sy; ++y ) {
                for( size_t x = 0; x < sx; ++x ) {
                    if( (!mask || !mask[y][x]) && predicate( array[y][x] ) ) tmp.insert( std::pair<size_t, T>( offset, filler( y, x ) ) );
                    ++offset;
                }
            }
            size_t cnt = 0;
            for( auto &it : tmp ) {
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

        template <typename T, typename U=uint8_t>
        void fillPixels( T** image, size_t sy, size_t sx, U** mask=nullptr ) {
            
            size_t offset(0);
            const size_t nPixels = sy*sx;

            std::map<size_t, T> values;
            size_t o;
            while( (o = offset++) < nPixels ) {
                if( !mask || mask[0][o%nPixels] ) {
                    size_t y = o/sy;
                    size_t x = o%sy;
                    values.insert( std::pair<size_t, T>( o, inverseDistanceWeight( image, sy, sx, y, x ) ) );
                    //values.insert( std::pair<size_t, T>( o, horizontalInterpolation( image, sy, sx, y, x ) ) );
                }
            }
            for( auto &it: values ) {
                image[0][it.first] = it.second;
            }

        }

        
        template <typename T, typename U=uint8_t>
        void fillPixels( T*** image, size_t nImages, size_t sy, size_t sx, U** mask=nullptr, unsigned int nThreads=std::thread::hardware_concurrency() ) {
            
            std::atomic<size_t> offset(0);
            const size_t nPixels = sy*sx;
            
            auto nextbad = [&](void) {
                size_t o;
                while( (o = offset.fetch_add(1)) < nPixels ) {
                    if( !mask || mask[0][o%nPixels]) return o;
                }
                return nPixels;
            };

            std::vector<std::thread> threads;
            for( unsigned int t=0; t<nThreads; ++t ) {
                threads.push_back( std::thread(
                    [&](){
                        size_t myOffset;
                        std::map<size_t, T> values;
                        while( (myOffset=nextbad()) < nPixels ) {
                            size_t y = myOffset/sy;
                            size_t x = myOffset%sy;
                            values.insert( std::pair<size_t, T>( myOffset, inverseDistanceWeight( image, sy, sx, y, x ) ) );
                        }

                        for( auto &it: values ) {
                            image[0][it.first] = it.second;
                        }
                    }));
            }
            for( auto& th : threads ) th.join();

        }

        template <typename T>
        void apodizeInPlace( redux::util::Array<T>& array, size_t blendRegion );
        template <typename T>
        redux::util::Array<T> apodize( const redux::util::Array<T>& array, size_t blendRegion );
        
        template <typename T>
        void normalizeIfMultiFrames( redux::image::Image<T>& img );
        
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
