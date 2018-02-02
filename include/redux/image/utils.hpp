#ifndef REDUX_IMAGE_UTIL_HPP
#define REDUX_IMAGE_UTIL_HPP

#include "redux/image/image.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/progresswatch.hpp"
#include "redux/types.hpp"

#ifdef RDX_WITH_OPENCV
#   include "redux/util/opencv.hpp"
#endif

#include <atomic>
#include <functional>
#include <map>
#include <queue>
#include <thread>

#include <boost/asio.hpp>
#include <boost/thread/thread.hpp>

namespace redux {

    namespace image {
        
        double makePupil( double** pupil, uint32_t nPoints, double outerRadius, double innerRadius=0.0 );

        void makeZernike( double** mode, int j, uint32_t nPoints, double radius, double angle=0 );

        template <typename T>
        void makeZernike( util::Array<T>& mode, int j, uint32_t nPoints, double radius, double angle=0) {
            mode.resize(nPoints,nPoints);
            T **ptr = redux::util::makePointers(mode.get(),nPoints,nPoints);
            makeZernike(ptr,nPoints,radius,angle); // FIXME: temporarily using MvN Z-maker for easier debugging.
            redux::util::delPointers(ptr);
        }

        
        template <typename T>
        redux::util::Array<T> fitPlane( const redux::util::Array<T>& in, bool subtract_mean=false, double* coeffs=nullptr );

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
            size_t nDims = img.nDimensions();
            if( nDims < 2 ) return;
            std::vector<size_t> finalClip;
            size_t nImgs = 1;
            if( nDims > 2 ) {
                for( size_t i=0; i<nDims-2; ++i ) {
                    size_t sz = img.dimSize(i);
                    nImgs *= sz;
                    finalClip.push_back(0);
                    finalClip.push_back(sz-1);
                }
            }
            std::vector<int16_t> tmpClip = clip;
            if( tmpClip.size() == 2 ) {        // apply same values to both dimensions.
                tmpClip.insert( tmpClip.end(), clip.begin(), clip.end() );
            }
            
            if( tmpClip.size() == 4 ) {
                bool flipX = false, flipY = false;
                // we have the y (row/slow) dimension first, momfbd cfg-files (and thus alignClip) has x first.
                if ( tmpClip[0] > tmpClip[1] ) {
                    std::swap( tmpClip[0], tmpClip[1] );
                    flipX = true;
                }
                if ( tmpClip[2] > tmpClip[3] ) {
                    std::swap( tmpClip[2], tmpClip[3] );
                    flipY = true;
                }
                for( auto & index : tmpClip )
                    --index;       // NOTE: momfbd cfg files uses 1-based indexes, internally we start with 0.
                size_t sy = tmpClip[3] - tmpClip[2] + 1;
                size_t sx = tmpClip[1] - tmpClip[0] + 1;
                if( symmetricClip ) {
                    const std::vector<size_t>& dims = img.dimensions();
                    int skewY = (dims[0] - sy) / 2  - tmpClip[2];
                    int skewX = (dims[1] - sx) / 2  - tmpClip[0];
                    tmpClip[0] += skewX;
                    tmpClip[1] += skewX;
                    tmpClip[2] += skewY;
                    tmpClip[3] += skewY;
                }
                
                finalClip.push_back(tmpClip[2]);
                finalClip.push_back(tmpClip[3]);
                finalClip.push_back(tmpClip[0]);
                finalClip.push_back(tmpClip[1]);
                img.setLimits( finalClip );
                img.trim();

                if( flipX || flipY ) {
                    std::shared_ptr<T*> arrayPtr = img.reshape(nImgs*sy, sx);
                    T** imgPtr = arrayPtr.get();
                    for( size_t i=0; i<nImgs; ++i) {
                        if (flipX) redux::util::reverseX(imgPtr, sy, sx);
                        if (flipY) redux::util::reverseY(imgPtr, sy, sx);
                        imgPtr += sy;
                    }
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
        void apodizeInPlace( T** data, size_t nRows, size_t nCols, size_t rowBlend, size_t colBlend, size_t rowMargin=0, size_t colMargin=0 );
        template <typename T>
        void apodizeInPlace( redux::util::Array<T>& array, size_t blendRegion, size_t margin=0 ) {
            size_t nRows = array.dimSize(0);
            size_t nCols = array.dimSize(1);
            std::shared_ptr<T*> tmp = array.reshape( nRows, nCols );
            apodizeInPlace( tmp.get(), nRows, nCols, blendRegion, blendRegion, margin, margin );
        }
        template <typename T>
        void apodizeInPlace( redux::util::Array<T>& array, size_t rowBlend, size_t colBlend, size_t rowMargin, size_t colMargin ) {
            size_t nRows = array.dimSize(0);
            size_t nCols = array.dimSize(1);
            std::shared_ptr<T*> tmp = array.reshape( nRows, nCols );
            apodizeInPlace( tmp.get(), nRows, nCols, rowBlend, colBlend, rowMargin, colMargin );
        }

        template <typename T>
        redux::util::Array<T> apodize( const redux::util::Array<T>& in, size_t blendRegion, size_t margin=0 ) {
            if( !blendRegion && !margin ) return in;       // nothing to do
            redux::util::Array<T> array;
            in.copy( array );
            apodizeInPlace( array, blendRegion, margin );
            return std::move(array);
        }
        template <typename T>
        redux::util::Array<T> apodize( const redux::util::Array<T>& in, size_t rowBlend, size_t colBlend, size_t rowMargin, size_t colMargin ) {
            if( !rowBlend && !colBlend && !rowMargin && !colMargin ) return in;       // nothing to do
            redux::util::Array<T> array;
            in.copy( array );
            apodizeInPlace( array, rowBlend, colBlend, rowMargin, colMargin );
            return std::move(array);
        }
        
        template <typename T, typename U>
        void insertPatch( T** img, size_t imgRows, size_t imgCols, const U** patch, size_t pRows, size_t pCols,
                                       size_t posY, size_t posX, double** weight=nullptr, double** globalSum=nullptr, bool transpose=false ) {
           
            if( transpose ) std::swap( pRows, pCols );
            if( posY + pRows >= imgRows ) pRows = imgRows-posY;
            if( posX + pCols >= imgCols ) pCols = imgCols-posX;
            
            double* gPtr = nullptr;
            if( transpose ) {
                for( size_t r=0; r<pRows; ++r ) {
                    T* imgRow = img[posY+r] + posX;
                    if( globalSum ) gPtr = globalSum[posY+r] + posX;
                    for( size_t c=0; c<pCols; ++c ) {
                        if( weight ) {
                            if( globalSum ) {
                                gPtr[c] += weight[c][r];
                            } else {
                                imgRow[c] *= 1.0 - weight[c][r];
                            }
                            imgRow[c] += weight[c][r]*static_cast<double>(patch[c][r]);
                        } else {
                            imgRow[c] = static_cast<double>(patch[c][r]);
                        }
                    }
                }
            } else {
                for( size_t r=0; r<pRows; ++r ) {
                    T* imgRow = img[posY+r] + posX;
                    if( globalSum ) gPtr = globalSum[posY+r] + posX;
                    if( weight ) {
                        for( size_t c=0; c<pCols; ++c ) {
                            if( globalSum ) {
                                gPtr[c] += weight[r][c];
                            } else {
                                imgRow[c] *= 1.0 - weight[r][c];
                            }
                            imgRow[c] += weight[r][c]*static_cast<double>(patch[r][c]);
                        }
                    } else {
                        for( size_t c=0; c<pCols; ++c ) {
                            imgRow[c] = static_cast<double>(patch[r][c]);
                        }
                    }
                }
            }

            
        }

        template <typename T, typename U>
        void mozaic( boost::asio::io_service& service, size_t nThreads, T** img, size_t imgRows, size_t imgCols, const U*** patches, size_t nPatches,
                     size_t pRows, size_t pCols, const int32_t* posY, const int32_t* posX, int32_t blend, int32_t margin, bool transpose=false ) {

            struct patch_info {
                patch_info() : pPtr(nullptr), pY(0), pX(0) {};
                patch_info(const U** pPtr_, int32_t pY_, int32_t pX_ ) : pPtr(pPtr_), pY(pY_), pX(pX_) {};
                const U** pPtr;
                int32_t pY;
                int32_t pX;
            };
            std::queue<patch_info> q;
            for( size_t i=0; i<nPatches; ++i ) q.push( patch_info( patches[i], posY[i], posX[i]) );

            size_t nPixels = imgRows*imgCols;
            double** sum = redux::util::newArray<double>( imgRows, imgCols );
            memset( *sum, 0, nPixels*sizeof(double) );
            double** tmpImg = redux::util::newArray<double>( imgRows, imgCols );
            memset( *tmpImg, 0, nPixels*sizeof(double) );
            
            double** weights = redux::util::newArray<double>( pRows, pCols );
            std::fill( *weights, *weights+pRows*pCols, 1.0 );
            apodizeInPlace( weights, pRows, pCols, blend, blend, margin, margin );
            nThreads = std::min<size_t>( nPatches, nThreads );
            std::mutex mtx;
            redux::util::ProgressWatch pw;
            pw.set( nPatches+nThreads );
            for( size_t i=0; i<nThreads; ++i ) {
                service.post([&](){
                    double** tsum = redux::util::newArray<double>( imgRows, imgCols );
                    memset( *tsum, 0, nPixels*sizeof(double) );
                    double** timg = redux::util::newArray<double>( imgRows, imgCols );
                    memset( *timg, 0, nPixels*sizeof(double) );
                    std::unique_lock<std::mutex> lock(mtx);
                    while( !q.empty() ) {
                        patch_info pi = q.front();
                        q.pop();
                        lock.unlock();
                        insertPatch( timg, imgRows, imgCols, pi.pPtr, pRows, pCols, pi.pY, pi.pX, weights, tsum, transpose );
                        lock.lock();
                        ++pw;
                    }
                    std::transform( *tmpImg, *tmpImg+nPixels, *timg, *tmpImg, std::plus<double>() );
                    std::transform( *sum, *sum+nPixels, *tsum, *sum, std::plus<double>() );
                    redux::util::delArray( timg );
                    redux::util::delArray( tsum );
                    ++pw;
                    lock.unlock();
                });
            }
            
            pw.wait();
            std::transform( *tmpImg, *tmpImg+nPixels, *sum, *img,
                            []( const double& a, const double& b ){
                                if( b > 1E-6 ) return a/b;
                                else return a;
                            });
            
            redux::util::delArray( weights );
            redux::util::delArray( tmpImg );
            redux::util::delArray( sum );
            
         }

        template <typename T, typename U>
        void mozaic( T** img, size_t imgRows, size_t imgCols, const U*** patches, size_t nPatches, size_t pRows, size_t pCols,
                     const int32_t* posY, const int32_t* posX, int32_t blend, int32_t margin, bool transpose=false ) {

            size_t nThreads = std::min<size_t>( nPatches, std::thread::hardware_concurrency() );
            boost::asio::io_service service;
            boost::thread_group pool;
            {
                boost::asio::io_service::work workLoop(service);
                for( size_t t = 0; t < nThreads; ++t ) {
                    pool.create_thread( boost::bind(&boost::asio::io_service::run, &service) );
                }
                mozaic( service, nThreads, img, imgRows, imgCols, patches, nPatches, pRows, pCols, posY, posX, blend, margin, transpose );
            }
            pool.join_all();
            
        }

        template <typename T>
        void img_trim( T**& img, size_t& imgRows, size_t& imgCols, float threshold=1E-6 ) {
            std::vector<T> rowSums( imgCols, 0 );
            std::vector<T> colSums( imgRows, 0 );
            for ( size_t c=0; c < imgCols; ++c ) {
                for ( size_t r=0; r < imgRows; ++r ) {
                    T val = img[r][c];
                    rowSums[c] += val;
                    colSums[r] += val;
                }
            }

            size_t firstCol(0), lastCol(imgCols-1), firstRow(0), lastRow(imgRows-1);
            while ( firstCol < lastCol && rowSums[firstCol] < threshold ) ++firstCol;
            while ( lastCol && rowSums[lastCol] < threshold  ) --lastCol;
            while ( firstRow < lastRow && colSums[firstRow] < threshold  ) ++firstRow;
            while ( lastRow && colSums[lastRow] < threshold  ) --lastRow;

            size_t imgRows2 = lastRow-firstRow+1;
            size_t imgCols2 = lastCol-firstCol+1;

            if( imgRows2 != imgRows || imgCols2 != imgCols ) {
                T** tmp = redux::util::newArray<T>( imgRows2, imgCols2 );
                for( size_t r=0; r<imgRows2; ++r ) {
                    memcpy( tmp[r], img[r+firstRow]+firstCol, imgCols2*sizeof(T) );
                }
                redux::util::delArray( img );
                img = tmp;
                imgRows = imgRows2;
                imgCols = imgCols2;
            }

        }

        template <typename T>
        void normalizeIfMultiFrames( redux::image::Image<T>& img );
        
        template <typename T, typename U>
        void make_mask( T* input, U* mask, size_t ySize, size_t xSize, double thres=0, int smooth=5, bool filterLarger=false, bool invert=false ) {
#ifdef RDX_WITH_OPENCV
            cv::Mat inMat( ySize, xSize, cv::cvType<T>(), input );
            cv::Mat maskMat( ySize, xSize, cv::cvType<U>(), mask );
            cv::make_mask( inMat, maskMat, thres, smooth, filterLarger, invert );
#else
            std::cerr << "make_mask is not yet implemented for non-OpenCV builds." << std::endl;
#endif            
        }
        
        template <typename T, typename U>
        void inpaint( T* img, U* mask, T* out, size_t ySize, size_t xSize, double radius, int flags=0 );
        
        template <typename T, typename U>
        void resize( T* in, size_t inSizeY, size_t inSizeX, U* out, size_t outSizeY, size_t outSizeX );

        
    }   // image

}   // redux


#endif  // REDUX_IMAGE_UTIL_HPP
