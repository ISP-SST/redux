#ifndef REDUX_IMAGE_FOURIERTRANSFORM_HPP
#define REDUX_IMAGE_FOURIERTRANSFORM_HPP

#include "redux/types.hpp"
#include "redux/util/array.hpp"
#include "redux/util/point.hpp"
#include "redux/util/stringutil.hpp"

#include <cassert>
#include <mutex>
#include <numeric>

namespace redux {

    namespace image {
        
        
        // TODO: optimize 
        // TODO: support >2 dimensions.

        enum FT_FLAGS { REORDER_FT=1,           //!< Re-order transform to have the 0 frequency at the center (FFTW expects it at (0,0)
                        NORMALIZE_FT,           //!< Normalize transform
                        FULLCOMPLEX=4,          //!< Force full complex format (default is to use "r2c half-complex" format for real input data)
                        REORDER_IMG=8,          //!< Reorder input data before performing the FFT.
                        FT_ALL=15
                      };

        class FourierTransform : public redux::util::Array<complex_t> {
            
            struct PlansContainer {
                PlansContainer() { fftw_init_threads(); };
                ~PlansContainer(){ fftw_cleanup_threads(); };
                std::mutex mtx;
            };

        public:
            struct Plan
#ifdef RDX_TRACE_MEM
            : public redux::util::TraceObject<Plan>
#endif
            {
                typedef std::shared_ptr<const Plan> Ptr;
                enum TYPE { R2C=1, C2C };
                struct Index {
                    Index(const std::vector<size_t>& dims, TYPE t, uint8_t nt);
                    Index( size_t sizeY, size_t sizeX, TYPE t, uint8_t nt);
                    bool operator<( const Index& rhs ) const;
                    TYPE tp;
                    uint8_t nThreads;
                    std::vector<size_t> sizes;
                } id;
                fftw_plan forward_plan, backward_plan;
                explicit Plan( const Index& );
                ~Plan();
                static Plan::Ptr get(const std::vector<size_t>& dims, Plan::TYPE tp, uint8_t nThreads=1);
                static Plan::Ptr get(size_t sizeY, size_t sizeX, Plan::TYPE tp, uint8_t nThreads=1);
                static void clear(void);
                static PlansContainer pc;
                void init( void );
                void forward( double* __restrict__ in, fftw_complex* __restrict__ out ) const ;
                void backward( fftw_complex* __restrict__ in, double* __restrict__ out ) const;

                bool operator<( const Plan& rhs ) const { return (id < rhs.id); };
            };

            
            ~FourierTransform() { redux::util::Array<complex_t>::clear(); resize(); };
            FourierTransform();
            FourierTransform( const FourierTransform& );
            FourierTransform( FourierTransform&& );
            FourierTransform( size_t ySize, size_t xSize, int flags=0, uint8_t nThreads=1 );
            template <typename T> FourierTransform( const Array<T>& rhs, int flags=0, uint8_t nT=1 ) :
                centered (false), normalized (false), currentFlags(flags), nThreads(nT),
                inputSize(0), ftSize(0), inPixels(0), ftPixels(0), currentBlockSize(0) {
                    const auto& dims = rhs.dimensions(true);
                    if( dims.size() != 2 ) {
                        throw std::logic_error ("FourierTransform only supports 2 dimensions at the moment: "
                            + redux::util::printArray(dimensions(), "dims"));
                    }
                    const T* __restrict__ tPtr = rhs.get();
                    std::shared_ptr<T> tmp;
                    if( !rhs.dense() ) {
                        tmp = redux::util::rdx_get_shared<T>( dims[0]*dims[1] );
                        tPtr = const_cast<const T*>( tmp.get() );
                        rhs.template copyTo<T>( tmp.get() );
                    }
                    init( tPtr, dims[0], dims[1], flags, nT );
            }
            template <typename T> FourierTransform( const T* in, size_t ySize, size_t xSize, int flags=0, uint8_t nT=1 ) :
                centered (false),
                normalized (false), currentFlags(flags), inputSize(0), ftSize(0), inPixels(0), ftPixels(0), currentBlockSize(0) {
                    init( in, ySize, xSize, flags, nT );
            }

            void center( bool force=false );
            void decenter( bool force=false );
            
            template <typename T>
            static void reorder( redux::util::Array<T>& in ) {
                reorder( in.get(), in.dimSize(0), in.dimSize(1) );
            }
            void reorder( void ) {
                if( currentFlags&FULLCOMPLEX ) {
                    reorderInto( ftPtr, inputSize.y , inputSize.x, tmpPtr );
                    std::swap( ftPtr, tmpPtr );
                    wrap();
                } else {
                    //  TODO: implement reorderHalf
                }
            }
            template <typename T, typename U>
            static void reorderInto( const T* __restrict__ inSmall, size_t inSizeY, size_t inSizeX,
                                           U* __restrict__ outBig, size_t outSizeY=0, size_t outSizeX=0,
                                           bool center=false ) {

                if( !inSmall || !outBig ) {
                    throw std::logic_error("FourierTransform::reorderInto: nullptr   (in=" +
                        redux::util::hexString(inSmall) + "  out=" + redux::util::hexString(outBig) + ")" );
                }
                
                if( inSmall == reinterpret_cast<T*>(outBig) ) {
                    throw std::logic_error("FourierTransform::reorderInto can not have the same in/out pointers, it will mangle data!");
                }
                
                if( !outSizeY && !outSizeX ) {          // no outSize -> use same as inSize
                    outSizeY = inSizeY;
                    outSizeX = inSizeX;
                }
                
                if( (inSizeY>outSizeY) || (inSizeX>outSizeX) ) {
                    throw std::logic_error("FourierTransform::reorderInto: Target has to be larger, or equal to, the source.");
                }
                
                if( (inSizeY==outSizeY) && (inSizeX==outSizeX) ) {      // same size -> center == nocenter
                    center = false;
                }
                
                const size_t inMidY = inSizeY/2;
                const size_t inMidX = inSizeX/2;
                const size_t oddY = inSizeY%2;      // for odd sizes, the lower-block is 1px larger, in both x & y
                const size_t oddX = inSizeX%2;
                const size_t diffY = (outSizeY-inSizeY);
                const size_t diffX = (outSizeX-inSizeX);
                const size_t diffYh = diffY/2;
                const size_t diffXh = diffX/2;
                
                const int wDx = (center?(diffXh):0);
                const int eDx = (center?(diffXh-diffX):0);
                const int sDy = (center?(diffYh*outSizeX):0);
                const int nDy = (center?((diffYh-diffY)*outSizeX):0);

                // define pointers to the 4 quadrants of input (southWestIn = inSmall)
                const T* const __restrict__ southEastIn = inSmall + inMidX + oddX;
                const T* const __restrict__ northWestIn = inSmall + (inMidY + oddY)*inSizeX;
                const T* const __restrict__ northEastIn = southEastIn + (inMidY + oddY)*inSizeX;
                
                // define pointers to the 4 quadrants of output
                U* const __restrict__ southWestOut = outBig + sDy+wDx;
                U* const __restrict__ southEastOut = outBig + outSizeX - inMidX - oddX + sDy+eDx;
                U* const __restrict__ northWestOut = outBig + (outSizeY-inMidY-oddY)*outSizeX + nDy+wDx;
                U* const __restrict__ northEastOut = outBig + (outSizeY-inMidY-oddY+1)*outSizeX - inMidX - oddX + nDy+eDx;

                for( size_t y = 0; y < inMidY; ++y ) {
                    std::copy_n( northEastIn+y*inSizeX, inMidX, southWestOut+y*outSizeX );
                    std::copy_n( northWestIn+y*inSizeX, inMidX+oddX, southEastOut+y*outSizeX );
                }
                for( int y(inMidY+oddY-1); y >= 0; --y ) {
                    std::copy_n( inSmall+y*inSizeX, inMidX+oddX, northEastOut+y*outSizeX );
                    std::copy_n( southEastIn+y*inSizeX, inMidX, northWestOut+y*outSizeX );
                }
            }
            
            
            template <typename T>
            static void reorder( T* __restrict__ in, size_t inSizeY, size_t inSizeX ) {
                const size_t N(inSizeX*inSizeY);
                std::unique_ptr<T[]> buf( new T[N] );
                reorderInto( in, inSizeY, inSizeX, buf.get() );
                std::copy_n( buf.get(), N, in );
            }


            void reorderHalf( void ) {
                // TODO
            }
            
            FourierTransform reordered (void) const;
            
            
            FourierTransform& operator*=( const FourierTransform& rhs );
            FourierTransform& operator=( const FourierTransform& rhs );

            template <typename T>
            const FourierTransform& operator=( const T& val ) { Array<complex_t>::operator=(val); return *this; }
            template <typename T>
            const FourierTransform& operator*=( const T& val ) { Array<complex_t>::operator*=(val); return *this; }
            
            template <typename T> int ft( const T* in, complex_t* out, int flags=0 ) const;
            template <typename T> void ft( const T* __restrict__ in, int flags=-1 ) {
                if( flags < 0 ) {   // use current settings
                    flags = currentFlags;
                } else {    // called with flags set -> re-initialize
                    if( flags&FULLCOMPLEX ) {
                        ftSize.x = inputSize.x;
                        ftPixels = inPixels;
                    } else {
                        ftSize.x = inputSize.x/2+1;
                        ftPixels = ftSize.y*ftSize.x;
                    }
                    wrap();
                }
                currentFlags = ft( in, ftPtr, flags );
                centered = currentFlags&REORDER_FT;
                normalized = currentFlags&NORMALIZE_FT;
            }
            template <typename T> int ift( const complex_t* in, T* out, int flags=0 ) const;
            template <typename T> void ift( T* __restrict__ out, int flags=0 ) const {
                flags |= currentFlags;
                if( centered ) flags |= REORDER_FT;         // might be manually centered/normalized.
                if( normalized ) flags &= ~NORMALIZE_FT;
                else flags |= NORMALIZE_FT;
                ift( ftPtr, out, flags );
            }
            template <typename T> void ift( redux::util::Array<T>& out, int flags=0 ) const {
                util::PointI outSize( out.dimSize(0), out.dimSize(1) );
                if( (outSize != inputSize) || !out.dense() ) {
                    out.resize( inputSize.y, inputSize.x );
                }
                ift( out.get(), flags );
            }

            void getIFTx( double* out );        //!> NOTE fftw:c2r is destructive, so only use these direct functions in one-shot contexts.
            void getIFTx( complex_t* out ) const;
            void getIFTn( complex_t* out ) const;
            
            void getIFT( double* out ) { getIFTx(out); };
            void getIFT( complex_t* out ) const { getIFTx(out); };
            
            // Conjugate FT
            template <typename T> static void conj( const complex_t* __restrict__ in, T* __restrict__ out, size_t N ) {
                std::transform( in, in+N, out, [](const complex_t&a){ return std::conj(a); } );
            }
            template <typename T> inline void conj( T* __restrict__ out ) const { conj( ftPtr, out, ftPixels ); }
            void conj( void ) { std::transform( ftPtr, ftPtr+ftPixels, ftPtr, [](const complex_t&a){ return std::conj(a); } ); };

            // Norm of FT (used for power-spectrum & auto-correlation)
            template <typename T> static void norm( const complex_t* __restrict__ in, T* __restrict__ out, size_t N ) {
                transform( in, in+N, out, [](const complex_t&a){ return std::norm(a); } );
            }
            template <typename T> inline void norm( T* __restrict__ out ) const { norm( ftPtr, out, ftPixels ); }
            void norm( void ) { transform( ftPtr, ftPtr+ftPixels, ftPtr, [](const complex_t&a){ return std::norm(a); } ); };
            void norm( double scale ) { transform( ftPtr, ftPtr+ftPixels, ftPtr, [scale](const complex_t&a){ return std::norm(a)*scale; } ); };

            template <typename T=float> redux::util::Array<T> power( bool center=true ) const;
            template <typename T> void power( T* __restrict__ out, bool center=true ) const {
                if( center && !centered ) {
                    double* tmpP = reinterpret_cast<double*>( const_cast<complex_t*>(tmpPtr) );
                    norm( ftPtr, tmpP, ftPixels );
                    reorderInto( tmpP, ftSize.y, ftSize.x, out );
                } else {
                    norm( ftPtr, out, ftPixels );
                }
                normalize( out, ftPixels );
            }
            
            double noise( int mask=0, double cutoff=0 ) const;
            
            static void normalize( FourierTransform& );
            template <typename T> void normalize( T* __restrict__ data, size_t N, bool force=false ) const {
                if( force || !normalized ) {
                    const double nrm = 1.0/inPixels;
                    std::transform( data, data+N, data, [nrm](const T& a){ return a*nrm; } );
                }
            }
            void normalize( bool force=false ) {
                normalize( ftPtr, ftPixels, force );
                normalized = true;
            }

            template <typename T, typename U>
            void autocorrelate( const T* in, U* out, bool center=false ) const {
                ft( in, tmpPtr2, FULLCOMPLEX );
                std::transform( tmpPtr2, tmpPtr2+inPixels, tmpPtr2, [](const complex_t&a){ return std::norm(a); } );
                ift( tmpPtr2, out, NORMALIZE_FT|FULLCOMPLEX );
                if( center ) reorder( out, inputSize.y, inputSize.x );
            }
            template <typename T>
            inline void autocorrelate( T* inout, bool center=false ) const {
                autocorrelate( inout, inout, center );
            }
            template <typename T>
            static void autocorrelate( redux::util::Array<T>& inout, bool center=false ) {
                FourierTransform ft( inout, FULLCOMPLEX );
                ft.norm();
                ft.ift( inout.get(), NORMALIZE_FT|FULLCOMPLEX );
                if( center ) reorder( inout.get(), ft.inputSize.y, ft.inputSize.x );
            }
            template <typename T, typename U>
            static void autocorrelate( const redux::util::Array<T>& in, redux::util::Array<T>& out, bool center=false ) {
                FourierTransform ft( in, FULLCOMPLEX );
                ft.autocorrelate( in.get(), out.get(), center );
            }
            template <typename T>
            static void autocorrelate( T* inout, size_t nY, size_t nX, bool center=false ) {
                FourierTransform ft( inout, nY, nX, FULLCOMPLEX );
                ft.autocorrelate( inout, inout, center );
            }

            template <typename T, typename U>
            void convolve( const T* in, U* out ) const {
                int flags(0);
                if( currentFlags&FULLCOMPLEX ) flags |= FULLCOMPLEX;
                int iflags = flags;
                if( !(currentFlags&REORDER_IMG) ) flags |= REORDER_IMG;
                if( centered ) flags |= REORDER_FT;
                if( !normalized ) flags |= NORMALIZE_FT;
                ft( in, tmpPtr2, flags );
                std::transform( tmpPtr2, tmpPtr2+ftPixels, ftPtr, tmpPtr2, std::multiplies<complex_t>() );
                if( centered ) {
                    reorderInto( tmpPtr2, inputSize.y, inputSize.x, tmpPtr );
                    std::swap( tmpPtr, tmpPtr2 ); // swap the temps, because ift will use tmpPtr internally
                }
                ift( tmpPtr2, out, iflags );
            }

            template <typename T, typename U>
            void convolve( const redux::util::Array<T>& in, redux::util::Array<U>& out ) const {
                if( !out.sameSizes(in) ) {
                    out.resize( in. dimensions() );
                }
                convolve( in.get(), out.get() );
            }
            
            template <typename T>
            redux::util::Array<T> convolve( const redux::util::Array<T>& in ) const {
                redux::util::Array<double> tmp( in.dimensions() );
                convolve( in.get(), tmp.get() );
                return tmp.copy<T>();
            }

            template <typename T>
            void convolveInPlace( redux::util::Array<T>& inout ) const {
                convolve( inout.get(), inout.get() );
            }

            template <typename T, typename U>
            void correlate( const T* in, U* out ) const {
                int flags(0);
                if( currentFlags&FULLCOMPLEX ) flags |= FULLCOMPLEX;
                int iflags = flags|REORDER_IMG;
                if( centered ) flags |= REORDER_FT;
                if( !normalized ) flags |= NORMALIZE_FT;
                ft( in, tmpPtr2, flags );
                std::transform( tmpPtr2, tmpPtr2+ftPixels, ftPtr, tmpPtr2, 
                                [](const complex_t& a, const complex_t& b) { return b*std::conj(a); } );
                if( centered ) {
                    reorderInto( tmpPtr2, inputSize.y, inputSize.x, tmpPtr );
                    std::swap( tmpPtr, tmpPtr2 ); // swap the temps, because ift will use tmpPtr internally
                }
                ift( tmpPtr2, out, iflags );
            }

            template <typename T, typename U>
            void correlate( const redux::util::Array<T>& in, redux::util::Array<U>& out ) const {
                if( !out.sameSizes(in) ) {
                    out.resize( in. dimensions() );
                }
                correlate( in.get(), out.get() );
            }
            
            template <typename T>
            redux::util::Array<T> correlate( const redux::util::Array<T>& in ) const {
                redux::util::Array<double> tmp( in.dimensions() );
                correlate( in.get(), tmp.get() );
                return tmp.copy<T>();
            }

            template <typename T>
            void correlateInPlace( redux::util::Array<T>& inout ) const {
                correlate( inout.get(), inout.get() );
            }

           
            // Accessors to members
            Plan::Ptr getFullPlan (void) const { return plan_full; };
            Plan::Ptr getHalfPlan (void) const { return plan_half; };
            const util::PointI& getInputSize(void) const { return inputSize; };
            const util::PointI& getFtSize(void) const { return ftSize; };
            const size_t& getInputPixels(void) const { return inPixels; };
            const size_t& getFtPixels(void) const { return ftPixels; };
            const int& getFlags(void) const { return currentFlags; };
            void setThreads(uint8_t nT) { nThreads = nT; };
            const uint8_t& getThreads(void) const { return nThreads; };
            bool isFullComplex( void ) const { return (currentFlags&FULLCOMPLEX); };
            void setCentered( bool c=true ) { centered = c; };
            bool isCentered( void ) const { return centered; };
            void setNormalized( bool n=true ) { normalized = n; };
            bool isNormalized( void ) const { return normalized; };
            
            void init( size_t ySize, size_t xSize, int flags=0, uint8_t nT=1 );
            template <typename T> void init( const T*, size_t, size_t, int flags=0, uint8_t nT=1 );
            template <typename T> void init( const Array<T>& in, int flags=0, uint8_t nT=1 ) {
                init( in.get(), in.dimSize(0), in.dimSize(1), flags, nT );
            }

            void resize( size_t blockSize=0 );
            
        private:
            
            void wrap( void );
            
            Plan::Ptr plan_full;
            Plan::Ptr plan_half;

            bool centered;                     // true if the FT is centered, output from FFTW is NOT centered
            bool normalized;
            int currentFlags;
            uint8_t nThreads;
            
            util::PointI inputSize;
            util::PointI ftSize;
            size_t inPixels;
            size_t ftPixels;
            size_t currentBlockSize;
            
            std::shared_ptr<complex_t> ftData;
            complex_t* ftPtr;
            mutable complex_t* tmpPtr;
            mutable complex_t* tmpPtr2;
            
        };
       

    }   // image

}   // redux


#endif  // REDUX_IMAGE_FOURIERTRANSFORM_HPP
