#include "redux/file/fileana.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/image.hpp"
#include "redux/image/grid.hpp"
#include "redux/image/utils.hpp"
#include "redux/math/functions.hpp"
#include "redux/util/stringutil.hpp"

#include <iostream>

//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

using namespace redux::file;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

using namespace boost::unit_test;
namespace tt = boost::test_tools;

namespace {
    
    const float tol = 1E-7;     // global tolerance, override locally if this is too strict.
    
    template <typename T> void fft_construct( size_t nY, size_t nX=0, int flags=0, uint8_t nThreads=1 ) {
        
        if( nY == 0 ) throw runtime_error("fft_construct can not be called with size=0");
        if( nX == 0 ) nX = nY;              // make square if only nX is given
        
        PointI nyquist_pos( 0, 0 );
        size_t ft_sizeX = nX/2+1;
        if( (flags&FULLCOMPLEX) ) {         // if full-complex
            ft_sizeX = nX;
            if( (flags&REORDER_FT) ) {      // if full-complex & re-ordered
                nyquist_pos.y = nY/2;
                nyquist_pos.x = nX/2;
            }
        }

        PointI inputSize( nY, nX );
        Array<T> input( nY, nX );           // Input to use
        input = 1;                          // Constant, will make a delta FT.
        
        Array<double> tmpD;       // Temporaries
        Array<complex_t> tmpC;
        
        {       // Default constructor, should make an empty FT
            FourierTransform empty;
            BOOST_TEST( empty.nDimensions() == 0 );
            BOOST_TEST( empty.nElements() == 0 );
            BOOST_TEST( empty.getInputSize() == PointI(0,0) );
            BOOST_TEST( empty.getFtSize() == PointI(0,0) );
            BOOST_TEST( empty.getInputPixels() == 0 );
            BOOST_TEST( empty.getFtPixels() == 0 );
            BOOST_TEST( empty.getThreads() == 1 );
            BOOST_TEST( empty.get() == nullptr);
        }

        {       // Construct by size, but without supplying data
            FourierTransform ft( nY, nX, flags, nThreads );
            BOOST_TEST( ft.nDimensions() == 2 );
            BOOST_TEST( ft.nElements() == nY*ft_sizeX );
            BOOST_TEST( ft.getInputSize() == PointI(nY,nX) );
            BOOST_TEST( ft.getFtSize() == PointI(nY,ft_sizeX) );
            BOOST_TEST( ft.getInputPixels() == nX*nY );
            BOOST_TEST( ft.getFtPixels() == nY*ft_sizeX );
            BOOST_TEST( ft.getThreads() == nThreads );
            BOOST_TEST( ft.get() != nullptr);
        }
        
        {       // Construct from input
            FourierTransform ft1( input, flags, nThreads );
            FourierTransform ft2( input.get(), nY, nX, flags, nThreads );   // These two constructors should yield identical results.
            PointI expectedFTsize = inputSize;
            if( !(ft1.getFlags()&FULLCOMPLEX) ) {
                expectedFTsize.x = expectedFTsize.x/2+1;
            }
            size_t expectedElements = expectedFTsize.y*expectedFTsize.x;
            BOOST_TEST( ft1.nDimensions() == 2 );
            BOOST_TEST( ft1.nElements() == expectedElements );
            BOOST_TEST( ft1.getInputSize() == inputSize );
            BOOST_TEST( ft1.getFtSize() == expectedFTsize );
            BOOST_TEST( ft1.getInputPixels() == nX*nY );
            BOOST_TEST( ft1.getFtPixels() == expectedElements );
            BOOST_TEST( ft1.getThreads() == nThreads );
            BOOST_TEST( ft1.get() != nullptr);
            
            BOOST_TEST( ft2.nDimensions() == 2 );
            BOOST_TEST( ft2.nElements() == expectedElements );
            BOOST_TEST( ft2.getInputSize() == inputSize );
            BOOST_TEST( ft2.getFtSize() == expectedFTsize );
            BOOST_TEST( ft2.getInputPixels() == nX*nY );
            BOOST_TEST( ft2.getFtPixels() == expectedElements );
            BOOST_TEST( ft2.getThreads() == nThreads );
            BOOST_TEST( ft2.get() != nullptr);

            tmpC.resize(nY,nX);
            tmpC.zero();                    // Set to expected results, i.e. a delta-function at the nyquist location
            double exptectedValue = (nX*nY);
            if( flags&NORMALIZE_FT ) {
                exptectedValue = 1.0;
            }
            PointI expectedPos( 0, 0 );
            if( (ft1.getFlags()&REORDER_FT) ) {      // if full-complex & re-ordered
                nyquist_pos.y = nY/2;
                nyquist_pos.x = nX/2;
            }
            tmpC( nyquist_pos.y, nyquist_pos.x ) = exptectedValue;
            for( size_t y(0); y<nY; ++y) { 
                for( size_t x(0); x<ft_sizeX; ++x ) {
                    BOOST_TEST( abs(ft1(y,x) - tmpC(y,x)) < tol );     // check that FT matches expectations
                    BOOST_TEST( abs(ft2(y,x) - tmpC(y,x)) < tol );     // check that FT matches expectations
                    //BOOST_TEST( ft1(y,x) == tmpC(y,x) );     // check that FT matches expectations
                    //BOOST_TEST( ft2(y,x) == tmpC(y,x) );     // check that FT matches expectations
                }
            }

            // This constructor should throw if input is not 2D.
            tmpD.resize(nX);
            BOOST_CHECK_THROW( FourierTransform ft3( tmpD, flags, nThreads ), logic_error );
            tmpD.resize(nX,nY,nX);
            BOOST_CHECK_THROW( FourierTransform ft3( tmpD, flags, nThreads ), logic_error );
        }

        {       // Copy-construct
            FourierTransform ft1( nY, nX, flags, nThreads );    // The sanity of ft1 was tested above.
            FourierTransform ft2( ft1 );
            BOOST_TEST( ft2.nDimensions() == ft1.nDimensions() );
            BOOST_TEST( ft2.nElements() == ft1.nElements() );
            BOOST_TEST( ft2.getInputSize() == ft1.getInputSize() );
            BOOST_TEST( ft2.getFtSize() == ft1.getFtSize() );
            BOOST_TEST( ft2.getInputPixels() == ft1.getInputPixels() );
            BOOST_TEST( ft2.getFtPixels() == ft1.getFtPixels() );
            BOOST_TEST( ft2.getThreads() == ft1.getThreads() );
            BOOST_TEST( ft2.get() != ft1.get() );       // should not share data
            BOOST_TEST( ft2.get() != nullptr );         // ...but should have it's own
        }

    }
    
    
    template <typename T> void fft_operators( size_t nX, size_t nY=0 ) {
        
        if( nX == 0 ) throw runtime_error("fft_operators can not be called with size=0");
        if( nY == 0 ) nY = nX;
        
        size_t halfSizeX = nX/2;
        Array<T> input( nY, nX );
        Array<double> tmpD( nY, nX );
        Array<complex_t> tmpC( nY, nX );

        {   // *= operator  (no Fourier-transforming yet, just basic operations on the FourierTransform class
            FourierTransform fullFT( nY, nX, FULLCOMPLEX );    // full-complex
            BOOST_TEST( fullFT.dimSize(0) == nY );
            BOOST_TEST( fullFT.dimSize(1) == nX );

            FourierTransform halfFT( nY, nX );                   // half-complex
            BOOST_TEST( halfFT.dimSize(0) == nY );
            BOOST_TEST( halfFT.dimSize(1) == halfSizeX+1 );

            complex_t salt(rand()%100,rand()%100);               // some noise to add to test-data.
            
            // test *= operator for full-complex to full-complex AND half-complex to half-complex
            for( unsigned int y=0; y<nY; ++y) { 
                for( unsigned int x=0; x<halfSizeX+1 ; ++x) {
                    complex_t val(y,x);
                    val += salt;
                    halfFT(y,x) = fullFT(y,x) = val;
                    if( x ) fullFT(y,nX-x) = val;
                }
            }
            fullFT *= fullFT;
            halfFT *= halfFT;
            for( unsigned int y=0; y<nY; ++y) { 
                for( unsigned int x=0; x<halfSizeX+1; ++x) {
                    complex_t val(y,x);
                    val += salt;
                    val *= val;
                    BOOST_TEST( norm(halfFT(y,x)-val) < tol );
                    BOOST_TEST( norm(fullFT(y,x)-val) < tol );
                    if( x ) BOOST_TEST( norm(fullFT(y,nX-x)-val) < tol );
                }
            }

            // test *= operator for full-complex to half-complex
            fullFT = salt;          // values in the high-x half should not matter, so set them to some garbage.
            for( unsigned int y=0; y<nY; ++y ) { 
                for( unsigned int x=0; x<halfSizeX+1; ++x ) {
                    complex_t val(y,x);
                    val += salt;
                    halfFT(y,x) = val;
                    fullFT(y,x) = val+salt*1.44;
                }
            }
            halfFT *= fullFT;
            for( unsigned int y=0; y<nY; ++y) { 
                for( unsigned int x=0; x<halfSizeX+1; ++x) {
                    complex_t val (y,x);
                    val += salt;
                    complex_t fullFTv = fullFT(y,x);
                    BOOST_TEST( norm(halfFT(y,x) - val*fullFTv ) < tol );
                }
            }

            // test *= operator for half-complex to full-complex
            for( unsigned int y=0; y<nY; ++y) { 
                for( unsigned int x=0; x<nX; ++x) {
                    complex_t val(y,x);
                    val += salt;
                    fullFT(y,x) = val;
                    if( x <= halfSizeX ) halfFT(y,x) = val+salt*0.13;
                }
            }
            fullFT *= halfFT;
            for( unsigned int y=0; y<nY; ++y) { 
                for( unsigned int x=0; x<halfSizeX; ++x) {
                    complex_t val(y,x);
                    val += salt;
                    complex_t halfFTv;
                    if( x <= halfSizeX ) halfFTv = halfFT(y,x);
                    else if( y ) halfFTv = ::conj(halfFT( nX-y, nX-x ));
                    else halfFTv = ::conj(halfFT(0,nX-x));
                    BOOST_TEST( norm(fullFT(y,x) - val*halfFTv ) < tol );
                }
            }

        }

    }
    
    
    template <typename T> void forward_FFT ( size_t nX, size_t nY=0 ) {
        
        if( nX == 0 ) throw runtime_error("forward_FFT can not be called with size=0");
        if( nY == 0 ) nY = nX;
        
        size_t halfSizeX = (nX+1)/2;
        Array<T> input( nY, nX );
        Array<double> tmpD( nY, nX );
        Array<complex_t> tmpC( nY, nX );

        {   // test some simple transforms
            
            
            // Check the transform of a constant function (should be a delta)
            input = 1;
            FourierTransform fullFT( input, FULLCOMPLEX );
            FourierTransform halfFT( input );
            tmpC.zero();                    // Set to expected results, i.e. a delta-function
            tmpC(0,0) = (nX*nY);            // DC component is stored in (0,0), and is not normalized by default, i.e. should be equal to (nX*nY)
            for( size_t y(0); y<nY; ++y) { 
                for( size_t x(0); x<nX; ++x ) {
                    BOOST_TEST( norm(fullFT(y,x)-tmpC(y,x)) < tol );                         // check that fullFT matches expectations
                    if( x < halfSizeX ) BOOST_TEST( norm(halfFT(y,x)-tmpC(y,x)) < tol );    // check that halfFT matches expectations
                }
            }
            
            // The same, but with "REORDER" flag
            fullFT.init( input, FULLCOMPLEX|REORDER_FT );
            tmpC.zero();                                // Set to expected results, i.e. a delta-function
            tmpC(nY/2,nX/2) = (nX*nY);             // DC component is stored in (nY/2,nX/2), and is not normalized by default, i.e. should be equal to (nX*nY)
            for( size_t y(0); y<nY; ++y) { 
                for( size_t x(0); x<nX; ++x ) {
                    BOOST_TEST( norm(fullFT(y,x)-tmpC(y,x)) < tol );                         // check that fullFT matches expectations
                }
            }

            // The same, but with "NORMALIZE" flag 
            fullFT.init( input, FULLCOMPLEX|NORMALIZE_FT );
            halfFT.init( input, NORMALIZE_FT );
            tmpC.zero();                                // Set to expected results, i.e. a delta-function
            tmpC(0,0) = 1;                              // DC component is stored in (nY/2,halfSizeX), when normalized it should be 1.0
            for( size_t y(0); y<nY; ++y) { 
                for( size_t x(0); x<nX; ++x ) {
                    BOOST_TEST( norm(fullFT(y,x)-tmpC(y,x)) < tol );                         // check that fullFT matches expectations
                    if( x < halfSizeX ) BOOST_TEST( norm(halfFT(y,x)-tmpC(y,x)) < tol );    // check that halfFT matches expectations
                }
            }

            // Check the transform of a delta-function (centered at the Nyquist location). It should be a constant (=1)
            input.zero();
            input( 0, 0 ) = 1;
            fullFT.init( input, FULLCOMPLEX );
            halfFT.init( input );
            tmpC = complex_t(1);
            for( size_t y(0); y<nY; ++y) { 
                for( size_t x(0); x<nX; ++x ) {
                    BOOST_TEST( norm(fullFT(y,x)-tmpC(y,x)) < tol );                         // check that fullFT matches expectations
                    if( x < halfSizeX ) BOOST_TEST( norm(halfFT(y,x)-tmpC(y,x)) < tol );     // check that halfFT matches expectations
                }
            }
            
            // The same, but with the delta-function centered in inout and using the "REORDER_IMG" flag instead 
            input.zero();
            input( nY/2, halfSizeX ) = 1;
            fullFT.init( input, FULLCOMPLEX|REORDER_IMG );
            halfFT.init( input, REORDER_IMG );
            for( size_t y(0); y<nY; ++y) { 
                for( size_t x(0); x<nX; ++x ) {
                    BOOST_TEST( norm(fullFT(y,x)-tmpC(y,x)) < tol );                         // check that fullFT matches expectations
                    if( x < halfSizeX ) BOOST_TEST( norm(halfFT(y,x)-tmpC(y,x)) < tol );    // check that halfFT matches expectations
                }
            }

            // Check the transform of a cosine-function with period = nY in the y-direction, and nX/2 in the x-direction
            input.zero();
            for( size_t y(0); y<nY; ++y) {
                double y_val = cos( y*2.0*M_PI/nY );
                for( size_t x(0); x<nX; ++x ) {
                    input(y,x) = y_val + cos( x*4.0*M_PI/nX );
                }
            }
            
            tmpC.zero();
            tmpC(nY/2+1,nX/2) = tmpC(nY/2-1,nX/2) = 0.5;
            tmpC(nY/2,nX/2+2) = tmpC(nY/2,nX/2-2) = 0.5;
            fullFT.init( input, FULLCOMPLEX|NORMALIZE_FT|REORDER_FT );
            for( size_t y(0); y<3; ++y) { 
                for( size_t x(0); x<3; ++x ) {
                    BOOST_TEST( norm(fullFT(y,x)-tmpC(y,x)) < 1E-4 );  // check that fullFT matches expectations
                }
            }

        }
    
    }
    
    
    template <typename T> void backward_FFT ( size_t nX, size_t nY=0 ) {
        
        if( nX == 0 ) throw runtime_error("backward_FFT can not be called with size=0");
        if( nY == 0 ) nY = nX;

        Array<T> input( nY, nX );
        Array<double> tmpD1( nY, nX );
        Array<double> tmpD2( nY, nX );

        {   // test some simple transforms
            
            // Arbitrary function, we'll use a gaussian with fwhm (nY/2,nX/4) centered at (mid+2,mid+4)
            redux::math::gauss( input.get(), nY, nX, nY/2, nX/4, nY/2+2, nX/2+4 );

            FourierTransform fullFT( input, FULLCOMPLEX );
            FourierTransform halfFT( input );
            tmpD1.zero();
            tmpD2.zero();
            fullFT.ift( tmpD1.get() );
            halfFT.ift( tmpD2.get() );
            for( size_t y(0); y<nY; ++y) { 
                for( size_t x(0); x<nX; ++x ) {
                    BOOST_TEST( abs(tmpD1(y,x)-input(y,x)) < tol );     // check that we recover input function
                    BOOST_TEST( abs(tmpD2(y,x)-input(y,x)) < tol );
                }
            }

            fullFT.init( input, FULLCOMPLEX|NORMALIZE_FT );             // same test, but initialize as normalized FT
            halfFT.init( input, NORMALIZE_FT );
            tmpD1.zero();
            tmpD2.zero();
            fullFT.ift( tmpD1.get() );
            halfFT.ift( tmpD2.get() );
            for( size_t y(0); y<nY; ++y) { 
                for( size_t x(0); x<nX; ++x ) {
                    BOOST_TEST( abs(tmpD1(y,x)-input(y,x)) < tol );     // check that we recover input function
                    BOOST_TEST( abs(tmpD2(y,x)-input(y,x)) < tol );
                }
            }

            fullFT.init( input, FULLCOMPLEX|REORDER_FT );               // same test, but reordered FT
            tmpD1.zero();
            fullFT.ift( tmpD1.get() );
            for( size_t y(0); y<nY; ++y) { 
                for( size_t x(0); x<nX; ++x ) {
                    BOOST_TEST( abs(tmpD1(y,x)-input(y,x)) < tol );     // check that we recover input function
                }
            }

            fullFT.init( input, FULLCOMPLEX, REORDER_IMG );             // same test, but reorder img before FT
            halfFT.init( input, REORDER_IMG );
            tmpD1.zero();
            tmpD2.zero();
            fullFT.ift( tmpD1.get() ); //, REORDER_IMG );
            halfFT.ift( tmpD2.get() ); //, REORDER_IMG );
            for( size_t y(0); y<nY; ++y) { 
                for( size_t x(0); x<nX; ++x ) {
                    BOOST_TEST( abs(tmpD1(y,x)-input(y,x)) < tol );     // check that we recover input function
                    BOOST_TEST( abs(tmpD2(y,x)-input(y,x)) < tol );
                }
            }
            
        }
    
    }
    
}


void fft_methods( void ) {
    
    // test some elementary methods/operations (reorder, conj, normalize, power)

    {   // test reorder
        
        const size_t N = 100;
        array<int,N> input;     // generate a list of integers [1,...,100]
        std::iota( input.begin(), input.end(), 1 );
        
        // test some simple/small cases of odd/even sizes
        
        // 4x4 input
        array<int,N> tmpA;     // temporary array
        vector<int> tmpV;
        FourierTransform::reorderInto( input.data(), 4, 4, tmpA.data() );
        vector<int> expected = { 11, 12,  9, 10,        // input:  1  2  3  4
                                 15, 16, 13, 14,        //         5  6  7  8
                                  3,  4,  1,  2,        //         9 10 11 12
                                  7,  8,  5,  6 };      //        13 14 15 16 
        tmpV.assign( tmpA.begin(), tmpA.begin()+16 );
        BOOST_TEST( tmpV == expected );

        // Again/reverse (only for even sizes, for odd sizes reordering is not it's own inverse!!)
        FourierTransform::reorderInto( tmpV.data(), 4, 4, tmpA.data() );
        for( size_t i(0); i<16; ++i ) {
            BOOST_TEST( tmpA[i] == input[i] );
        }
        
        // 5x5 input
        FourierTransform::reorderInto( input.data(), 5, 5, tmpA.data() );
        expected = { 19, 20, 16, 17, 18,        // input:  1  2  3  4  5
                     24, 25, 21, 22, 23,        //         6  7  8  9 10
                      4,  5,  1,  2,  3,        //        11 12 13 14 15
                      9, 10,  6,  7,  8,        //        16 17 18 19 20
                     14, 15, 11, 12, 13 };      //        21 22 23 24 25
        tmpV.assign( tmpA.begin(), tmpA.begin()+25 );
        BOOST_TEST( tmpV == expected );
        
        // test reordering into a larger target
        tmpA.fill(0);
        FourierTransform::reorderInto( input.data(), 4, 4, tmpA.data(), 5, 6 );
        expected = { 11, 12,  0,  0,  9, 10,        // input:  1  2  3  4
                     15, 16,  0,  0, 13, 14,        //         5  6  7  8
                      0,  0,  0,  0,  0,  0,        //         9 10 11 12
                      3,  4,  0,  0,  1,  2,        //        13 14 15 16
                      7,  8,  0,  0,  5,  6 };
        tmpV.assign( tmpA.begin(), tmpA.begin()+30 );
        BOOST_TEST( tmpV == expected );
        
        // same, but wih center-option
        tmpA.fill(0);
        FourierTransform::reorderInto( input.data(), 4, 4, tmpA.data(), 5, 6, true );
        expected = { 0, 11, 12,  9, 10,  0,        // input:  1  2  3  4
                     0, 15, 16, 13, 14,  0,        //         5  6  7  8
                     0,  3,  4,  1,  2,  0,        //         9 10 11 12
                     0,  7,  8,  5,  6,  0,        //        13 14 15 16
                     0,  0,  0,  0,  0,  0 };
        tmpV.assign( tmpA.begin(), tmpA.begin()+30 );
        BOOST_TEST( tmpV == expected );
        
        // reorderInto should throw when called with nullptr
        BOOST_CHECK_THROW( FourierTransform::reorderInto( input.data(), 4, 4, (int*)nullptr ), logic_error );
        BOOST_CHECK_THROW( FourierTransform::reorderInto( (int*)nullptr, 4, 4, tmpA.data() ), logic_error );

        // reorderInto should throw when inSize > outSize
        BOOST_CHECK_THROW( FourierTransform::reorderInto( input.data(), 4, 4, tmpA.data(), 3, 4 ), logic_error );
        BOOST_CHECK_THROW( FourierTransform::reorderInto( input.data(), 4, 4, tmpA.data(), 4, 3 ), logic_error );

        // reorderInto should throw when inPtr == outPtr
        BOOST_CHECK_THROW( FourierTransform::reorderInto( tmpA.data(), 4, 4, tmpA.data(), 5, 5 ), logic_error );
        
    }

    {   // test center/decenter
        
        size_t nX = 64;
        size_t nY = 64;
        
        Array<double> peak( nY, nX );
        peak.zero();
        
        // Arbitrary function, we'll use a gaussian with fwhm (nY/2,nX/4) centered at (mid+2,mid+4)
        redux::math::gauss( peak.get(), nY, nX, nY/2, nX/4, nY/2+2, nX/2+4 );
        FourierTransform ft( peak, FULLCOMPLEX );
        
        FourierTransform ft_uncentered(ft);       // make copies, for comparisons
        FourierTransform ft_centered(ft);
        FourierTransform::reorder( ft_centered );
        
        BOOST_TEST( (ft == ft_uncentered) );    // trivial check, since we just made a copy
        
        // elementary checks, FT should be un-centered and full-complex.
        BOOST_TEST( ft.isFullComplex() );
        BOOST_TEST( !ft.isCentered() );
        ft.decenter();                          // should do nothing (not centered)
        BOOST_TEST( !ft.isCentered() );
        ft.center();                            // center/reorder the FT
        BOOST_TEST( ft.isCentered() );
        BOOST_TEST( (ft == ft_centered) );      // compare data
        ft.center();                            // should do nothing (already centered)
        BOOST_TEST( ft.isCentered() );
        ft.decenter();                          // decenter/reorder the FT
        BOOST_TEST( !ft.isCentered() );

        // do similar tests, but initalize FT to centered.
        FourierTransform ft2( peak, FULLCOMPLEX|REORDER_FT );
        ft_uncentered = FourierTransform(ft2);  // make copies, for comparisons
        ft_centered = FourierTransform(ft2);
        FourierTransform::reorder( ft_uncentered );
        
        BOOST_TEST( (ft2 == ft_centered) );    // trivial check, since we just made a copy
        
        // elementary checks, FT should be centered and full-complex.
        BOOST_TEST( ft2.isFullComplex() );
        BOOST_TEST( ft2.isCentered() );
        
        ft2.center();                           // should do nothing (already centered)
        BOOST_TEST( ft2.isCentered() );
        ft2.decenter();                         // decenter/reorder the FT
        BOOST_TEST( !ft2.isCentered() );
        BOOST_TEST( (ft == ft_uncentered) );    // compare data
        ft2.decenter();                         // should do nothing (not centered)
        BOOST_TEST( !ft2.isCentered() );
        ft2.center();                           // center/reorder the FT
        BOOST_TEST( ft2.isCentered() );

        Array<double> bla( nY, nX );
        bla.zero();
        ft.init( peak, FULLCOMPLEX|REORDER_FT );
        
        ft.ift( bla );
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                //BOOST_TEST( bla(y,x)-peak(y,x) < tol );
            }
        }
        
    }

        
    {   // test conjugate & norm
        
        size_t nX = 64;
        size_t nY = 64;
        
        Array<double> tmp( nY, nX );
        tmp.zero();
        
        // Arbitrary function, we'll use a gaussian with fwhm (nY/2,nX/4) centered at (mid+2,mid+4)
        redux::math::gauss( tmp.get(), nY, nX, nY/2, nX/4, nY/2+2, nX/2+4 );

        FourierTransform fullFT( tmp, FULLCOMPLEX );
        FourierTransform halfFT( tmp );

        FourierTransform ft_ref(fullFT);       // make reference-copy
        
        fullFT.conj();      // conjugate FT's
        halfFT.conj();
        
        for( size_t y(0); y<nY; ++y ) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(fullFT(y,x) - std::conj(ft_ref(y,x))) < tol );
                BOOST_TEST( abs(norm(fullFT(y,x)) - norm(ft_ref(y,x))) < tol );
                if( x <= nX/2 ) {
                    BOOST_TEST( abs(halfFT(y,x) - std::conj(ft_ref(y,x))) < tol );
                    BOOST_TEST( abs(norm(halfFT(y,x)) - norm(ft_ref(y,x))) < tol );
                }
            }
        }
        
    }
        
    {   // test power
        
        size_t nX = 64;
        size_t nY = 64;
        
        Array<double> tmp( nY, nX );
        tmp.zero();
        
        // Arbitrary function, we'll use a gaussian with fwhm (nY/2,nX/4) centered at (mid+2,mid+4)
        redux::math::gauss( tmp.get(), nY, nX, nY/2, nX/4, nY/2+2, nX/2+4 );

        FourierTransform fullFT( tmp, FULLCOMPLEX );
        
        Array<float> pwr = fullFT.power();

        // TODO: test for expected result
//         Ana::write( "data.f0", tmp );
//         Ana::write( "ft.f0", fullFT );
//         Ana::write( "pwr.f0", pwr );

        
    }
    
    
    {   // test convolution
        
        size_t nX = 64;
        size_t nY = 64;
        size_t midY = nY/2;
        size_t midX = nX/2;
        
        Array<double> data( nY, nX );       // almost-delta-function
        data.zero();
        data( midY, midX ) = data( midY, midX+1 ) = 0.5;
                           
        Array<double> psf( nY, nX );        // make a cross
        psf.zero();
        psf( midY, midX ) = 0.5;
        psf( midY+1, midX ) = psf( midY-1, midX )
                            = psf( midY, midX+1 ) = psf( midY, midX-1 ) = 0.125;
        
        Array<double> expected( nY, nX );
        expected.zero();
        Array<double> eview( expected, midY-2, midY+2, midX-2, midX+2 );
        Array<double> pview( psf, midY-2, midY+2, midX-2, midX+2 );
        eview.assign(pview);
        eview.shift(1,1);
        eview += pview;
        expected *= 0.5;

        FourierTransform fullFT( data, FULLCOMPLEX );   // test full-complex
        Array<double> tmp = psf.copy<double>();         // temporaries
        Array<double> tmp2 = data.copy<double>();
        fullFT.convolve( tmp.get(), tmp2.get() );       // test version with naked data-pointers
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }

        fullFT.convolve( tmp, tmp2 );      // test version with Array<T> for in/out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        tmp2 = fullFT.convolve( tmp );      // test version with Array<T> for in, and returning out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.convolveInPlace( tmp );      // test in-place version with Array<T>
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        
        fullFT.init( data, FULLCOMPLEX|REORDER_FT|NORMALIZE_FT );    // test full-complex + centered + normalized

        tmp = psf.copy<double>();         // temporaries
        tmp2 = data.copy<double>();
        fullFT.convolve( tmp.get(), tmp2.get() );       // test version with naked data-pointers
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.convolve( tmp, tmp2 );      // test version with Array<T> for in/out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        tmp2 = fullFT.convolve( tmp );      // test version with Array<T> for in, and returning out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.convolveInPlace( tmp );      // test in-place version with Array<T>
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        
        fullFT.init( data, FULLCOMPLEX|REORDER_FT|REORDER_IMG );    // test full-complex + centered + re-ordered input

        tmp = psf.copy<double>();         // temporaries
        tmp2 = data.copy<double>();
        fullFT.convolve( tmp.get(), tmp2.get() );       // test version with naked data-pointers
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.convolve( tmp, tmp2 );      // test version with Array<T> for in/out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        tmp2 = fullFT.convolve( tmp );      // test version with Array<T> for in, and returning out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.convolveInPlace( tmp );      // test in-place version with Array<T>
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        
        FourierTransform halfFT( data );            // test half-complex
        tmp = psf.copy<double>();                   // reset temporaries
        tmp2 = data.copy<double>();
        fullFT.convolve( tmp.get(), tmp2.get() );      // test version with naked data-pointers
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.convolve( tmp, tmp2 );      // test version with Array<T> for in/out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        tmp2 = fullFT.convolve( tmp );      // test version with Array<T> for in, and returning out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.convolveInPlace( tmp );      // test in-place version with Array<T>
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        
        halfFT.init( data, NORMALIZE_FT|REORDER_IMG );            // test half-complex + normalized +  re-ordered input
        tmp = psf.copy<double>();                   // reset temporaries
        tmp2 = data.copy<double>();
        fullFT.convolve( tmp.get(), tmp2.get() );      // test version with naked data-pointers
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.convolve( tmp, tmp2 );      // test version with Array<T> for in/out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        tmp2 = fullFT.convolve( tmp );      // test version with Array<T> for in, and returning out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.convolveInPlace( tmp );      // test in-place version with Array<T>
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
    }
    
    
    {   // test correlation
        
        size_t nX = 64;
        size_t nY = 64;
        size_t midY = nY/2;
        size_t midX = nX/2;
        
        Array<double> img1( nY, nX );       // almost-delta-function
        img1.zero();
        img1 ( midY, midX ) = img1( midY, midX+1 ) = 0.5;
                           
        Array<double> img2( nY, nX );        // make a cross
        img2.zero();
        img2( midY, midX ) = 0.5;
        img2( midY+1, midX ) = img2( midY-1, midX )
                             = img2( midY, midX+1 ) = img2( midY, midX-1 ) = 0.125;
        
        Array<double> expected( nY, nX );
        expected.zero();
        Array<double> eview( expected, midY-2, midY+2, midX-2, midX+2 );
        Array<double> view2 ( img2, midY-2, midY+2, midX-2, midX+2 );
        eview.assign( view2 );
        eview.shift(1,1);
        eview += view2;
        expected *= 0.5;
        
        FourierTransform fullFT( img1, FULLCOMPLEX );   // test full-complex
        Array<double> tmp = img2.copy<double>();         // temporaries
        Array<double> tmp2 = img1.copy<double>();
        fullFT.correlate( tmp.get(), tmp2.get() );       // test version with naked data-pointers
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.correlate( tmp, tmp2 );      // test version with Array<T> for in/out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        tmp2 = fullFT.correlate( tmp );      // test version with Array<T> for in, and returning out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.correlateInPlace( tmp );      // test in-place version with Array<T>
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        
        FourierTransform halfFT( img1 );             // test half-complex
        tmp = img2.copy<double>();                   // reset temporaries
        tmp2 = img1.copy<double>();
        fullFT.correlate( tmp.get(), tmp2.get() );      // test version with naked data-pointers
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.correlate( tmp, tmp2 );      // test version with Array<T> for in/out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        tmp2 = fullFT.correlate( tmp );      // test version with Array<T> for in, and returning out
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp2(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.correlateInPlace( tmp );      // test in-place version with Array<T>
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        
    }
    
    
    {   // test auto-correlation
        
        size_t nX = 64;
        size_t nY = 64;
        
        FourierTransform fullFT( nY, nX, FULLCOMPLEX );
        Array<double> input( nY, nX );
        input.zero();
        input( 5, 5 ) = input( 6, 6 ) = input( 6, 5 ) = input( 5, 6 ) = 3;
        Array<double> expected( nY, nX );
        expected.zero();
        expected( 0, 0 ) = 36;
        expected( 0, 1 ) = expected( 1, 0 ) = expected( nY-1, 0 ) = expected( 0, nX-1 ) = 18;
        expected( 1, 1 ) = expected( nY-1, nX-1 ) = expected( nY-1, 1 ) = expected( 1, nX-1 ) = 9;
       
        Array<double> tmp = input.copy<double>();       // just a copy we can destroy
        fullFT.autocorrelate( tmp.get() );              // test version with naked data-pointer
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        tmp = input.copy<double>();
        fullFT.autocorrelate( tmp );                    // test version with Array<T>
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.autocorrelate( input.get(), tmp.get(), false );  // test version with separate in/out pointers
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        
        FourierTransform::reorder( expected.get(), nY, nX );        // test with center-flag
        
        tmp = input.copy<double>();
        fullFT.autocorrelate( tmp.get(), true );              // test version with naked data-pointer
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        tmp = input.copy<double>();
        fullFT.autocorrelate( tmp, true );                    // test version with Array<T>
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        fullFT.autocorrelate( input.get(), tmp.get(), true );  // test version with separate in/out pointers
        for( size_t y(0); y<nY; ++y) { 
            for( size_t x(0); x<nX; ++x ) {
                BOOST_TEST( abs(tmp(y,x)-expected(y,x)) < tol );
            }
        }
        
//         Ana::write( "input.f0", input );
//         Ana::write( "tmp.f0", tmp );
//         Ana::write( "expected.f0", expected );
//         return;
        
        
        
    }
        
}



void fft_constructors( void ) {
    
    for( int f(0); f<9; ++f ) {
        for( int t(1); t<5; t+=2 ) {
            fft_construct<double>( 64, 0, f, t );
            fft_construct<float>( 53, 0, f, t );
            fft_construct<complex_t>( 71, 53, f, t );
            fft_construct<int32_t>( 53, 88, f, t );
            //fft_construct<int16_t>( 53, 88, f, t );
            //fft_construct<int8_t>( 53, 88, f, t );
            //fft_construct<double>( 64, 53, f, t );
            //fft_construct<float>( 53, 88, f, t );
            //fft_construct<double>( 64, 53, f, t );
        }
    }
}


void fft_manipulations( void ) {
    

//    fft_operators<double>( 64 );
//    fft_operators<float>( 53 );
    
}


void forward_fft( void ) {
    
    forward_FFT<double>( 64, 64 );
    forward_FFT<double>( 63, 64 );
    forward_FFT<float>( 53, 64 );
    forward_FFT<int32_t>( 53, 64 );
    
/*
    Array<double> peak(nPoints,nPoints);
    peak.zero();
    redux::math::gauss( peak.get(), nPoints, nPoints, nPoints/20, nPoints/40, nPoints/2, nPoints/2);    // centered gaussian with fwhm (20,20)
//     for( int i = 0; i < nPoints; ++i ) {
//         double y = i-nPoints/2;
//         for( int j = 0; j < nPoints; ++j ) {
//             double x = j-nPoints/2;
//             if(x*x+y*y < nPoints/2) {
//                 peak(i,j) = 1;
//             }
// //             peak(i,j) *= 7;
// //             peak(i,j) += (y+0.5*x)/800;
// //             peak(i,j) += sin(x/50)*cos(y/200);
//         }
//     }


Array<double> peak2(peak,nPoints,nPoints/2);

    FourierTransform ft(peak,NORMALIZE_FT);
    //ft.reorder();
    FourierTransform ft2(peak,FULLCOMPLEX|NORMALIZE_FT);
    //ft2.reorder();
//    ft.convolveInPlace(peak);
//    ft.inv(peak);
//    ft2.convolveInPlace(peak);
//    ft2.inv(peak);

    //auto bla1 = ft.convolve(peak);
    
    redux::util::Array<double> tmp;
    peak.copy(tmp);
    //ft.convolveInPlace( tmp );
    //FourierTransform ft3(ft2,ft.dimensions());
    //Array<complex_t> ft3( reinterpret_cast<Array<complex_t>&>(ft2), ft.dimensions() );
    //FourierTransform ft4 = ft;
    //FourierTransform ft5 = ft2;
    //int bla=1;
    //for(auto & it: ft) it = bla++;
    //ft *= complex_t(3,0);
    //Array<complex_t> ft3;
    //ft.copy(ft3);
    //ft3 = complex_t(0.1,0);

    //for(auto & it: ft3) it = bla++;
    //ft3 *= ft;
    //ft *= ft4;
    //redux::file::Ana::write( "fthh.f0", ft );
    //ft2 *= ft5;
    //redux::file::Ana::write( "ftff.f0", ft2 );
    //ft = ft4;
    //ft2 = ft5;
    //ft4.copy(ft);
    //ft5.copy(ft2);
    ft *= ft2;

    //ft2 *= ft4;
    //redux::file::Ana::write( "fthf.f0", ft2 );


return;
    auto bla2 = ft2.convolve(peak);
    //ft.autocorrelate();
    //ft2.autocorrelate();
    
return;

    ft.ift(peak);
    ft2.ift(peak);
   */

}

void backward_fft( void ) {
    
    backward_FFT<double>( 64 );
    //backward_FFT<float>( 53, 64 );

    
}

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
                BOOST_TEST( abs(distPtr[y][x] - (double)sqrt(yDist*yDist+xDist*xDist)) < 0.0001 );
                if( yDist || xDist ) {
                    BOOST_TEST( abs(anglePtr[y][x] - (double)atan2(yDist, xDist)) < 0.0001 );
                } else {
                    BOOST_TEST( anglePtr[y][x] == 0.0 );
                }
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




namespace testsuite {

    namespace image {

        void add_tests( test_suite* ts ) {

            srand(time(NULL));
            
            // FourierTransform
            ts->add( BOOST_TEST_CASE( &fft_methods ) );
            ts->add( BOOST_TEST_CASE( &fft_constructors ) );
            ts->add( BOOST_TEST_CASE( &fft_manipulations ) );
            ts->add( BOOST_TEST_CASE( &forward_fft ) );
            ts->add( BOOST_TEST_CASE( &backward_fft ) );
            
            // Grid, Pupil & Modes
            ts->add( BOOST_TEST_CASE( &grid_test ) );

        }
    }
}
