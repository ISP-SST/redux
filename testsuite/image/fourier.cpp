#include "redux/image/fouriertransform.hpp"
#include "redux/math/functions.hpp"

#include "catch.hpp"

#include <array>
#include <cstdlib>
#include <ctime>
#include <numeric>

using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;


namespace testsuite {

    namespace image {

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
                CHECK( empty.nDimensions() == 0 );
                CHECK( empty.nElements() == 0 );
                CHECK( empty.getInputSize() == PointI(0,0) );
                CHECK( empty.getFtSize() == PointI(0,0) );
                CHECK( empty.getInputPixels() == 0 );
                CHECK( empty.getFtPixels() == 0 );
                CHECK( empty.getThreads() == 1 );
                CHECK( empty.get() == nullptr);
            }

            {       // Construct by size, but without supplying data
                FourierTransform ft( nY, nX, flags, nThreads );
                CHECK( ft.nDimensions() == 2 );
                CHECK( ft.nElements() == nY*ft_sizeX );
                CHECK( ft.getInputSize() == PointI(nY,nX) );
                CHECK( ft.getFtSize() == PointI(nY,ft_sizeX) );
                CHECK( ft.getInputPixels() == nX*nY );
                CHECK( ft.getFtPixels() == nY*ft_sizeX );
                CHECK( ft.getThreads() == nThreads );
                CHECK( ft.get() != nullptr);
            }

            {       // Construct from input
                FourierTransform ft1( input, flags, nThreads );
                FourierTransform ft2( input.get(), nY, nX, flags, nThreads );   // These two constructors should yield identical results.
                PointI expectedFTsize = inputSize;
                if( !(ft1.getFlags()&FULLCOMPLEX) ) {
                    expectedFTsize.x = expectedFTsize.x/2+1;
                }
                size_t expectedElements = expectedFTsize.y*expectedFTsize.x;
                CHECK( ft1.nDimensions() == 2 );
                CHECK( ft1.nElements() == expectedElements );
                CHECK( ft1.getInputSize() == inputSize );
                CHECK( ft1.getFtSize() == expectedFTsize );
                CHECK( ft1.getInputPixels() == nX*nY );
                CHECK( ft1.getFtPixels() == expectedElements );
                CHECK( ft1.getThreads() == nThreads );
                CHECK( ft1.get() != nullptr);

                CHECK( ft2.nDimensions() == 2 );
                CHECK( ft2.nElements() == expectedElements );
                CHECK( ft2.getInputSize() == inputSize );
                CHECK( ft2.getFtSize() == expectedFTsize );
                CHECK( ft2.getInputPixels() == nX*nY );
                CHECK( ft2.getFtPixels() == expectedElements );
                CHECK( ft2.getThreads() == nThreads );
                CHECK( ft2.get() != nullptr);

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
                        CHECK( abs(ft1(y,x) - tmpC(y,x)) < tol );     // check that FT matches expectations
                        CHECK( abs(ft2(y,x) - tmpC(y,x)) < tol );     // check that FT matches expectations
                        //CHECK( (ft1(y,x)) == (tmpC(y,x)) );     // check that FT matches expectations
                        //CHECK( (ft2(y,x)) == (tmpC(y,x)) );     // check that FT matches expectations
                    }
                }

                // This constructor should throw if input is not 2D.
                tmpD.resize(nX);
                CHECK_THROWS_AS( FourierTransform( tmpD, flags, nThreads ), logic_error );
                tmpD.resize(nX,nY,nX);
                CHECK_THROWS_AS( FourierTransform( tmpD, flags, nThreads ), logic_error );
            }

            {       // Copy-construct
                FourierTransform ft1( nY, nX, flags, nThreads );    // The sanity of ft1 was tested above.
                FourierTransform ft2( ft1 );
                CHECK( ft2.nDimensions() == ft1.nDimensions() );
                CHECK( ft2.nElements() == ft1.nElements() );
                CHECK( ft2.getInputSize() == ft1.getInputSize() );
                CHECK( ft2.getFtSize() == ft1.getFtSize() );
                CHECK( ft2.getInputPixels() == ft1.getInputPixels() );
                CHECK( ft2.getFtPixels() == ft1.getFtPixels() );
                CHECK( ft2.getThreads() == ft1.getThreads() );
                CHECK( ft2.get() != ft1.get() );       // should not share data
                CHECK( ft2.get() != nullptr );         // ...but should have it's own
            }

        }

        template <typename T> void fft_operators( size_t nX, size_t nY=0 ) {

            if( nX == 0 ) throw runtime_error("fft_operators can not be called with size=0");
            if( nY == 0 ) nY = nX;

            srand(time(NULL));

            size_t halfSizeX = nX/2;
            Array<T> input( nY, nX );
            Array<double> tmpD( nY, nX );
            Array<complex_t> tmpC( nY, nX );

            {   // *= operator  (no Fourier-transforming yet, just basic operations on the FourierTransform class
                FourierTransform fullFT( nY, nX, FULLCOMPLEX );    // full-complex
                CHECK( fullFT.dimSize(0) == nY );
                CHECK( fullFT.dimSize(1) == nX );

                FourierTransform halfFT( nY, nX );                   // half-complex
                CHECK( halfFT.dimSize(0) == nY );
                CHECK( halfFT.dimSize(1) == halfSizeX+1 );

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
                        CHECK( norm(halfFT(y,x)-val) < tol );
                        CHECK( norm(fullFT(y,x)-val) < tol );
                        if( x ) CHECK( norm(fullFT(y,nX-x)-val) < tol );
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
                        CHECK( norm(halfFT(y,x) - val*fullFTv ) < tol );
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
                        CHECK( norm(fullFT(y,x) - val*halfFTv ) < tol );
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
                        CHECK( norm(fullFT(y,x)-tmpC(y,x)) < tol );                         // check that fullFT matches expectations
                        if( x < halfSizeX ) CHECK( norm(halfFT(y,x)-tmpC(y,x)) < tol );    // check that halfFT matches expectations
                    }
                }

                // The same, but with "REORDER" flag
                fullFT.init( input, FULLCOMPLEX|REORDER_FT );
                tmpC.zero();                                // Set to expected results, i.e. a delta-function
                tmpC(nY/2,nX/2) = (nX*nY);             // DC component is stored in (nY/2,nX/2), and is not normalized by default, i.e. should be equal to (nX*nY)
                for( size_t y(0); y<nY; ++y) {
                    for( size_t x(0); x<nX; ++x ) {
                        CHECK( norm(fullFT(y,x)-tmpC(y,x)) < tol );                         // check that fullFT matches expectations
                    }
                }

                // The same, but with "NORMALIZE" flag
                fullFT.init( input, FULLCOMPLEX|NORMALIZE_FT );
                halfFT.init( input, NORMALIZE_FT );
                tmpC.zero();                                // Set to expected results, i.e. a delta-function
                tmpC(0,0) = 1;                              // DC component is stored in (nY/2,halfSizeX), when normalized it should be 1.0
                for( size_t y(0); y<nY; ++y) {
                    for( size_t x(0); x<nX; ++x ) {
                        CHECK( norm(fullFT(y,x)-tmpC(y,x)) < tol );                         // check that fullFT matches expectations
                        if( x < halfSizeX ) CHECK( norm(halfFT(y,x)-tmpC(y,x)) < tol );    // check that halfFT matches expectations
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
                        CHECK( norm(fullFT(y,x)-tmpC(y,x)) < tol );                         // check that fullFT matches expectations
                        if( x < halfSizeX ) CHECK( norm(halfFT(y,x)-tmpC(y,x)) < tol );     // check that halfFT matches expectations
                    }
                }

                // The same, but with the delta-function centered in inout and using the "REORDER_IMG" flag instead
                input.zero();
                input( nY/2, halfSizeX ) = 1;
                fullFT.init( input, FULLCOMPLEX|REORDER_IMG );
                halfFT.init( input, REORDER_IMG );
                for( size_t y(0); y<nY; ++y) {
                    for( size_t x(0); x<nX; ++x ) {
                        CHECK( norm(fullFT(y,x)-tmpC(y,x)) < tol );                         // check that fullFT matches expectations
                        if( x < halfSizeX ) CHECK( norm(halfFT(y,x)-tmpC(y,x)) < tol );    // check that halfFT matches expectations
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
                        CHECK( norm(fullFT(y,x)-tmpC(y,x)) < 1E-4 );  // check that fullFT matches expectations
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
                        CHECK( abs(tmpD1(y,x)-input(y,x)) < tol );     // check that we recover input function
                        CHECK( abs(tmpD2(y,x)-input(y,x)) < tol );
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
                        CHECK( abs(tmpD1(y,x)-input(y,x)) < tol );     // check that we recover input function
                        CHECK( abs(tmpD2(y,x)-input(y,x)) < tol );
                    }
                }

                fullFT.init( input, FULLCOMPLEX|REORDER_FT );               // same test, but reordered FT
                tmpD1.zero();
                fullFT.ift( tmpD1.get() );
                for( size_t y(0); y<nY; ++y) {
                    for( size_t x(0); x<nX; ++x ) {
                        CHECK( abs(tmpD1(y,x)-input(y,x)) < tol );     // check that we recover input function
                    }
                }

                fullFT.init( input, FULLCOMPLEX, REORDER_IMG );             // same test, but reorder img before FT
                halfFT.init( input, REORDER_IMG );
                tmpD1.zero();
                tmpD2.zero();
                fullFT.ift( tmpD1.get() );
                halfFT.ift( tmpD2.get() );
                for( size_t y(0); y<nY; ++y) {
                    for( size_t x(0); x<nX; ++x ) {
                        CHECK( abs(tmpD1(y,x)-input(y,x)) < tol );     // check that we recover input function
                        CHECK( abs(tmpD2(y,x)-input(y,x)) < tol );
                    }
                }

            }

        }

        void checkClose2D( const Array<double>& actual, const Array<double>& expected, double tolerance ) {
            size_t nY = actual.dimSize(0);
            size_t nX = actual.dimSize(1);
            for( size_t y(0); y<nY; ++y ) {
                for( size_t x(0); x<nX; ++x ) {
                    CHECK( abs(actual(y,x)-expected(y,x)) < tolerance );
                }
            }
        }

        void checkConvolve( FourierTransform& fft, const Array<double>& psf, const Array<double>& data, const Array<double>& expected ) {
            Array<double> tmp = psf.copy<double>();
            Array<double> tmp2 = data.copy<double>();
            fft.convolve( tmp.get(), tmp2.get() );       // naked data-pointers
            checkClose2D( tmp2, expected, tol );

            tmp = psf.copy<double>();
            tmp2 = data.copy<double>();
            fft.convolve( tmp, tmp2 );                   // Array<T> for in/out
            checkClose2D( tmp2, expected, tol );

            tmp = psf.copy<double>();
            tmp2 = fft.convolve( tmp );                  // Array<T> for in, returning out
            checkClose2D( tmp2, expected, tol );

            tmp = psf.copy<double>();
            fft.convolveInPlace( tmp );                  // in-place, Array<T>
            checkClose2D( tmp, expected, tol );
        }

        // Same 4 call-styles, for correlate()/correlateInPlace().
        void checkCorrelate( FourierTransform& fft, const Array<double>& img2, const Array<double>& img1, const Array<double>& expected ) {
            Array<double> tmp = img2.copy<double>();
            Array<double> tmp2 = img1.copy<double>();
            fft.correlate( tmp.get(), tmp2.get() );
            checkClose2D( tmp2, expected, tol );

            tmp = img2.copy<double>();
            tmp2 = img1.copy<double>();
            fft.correlate( tmp, tmp2 );
            checkClose2D( tmp2, expected, tol );

            tmp = img2.copy<double>();
            tmp2 = fft.correlate( tmp );
            checkClose2D( tmp2, expected, tol );

            tmp = img2.copy<double>();
            fft.correlateInPlace( tmp );
            checkClose2D( tmp, expected, tol );
        }

        void checkAutocorrelate( FourierTransform& fft, const Array<double>& input, const Array<double>& expected, bool center ) {
            Array<double> tmp = input.copy<double>();
            fft.autocorrelate( tmp.get(), center );
            checkClose2D( tmp, expected, tol );

            tmp = input.copy<double>();
            fft.autocorrelate( tmp, center );
            checkClose2D( tmp, expected, tol );

            tmp = input.copy<double>();
            fft.autocorrelate( input.get(), tmp.get(), center );
            checkClose2D( tmp, expected, tol );
        }

        }  // namespace


        TEST_CASE( "FT Constructors", "[image]" ) {
            for( int f(0); f<9; ++f ) {
                for( int t(1); t<5; t+=2 ) {
                    fft_construct<double>( 64, 0, f, t );
                    fft_construct<float>( 53, 0, f, t );
                    fft_construct<complex_t>( 71, 53, f, t );
                    fft_construct<int32_t>( 53, 88, f, t );
                }
            }
        }

        TEST_CASE( "FT reorder", "[image]" ) {

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
            CHECK( tmpV == expected );

            FourierTransform::reorderInto( tmpV.data(), 4, 4, tmpA.data() );
            for( size_t i(0); i<16; ++i ) {
                CHECK( tmpA[i] == input[i] );
            }

            // 5x5 input
            FourierTransform::reorderInto( input.data(), 5, 5, tmpA.data() );
            expected = { 19, 20, 16, 17, 18,        // input:  1  2  3  4  5
                        24, 25, 21, 22, 23,        //         6  7  8  9 10
                        4,  5,  1,  2,  3,        //        11 12 13 14 15
                        9, 10,  6,  7,  8,        //        16 17 18 19 20
                        14, 15, 11, 12, 13 };      //        21 22 23 24 25
            tmpV.assign( tmpA.begin(), tmpA.begin()+25 );
            CHECK( tmpV == expected );

            // test reordering into a larger target
            tmpA.fill(0);
            FourierTransform::reorderInto( input.data(), 4, 4, tmpA.data(), 5, 6 );
            expected = { 11, 12,  0,  0,  9, 10,        // input:  1  2  3  4
                        15, 16,  0,  0, 13, 14,        //         5  6  7  8
                        0,  0,  0,  0,  0,  0,        //         9 10 11 12
                        3,  4,  0,  0,  1,  2,        //        13 14 15 16
                        7,  8,  0,  0,  5,  6 };
            tmpV.assign( tmpA.begin(), tmpA.begin()+30 );
            CHECK( tmpV == expected );

            // same, but wih center-option
            tmpA.fill(0);
            FourierTransform::reorderInto( input.data(), 4, 4, tmpA.data(), 5, 6, true );
            expected = { 0, 11, 12,  9, 10,  0,        // input:  1  2  3  4
                        0, 15, 16, 13, 14,  0,        //         5  6  7  8
                        0,  3,  4,  1,  2,  0,        //         9 10 11 12
                        0,  7,  8,  5,  6,  0,        //        13 14 15 16
                        0,  0,  0,  0,  0,  0 };
            tmpV.assign( tmpA.begin(), tmpA.begin()+30 );
            CHECK( tmpV == expected );

            // reorderInto should throw when called with nullptr
            CHECK_THROWS_AS( FourierTransform::reorderInto( input.data(), 4, 4, (int*)nullptr ), logic_error );
            CHECK_THROWS_AS( FourierTransform::reorderInto( (int*)nullptr, 4, 4, tmpA.data() ), logic_error );

            // reorderInto should throw when inSize > outSize
            CHECK_THROWS_AS( FourierTransform::reorderInto( input.data(), 4, 4, tmpA.data(), 3, 4 ), logic_error );
            CHECK_THROWS_AS( FourierTransform::reorderInto( input.data(), 4, 4, tmpA.data(), 4, 3 ), logic_error );

            // reorderInto should throw when inPtr == outPtr
            CHECK_THROWS_AS( FourierTransform::reorderInto( tmpA.data(), 4, 4, tmpA.data(), 5, 5 ), logic_error );

        }

        TEST_CASE( "FT center/decenter", "[image]" ) {

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

            CHECK( (ft == ft_uncentered) );    // trivial check, since we just made a copy

            // elementary checks, FT should be un-centered and full-complex.
            CHECK( ft.isFullComplex() );
            CHECK( !ft.isCentered() );
            ft.decenter();                          // should do nothing (not centered)
            CHECK( !ft.isCentered() );
            ft.center();                            // center/reorder the FT
            CHECK( ft.isCentered() );
            CHECK( (ft == ft_centered) );      // compare data
            ft.center();                            // should do nothing (already centered)
            CHECK( ft.isCentered() );
            ft.decenter();                          // decenter/reorder the FT
            CHECK( !ft.isCentered() );

            // do similar tests, but initalize FT to centered.
            FourierTransform ft2( peak, FULLCOMPLEX|REORDER_FT );
            ft_uncentered = FourierTransform(ft2);  // make copies, for comparisons
            ft_centered = FourierTransform(ft2);
            FourierTransform::reorder( ft_uncentered );

            CHECK( (ft2 == ft_centered) );    // trivial check, since we just made a copy

            // elementary checks, FT should be centered and full-complex.
            CHECK( ft2.isFullComplex() );
            CHECK( ft2.isCentered() );

            ft2.center();                           // should do nothing (already centered)
            CHECK( ft2.isCentered() );
            ft2.decenter();                         // decenter/reorder the FT
            CHECK( !ft2.isCentered() );
            CHECK( (ft == ft_uncentered) );    // compare data
            ft2.decenter();                         // should do nothing (not centered)
            CHECK( !ft2.isCentered() );
            ft2.center();                           // center/reorder the FT
            CHECK( ft2.isCentered() );

            Array<double> bla( nY, nX );
            bla.zero();
            ft.init( peak, FULLCOMPLEX|REORDER_FT );

            ft.ift( bla );
            for( size_t y(0); y<nY; ++y) {
                for( size_t x(0); x<nX; ++x ) {
                    //CHECK( (bla(y,x)-peak(y,x)) < (tol) );
                }
            }

        }

        TEST_CASE( "FT conjugate and norm", "[image]" ) {

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
                    CHECK( abs(fullFT(y,x) - std::conj(ft_ref(y,x))) < tol );
                    CHECK( abs(norm(fullFT(y,x)) - norm(ft_ref(y,x))) < tol );
                    if( x <= nX/2 ) {
                        CHECK( abs(halfFT(y,x) - std::conj(ft_ref(y,x))) < tol );
                        CHECK( abs(norm(halfFT(y,x)) - norm(ft_ref(y,x))) < tol );
                    }
                }
            }

        }

        TEST_CASE( "FT power", "[image]" ) {

            size_t nX = 64;
            size_t nY = 64;

            Array<double> tmp( nY, nX );
            tmp.zero();

            // Arbitrary function, we'll use a gaussian with fwhm (nY/2,nX/4) centered at (mid+2,mid+4)
            redux::math::gauss( tmp.get(), nY, nX, nY/2, nX/4, nY/2+2, nX/2+4 );

            FourierTransform fullFT( tmp, FULLCOMPLEX );

            Array<float> pwr = fullFT.power();

            // TODO: test for expected result (never implemented in the original Boost test either)

        }

        TEST_CASE( "FT convolution", "[image]" ) {

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
            checkConvolve( fullFT, psf, data, expected );

            fullFT.init( data, FULLCOMPLEX|REORDER_FT|NORMALIZE_FT );    // test full-complex + centered + normalized
            checkConvolve( fullFT, psf, data, expected );

            fullFT.init( data, FULLCOMPLEX|REORDER_FT|REORDER_IMG );    // test full-complex + centered + re-ordered input
            checkConvolve( fullFT, psf, data, expected );

            FourierTransform halfFT( data );            // test half-complex
            checkConvolve( fullFT, psf, data, expected );

            halfFT.init( data, NORMALIZE_FT|REORDER_IMG );            // test half-complex + normalized +  re-ordered input
            checkConvolve( fullFT, psf, data, expected );

        }

        TEST_CASE( "FT correlation", "[image]" ) {

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
            checkCorrelate( fullFT, img2, img1, expected );

            FourierTransform halfFT( img1 );             // test half-complex
            checkCorrelate( fullFT, img2, img1, expected );

        }

        TEST_CASE( "FT autocorrelation", "[image]" ) {

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

            checkAutocorrelate( fullFT, input, expected, false );

            FourierTransform::reorder( expected.get(), nY, nX );        // test with center-flag
            checkAutocorrelate( fullFT, input, expected, true );

        }

        TEST_CASE( "FT Forward", "[image]" ) {
            forward_FFT<double>( 64, 64 );
            forward_FFT<double>( 63, 64 );
            forward_FFT<float>( 53, 64 );
            forward_FFT<int32_t>( 53, 64 );
        }

        TEST_CASE( "FT Reverse", "[image]" ) {
            backward_FFT<double>( 64 );
            // TODO: backward_FFT<float>(53, 64) was already commented out in the original
        }

        TEST_CASE( "FT Manipulations", "[image]" ) {
        }

    }

}
