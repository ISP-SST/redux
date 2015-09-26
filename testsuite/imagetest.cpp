#include <boost/test/unit_test.hpp>


#include "redux/file/fileana.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/image.hpp"
#include "redux/image/grid.hpp"
#include "redux/image/utils.hpp"
#include "redux/math/functions.hpp"
#include "redux/util/stringutil.hpp"

#include <iostream>


using namespace redux::file;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

using namespace boost::unit_test_framework;

namespace {
     

}


void fourierTest( void ) {
    
    size_t nPoints = 10;
    size_t halfSize = nPoints/2;
    Array<double> input(nPoints,nPoints);
    Array<double> result(nPoints,nPoints);
    Array<double> expected(nPoints,nPoints);
    double tolerance = 1E-9;
    srand(time(NULL));
    {   // *= operator
        FourierTransform fullFT(nPoints,nPoints,FT_FULLCOMPLEX);    // full-complex
        BOOST_CHECK( fullFT.dimSize(0) == nPoints );
        BOOST_CHECK( fullFT.dimSize(1) == nPoints );

        FourierTransform halfFT(nPoints,nPoints);                   // half-complex
        BOOST_CHECK( halfFT.dimSize(0) == nPoints );
        BOOST_CHECK( halfFT.dimSize(1) == halfSize+1 );

        complex_t salt(rand()%100,rand()%100);
        // test *= operator for full-complex to full-complex AND half-complex to half-complex
        for( uint y=0; y<nPoints; ++y) { 
            for( uint x=0; x<halfSize+1 ; ++x) {
                complex_t val(y,x);
                val += salt;
                halfFT(y,x) = fullFT(y,x) = val;
                if( x ) fullFT(y,nPoints-x) = val;
            }
        }
        fullFT *= fullFT;
        halfFT *= halfFT;
        for( uint y=0; y<nPoints; ++y) { 
            for( uint x=0; x<halfSize+1; ++x) {
                complex_t val(y,x);
                val += salt;
                val *= val;
                BOOST_CHECK_SMALL( norm(halfFT(y,x)-val), tolerance );
                BOOST_CHECK_SMALL( norm(fullFT(y,x)-val), tolerance );
                if( x ) BOOST_CHECK_SMALL( norm(fullFT(y,nPoints-x)-val), tolerance );
            }
        }

        // test *= operator for full-complex to half-complex
        for( uint y=0; y<nPoints; ++y) { 
            for( uint x=0; x<halfSize+1 ; ++x) {
                complex_t val(y,x);
                val += salt;
                halfFT(y,x) = fullFT(y,x) = val;
                if( x ) fullFT(y,nPoints-x) = val;
            }
        }
        halfFT *= fullFT;
        for( uint y=0; y<nPoints; ++y) { 
            for( uint x=0; x<halfSize+1; ++x) {
                complex_t val(y,x);
                val += salt;
                val *= val;
                BOOST_CHECK_SMALL( norm(halfFT(y,x)-val), tolerance );
            }
        }

        // test *= operator for half-complex to full-complex
        for( uint y=0; y<nPoints; ++y) { 
            for( uint x=0; x<halfSize+1 ; ++x) {
                complex_t val(y,x);
                val += salt;
                halfFT(y,x) = fullFT(y,x) = val;
                if( x ) fullFT(y,nPoints-x) = val;
            }
        }
        fullFT *= halfFT;
        for( uint y=0; y<nPoints; ++y) { 
            for( uint x=0; x<halfSize+1; ++x) {
                complex_t val(y,x);
                val += salt;
                val *= val;
                BOOST_CHECK_SMALL( norm(fullFT(y,x)-val), tolerance );
                if( x ) BOOST_CHECK_SMALL( norm(fullFT(y,nPoints-x)-val), tolerance );
            }
        }

    }
    
    {   // test some simple transformas
        double salt(rand()%100+1);
        //input.zero();
        //input(0,0) = 1;      // a delta function
        //input(halfSize,halfSize) = 1;      // a delta function
        input = salt;             //  Set the array to some arbitrary constant value
        {                       //  The transform of a constant will only have 1 non-zero value, located at (0,0)
            FourierTransform fullFT(input,FT_FULLCOMPLEX|FT_NORMALIZE);
            FourierTransform halfFT(input,FT_NORMALIZE);
            complex_t val(salt,0);                  // expected value of zero-frequency component
            BOOST_CHECK_SMALL( norm(fullFT(0,0)-val), tolerance );
            BOOST_CHECK_SMALL( norm(halfFT(0,0)-val), tolerance );
            for( uint y=0; y<nPoints; ++y) { 
                for( uint x=0; x<halfSize+1; ++x) {
                    if( x || y ) BOOST_CHECK_SMALL( norm(halfFT(y,x)), tolerance );         // check all except (0,0)
                    if( x || y ) BOOST_CHECK_SMALL( norm(fullFT(y,x)), tolerance );         // check all except (0,0)
                    if( x ) BOOST_CHECK_SMALL( norm(fullFT(y,nPoints-x)), tolerance );
                }
            }
            fullFT.reorder();                       // check that the zero-frequency value is reordered to the right location
            for( uint y=0; y<nPoints; ++y) { 
                for( uint x=0; x<nPoints; ++x) {
                    if( (x == halfSize) && (y == halfSize) ) BOOST_CHECK_SMALL( norm(fullFT(y,x)-val), tolerance );
                    else BOOST_CHECK_SMALL( norm(fullFT(y,x)), tolerance );
                }
            }
            fullFT.reorder();
            
            // get inverse and compare
            halfFT.inv(result);
            BOOST_CHECK( result == input );
            // get inverse (of FT) and compare
            fullFT.inv(result);
            BOOST_CHECK( result == input );

            input.zero();
            input(0,0) = 3;      // a delta function
            expected.zero();
            expected(halfSize,halfSize) = 9;        // centered delta with expected value
            fullFT.reset(input,FT_FULLCOMPLEX);
            halfFT.reset(input);
            //result(nPoints/2,nPoints/2) = 1;      // a delta function
            input.copy(result);     // convolveInPlace is destructive, use temporary
            fullFT.convolveInPlace(result);
            BOOST_CHECK( result == expected );
            input.copy(result);     // convolveInPlace is destructive, use temporary
            halfFT.convolveInPlace(result);
            //BOOST_CHECK( result == expected );
            for( uint y=0; y<nPoints; ++y) { 
                for( uint x=0; x<nPoints; ++x) {
                    BOOST_CHECK_SMALL( result(y,x)-expected(y,x), tolerance );
                }
            }
           
            input.zero();
            input(0,0) = input(0,1) = input(1,0) = input(1,1) = 1;      // 2x2
            input.copy(result);     // autocorelate is destructive, use temporary
            expected.zero();
            expected(halfSize,halfSize) = 4;        // create expected result from autocorrelation
            expected(halfSize+1,halfSize) = expected(halfSize,halfSize+1) = expected(halfSize-1,halfSize) = expected(halfSize,halfSize-1) = 2;
            expected(halfSize+1,halfSize+1) = expected(halfSize-1,halfSize+1) = expected(halfSize-1,halfSize-1) = expected(halfSize+1,halfSize-1) = 1;
            FourierTransform::autocorrelate(result);
            for( uint y=0; y<nPoints; ++y) { 
                for( uint x=0; x<nPoints; ++x) {
                    BOOST_CHECK_SMALL( result(y,x)-expected(y,x), tolerance );
                }
            }

/*
 
            fullFT.convolveInPlace(tmp2,FT_FULLCOMPLEX|FT_NORMALIZE);
            redux::file::Ana::write( "tmp.f0", tmp );
            redux::file::Ana::write( "tmp2.f0", tmp2 );
            BOOST_CHECK( tmp2 == tmp );
            */
        }
        
        return;
        
        {                       //  The same but as normalized FTs
            FourierTransform fullFT(input,FT_FULLCOMPLEX|FT_NORMALIZE);
            FourierTransform halfFT(input,FT_NORMALIZE);
            complex_t val(salt,0);          // geometry factor normalized out, only the mean-value remains
            BOOST_CHECK_EQUAL( fullFT(0,0), val );
            BOOST_CHECK_EQUAL( halfFT(0,0), val );
            for( uint y=0; y<nPoints; ++y) { 
                for( uint x=0; x<nPoints/2+1; ++x) {
                    if( x || y ) BOOST_CHECK_EQUAL( fullFT(y,x), complex_t(0,0) );         // check all except (0,0)
                    if( x ) BOOST_CHECK_EQUAL( fullFT(y,nPoints-x), complex_t(0,0) );
                }
            }
            fullFT.reorder();
            for( uint y=0; y<nPoints; ++y) { 
                for( uint x=0; x<nPoints; ++x) {
                    if( (x == nPoints/2) && (y == nPoints/2) ) BOOST_CHECK_EQUAL( fullFT(y,x), val );
                    else  BOOST_CHECK_EQUAL( fullFT(y,x), complex_t(0,0) );
                }
            }
            Array<double> tmp2;
            fullFT.inv(tmp2);
            BOOST_CHECK( tmp2 == input );
            halfFT.inv(tmp2);
            BOOST_CHECK( tmp2 == input );
        }
        {
//             Array<double> tmp2(nPoints,nPoints);
//             tmp2.zero();
//             tmp2(0,0) = 1;      // a delta function
//             
        }
    }
return;
    
    Array<double> peak(nPoints,nPoints);
    peak.zero();
    redux::math::gauss( peak.get(), nPoints, nPoints, nPoints/20, nPoints/40, nPoints/2, nPoints/2);    // centered gaussian with fwhm (20,20)
/*    for( int i = 0; i < nPoints; ++i ) {
        double y = i-nPoints/2;
        for( int j = 0; j < nPoints; ++j ) {
            double x = j-nPoints/2;
            if(x*x+y*y < nPoints/2) {
                peak(i,j) = 1;
            }
//             peak(i,j) *= 7;
//             peak(i,j) += (y+0.5*x)/800;
//             peak(i,j) += sin(x/50)*cos(y/200);
        }
    }
*/

Array<double> peak2(peak,nPoints,nPoints/2);

    FourierTransform ft(peak,FT_NORMALIZE);
    //ft.reorder();
    FourierTransform ft2(peak,FT_FULLCOMPLEX|FT_NORMALIZE);
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

    ft.inv(peak);
    ft2.inv(peak);
   

}


void gridTest( void ) {
    
    // Grids
    // TODO: Better testcase, this is basically the same math as int the Grid-constructor.

    int64_t nPoints = 50;
    float originX = 11;
    float originY = 11.5;
    Grid grid = Grid::get(nPoints,originY,originX);
    float** distPtr = grid.distance.get();
    float** anglePtr = grid.angle.get();
    for( int i = 0; i < nPoints; ++i ) {
        double yDist = i - originY;
        for( int j = 0; j < nPoints; ++j ) {
            double xDist = j - originX;
            BOOST_CHECK_CLOSE( distPtr[i][j], (float)sqrt(yDist*yDist+xDist*xDist), 1E-10);
            if(yDist || xDist) {
                BOOST_CHECK_CLOSE( anglePtr[i][j], (float)atan2(yDist, xDist), 1E-5 );
            } else  BOOST_CHECK_EQUAL( anglePtr[i][j], 0.0 );
        }
    }

    
}




namespace testsuite {

    namespace image {

        void imageTest( void ) {

            test_suite* ts = BOOST_TEST_SUITE( "IMAGE" );

            ts->add( BOOST_TEST_CASE( &fourierTest ) );
            ts->add( BOOST_TEST_CASE( &gridTest ) );

            framework::master_test_suite().add( ts );

        }

    }

}
