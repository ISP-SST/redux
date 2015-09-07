#include <boost/test/unit_test.hpp>


#include "redux/file/fileana.hpp"
#include "redux/image/image.hpp"
#include "redux/util/arraystats.hpp"
#include "redux/image/utils.hpp"
#include "redux/math/functions.hpp"
#include "redux/math/helpers.hpp"
#include "redux/util/stringutil.hpp"

#include <iostream>


using namespace redux::file;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

using namespace boost::unit_test_framework;

namespace {
        
    double std_dev(const Array<double>& pic,int np)  //   standard deviation
    {
      double s=0.0,ss=0.0;
      for(int x=0;x<np;++x)
        for(int y=0;y<np;++y){
            double tmp = pic(x,y);
          ss += tmp*tmp;
          s += tmp;
        }
      s/=(double)(np*np);
      ss/=(double)(np*np);
      return sqrt(ss-s*s);
    }
    
    template <typename T, typename ...S>
    void testArray(S ...s) {
        
        std::vector<size_t> dims({static_cast<size_t>(s)...});
        size_t nDims = dims.size();
        if ( nDims == 0 ) return;
        
        std::vector<size_t> subFirst=dims, subLast=dims;    // dimensions of a centered & half-sized subarray
        for( size_t i=0; i<nDims; ++i) {
            subFirst[i] *= 0.25;
            subLast[i]  *= 0.75;
        }
        
        Array<T> array( 2, s... );      // test if constructor works if a new dimension is prepended
        BOOST_CHECK_EQUAL( array.nDimensions(), nDims+1 );
        BOOST_CHECK_EQUAL( array.dimSize(0), 2 );
        for( size_t i=0; i<nDims; ++i) {
            BOOST_CHECK_EQUAL( array.dimSize(i+1), dims[i] );
        }
        
        array.resize( dims );
        size_t nElements = array.nElements();
        BOOST_CHECK_EQUAL( array.nDimensions(), nDims );
        for( size_t i=0; i<nDims; ++i) {
            BOOST_CHECK_EQUAL( array.dimSize(i), dims[i] );
        }
        
        const std::vector<size_t>& strides = array.strides();
        size_t stride = 1;
        for( size_t i=nDims-1; i>0; --i) {
            BOOST_CHECK_EQUAL( strides[i], stride );
            stride *= dims[i];
        }
        
        {
            // test comparison and assignment
            Array<T> array2(s...);
            array2 = 3;      // set all values to 3
            array =  3;      // set all values to 3
            BOOST_CHECK_EQUAL( array2.nDimensions(), nDims );
            for( size_t i=0; i<nDims; ++i) {
                BOOST_CHECK_EQUAL( array.dimSize(i), dims[i] );
            }
            BOOST_CHECK( array2 == array );
            *array2.ptr() = 123;      // 1 wrong value
            BOOST_CHECK( !( array2 == array ) );

            // test that arrays of different dimensions compare as different.
            array2.resize(2,s...);
            BOOST_CHECK_EQUAL( array2.nDimensions(), nDims+1 );
            array2 = 3;      // all values equal
            BOOST_CHECK( !( array2 == array ) );
        }
        
        // test operators with scalars
        array = 1;
        for( auto it : array ) { BOOST_CHECK_EQUAL( it, 1 ); }
        array += 10;
        for( auto it : array ) { BOOST_CHECK_EQUAL( it, 11 ); }
        array -= 1;
        for( auto it : array ) { BOOST_CHECK_EQUAL( it, 10 ); }
        array *= 2;
        for( auto it : array ) { BOOST_CHECK_EQUAL( it, 20 ); }
        array /= 5;
        for( auto it : array ) { BOOST_CHECK_EQUAL( it, 4 ); }
        array.zero();
        for( auto it : array ) { BOOST_CHECK_EQUAL( it, 0 ); }

        /***** test iterators *****/
        {

            Array<T> subarray( array, subFirst, subLast );
            BOOST_CHECK_EQUAL( subarray.nDimensions(), nDims );
            for( size_t i=0; i<nDims; ++i) {
                BOOST_CHECK_EQUAL( array.dimSize(i), dims[i] );
            }
            T* rawPtr = array.get();
            for( uint x = 0; x < nElements; ++x ) {
                rawPtr[x] = x;          // set values equal to real offsets: 0,1,2,3...
            }

            typename Array<T>::iterator it = array.begin();
            typename Array<T>::const_iterator cit = array.begin();
            typename Array<T>::iterator sit = subarray.begin();
            
            // check that subarray spans from subFirst to subLast
            BOOST_CHECK_EQUAL( array.pos(subFirst).pos(), sit.pos() );
            BOOST_CHECK_EQUAL( array.pos(subLast).pos(), (--subarray.end()).pos() );
            
            // check that whole datablok is accessing the right element
            for( uint i = 0; i < nElements; ++i ) {
                BOOST_CHECK_EQUAL( it.pos(), i );
                BOOST_CHECK_EQUAL( cit.pos(), i );
                BOOST_CHECK_EQUAL( *it, i );
                BOOST_CHECK_EQUAL( *cit, i );
                it++,cit++;
            }
           
            it = array.begin();
            cit = array.begin();
            
            // postfix ++
            for( uint i = 0; i < nElements; ++i ) {
                BOOST_CHECK_EQUAL( *it++, i );
                BOOST_CHECK_EQUAL( *cit++, i );
            }

            for( uint i = 0; i<subarray.nElements(); ++i ) {
                sit++;
            }

            // at end()
            BOOST_CHECK( it == array.end() );
            BOOST_CHECK( cit == array.end() );
            BOOST_CHECK( sit == subarray.end() );
            
            // prefix --
            for( int i = nElements-1; i >= 0; --i ) {
                BOOST_CHECK_EQUAL( *--it, i );
                BOOST_CHECK_EQUAL( *--cit, i );
            }
            for( uint i = 0; i<subarray.nElements(); ++i ) {
                --sit;
            }

            // at begin()
            BOOST_CHECK( it == array.begin() );
            BOOST_CHECK( cit == array.begin() );
            BOOST_CHECK( sit == subarray.begin() );

            // prefix ++
            for( uint i = 0; i < nElements; ++i ) {
                BOOST_CHECK_EQUAL( *it, i );
                BOOST_CHECK_EQUAL( *cit, i );
                ++it;
                ++cit;
            }
            for( int i = 0; i<subarray.nElements(); ++i ) {
                ++sit;
            }


            // at end()
            BOOST_CHECK( it == array.end() );
            BOOST_CHECK( cit == array.end() );
            BOOST_CHECK( sit == subarray.end() );


            // postfix --
            for( int i = nElements-1; i >= 0; --i ) {
                it--;
                cit--;
                BOOST_CHECK_EQUAL( *it, i );
                BOOST_CHECK_EQUAL( *cit, i );
            }
            for( int i = 0; i<subarray.nElements(); ++i ) {
                sit--;
            }

            // at begin()
            BOOST_CHECK( it == array.begin() );
            BOOST_CHECK( cit == array.begin() );
            BOOST_CHECK( sit == subarray.begin() );

            // set/check values via auto.
            int cnt = 0;
            for( auto & ait : array ) {
                ait = ++cnt;
            }

            cnt = 0;
            for( auto ait : array ) {
                BOOST_CHECK_EQUAL( ait, ++cnt );
            }

        }
        /**************************/
        
        {
            Array<T> subarray( array, subFirst, subLast );

            // compare pointers to verify that the data is shared
            BOOST_CHECK( subarray.get() == array.get() );
            
            // check if subarrays are dense
            BOOST_CHECK( array.dense() );
            if(nDims == 1) {
                BOOST_CHECK( subarray.dense() );        // 1D arrays are always dense
            } else {
                BOOST_CHECK( !subarray.dense() );
            }

            // scalar assign to subarray
            array = 2;
            for( auto it : array ) {
                BOOST_CHECK_EQUAL( it, 2 );
            }
            subarray = 13;
            for( auto it : subarray ) {
                //cout << "  it = " << it << endl;
                BOOST_CHECK_EQUAL( (int)it, 13 );
            }
            
            typename Array<T>::iterator ait = array.begin();
            typename Array<T>::iterator sit = subarray.begin();
            size_t count1(0),count2(0);
            while( ait < array.end() ) {
                count1++;
                if( ait.pos() < sit.pos() || sit == subarray.end()) {   // outside subarray, should have value = 2
                    BOOST_CHECK_EQUAL( (int)*ait++, 2 );
                } else {                        // inside subarray, should have value = 13
                    BOOST_CHECK_EQUAL( (int)*ait++, 13 );
                    if(sit < subarray.end()) BOOST_CHECK_EQUAL( *sit++, 13 );
                    count2++;
                }
            }
            BOOST_CHECK_EQUAL( count1, nElements );
            BOOST_CHECK_EQUAL( count2, subarray.nElements() );
            
            // elementwise modification
            array = 2;
            for( auto it : array ) {
                BOOST_CHECK_EQUAL( it, 2 );
            }

            for( auto it : subarray ) {
                BOOST_CHECK_EQUAL( it, 2 );
            }
            
            ait = array.begin();
            sit = subarray.begin();
            while( sit < subarray.end() ) {
                *sit *= 2;
                BOOST_CHECK_EQUAL( *(ait+sit.pos()), 4 );   // check value changed
                sit++;
            }
            
            // compare pointers to verify that the data is still shared
            BOOST_CHECK( subarray.get() == array.get() );

            // shift subarray
            vector<int64_t> tmpV(nDims);
            for( int i = 0; i < nDims; ++i ) {
                tmpV[i] = subFirst[i];
            }
            for( int i = 0; i < nDims; ++i ) {
                Array<T> subarray2 = subarray;
                vector<int64_t> tmp = tmpV;
                BOOST_CHECK_EQUAL( subarray2.shift(i,-100000), -tmp[i]);                // can only shift to the edge
                tmp[i] = 0;
                BOOST_CHECK_EQUAL( subarray2.begin().pos(), array.pos(tmp).pos() );
            }
            for( int i = 0; i < nDims; ++i ) {
                tmpV[i] = subLast[i];
            }
            for( int i = 0; i < nDims; ++i ) {
                Array<T> subarray2 = subarray;
                vector<int64_t> tmp = tmpV;
                BOOST_CHECK_EQUAL( subarray2.shift(i,100000), dims[i]-tmp[i]-1);        // can only shift to the edge
                tmp[i] = dims[i]-1;
                BOOST_CHECK_EQUAL( (--subarray2.end()).pos(), array.pos(tmp).pos() );
            }
          
        }
        
        {   // test pack/unpack
            
            auto buf = sharedArray<char>( array.size() );
            char* ptr = buf.get();
            uint64_t count = array.pack(ptr);
            BOOST_CHECK_EQUAL( count, array.size() );
            
            Array<T> tmp;
            count = tmp.unpack(ptr,false);
            BOOST_CHECK_EQUAL( count, array.size() );
            BOOST_CHECK_EQUAL( count, tmp.size() );
            BOOST_CHECK( tmp == array );
            
            Array<T> tmp2(array,subFirst,subLast);          // pack/unpack subarray;
            buf = sharedArray<char>( tmp2.size() );
            ptr = buf.get();
            count = tmp2.pack(ptr);
            BOOST_CHECK_EQUAL( count, tmp2.size() );
            
            count = tmp.unpack(ptr,false);
            BOOST_CHECK_EQUAL( count, tmp2.size() );
            BOOST_CHECK_EQUAL( count, tmp.size() );
            //BOOST_CHECK( tmp == tmp2 );  FIXME
            
            tmp2.resize();                                  // pack/unpack empty array;
            buf = sharedArray<char>( tmp2.size() );
            ptr = buf.get();
            count = tmp2.pack(ptr);
            BOOST_CHECK_EQUAL( count, tmp2.size() );
            
            count = tmp.unpack(ptr,false);
            BOOST_CHECK_EQUAL( count, tmp2.size() );
            BOOST_CHECK_EQUAL( count, tmp.size() );
            BOOST_CHECK( tmp == tmp2 );

            Array<T> tmp3;
            buf = sharedArray<char>( tmp3.size() );
            ptr = buf.get();
            count = tmp3.pack(ptr);
            BOOST_CHECK_EQUAL( count, tmp3.size() );
            
            count = tmp.unpack(ptr,false);
            BOOST_CHECK_EQUAL( count, tmp3.size() );
            BOOST_CHECK_EQUAL( count, tmp.size() );
            BOOST_CHECK( tmp == tmp3 );
        }
        
    }
    
    
    template <typename T, typename ...S>
    void testStat(S ...s) {
        
        std::vector<size_t> dims({static_cast<size_t>(s)...});
        size_t nDims = dims.size();
        if ( nDims == 0 ) return;
        
        std::vector<size_t> subFirst=dims, subLast=dims;
        for( size_t i=0; i<nDims; ++i) {
            subFirst[i] *= 0.25;
            subLast[i]  *= 0.75;
        }
        
        Array<T> array( s... );
        T* rawPtr = array.get();
        for( int x = 0; x < array.nElements(); ++x ) {
            rawPtr[x] = x;          // set values equal to real offsets: 0,1,2,3...
        }
        
        ArrayStats stats;
        stats.getMinMaxMean(array);
        T mn,mx;
        double avg;
        redux::math::minMaxMean(rawPtr, array.nElements(), mn, mx, avg);
        BOOST_CHECK_CLOSE( stats.min, mn, 1E-20);
        BOOST_CHECK_CLOSE( stats.max, mx, 1E-20);
        BOOST_CHECK_CLOSE( stats.mean, avg, 1E-20);
        
        stats.getRmsStddev(array);
        double rms,stddev;
        redux::math::rmsStddev(rawPtr,array.nElements(),rms,stddev);
        BOOST_CHECK_CLOSE( stats.rms, rms, 1E-20);
        BOOST_CHECK_CLOSE( stats.stddev, stddev, 1E-20);
        
    }

}


void arrayTests( void ) {
    
    // test different dimensions and datatypes
    testArray<int8_t>(100);
    testArray<uint16_t>(44,55);
    testArray<int>(4,14,9);
    testArray<size_t>(4,7,6,5);
    testArray<double>(4,5,6,7,8);
    testArray<float>(100,100,100);
    testArray<uint32_t>(13,1,9);
    
//    testStat<int>(100,103); FIXME
    
}


void statTest( void ) {
    

    int64_t nPoints = 800;
    Array<double> peak(nPoints,nPoints);
    redux::math::gauss( peak.get(), nPoints, nPoints, 200, 500, nPoints/2, nPoints/2);    // centered gaussian with fwhm (50,50)
    for( int i = 0; i < nPoints; ++i ) {
        double y = i;
        for( int j = 0; j < nPoints; ++j ) {
            double x = j;
            peak(i,j) *= 7;
            peak(i,j) += (y+0.5*x)/800;
            peak(i,j) += sin(x/50)*cos(y/200);
        }
    }
    
    using namespace std;
    ArrayStats stats;
    stats.getMinMaxMean(peak);
    //cout << std::setprecision(9) << "Min = " << stats.min << "  Max = " << stats.max << "   Mean = " << stats.mean << endl;
    stats.getRmsStddev(peak);
    //cout << "StdDev = " << stats.stddev << "  RMS = " << stats.rms << endl;
    
    int mask = 10; //10;
    double limit = -1; //130;
    
    FourierTransform ft1(peak.copy<complex_t>());
    double noise1 = ft1.noise(mask,limit);
    FourierTransform ft2(peak);
    double noise2 = ft2.noise(mask,limit);
    //cout << "Noise1 = " << noise1 << "  Noise2 = " << noise2 << endl;
    
    //peak -= stats.mean;
    double std_mvn = std_dev(peak,nPoints);
    //cout << "StdDev_mvn = " << std_mvn << endl;


}



void gridTest( void ) {
    
    // Grids
    // TODO: Better testcase, this is basically the same math as int the Grid-constructor.

    int64_t nPoints = 50;
    float originX = 11;
    float originY = 11.5;
    Grid grid(nPoints,originY,originX);
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
//             redux::file::Ana::write( "input.f0", input );
//             redux::file::Ana::write( "result.f0", result );
//             redux::file::Ana::write( "expected.f0", expected );
/*
datadir='/home/tomas/build/redux/'
fullft=f0(datadir+'fullft.f0')
halfft=f0(datadir+'halfft.f0')
input=f0(datadir+'input.f0')
result=f0(datadir+'result.f0')
expected=f0(datadir+'expected.f0')
;idlco=convol(tmp0,tmp0)
;diff = result-expected
;print,min(diff),max(diff),mean(diff)
print,input
print,result
print,expected

inft=f0('inft.f0')
tmp0=f0('tmp0.f0')
idlft=fft(tmp0)

 
            fullFT.convolveInPlace(tmp2,FT_FULLCOMPLEX|FT_NORMALIZE);
            redux::file::Ana::write( "tmp.f0", tmp );
            redux::file::Ana::write( "tmp2.f0", tmp2 );
            BOOST_CHECK( tmp2 == tmp );*/
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

/*

halfft=f0('halfft.f0')
fullft=f0('fullft.f0')
tmp=f0('tmp.f0')
tmp2=f0('tmp2.f0')
inft=f0('inft.f0')

tvscl,inft
tvscl,tmp2,105,0
;tvscl,prod,310,0

print,min(tmp),max(tmp),mean(tmp)
print,min(tmp2),max(tmp2),mean(tmp2)

;print,real_part(halfft)
;print,real_part(fullft)
;print,real_part(prod)

peak=f0('peak.f0')
ft=f0('ft.f0')
fc=f0('fc.f0')
ft3=f0('ft3.f0')
ft4=f0('ft4.f0')
fthh=f0('fthh.f0')
ftff=f0('ftff.f0')
fthf=f0('fthf.f0')
ftfh=f0('ftfh.f0')
idlft = FFT(peak)  
ftinv=f0('ftinv.f0')
fcinv=f0('fcinv.f0')
ftci=f0('ftci.f0')
fcci=f0('fcci.f0')
ftac=f0('ftac.f0')
fcac=f0('fcac.f0')
ftacinv=f0('ftacinv.f0')
fcacinv=f0('fcacinv.f0')

;tvscl,peak
tvscl,ft*ft
tvscl,fc*fc,105,0
tvscl,fthh,310,0
tvscl,ftff,415,0
tvscl,ftfh,620,0
tvscl,fthf,825,0
print,min(ft*ft),max(ft*ft),mean(ft*ft)
print,min(ftfh),max(ftfh),mean(ftfh)

print,min(ftff),max(ftff),mean(ftff)
print,min(fthf),max(fthf),mean(fthf)
print,ft(99:100,0:3)
print,ft3(99:103,0:3)
print,fc(99:103,0:3)

tvscl,ftinv,0,200
tvscl,fcinv,200,200
tvscl,ftci,400,200
tvscl,fcci,600,200

tvscl,ftac,0,400
tvscl,fcac,200,400
tvscl,ftacinv,400,400
tvscl,fcacinv,600,400

print,min(peak),max(peak)
print,min(ft),max(ft)
print,min(fc),max(fc)
print,min(ft-fc(0:100,*)),max(ft-fc(0:100,*))
print,min(ftac),max(ftac)
print,min(fcac),max(fcac)
print,min(ftac-fcac(0:100,*)),max(ftac-fcac(0:100,*))
print,min(fcci-fcac(0:100,*)),max(fcci-fcac(0:100,*))

   rhs #100 = (-2.70823e-08,-5.55051e-16)   oValue #100 = (-2.70823e-08,-5.55051e-16)   nValue #100 = (7.33452e-16,3.00641e-23)

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
    
    
    return;
    testArray<int8_t>(100);
    testArray<uint16_t>(44,55);
    testArray<int>(4,14,9);
    testArray<size_t>(4,7,6,5);
    testArray<double>(4,5,6,7,8);
    testArray<float>(100,100,100);
    testArray<uint32_t>(13,1,9);
    
    testStat<int>(100,103);

return;    
    // Grids
    // TODO: Better testcase, this is basically the same math as int the Grid-constructor.
    {
        int64_t nPoints = 50;
        float originX = 11;
        float originY = 11.5;
        Grid grid(nPoints,originY,originX);
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

    // Stats
    {
        int64_t nPoints = 800;
        Array<double> peak(nPoints,nPoints);
        redux::math::gauss( peak.get(), nPoints, nPoints, 200, 500, nPoints/2, nPoints/2);    // centered gaussian with fwhm (50,50)
        for( int i = 0; i < nPoints; ++i ) {
            double y = i;
            for( int j = 0; j < nPoints; ++j ) {
                double x = j;
                peak(i,j) *= 7;
                peak(i,j) += (y+0.5*x)/800;
                peak(i,j) += sin(x/50)*cos(y/200);
            }
        }
        
        using namespace std;
        ArrayStats stats;
        stats.getMinMaxMean(peak);
        //cout << std::setprecision(9) << "Min = " << stats.min << "  Max = " << stats.max << "   Mean = " << stats.mean << endl;
        stats.getRmsStddev(peak);
        //cout << "StdDev = " << stats.stddev << "  RMS = " << stats.rms << endl;
        
        int mask = 10; //10;
        double limit = -1; //130;
        
        FourierTransform ft1(peak.copy<complex_t>());
        double noise1 = ft1.noise(mask,limit);
        FourierTransform ft2(peak);
        double noise2 = ft2.noise(mask,limit);
        //cout << "Noise1 = " << noise1 << "  Noise2 = " << noise2 << endl;
        
        //peak -= stats.mean;
        double std_mvn = std_dev(peak,nPoints);
        //cout << "StdDev_mvn = " << std_mvn << endl;
    }

}


void transformTest( void ) {
    
    int64_t sizeX = 5;
    int64_t sizeY = 4;
    int64_t sizeZ = 3;
    Array<int> array( sizeZ, sizeY, sizeX );
    
    int cnt = 0;
    for(auto& it: array) it = ++cnt;
  
    // check values
    for( int i = 0; i < sizeZ; ++i ) {
        for( int j = 0; j < sizeY; ++j ) {
            for( int k = 0; k < sizeX; ++k ) {
                BOOST_CHECK_EQUAL( array( i, j, k ), sizeX*sizeY*i + j*sizeX + k + 1 );
            }
        }
    }
    
    Array<int> array2 = array.copy();           // make a copy for comparison
    std::shared_ptr<int**> sharedPtr = array.reshape(sizeZ,sizeY,sizeX);
    int*** rawPtr = sharedPtr.get();

    // check reverseX
    for( int i = 0; i < sizeZ; ++i ) {
        reverseX(rawPtr[i],sizeY,sizeX);
        for( int j = 0; j < sizeY; ++j ) {
            for( int k = 0; k < sizeX; ++k ) {
                BOOST_CHECK_EQUAL( array( i, j, k ), array2( i, j, sizeX-k-1) );
            }
        }
        reverseX(rawPtr[i],sizeY,sizeX);        // flip back
    }

    // check reverseY
    for( int i = 0; i < sizeZ; ++i ) {
        reverseY(rawPtr[i],sizeY,sizeX);
        for( int j = 0; j < sizeY; ++j ) {
            for( int k = 0; k < sizeX; ++k ) {
                BOOST_CHECK_EQUAL( array( i, j, k ), array2( i, sizeY-j-1, k) );
            }
        }
        reverseY(rawPtr[i],sizeY,sizeX);        // flip back
    }

    // check transpose
    array.permuteDimensions(1,2);
    for( int i = 0; i < sizeZ; ++i ) {
        transpose(*rawPtr[i],sizeY,sizeX);
        for( int j = 0; j < sizeY; ++j ) {
            for( int k = 0; k < sizeX; ++k ) {
                BOOST_CHECK_EQUAL( array(i,k,j), array2(i,j,k) );
            }
        }
        transpose(*rawPtr[i],sizeX,sizeY);        // transpose back
    }
    
    BOOST_CHECK( array == array2 );             // Just to verify they are still the same.
    
}



namespace testsuite {

    namespace image {

        void imageTest( void ) {

            test_suite* ts = BOOST_TEST_SUITE( "IMAGE" );

            ts->add( BOOST_TEST_CASE( &arrayTests ) );
            ts->add( BOOST_TEST_CASE( &statTest ) );
            ts->add( BOOST_TEST_CASE( &fourierTest ) );
            //ts->add( BOOST_TEST_CASE( &transformTest ) );

            framework::master_test_suite().add( ts );

        }

    }

}
