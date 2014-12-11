#include <boost/test/unit_test.hpp>


#include "redux/file/fileana.hpp"
#include "redux/image/image.hpp"
#include "redux/image/statistics.hpp"
#include "redux/image/utils.hpp"
#include "redux/math/functions.hpp"
#include "redux/math/helpers.hpp"
#include "redux/util/stringutil.hpp"

#include <iostream>

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
        
        std::vector<size_t> subFirst=dims, subLast=dims;
        for( size_t i=0; i<nDims; ++i) {
            subFirst[i] *= 0.25;
            subLast[i]  *= 0.75;
        }
        
        Array<T> array( 2, s... );
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
                BOOST_CHECK_EQUAL( it, 13 );
            }
            
            typename Array<T>::iterator ait = array.begin();
            typename Array<T>::iterator sit = subarray.begin();
            size_t count1(0),count2(0);
            while( ait < array.end() ) {
                count1++;
                if( ait.pos() < sit.pos() || sit == subarray.end()) {   // outside subarray, should have value = 2
                    BOOST_CHECK_EQUAL( *ait++, 2 );
                } else {                        // inside subarray, should have value = 13
                    BOOST_CHECK_EQUAL( *ait++, 13 );
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
            BOOST_CHECK( tmp == tmp2 );
            
            tmp2.resize();                                  // pack/unpack empty array;
            buf = sharedArray<char>( tmp2.size() );
            ptr = buf.get();
            count = tmp2.pack(ptr);
            BOOST_CHECK_EQUAL( count, tmp2.size() );
            
            count = tmp.unpack(ptr,false);
            BOOST_CHECK_EQUAL( count, tmp2.size() );
            BOOST_CHECK_EQUAL( count, tmp.size() );
            BOOST_CHECK( tmp == tmp2 );

            
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
        
        Statistics stats;
        stats.getMinMaxMean(array);
        T mn,mx;
        double avg;
        redux::math::minMaxMean(rawPtr, array.nElements(), mn, mx, avg);
        BOOST_CHECK_CLOSE( stats.min, mn, 1E-20);
        BOOST_CHECK_CLOSE( stats.max, mx, 1E-20);
        BOOST_CHECK_CLOSE( stats.mean, avg, 1E-20);
        
        stats.getRmsStddev(array,avg);
        double rms,stddev;
        redux::math::rmsStddev(rawPtr,array.nElements(),rms,stddev);
        BOOST_CHECK_CLOSE( stats.rms, rms, 1E-20);
        BOOST_CHECK_CLOSE( stats.stddev, stddev, 1E-20);
        
    }

}

void arrayTests( void ) {
    
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

        redux::file::Ana::write( "stats.f0", peak );
        
        using namespace std;
        Statistics stats;
        stats.getMinMaxMean(peak);
        cout << std::setprecision(9) << "Min = " << stats.min << "  Max = " << stats.max << "   Mean = " << stats.mean << endl;
        stats.getRmsStddev(peak,stats.mean);
        cout << "StdDev = " << stats.stddev << "  RMS = " << stats.rms << endl;
        
        int mask = 10; //10;
        double limit = -1; //130;
        
        FourierTransform ft1(peak.copy<complex_t>());
        redux::file::Ana::write( "pow1.f0", ft1.power() );
        double noise1 = ft1.noise(mask,limit);
        FourierTransform ft2(peak);
        redux::file::Ana::write( "pow2.f0", ft2.power() );
        double noise2 = ft2.noise(mask,limit);
        cout << "Noise1 = " << noise1 << "  Noise2 = " << noise2 << endl;
        
        //peak -= stats.mean;
        double std_mvn = std_dev(peak,nPoints);
        cout << "StdDev_mvn = " << std_mvn << endl;
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
    std::shared_ptr<int**> sharedPtr = array.get(sizeZ,sizeY,sizeX);
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
            ts->add( BOOST_TEST_CASE( &transformTest ) );

            framework::master_test_suite().add( ts );

        }

    }

}
