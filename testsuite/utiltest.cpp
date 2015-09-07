#include <boost/test/unit_test.hpp>


#include "redux/constants.hpp"
#include "redux/types.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/math/helpers.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arraystats.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/bitoperations.hpp"
#include "redux/util/boundvalue.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <iostream>
#include <memory>

using namespace redux::image;
using namespace redux::math;
using namespace redux::util;
using namespace redux;
using namespace std;
using namespace boost::unit_test_framework;

namespace {

    uint64_t naiveCount( uint64_t v ) {
        uint64_t cnt( 0 );
        do {
            if( v & 1 ) cnt++;
        }
        while( v >>= 1 );
        return cnt;
    }

       
    double std_dev(const Array<double>& pic,int np) { //   standard deviation
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
            for( uint i = 0; i<subarray.nElements(); ++i ) {
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
            for( uint i = 0; i<subarray.nElements(); ++i ) {
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
            for( uint i = 0; i < nDims; ++i ) {
                tmpV[i] = subFirst[i];
            }
            for( uint i = 0; i < nDims; ++i ) {
                Array<T> subarray2 = subarray;
                vector<int64_t> tmp = tmpV;
                BOOST_CHECK_EQUAL( subarray2.shift(i,-100000), -tmp[i]);                // can only shift to the edge
                tmp[i] = 0;
                BOOST_CHECK_EQUAL( subarray2.begin().pos(), array.pos(tmp).pos() );
            }
            for( uint i = 0; i < nDims; ++i ) {
                tmpV[i] = subLast[i];
            }
            for( uint i = 0; i < nDims; ++i ) {
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
    
    template <typename T>
    void transformTest( void ) {
        
        int64_t sizeX = 5;
        int64_t sizeY = 4;
        int64_t sizeZ = 3;
        Array<T> array( sizeZ, sizeY, sizeX );
        
        T cnt(0);
        for(auto& it: array) it = (cnt+=1);
    
        // check values
        for( int i = 0; i < sizeZ; ++i ) {
            for( int j = 0; j < sizeY; ++j ) {
                for( int k = 0; k < sizeX; ++k ) {
                    BOOST_CHECK_EQUAL( array( i, j, k ), sizeX*sizeY*i + j*sizeX + k + 1 );
                }
            }
        }
        
        Array<T> array2 = array.copy();           // make a copy for comparison
        std::shared_ptr<T**> sharedPtr = array.reshape(sizeZ,sizeY,sizeX);
        T*** rawPtr = sharedPtr.get();

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


void arrayTest( void ) {
    
    Array<int> array4x5( 4, 5 );
    Array<int> array3x3;

    // check that assign throws for nValues > nElements;
    BOOST_CHECK_THROW( array4x5.assign( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 ), logic_error );

    array4x5.assign( 1,  2,  3,  4,  5.0,   // all values will be cast to int in assign
                     6,  7,  8,  9, 10,
                    11, 12, 13, 14, 15,
                    16, 17, 18, 19, 20 );

    // check the iterator is accessing the right elements.
    Array<int>::iterator it = array4x5.begin();
    BOOST_CHECK_EQUAL( it.pos(), 0 );
    it = array4x5.end();
    BOOST_CHECK_EQUAL( it.pos(), 20 );
    BOOST_CHECK_EQUAL( *--it, 20 );
    it = array4x5.pos(1,1);
    BOOST_CHECK_EQUAL( it.pos(), 6 );
    BOOST_CHECK_EQUAL( it.step(1).pos(), 7 );
    BOOST_CHECK_EQUAL( it.step(-1).pos(), 6 );
    BOOST_CHECK_EQUAL( it.step(1,1).pos(), 12 );
    BOOST_CHECK_EQUAL( it.step(-1,-1).pos(), 6 );
    BOOST_CHECK_EQUAL( it.step(1,-1).pos(), 10 );
    BOOST_CHECK_EQUAL( it.step(-1,1).pos(), 6 );
   
    // check index out of bounds when using at()
    BOOST_CHECK_THROW( array4x5.at( 1, 5 ), out_of_range );
    BOOST_CHECK_THROW( array4x5.at( 5, 1 ), out_of_range );
    BOOST_CHECK_THROW( array4x5.at( -1, 1 ), out_of_range );
    BOOST_CHECK_THROW( array4x5.at( 1, -1 ), out_of_range );
    
    // test iterator returned by Array::pos() for some values
    BOOST_CHECK_EQUAL( *array4x5.pos( 0, 0 ), 1 );
    BOOST_CHECK_EQUAL( *array4x5.pos( 1, 1 ), 7 );
    BOOST_CHECK_EQUAL( *array4x5.pos( 3, 4 ), 20 );

   
    // check values
    for( int i = 0; i < 4; ++i ) {
        for( int j = 0; j < 5; ++j ) {
            BOOST_CHECK_EQUAL( array4x5( i, j ), 5 * i + j + 1 );
        }
    }

    // test iterators
    it = array4x5.begin();
    Array<int>::iterator endit = array4x5.end();
    Array<int>::const_iterator cit = array4x5.begin();
    Array<int>::const_iterator cendit = array4x5.end();

    // prefix ++
    for( int i = 1; i <= 20; ++i ) {
        BOOST_CHECK_EQUAL( *it, i );
        BOOST_CHECK_EQUAL( *cit, i );
        ++it;
        ++cit;
    }

    // at end()
    BOOST_CHECK( it == endit );
    BOOST_CHECK( cit == cendit );


    // prefix --
    for( int i = 20; i > 0; --i ) {
        BOOST_CHECK_EQUAL( *--it, i );
        BOOST_CHECK_EQUAL( *--cit, i );
    }

    // postfix ++
    for( int i = 1; i <= 20; ++i ) {
        BOOST_CHECK_EQUAL( *it++, i );
        BOOST_CHECK_EQUAL( *cit++, i );
    }

    // at end()
    BOOST_CHECK( it == endit );
    BOOST_CHECK( cit == cendit );

    // postfix --
    for( int i = 20; i > 0; --i ) {
        it--;
        cit--;
        BOOST_CHECK_EQUAL( *it, i );
        BOOST_CHECK_EQUAL( *cit, i );
    }

    // test with std auto iterator (by reference)
    int i = 1;
    for( auto & ait : array4x5 ) {
        BOOST_CHECK_EQUAL( ait, i++ );
    }

    // shallow copy (shared data)
    array3x3 = array4x5;
    BOOST_CHECK( array3x3 == array4x5 );
    // verify that data is shared
    BOOST_CHECK( array3x3.ptr() == array4x5.ptr() );

    // deep copy
    array3x3 = array4x5.copy();
    BOOST_CHECK( array3x3 == array4x5 );
    // verify that data is not shared
    BOOST_CHECK( array3x3.ptr() != array4x5.ptr() );

    // test sub-array constructor
    Array<int> subarray( array4x5, 1, 3, 1, 3 ); // 3x3 sub-array, with offset (1,1)
    
    BOOST_CHECK( array4x5.dense() );
    BOOST_CHECK( !subarray.dense() );

    // verify that data is shared
    BOOST_CHECK( subarray.ptr() == array4x5.ptr() );
    // define an array equal to the center 3x3 elements of the 5x5 array above for reference
    array3x3.resize( 3, 3 );
    array3x3.assign( 7,  8,  9,
                    12, 13, 14,
                    17, 18, 19 );
    // verify that data matches
    BOOST_CHECK( subarray == array3x3 );


    // test subiterator::step() to step along any dimension
    BOOST_CHECK_EQUAL( *(subarray.pos( 1, 1 ).step()), 14 );
    BOOST_CHECK_EQUAL( *(subarray.pos( 1, 1 ).step(1)), 14 );
    BOOST_CHECK_EQUAL( *(subarray.pos( 1, 1 ).step(1,-1)), 17 );
    BOOST_CHECK_EQUAL( *(subarray.pos( 1, 1 ).step(-1,-1)), 7 );
    BOOST_CHECK_EQUAL( *(subarray.pos( 1, 1 ).step(-1,1)), 9 );

    // test assigning to sub-array
    array3x3.assign( 70,  80,  90,
                    120, 130, 140,
                    170, 180, 190 );
    subarray.assign( array3x3 );

    // array4x5 should now be:
    //     1,   2,   3,   4,  5,
    //     6,  70,  80,  90, 10,
    //    11, 120, 130, 140, 15,
    //    16, 170, 180, 190, 20 );
    //
    // check values
    for( int i = 0; i < 4; ++i ) {
        for( int j = 0; j < 5; ++j ) {
            int val = 5 * i + j + 1;
            if( i >= 1 && i <= 3 && j >= 1 && j <= 3 ) val *= 10;
            BOOST_CHECK_EQUAL( array4x5( i, j ), val );
        }
    }

    Array<int> array4x5x6( 4, 5, 6 );
    size_t cnt( 0 );
    for( auto & it : array4x5x6 ) {
        it = ++cnt;
    }
    BOOST_CHECK_EQUAL( cnt, 120 );
    
    {
        // arbitrary subarray.
        Array<int> subarray( array4x5x6, 1, 2, 2, 3, 3, 4 );
        
        // check that the subarray is accessing the right elements.
        for( size_t i = 0; i < 2; ++i ) {
            for( size_t j = 0; j < 2; ++j ) {
                for( size_t k = 0; k < 2; ++k ) {
                    BOOST_CHECK_EQUAL( array4x5x6( i + 1, j + 2, k + 3 ), subarray( i, j, k ) );
                }
            }
        }

        // assign to subarray
        subarray = 999;
        cnt = 0;
        for( auto cit=subarray.begin(); cit != subarray.end(); ++cit ) {
            BOOST_CHECK_EQUAL( *cit, 999 );
            cnt++;
        }
        BOOST_CHECK_EQUAL( cnt, 8 );
        
        cnt = 0;
        for( auto cit=subarray.end(); cit != subarray.begin();  ) {
            BOOST_CHECK_EQUAL( *--cit, 999 );
            cnt++;
        }
        BOOST_CHECK_EQUAL( cnt, 8 );

    }
    
    {
        // test access as raw multidimensional array
        cnt=0;
        for( auto & it : array4x5x6 ) {
            it = ++cnt;
        }
        shared_ptr<int**> raiiArray = array4x5x6.reshape(4,5,6);
        int*** rawArray = raiiArray.get();
        cnt = 0;
        for( size_t i=0; i<4; ++i) {
            for( size_t j=0; j<5; ++j) {
                for( size_t k=0; k<6; ++k) {
                    BOOST_CHECK_EQUAL( rawArray[i][j][k], ++cnt );
                }
            }
        }
    }
    
    // test different dimensions and datatypes
    testArray<int8_t>(100);
    testArray<uint16_t>(44,55);
    testArray<int>(4,14,9);
    testArray<size_t>(4,7,6,5);
    testArray<double>(4,5,6,7,8);
    testArray<float>(100,100,100);
    testArray<uint32_t>(13,1,9);
    
    transformTest<int8_t>();
    transformTest<int>();
    transformTest<double>();
    transformTest<size_t>();
    
}



void bitTest( void ) {

    // deBruijn based log2
    for( uint8_t i( 1 ); i < 32; ++i ) {
        BOOST_CHECK_EQUAL( redux::util::log2( ( 1 << i ) + 1 ), i );
    }

    // deBruijn LSB detector
    for( uint8_t i( 0 ); i < 32; ++i ) {
        BOOST_CHECK_EQUAL( findLSB( 1 << i ), i );
    }

    // deBruijn LSB detector, 1-based index
    for( uint8_t i( 0 ); i < 32; ++i ) {
        BOOST_CHECK_EQUAL( findLSB1( 1 << i ), i + 1 );
    }

    // 64-bit deBruijn LSB detector, 1-based index
    uint64_t tmp64( 1 );
    for( uint8_t i( 0 ); i < 64; ++i ) {
        BOOST_CHECK_EQUAL( findLSB64( tmp64 << i ), i + 1 );
    }

    // next power of two
    for( uint8_t i( 0 ); i < 31; ++i ) {
        BOOST_CHECK_EQUAL( nextPowerOfTwo( ( 1 << i ) + 1 ), ( 1 << ( i + 1 ) ) );
    }

    uint32_t mask( ( 1 << 16 ) - 1 );
    // swap first 2 bytes of a.i[0] with the first 2 bytes of a.i[1]
    union {
        uint64_t l;
        uint32_t i[2];
        uint16_t s[4];
        uint8_t b[8];
    } a, b;
    for( uint16_t i( 0 ); i < 10; i++ ) {
        a.i[0] = b.i[0] = rand();
        a.i[1] = b.i[1] = rand();
        swapBits( a.i[0], a.i[1], mask );
        BOOST_CHECK_EQUAL( a.s[0], b.s[2] ); // verify
        BOOST_CHECK_EQUAL( a.s[2], b.s[0] );
        BOOST_CHECK_EQUAL( a.s[1], b.s[1] ); // check that the unmasked bits are unchanged.
        BOOST_CHECK_EQUAL( a.s[3], b.s[3] );
    }

    // swap first 2 bytes of each integer in array with matching bytes in another array
    uint32_t aa[10], bb[10], cc[10];
    memset( bb, 0, 10 * sizeof( uint32_t ) );
    for( uint16_t i( 0 ); i < 10; i++ ) {
        aa[i] = cc[i] = rand();
    }
    swapBits( aa, bb, mask, 10 );
    for( uint16_t i( 0 ); i < 10; i++ ) {
        BOOST_CHECK_EQUAL( cc[i] & ( mask << 16 ), aa[i] );
        BOOST_CHECK_EQUAL( cc[i]&mask, bb[i] );
    }

    memcpy( aa, cc, 10 * sizeof( uint32_t ) );
    // flip the first 16 bits in aa[10]
    flipBits( aa, mask, 10 );
    for( uint16_t i( 0 ); i < 10; i++ ) {
        BOOST_CHECK_EQUAL( cc[i] ^ mask, aa[i] );
    }

    for( int i( 0 ); i < 15; ++i ) {
        a.i[0] = rand();
        a.i[1] = rand();
        //BOOST_CHECK_EQUAL( countBits(a.b[0]), naiveCount(b.s[0]) );
        //BOOST_CHECK_EQUAL( countBits(a.s[0]), naiveCount(a.s[0]) );
        BOOST_CHECK_EQUAL( countBits( a.i[0] ), naiveCount( a.i[0] ) );
        BOOST_CHECK_EQUAL( countBits( a.l ), naiveCount( a.l ) );
    }


}


void boundValueTest( void ) {

    typedef BoundValue<int, detail::UNDEFINEDTRIM> invalidBV;
    BOOST_CHECK_THROW( invalidBV(), invalid_argument );

    BoundValue<int> intVal( 3, -9, 9 );                 // ordinary construction / assignment
    BOOST_CHECK_EQUAL( intVal, 3 );
    BOOST_CHECK_EQUAL( intVal.min(), -9 );
    BOOST_CHECK_EQUAL( intVal.max(),  9 );
    intVal = 14;
    BOOST_CHECK_EQUAL( intVal, 9 );
    intVal = -14.0;                                     // test implicit arithmetic cast.
    BOOST_CHECK_EQUAL( intVal, -9 );

    BoundValue<double> dblVal( 3, 9, -9 );              // reverse min/max order
    BOOST_CHECK_EQUAL( dblVal, 3.0 );
    BOOST_CHECK_EQUAL( dblVal.min(), -9.0 );
    BOOST_CHECK_EQUAL( dblVal.max(),  9.0 );
    dblVal = 14;
    BOOST_CHECK_EQUAL( dblVal, 9.0 );
    dblVal = ( int ) - 14;                              // test implicit arithmetic cast.
    BOOST_CHECK_EQUAL( dblVal, -9.0 );

    BoundValue<double> dblVal2( dblVal );               // copy construction of the same type
    BOOST_CHECK_EQUAL( dblVal2, -9.0 );
    BOOST_CHECK_EQUAL( dblVal2.min(), -9.0 );
    BOOST_CHECK_EQUAL( dblVal2.max(),  9.0 );

    BoundValue<int> intVal2( dblVal );                  // copy construction from different type
    BOOST_CHECK_EQUAL( intVal2, -9 );
    BOOST_CHECK_EQUAL( intVal2.min(), -9 );
    BOOST_CHECK_EQUAL( intVal2.max(),  9 );


    BoundValue<int, detail::WRAP> wrappedVal( 0, 0, 5 ); // range [0,5)  ->  equivalent to modulo 5
    wrappedVal = 6;
    BOOST_CHECK_EQUAL( wrappedVal, 1 );
    wrappedVal = -1;
    BOOST_CHECK_EQUAL( wrappedVal, 4 );
    wrappedVal = -126;
    BOOST_CHECK_EQUAL( wrappedVal, 4 );
    wrappedVal = 126;
    BOOST_CHECK_EQUAL( wrappedVal, 1 );


    using redux::PI;
    double delta( 0.05 );
    double epsilon = 1E-9;

    BoundValue<double, detail::WRAP> periodicVal( 0, 0, 2 * PI );
    periodicVal = 3 * PI + delta;
    BOOST_CHECK_CLOSE( ( double )periodicVal, PI + delta, epsilon );
    periodicVal = -3 * PI / 2 + delta;
    BOOST_CHECK_CLOSE( ( double )periodicVal, PI / 2 + delta, epsilon );


    BoundValue<double, detail::REFLECT> inclVal( 0, 0, PI ); // e.g. inclination = PI/2+\delta -> PI/2 - \delta
    inclVal = PI + delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, PI - delta, epsilon );
    inclVal = 2 * PI + delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, delta, epsilon );
    inclVal = 3 * PI + delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, PI - delta, epsilon );
    inclVal = -delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, delta, epsilon );
    inclVal = -PI - delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, PI - delta, epsilon );
    inclVal = -2 * PI - delta;
    BOOST_CHECK_CLOSE( ( double )inclVal, delta, epsilon );


    BoundValue<int, detail::REFLECT> reflectedVal( 0, 0, 4 );
    reflectedVal = 0;
    BOOST_CHECK_EQUAL( reflectedVal, 0 );
    reflectedVal = 5;
    BOOST_CHECK_EQUAL( reflectedVal, 3 );
    reflectedVal = 9;
    BOOST_CHECK_EQUAL( reflectedVal, 1 );
    reflectedVal = 11;
    BOOST_CHECK_EQUAL( reflectedVal, 3 );
    reflectedVal.setMax( 5 );
    reflectedVal = -3;
    BOOST_CHECK_EQUAL( reflectedVal, 3 );
    reflectedVal = -5;
    BOOST_CHECK_EQUAL( reflectedVal, 5 );
    reflectedVal = -9;
    BOOST_CHECK_EQUAL( reflectedVal, 1 );
    reflectedVal = -11;
    BOOST_CHECK_EQUAL( reflectedVal, 1 );

}


void dataTest( void ) {

    uint32_t a[10], b[10], aa[10], bb[10];
    for( uint16_t i( 0 ); i < 10; i++ ) {
        a[i] = aa[i] = rand();
        b[i] = bb[i] = rand();
    }

    // check the array wrapper for std::swap
    swap( a, b, 10 );
    for( uint16_t i( 0 ); i < 10; i++ ) {
        BOOST_CHECK_EQUAL( b[i], aa[i] );
        BOOST_CHECK_EQUAL( a[i], bb[i] );
    }

}


void endianTest( void ) {

    union {
        uint32_t i;
        uint8_t c[4];
        bool b[4];
    } a, aa;

    a.i = 1;
    BOOST_CHECK( a.b[0] == ( REDUX_BYTE_ORDER == REDUX_LITTLE_ENDIAN ) );

    for( int i( 0 ); i < 4; i++ ) {
        a.c[i] = aa.c[i] = i;
    }
    swapEndian( a.i );
    for( int i( 0 ); i < 4; i++ ) {
        BOOST_CHECK_EQUAL( a.c[i], 3 - i ); // verify
    }

    uint32_t b[] = { aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i, aa.i };
    swapEndian( b, 10 );
    for( int i( 0 ); i < 4; i++ ) {
        BOOST_CHECK_EQUAL( b[i], a.i ); // verify
    }

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


void stringTest( void ) {
    
    // checks
    BOOST_CHECK( onlyDigits( "1234567890" ) );
    BOOST_CHECK( !onlyDigits( "1234d567890" ) );
    BOOST_CHECK( onlyAlpha( "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
    BOOST_CHECK( !onlyAlpha( "abcdefghijklmnopqrstuvwxyz7ABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
    BOOST_CHECK( onlyAlnum( "1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
    BOOST_CHECK( !onlyAlnum( "1234567890abcdefghijklmnopqrstuvw@xyzABCDEFGHIJKLMNOPQRSTUVWXYZ" ) );
    BOOST_CHECK( onlyHex( "1234567890abcdefABCDEF" ) );
    BOOST_CHECK( !onlyHex( "1234567890aQbcdefABCDEF" ) );
    BOOST_CHECK( isInteger( "1234567890" ) );
    BOOST_CHECK( isInteger( "0x123" ) );
    BOOST_CHECK( isInteger( "-19" ) );
    BOOST_CHECK( !isInteger( "-19.4" ) );
    BOOST_CHECK( !isInteger( "" ) );


    // align
    BOOST_CHECK_EQUAL( alignLeft( "string", 10, '-' ),   "string----" );
    BOOST_CHECK_EQUAL( alignRight( "string", 10, '-' ),  "----string" );
    BOOST_CHECK_EQUAL( alignCenter( "string", 10, '|' ), "||string||" );
    BOOST_CHECK_EQUAL( alignCenter( "string", 11, '|' ), "|||string||" ); // for odd padding, the extra char goes on the left side


    // hexString
    BOOST_CHECK_EQUAL( hexString( uint8_t( 0xc ) ), "0xc" );
    BOOST_CHECK_EQUAL( hexString( int16_t( 0x301 ) ), "0x301" );
    BOOST_CHECK_EQUAL( hexString( int32_t( 0xf5a3b1 ) ), "0xf5a3b1" );
    BOOST_CHECK_EQUAL( hexString( size_t( 0xff0af5a3b1 ) ), "0xff0af5a3b1" );
    BOOST_CHECK_EQUAL( hexString( uint8_t( 0xc ), false ), "c" );
    BOOST_CHECK_EQUAL( hexString( int16_t( 0x301 ), false ), "301" );
    BOOST_CHECK_EQUAL( hexString( int32_t( 0xf5a3b1 ), false ), "f5a3b1" );
    BOOST_CHECK_EQUAL( hexString( size_t( 0xff0af5a3b1 ), false ), "ff0af5a3b1" );


    // bitString
    BOOST_CHECK_EQUAL( bitString( uint8_t( 0xc ) ), "00001100" );
    BOOST_CHECK_EQUAL( bitString( int16_t( 0x301 ) ), "00000011 00000001" );
    uint32_t tmp = 0xf5a3b1;
    BOOST_CHECK_EQUAL( bitString( tmp ), "00000000 11110101 10100011 10110001" );
    BOOST_CHECK_EQUAL( bitString( ( uint8_t* )&tmp, 4 ), "00000000 11110101 10100011 10110001" );


}

namespace testsuite {

    namespace util {

        void utilTest( void ) {

            test_suite* ts = BOOST_TEST_SUITE( "UTIL" );

            ts->add( BOOST_TEST_CASE( &arrayTest ) );
            ts->add( BOOST_TEST_CASE( &bitTest ) );
            ts->add( BOOST_TEST_CASE( &boundValueTest ) );
            ts->add( BOOST_TEST_CASE( &dataTest ) );
            ts->add( BOOST_TEST_CASE( &endianTest ) );
            ts->add( BOOST_TEST_CASE( &statTest ) );
            ts->add( BOOST_TEST_CASE( &stringTest ) );

            framework::master_test_suite().add( ts );

        }

    }

}

