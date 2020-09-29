
#include "redux/image/fouriertransform.hpp"
#include "redux/math/helpers.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arraystats.hpp"

#include <boost/test/unit_test.hpp>

using namespace redux::image;
using namespace redux::math;
using namespace redux::util;
using namespace redux;
using namespace std;

using namespace boost::unit_test_framework;

namespace testsuite {

    namespace util {

            
        template <typename T, typename ...S>
        void testArray(S ...s) {
            
            std::vector<size_t> dims({static_cast<size_t>(s)...});
            size_t nDims = dims.size();
            if ( nDims == 0 ) return;
            
            std::vector<size_t> subFirst=dims, subLast=dims;    // dimensions of a centered & half-sized subarray
            for( size_t i(0); i<nDims; ++i) {
                subFirst[i] *= 0.25;
                subLast[i]  *= 0.75;
            }
            
            Array<T> array( 2, s... );      // test if constructor works if a new dimension is prepended
            BOOST_CHECK_EQUAL( array.nDimensions(), nDims+1 );
            BOOST_CHECK_EQUAL( array.dimSize(0), 2 );
            for( size_t i(0); i<nDims; ++i) {
                BOOST_CHECK_EQUAL( array.dimSize(i+1), dims[i] );
            }
            
            array.resize( dims );
            size_t nElements = array.nElements();
            BOOST_CHECK_EQUAL( array.nDimensions(), nDims );
            for( size_t i(0); i<nDims; ++i) {
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
                for( size_t i(0); i<nDims; ++i) {
                    BOOST_CHECK_EQUAL( array.dimSize(i), dims[i] );
                }
                BOOST_CHECK( array2.sameSizes(array) );
                BOOST_CHECK( array.sameSizes(array2) );
                BOOST_CHECK( array2 == array );
                BOOST_CHECK( !(array2 != array) );
                *array2.ptr() = 123;      // 1 wrong value
                BOOST_CHECK( !( array2 == array ) );
                BOOST_CHECK( array2 != array );

                // test that arrays of different dimensions compare as different.
                array2.resize( 2, s... );
                BOOST_CHECK( !array2.sameSizes(array) );
                BOOST_CHECK( !array.sameSizes(array2) );
                BOOST_CHECK_EQUAL( array2.nDimensions(), nDims+1 );
                array2 = 3;      // all values equal
                BOOST_CHECK( !( array2 == array ) );
                BOOST_CHECK( array2 != array );
                
                vector<int64_t> sf = array2.first();
                vector<int64_t> sl = array2.last();
                sl[0] = sf[0];
                Array<T> subarray( array2, sf, sl );        // subarray of array2 of same size as array
                BOOST_CHECK( array.sameSizes(subarray) );   // size should match since only trivial dimensions differ
                BOOST_CHECK( array == subarray );           // ...as should the data for these arrays
                
            }
            
            // test operators with scalars
            array = 1;
            T* abeg = array.ptr();
            T* aend = abeg + array.nElements();
            std::for_each( abeg, aend, []( const T& t ){ BOOST_CHECK_EQUAL( t, 1 ); } );
            array += 10;
            std::for_each( abeg, aend, []( const T& t ){ BOOST_CHECK_EQUAL( t, 11 ); } );
            array -= 1;
            std::for_each( abeg, aend, []( const T& t ){ BOOST_CHECK_EQUAL( t, 10 ); } );
            array *= 2;
            std::for_each( abeg, aend, []( const T& t ){ BOOST_CHECK_EQUAL( t, 20 ); } );
            array /= 5;
            std::for_each( abeg, aend, []( const T& t ){ BOOST_CHECK_EQUAL( t, 4 ); } );
            array.zero();
            std::for_each( abeg, aend, []( const T& t ){ BOOST_CHECK_EQUAL( t, 0 ); } );

            /***** test iterators *****/
            {

                Array<T> subarray( array, subFirst, subLast );
                BOOST_CHECK_EQUAL( subarray.nDimensions(), nDims );
                for( size_t i(0); i<nDims; ++i) {
                    BOOST_CHECK_EQUAL( array.dimSize(i), dims[i] );
                }
                T* rawPtr = array.get();
                for( size_t x(0); x < nElements; ++x ) {
                    rawPtr[x] = x;          // set values equal to real offsets: 0,1,2,3...
                }

                typename Array<T>::iterator it = array.begin();
                typename Array<T>::const_iterator cit = array.begin();
                typename Array<T>::iterator sit = subarray.begin();
                
                // check that subarray spans from subFirst to subLast
                BOOST_CHECK_EQUAL( array.pos(subFirst).pos(), sit.pos() );
                BOOST_CHECK_EQUAL( array.pos(subLast).pos(), (--subarray.end()).pos() );
                
                // check that whole datablok is accessing the right element
                for( size_t i(0); i < nElements; ++i ) {
                    BOOST_CHECK_EQUAL( it.pos(), i );
                    BOOST_CHECK_EQUAL( cit.pos(), i );
                    BOOST_CHECK_EQUAL( *it, i );
                    BOOST_CHECK_EQUAL( *cit, i );
                    it++,cit++;
                }
            
                it = array.begin();
                cit = array.begin();
                
                // postfix ++
                for( size_t i(0); i < nElements; ++i ) {
                    BOOST_CHECK_EQUAL( *it++, i );
                    BOOST_CHECK_EQUAL( *cit++, i );
                }

                for( size_t i(0); i<subarray.nElements(); ++i ) {
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
                for( size_t i(0); i<subarray.nElements(); ++i ) {
                    --sit;
                }

                // at begin()
                BOOST_CHECK( it == array.begin() );
                BOOST_CHECK( cit == array.begin() );
                BOOST_CHECK( sit == subarray.begin() );

                // prefix ++
                for( size_t i(0); i < nElements; ++i ) {
                    BOOST_CHECK_EQUAL( *it, i );
                    BOOST_CHECK_EQUAL( *cit, i );
                    ++it;
                    ++cit;
                }
                for( size_t i(0); i<subarray.nElements(); ++i ) {
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
                for( size_t i(0); i<subarray.nElements(); ++i ) {
                    sit--;
                }

                // at begin()
                BOOST_CHECK( it == array.begin() );
                BOOST_CHECK( cit == array.begin() );
                BOOST_CHECK( sit == subarray.begin() );

                // set/check values via auto.
                int cnt = 0;
                for( auto & value : array ) {
                    value = ++cnt;
                }

                cnt = 0;
                for( auto& value : array ) {
                    BOOST_CHECK_EQUAL( value, ++cnt );
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
                for( auto& value : array ) {
                    BOOST_CHECK_EQUAL( value, 2 );
                }
                subarray = 13;
                for( auto& value : subarray ) {
                    //cout << "  it = " << it << endl;
                    BOOST_CHECK_EQUAL( (int)value, 13 );
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
                for( auto& value : array ) {
                    BOOST_CHECK_EQUAL( value, 2 );
                }

                for( auto& value : subarray ) {
                    BOOST_CHECK_EQUAL( value, 2 );
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
                for( unsigned int i = 0; i < nDims; ++i ) {
                    tmpV[i] = subFirst[i];
                }
                for( unsigned int i = 0; i < nDims; ++i ) {
                    Array<T> subarray2 = subarray;
                    vector<int64_t> tmp = tmpV;
                    BOOST_CHECK_EQUAL( subarray2.shift(i,-100000), -tmp[i]);                // can only shift to the edge
                    tmp[i] = 0;
                    BOOST_CHECK_EQUAL( subarray2.begin().pos(), array.pos(tmp).pos() );
                }
                for( unsigned int i = 0; i < nDims; ++i ) {
                    tmpV[i] = subLast[i];
                }
                for( unsigned int i = 0; i < nDims; ++i ) {
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
            for(auto& value: array) value = (cnt+=1);
        
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
        void testStat( S ...s ) {
            
            const size_t range(20);
            const size_t offset(100-range/2);
            
            std::vector<size_t> dims({static_cast<size_t>(s)...});
            size_t nDims = dims.size();
            if ( nDims == 0 ) return;                   //  Skip pathological test cases
            
            std::vector<size_t> subFirst=dims, subLast=dims;
            for( size_t i(0); i<nDims; ++i) {
                subFirst[i] *= 0.25;
                subLast[i]  *= 0.75;
            }
            
            Array<T> array( s... );
            T* rawPtr = array.get();
            size_t nEl = array.nElements();
            if( !nEl ) return;                          //  Skip pathological test cases
            
            T mn = std::numeric_limits<T>::max();
            T mx = std::numeric_limits<T>::min();
            double avg(0);
            double rms(0);
            for( size_t x(0); x < nEl; ++x ) {
                rawPtr[x] = static_cast<T>(offset+rand()%(range/5)+x%(range+1));          // generate some data
                double tmp = rawPtr[x];
                if( rawPtr[x] < mn ) mn = rawPtr[x];
                if( rawPtr[x] > mx ) mx = rawPtr[x];
                avg += tmp;
                rms += tmp*tmp;
            }
            avg /= nEl;
            rms /= nEl;
            double stddev = sqrt(rms - avg*avg);
            rms = sqrt(rms);
            
            ArrayStats stats;
            stats.getMinMaxMean( array );
            BOOST_CHECK_CLOSE( stats.min, mn, 1E-20 );
            BOOST_CHECK_CLOSE( stats.max, mx, 1E-20 );
            BOOST_CHECK_CLOSE( stats.mean, avg, 1E-20 );
            
            stats.getRmsStddev( array );
            BOOST_CHECK_CLOSE( stats.rms, rms, 1E-20 );
            BOOST_CHECK_CLOSE( stats.stddev, stddev, 1E-20 );

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
            for( int i(0); i < 4; ++i ) {
                for( int j(0); j < 5; ++j ) {
                    BOOST_CHECK_EQUAL( array4x5( i, j ), 5 * i + j + 1 );
                }
            }

            // test iterators
            {
                it = array4x5.begin();
                Array<int>::iterator endit = array4x5.end();
                Array<int>::const_iterator cit = array4x5.begin();
                Array<int>::const_iterator cendit = array4x5.end();

                // prefix ++
                for( int i(1); i <= 20; ++i ) {
                    BOOST_CHECK_EQUAL( *it, i );
                    BOOST_CHECK_EQUAL( *cit, i );
                    ++it;
                    ++cit;
                }

                // at end()
                BOOST_CHECK( it == endit );
                BOOST_CHECK( cit == cendit );


                // prefix --
                for( int i(20); i > 0; --i ) {
                    BOOST_CHECK_EQUAL( *--it, i );
                    BOOST_CHECK_EQUAL( *--cit, i );
                }

                // postfix ++
                for( int i(1); i <= 20; ++i ) {
                    BOOST_CHECK_EQUAL( *it++, i );
                    BOOST_CHECK_EQUAL( *cit++, i );
                }

                // at end()
                BOOST_CHECK( it == endit );
                BOOST_CHECK( cit == cendit );

                // postfix --
                for( int i(20); i > 0; --i ) {
                    it--;
                    cit--;
                    BOOST_CHECK_EQUAL( *it, i );
                    BOOST_CHECK_EQUAL( *cit, i );
                }

                // test with std auto iterator (by reference)
                {
                    int i(1);
                    for( auto& value : array4x5 ) {
                        BOOST_CHECK_EQUAL( value, i++ );
                    }
                }
            }
            
            // shallow copy (shared data)
            array3x3 = array4x5;
            BOOST_CHECK( array3x3 == array4x5 );
            // verify that data is shared
            BOOST_CHECK( array3x3.get() == array4x5.get() );

            // deep copy
            array3x3 = array4x5.copy();
            BOOST_CHECK( array3x3 == array4x5 );
            // verify that data is not shared
            BOOST_CHECK( array3x3.get() != array4x5.get() );

            BOOST_CHECK( array4x5.dense() );

            // test sub-array
            {
                Array<int> subarray( array4x5, 0, 0, 0, 4 ); // 1x5 sub-array (slice)
                for( unsigned int i(0); i<4; ++i ) {
                    BOOST_CHECK( subarray.get() == array4x5.get() );
                    BOOST_CHECK( subarray.dense() );
                    subarray.shift(0,1);
                }

                subarray.wrap( array4x5, 1, 3, 1, 3 ); // 3x3 sub-array, with offset (1,1)
                BOOST_CHECK( !subarray.dense() );

                // verify that data is shared
                BOOST_CHECK( subarray.get() == array4x5.get() );
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
                for( int i(0); i < 4; ++i ) {
                    for( int j(0); j < 5; ++j ) {
                        int val = 5*i + j + 1;
                        if( i >= 1 && i <= 3 && j >= 1 && j <= 3 ) val *= 10;
                        BOOST_CHECK_EQUAL( array4x5( i, j ), val );
                    }
                }
            }
            
            Array<int> array4x5x6( 4, 5, 6 );
            size_t cnt(0);
            for( auto& value : array4x5x6 ) {
                value = ++cnt;
            }
            BOOST_CHECK_EQUAL( cnt, 120 );

            {
                // arbitrary subarray.
                Array<int> subarray( array4x5x6, 1, 2, 2, 3, 3, 4 );
        
                // check that the subarray is accessing the right elements.
                for( size_t i(0); i < 2; ++i ) {
                    for( size_t j(0); j < 2; ++j ) {
                        for( size_t k = 0; k < 2; ++k ) {
                            BOOST_CHECK_EQUAL( array4x5x6( i+1, j+2, k+3 ), subarray( i, j, k ) );
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
                cnt = 0;
                for( auto & value : array4x5x6 ) {
                    value = ++cnt;
                }
                shared_ptr<int**> raiiArray = array4x5x6.reshape(4,5,6);
                int*** rawArray = raiiArray.get();
                cnt = 0;
                for( size_t i(0); i<4; ++i) {
                    for( size_t j(0); j<5; ++j) {
                        for( size_t k(0); k<6; ++k) {
                            BOOST_CHECK_EQUAL( rawArray[i][j][k], ++cnt );
                        }
                    }
                }
            }
            
            // test different dimensions and datatypes
            testArray<int8_t>(100);
            testArray<uint16_t>(20,30);
            testArray<int>(4,1,9);
            testArray<size_t>(1,5,7,1,5);
            testArray<double>(4,5,6,7);
            testArray<float>(100,100);
            testArray<uint32_t>(13,9);
            
            transformTest<int8_t>();
            transformTest<int>();
            transformTest<double>();
            transformTest<size_t>();
            
        }


        void arrayStatTest( void ) {
            
            // Numerical verification for various types.
            testStat<int8_t>( 100, 100 );
            testStat<int16_t>( 2, 100, 100 );
            testStat<int32_t>( 2, 3, 100, 100 );
            testStat<int64_t>( 100, 100 );
            testStat<uint8_t>( 100, 100 );
            testStat<uint16_t>( 100, 100 );
            testStat<uint32_t>( 100, 100 );
            testStat<uint64_t>( 100, 100 );
            testStat<float>( 4, 100, 100 );
            testStat<double>( 2, 3, 100, 100 );
            
            // TODO
            // Stats for sub-array
            // Stats with clip
            // Stats for vectors and/or generic iterators?
            

        }


        void add_array_tests( test_suite* ts ) {

            ts->add( BOOST_TEST_CASE_NAME( &arrayTest, "Array manipulations" ) );
            ts->add( BOOST_TEST_CASE_NAME( &arrayStatTest, "Statistics and numerical tools"  ) );

        }

    }

}

