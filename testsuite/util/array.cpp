
#include "redux/image/fouriertransform.hpp"
#include "redux/math/helpers.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arraystats.hpp"

#include "catch.hpp"

#include <boost/filesystem.hpp>

using namespace redux::image;
using namespace redux::math;
using namespace redux::util;
using namespace redux;
using namespace std;

namespace bfs = boost::filesystem;

namespace testsuite {

    namespace util {

        namespace {

            
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
            CHECK( (array.nDimensions()) == (nDims+1) );
            CHECK( (array.dimSize(0)) == (2) );
            for( size_t i(0); i<nDims; ++i) {
                CHECK( (array.dimSize(i+1)) == (dims[i]) );
            }
            
            array.resize( dims );
            size_t nElements = array.nElements();
            CHECK( (array.nDimensions()) == (nDims) );
            for( size_t i(0); i<nDims; ++i) {
                CHECK( (array.dimSize(i)) == (dims[i]) );
            }
            
            const std::vector<size_t>& strides = array.strides();
            size_t stride = 1;
            for( size_t i=nDims-1; i>0; --i) {
                CHECK( (strides[i]) == (stride) );
                stride *= dims[i];
            }
            
            {
                // test comparison and assignment
                Array<T> array2(s...);
                array2 = 3;      // set all values to 3
                array =  3;      // set all values to 3
                CHECK( (array2.nDimensions()) == (nDims) );
                for( size_t i(0); i<nDims; ++i) {
                    CHECK( (array.dimSize(i)) == (dims[i]) );
                }
                CHECK(  array2.sameSizes(array)  );
                CHECK(  array.sameSizes(array2)  );
                CHECK(  array2 == array  );
                CHECK(  !(array2 != array)  );
                *array2.ptr() = 123;      // 1 wrong value
                CHECK(  !( array2 == array )  );
                CHECK(  array2 != array  );

                // test that arrays of different dimensions compare as different.
                array2.resize( 2, s... );
                CHECK(  !array2.sameSizes(array)  );
                CHECK(  !array.sameSizes(array2)  );
                CHECK( (array2.nDimensions()) == (nDims+1) );
                array2 = 3;      // all values equal
                CHECK(  !( array2 == array )  );
                CHECK(  array2 != array  );
                
                vector<int64_t> sf = array2.first();
                vector<int64_t> sl = array2.last();
                sl[0] = sf[0];
                Array<T> subarray( array2, sf, sl );        // subarray of array2 of same size as array
                CHECK(  array.sameSizes(subarray)  );   // size should match since only trivial dimensions differ
                CHECK(  array == subarray  );           // ...as should the data for these arrays
                
            }
            
            // test operators with scalars
            array = 1;
            T* abeg = array.ptr();
            T* aend = abeg + array.nElements();
            std::for_each( abeg, aend, []( const T& t ){ CHECK( (t) == (1) ); } );
            array += 10;
            std::for_each( abeg, aend, []( const T& t ){ CHECK( (t) == (11) ); } );
            array -= 1;
            std::for_each( abeg, aend, []( const T& t ){ CHECK( (t) == (10) ); } );
            array *= 2;
            std::for_each( abeg, aend, []( const T& t ){ CHECK( (t) == (20) ); } );
            array /= 5;
            std::for_each( abeg, aend, []( const T& t ){ CHECK( (t) == (4) ); } );
            array.zero();
            std::for_each( abeg, aend, []( const T& t ){ CHECK( (t) == (0) ); } );

            /***** test iterators *****/
            {

                Array<T> subarray( array, subFirst, subLast );
                CHECK( (subarray.nDimensions()) == (nDims) );
                for( size_t i(0); i<nDims; ++i) {
                    CHECK( (array.dimSize(i)) == (dims[i]) );
                }
                T* rawPtr = array.get();
                for( size_t x(0); x < nElements; ++x ) {
                    rawPtr[x] = x;          // set values equal to real offsets: 0,1,2,3...
                }

                typename Array<T>::iterator it = array.begin();
                typename Array<T>::const_iterator cit = array.begin();
                typename Array<T>::iterator sit = subarray.begin();
                
                // check that subarray spans from subFirst to subLast
                CHECK( (array.pos(subFirst).pos()) == (sit.pos()) );
                CHECK( (array.pos(subLast).pos()) == ((--subarray.end()).pos()) );
                
                // check that whole datablok is accessing the right element
                for( size_t i(0); i < nElements; ++i ) {
                    CHECK( (it.pos()) == (i) );
                    CHECK( (cit.pos()) == (i) );
                    CHECK( (*it) == (i) );
                    CHECK( (*cit) == (i) );
                    it++,cit++;
                }
            
                it = array.begin();
                cit = array.begin();
                
                // postfix ++
                for( size_t i(0); i < nElements; ++i ) {
                    CHECK( (*it++) == (i) );
                    CHECK( (*cit++) == (i) );
                }

                for( size_t i(0); i<subarray.nElements(); ++i ) {
                    sit++;
                }

                // at end()
                CHECK(  it == array.end()  );
                CHECK(  cit == array.end()  );
                CHECK(  sit == subarray.end()  );
                
                // prefix --
                for( int i = nElements-1; i >= 0; --i ) {
                    CHECK( (*--it) == (i) );
                    CHECK( (*--cit) == (i) );
                }
                for( size_t i(0); i<subarray.nElements(); ++i ) {
                    --sit;
                }

                // at begin()
                CHECK(  it == array.begin()  );
                CHECK(  cit == array.begin()  );
                CHECK(  sit == subarray.begin()  );

                // prefix ++
                for( size_t i(0); i < nElements; ++i ) {
                    CHECK( (*it) == (i) );
                    CHECK( (*cit) == (i) );
                    ++it;
                    ++cit;
                }
                for( size_t i(0); i<subarray.nElements(); ++i ) {
                    ++sit;
                }


                // at end()
                CHECK(  it == array.end()  );
                CHECK(  cit == array.end()  );
                CHECK(  sit == subarray.end()  );


                // postfix --
                for( int i = nElements-1; i >= 0; --i ) {
                    it--;
                    cit--;
                    CHECK( (*it) == (i) );
                    CHECK( (*cit) == (i) );
                }
                for( size_t i(0); i<subarray.nElements(); ++i ) {
                    sit--;
                }

                // at begin()
                CHECK(  it == array.begin()  );
                CHECK(  cit == array.begin()  );
                CHECK(  sit == subarray.begin()  );

                // set/check values via auto.
                int cnt = 0;
                for( auto & value : array ) {
                    value = ++cnt;
                }

                cnt = 0;
                for( auto& value : array ) {
                    CHECK( (value) == (++cnt) );
                }

            }
            /**************************/
            
            {
                Array<T> subarray( array, subFirst, subLast );

                // compare pointers to verify that the data is shared
                CHECK(  subarray.get() == array.get()  );
                
                // check if subarrays are dense
                CHECK(  array.dense()  );
                if(nDims == 1) {
                    CHECK(  subarray.dense()  );        // 1D arrays are always dense
                } else {
                    CHECK(  !subarray.dense()  );
                }

                // scalar assign to subarray
                array = 2;
                for( auto& value : array ) {
                    CHECK( (value) == (2) );
                }
                subarray = 13;
                for( auto& value : subarray ) {
                    //cout << "  it = " << it << endl;
                    CHECK( ((int)value) == (13) );
                }
                
                typename Array<T>::iterator ait = array.begin();
                typename Array<T>::iterator sit = subarray.begin();
                size_t count1(0),count2(0);
                while( ait < array.end() ) {
                    count1++;
                    if( ait.pos() < sit.pos() || sit == subarray.end()) {   // outside subarray, should have value = 2
                        CHECK( ((int)*ait++) == (2) );
                    } else {                        // inside subarray, should have value = 13
                        CHECK( ((int)*ait++) == (13) );
                        if(sit < subarray.end()) CHECK( (*sit++) == (13) );
                        count2++;
                    }
                }
                CHECK( (count1) == (nElements) );
                CHECK( (count2) == (subarray.nElements()) );
                
                // elementwise modification
                array = 2;
                for( auto& value : array ) {
                    CHECK( (value) == (2) );
                }

                for( auto& value : subarray ) {
                    CHECK( (value) == (2) );
                }
                
                ait = array.begin();
                sit = subarray.begin();
                while( sit < subarray.end() ) {
                    *sit *= 2;
                    CHECK( (*(ait+sit.pos())) == (4) );   // check value changed
                    sit++;
                }
                
                // compare pointers to verify that the data is still shared
                CHECK(  subarray.get() == array.get()  );

                // shift subarray
                vector<int64_t> tmpV(nDims);
                for( unsigned int i = 0; i < nDims; ++i ) {
                    tmpV[i] = subFirst[i];
                }
                for( unsigned int i = 0; i < nDims; ++i ) {
                    Array<T> subarray2 = subarray;
                    vector<int64_t> tmp = tmpV;
                    CHECK( (subarray2.shift(i,-100000)) == (-tmp[i]) );                // can only shift to the edge
                    tmp[i] = 0;
                    CHECK( (subarray2.begin().pos()) == (array.pos(tmp).pos()) );
                }
                for( unsigned int i = 0; i < nDims; ++i ) {
                    tmpV[i] = subLast[i];
                }
                for( unsigned int i = 0; i < nDims; ++i ) {
                    Array<T> subarray2 = subarray;
                    vector<int64_t> tmp = tmpV;
                    CHECK( (subarray2.shift(i,100000)) == (dims[i]-tmp[i]-1) );        // can only shift to the edge
                    tmp[i] = dims[i]-1;
                    CHECK( ((--subarray2.end()).pos()) == (array.pos(tmp).pos()) );
                }
            
            }
            
            {   // test pack/unpack
                
                auto buf = sharedArray<char>( array.size() );
                char* ptr = buf.get();
                uint64_t count = array.pack(ptr);
                CHECK( (count) == (array.size()) );
                
                Array<T> tmp;
                count = tmp.unpack(ptr,false);
                CHECK( (count) == (array.size()) );
                CHECK( (count) == (tmp.size()) );
                CHECK(  tmp == array  );
                
                Array<T> tmp2(array,subFirst,subLast);          // pack/unpack subarray;
                buf = sharedArray<char>( tmp2.size() );
                ptr = buf.get();
                count = tmp2.pack(ptr);
                CHECK( (count) == (tmp2.size()) );
                
                count = tmp.unpack(ptr,false);
                CHECK( (count) == (tmp2.size()) );
                CHECK( (count) == (tmp.size()) );
                CHECK(  tmp == tmp2  );

                tmp2.resize();                                  // pack/unpack empty array;
                buf = sharedArray<char>( tmp2.size() );
                ptr = buf.get();
                count = tmp2.pack(ptr);
                CHECK( (count) == (tmp2.size()) );
                
                count = tmp.unpack(ptr,false);
                CHECK( (count) == (tmp2.size()) );
                CHECK( (count) == (tmp.size()) );
                CHECK(  tmp == tmp2  );

                Array<T> tmp3;
                buf = sharedArray<char>( tmp3.size() );
                ptr = buf.get();
                count = tmp3.pack(ptr);
                CHECK( (count) == (tmp3.size()) );
                
                count = tmp.unpack(ptr,false);
                CHECK( (count) == (tmp3.size()) );
                CHECK( (count) == (tmp.size()) );
                CHECK(  tmp == tmp3  );
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
                        CHECK( (array( i, j, k )) == (sizeX*sizeY*i + j*sizeX + k + 1) );
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
                        CHECK( (array( i, j, k )) == (array2( i, j, sizeX-k-1)) );
                    }
                }
                reverseX(rawPtr[i],sizeY,sizeX);        // flip back
            }

            // check reverseY
            for( int i = 0; i < sizeZ; ++i ) {
                reverseY(rawPtr[i],sizeY,sizeX);
                for( int j = 0; j < sizeY; ++j ) {
                    for( int k = 0; k < sizeX; ++k ) {
                        CHECK( (array( i, j, k )) == (array2( i, sizeY-j-1, k)) );
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
                        CHECK( (array(i,k,j)) == (array2(i,j,k)) );
                    }
                }
                transpose(*rawPtr[i],sizeX,sizeY);        // transpose back
            }
            
            CHECK(  array == array2  );             // Just to verify they are still the same.
            
        }

        template <typename T, typename ...S>
        void testStat( S ...s ) {
            
            const size_t range(20);
            const size_t offset(100-range/2);
            const double tiny(1E-15);
            
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
            CHECK( (stats.min) == Approx(double(mn)).epsilon((tiny)/100.0) );
            CHECK( (stats.max) == Approx(double(mx)).epsilon((tiny)/100.0) );
            CHECK( (stats.mean) == Approx(avg).epsilon((tiny)/100.0) );
            
            stats.getRmsStddev( array );
            CHECK( (stats.rms) == Approx(rms).epsilon((tiny)/100.0) );
            CHECK( (stats.stddev) == Approx(stddev).epsilon((tiny)/100.0) );

            // Same data via a plain std::vector's buffer
            std::vector<T> vec( rawPtr, rawPtr+nEl );
            ArrayStats vecStats;
            vecStats.getMinMaxMean( vec.data(), vec.size() );
            CHECK( (vecStats.min) == Approx(double(mn)).epsilon((tiny)/100.0) );
            CHECK( (vecStats.max) == Approx(double(mx)).epsilon((tiny)/100.0) );
            CHECK( (vecStats.mean) == Approx(avg).epsilon((tiny)/100.0) );

            vecStats.getRmsStddev( vec.data(), vec.size() );
            CHECK( (vecStats.rms) == Approx(rms).epsilon((tiny)/100.0) );
            CHECK( (vecStats.stddev) == Approx(stddev).epsilon((tiny)/100.0) );

            // Stats on a (possibly dense, non-zero-offset) sub-array view
            Array<T> sub( array, subFirst, subLast );
            size_t subN = sub.nElements();
            if( subN ) {
                T subMn = std::numeric_limits<T>::max();
                T subMx = std::numeric_limits<T>::min();
                double subAvg(0), subRms(0);
                for( const T& v: sub ) {
                    double tmp = v;
                    if( v < subMn ) subMn = v;
                    if( v > subMx ) subMx = v;
                    subAvg += tmp;
                    subRms += tmp*tmp;
                }
                subAvg /= subN;
                subRms /= subN;
                double subStddev = sqrt(subRms - subAvg*subAvg);
                subRms = sqrt(subRms);

                ArrayStats subStats;
                subStats.getMinMaxMean( sub );
                CHECK( (subStats.min) == Approx(double(subMn)).epsilon((tiny)/100.0) );
                CHECK( (subStats.max) == Approx(double(subMx)).epsilon((tiny)/100.0) );
                CHECK( (subStats.mean) == Approx(subAvg).epsilon((tiny)/100.0) );

                subStats.getRmsStddev( sub );
                CHECK( (subStats.rms) == Approx(subRms).epsilon((tiny)/100.0) );
                CHECK( (subStats.stddev) == Approx(subStddev).epsilon((tiny)/100.0) );
            }

        }

        // ArrayStats::getStats(borderClip, ...) is only instantiated for these 4 types.
        template <typename T, typename ...S>
        void testStatBorderClip( uint32_t borderClip, S ...s ) {

            const size_t range(20);
            const size_t offset(100-range/2);
            const double tiny(1E-15);

            std::vector<size_t> dims({static_cast<size_t>(s)...});
            size_t nDims = dims.size();
            if( nDims == 0 ) return;
            for( size_t d: dims ) {
                if( d > 1 && 2*borderClip >= d ) return;      // would clip everything away
            }

            Array<T> array( s... );
            T* rawPtr = array.get();
            size_t nEl = array.nElements();
            if( !nEl ) return;

            for( size_t x(0); x < nEl; ++x ) {
                rawPtr[x] = static_cast<T>(offset+rand()%(range/5)+x%(range+1));
            }

            std::vector<int64_t> first = array.first();
            std::vector<int64_t> last = array.last();
            for( size_t i(0); i<nDims; ++i ) {
                if( dims[i] > 1 ) {
                    first[i] += borderClip;
                    last[i]  -= borderClip;
                }
            }
            Array<T> clipped( array, first, last );

            size_t nClipped = clipped.nElements();
            T mn = std::numeric_limits<T>::max();
            T mx = std::numeric_limits<T>::min();
            double avg(0), rms(0);
            for( const T& v: clipped ) {
                double tmp = v;
                if( v < mn ) mn = v;
                if( v > mx ) mx = v;
                avg += tmp;
                rms += tmp*tmp;
            }
            avg /= nClipped;
            rms /= nClipped;
            double stddev = sqrt(rms - avg*avg);
            rms = sqrt(rms);

            ArrayStats stats;
            stats.getStats( borderClip, array, ST_ALL );
            CHECK( (stats.min) == Approx(double(mn)).epsilon((tiny)/100.0) );
            CHECK( (stats.max) == Approx(double(mx)).epsilon((tiny)/100.0) );
            CHECK( (stats.mean) == Approx(avg).epsilon((tiny)/100.0) );
            CHECK( (stats.rms) == Approx(rms).epsilon((tiny)/100.0) );
            CHECK( (stats.stddev) == Approx(stddev).epsilon((tiny)/100.0) );

        }

        }

        TEST_CASE( "Array manipulations", "[util]" ) {

            Array<int> array4x5( 4, 5 );
            Array<int> array3x3;

            // check that assign throws for nValues > nElements;
            CHECK_THROWS_AS( array4x5.assign( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 ), logic_error );

            array4x5.assign( 1,  2,  3,  4,  5.0,   // all values will be cast to int in assign
                            6,  7,  8,  9, 10,
                            11, 12, 13, 14, 15,
                            16, 17, 18, 19, 20 );

            // check the iterator is accessing the right elements.
            Array<int>::iterator it = array4x5.begin();
            CHECK( (it.pos()) == (0) );
            it = array4x5.end();
            CHECK( (it.pos()) == (20) );
            CHECK( (*--it) == (20) );
            it = array4x5.pos(1,1);
            CHECK( (it.pos()) == (6) );
            CHECK( (it.step(1).pos()) == (7) );
            CHECK( (it.step(-1).pos()) == (6) );
            CHECK( (it.step(1,1).pos()) == (12) );
            CHECK( (it.step(-1,-1).pos()) == (6) );
            CHECK( (it.step(1,-1).pos()) == (10) );
            CHECK( (it.step(-1,1).pos()) == (6) );
        
            // check index out of bounds when using at()
            CHECK_THROWS_AS( array4x5.at( 1, 5 ), out_of_range );
            CHECK_THROWS_AS( array4x5.at( 5, 1 ), out_of_range );
            CHECK_THROWS_AS( array4x5.at( -1, 1 ), out_of_range );
            CHECK_THROWS_AS( array4x5.at( 1, -1 ), out_of_range );
            
            // test iterator returned by Array::pos() for some values
            CHECK( (*array4x5.pos( 0, 0 )) == (1) );
            CHECK( (*array4x5.pos( 1, 1 )) == (7) );
            CHECK( (*array4x5.pos( 3, 4 )) == (20) );

        
            // check values
            for( int i(0); i < 4; ++i ) {
                for( int j(0); j < 5; ++j ) {
                    CHECK( (array4x5( i, j )) == (5 * i + j + 1) );
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
                    CHECK( (*it) == (i) );
                    CHECK( (*cit) == (i) );
                    ++it;
                    ++cit;
                }

                // at end()
                CHECK(  it == endit  );
                CHECK(  cit == cendit  );


                // prefix --
                for( int i(20); i > 0; --i ) {
                    CHECK( (*--it) == (i) );
                    CHECK( (*--cit) == (i) );
                }

                // postfix ++
                for( int i(1); i <= 20; ++i ) {
                    CHECK( (*it++) == (i) );
                    CHECK( (*cit++) == (i) );
                }

                // at end()
                CHECK(  it == endit  );
                CHECK(  cit == cendit  );

                // postfix --
                for( int i(20); i > 0; --i ) {
                    it--;
                    cit--;
                    CHECK( (*it) == (i) );
                    CHECK( (*cit) == (i) );
                }

                // test with std auto iterator (by reference)
                {
                    int i(1);
                    for( auto& value : array4x5 ) {
                        CHECK( (value) == (i++) );
                    }
                }
            }
            
            // shallow copy (shared data)
            array3x3 = array4x5;
            CHECK(  array3x3 == array4x5  );
            // verify that data is shared
            CHECK(  array3x3.get() == array4x5.get()  );

            // deep copy
            array3x3 = array4x5.copy();
            CHECK(  array3x3 == array4x5  );
            // verify that data is not shared
            CHECK(  array3x3.get() != array4x5.get()  );

            CHECK(  array4x5.dense()  );

            // test sub-array
            {
                Array<int> subarray( array4x5, 0, 0, 0, 4 ); // 1x5 sub-array (slice)
                for( unsigned int i(0); i<4; ++i ) {
                    CHECK(  subarray.get() == array4x5.get()  );
                    CHECK(  subarray.dense()  );
                    subarray.shift(0,1);
                }

                subarray.wrap( array4x5, 1, 3, 1, 3 ); // 3x3 sub-array, with offset (1,1)
                CHECK(  !subarray.dense()  );

                // verify that data is shared
                CHECK(  subarray.get() == array4x5.get()  );
                // define an array equal to the center 3x3 elements of the 5x5 array above for reference
                array3x3.resize( 3, 3 );
                array3x3.assign( 7,  8,  9,
                                12, 13, 14,
                                17, 18, 19 );
                // verify that data matches
                CHECK(  subarray == array3x3  );


                // test subiterator::step() to step along any dimension
                CHECK( (*(subarray.pos( 1, 1 ).step())) == (14) );
                CHECK( (*(subarray.pos( 1, 1 ).step(1))) == (14) );
                CHECK( (*(subarray.pos( 1, 1 ).step(1,-1))) == (17) );
                CHECK( (*(subarray.pos( 1, 1 ).step(-1,-1))) == (7) );
                CHECK( (*(subarray.pos( 1, 1 ).step(-1,1))) == (9) );

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
                        CHECK( (array4x5( i, j )) == (val) );
                    }
                }
            }
            
            Array<int> array4x5x6( 4, 5, 6 );
            size_t cnt(0);
            for( auto& value : array4x5x6 ) {
                value = ++cnt;
            }
            CHECK( (cnt) == (120) );

            {
                // arbitrary subarray.
                Array<int> subarray( array4x5x6, 1, 2, 2, 3, 3, 4 );
        
                // check that the subarray is accessing the right elements.
                for( size_t i(0); i < 2; ++i ) {
                    for( size_t j(0); j < 2; ++j ) {
                        for( size_t k = 0; k < 2; ++k ) {
                            CHECK( (array4x5x6( i+1, j+2, k+3 )) == (subarray( i, j, k )) );
                        }
                    }
                }

                // assign to subarray
                subarray = 999;
                cnt = 0;
                for( auto cit=subarray.begin(); cit != subarray.end(); ++cit ) {
                    CHECK( (*cit) == (999) );
                    cnt++;
                }
                CHECK( (cnt) == (8) );
                
                cnt = 0;
                for( auto cit=subarray.end(); cit != subarray.begin();  ) {
                    CHECK( (*--cit) == (999) );
                    cnt++;
                }
                CHECK( (cnt) == (8) );

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
                            CHECK( (rawArray[i][j][k]) == (++cnt) );
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
            
        }

        TEMPLATE_TEST_CASE( "Array value-transform helpers (reverseX/reverseY/transpose)", "[util]",
                             int8_t, int, double, size_t ) {
            transformTest<TestType>();
        }

        TEST_CASE( "Array statistics across types/dimensionalities", "[util]" ) {
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

            testStatBorderClip<int16_t>( 3, 100, 100 );
            testStatBorderClip<int32_t>( 3, 2, 100, 100 );
            testStatBorderClip<float>( 3, 4, 100, 100 );
            testStatBorderClip<double>( 3, 2, 3, 100, 100 );
        }

        TEST_CASE( "Array wrap() over raw external memory", "[util]" ) {

            std::vector<int> external(20);
            for( size_t i=0; i<external.size(); ++i ) external[i] = (int)i;

            {
                Array<int> view;
                view.wrap( external.data(), 4, 5 );
                CHECK( view.nDimensions() == 2 );
                CHECK( view.dimSize(0) == 4 );
                CHECK( view.dimSize(1) == 5 );
                CHECK( view.dense() );
                CHECK( view.get() == external.data() );        // shares the external buffer, no copy

                for( int i=0; i<4; ++i ) {
                    for( int j=0; j<5; ++j ) {
                        CHECK( view(i,j) == external[i*5+j] );
                    }
                }

                // writing through the view modifies the external buffer directly
                view(2,3) = 999;
                CHECK( external[2*5+3] == 999 );
            }
            // view destructed, the external buffer should still be intact (not freed).
            CHECK( external[2*5+3] == 999 );

        }

        TEST_CASE( "Array resetLimits() restores the full underlying view", "[util]" ) {

            Array<int> array(4,5);
            int cnt(0);
            for( auto& v: array ) v = ++cnt;

            Array<int> sub( array, 1,2,1,3 );
            CHECK( sub.nElements() == 2*3 );

            sub.resetLimits();
            CHECK( sub.dimSize(0) == 4 );
            CHECK( sub.dimSize(1) == 5 );
            CHECK( sub.nElements() == 20 );
            CHECK( sub.dense() );
            CHECK( sub == array );

        }

        TEST_CASE( "Array setLimits() clamps out-of-range and swaps first>last", "[util]" ) {

            Array<int> array(10);
            int cnt(0);
            for( auto& v: array ) v = ++cnt;        // values 1..10

            Array<int> sub( array, std::vector<int64_t>{2}, std::vector<int64_t>{20} );
            CHECK( sub.dimSize(0) == 8 );            // indices 2..9
            CHECK( sub(0) == 3 );
            CHECK( sub(7) == 10 );

            // first > last -> swapped, not an empty/inverted range
            Array<int> sub2( array, std::vector<int64_t>{7}, std::vector<int64_t>{3} );
            CHECK( sub2.dimSize(0) == 5 );           // indices 3..7
            CHECK( sub2(0) == 4 );
            CHECK( sub2(4) == 8 );

        }

        TEST_CASE( "Array shift() moves an in-bounds view by the exact requested delta", "[util]" ) {

            Array<int> array(10);
            int cnt(0);
            for( auto& v: array ) v = ++cnt;        // values 1..10

            Array<int> sub( array, std::vector<int64_t>{2}, std::vector<int64_t>{4} );   // indices 2..4 -> 3,4,5
            CHECK( sub(0) == 3 );
            CHECK( sub(2) == 5 );

            int actual = sub.shift( 0, 2 );          // well within bounds
            CHECK( actual == 2 );                    // unclamped: full requested shift performed
            CHECK( sub(0) == 5 );                    // now indices 4..6 -> 5,6,7
            CHECK( sub(2) == 7 );

            actual = sub.shift( 0, -1 );
            CHECK( actual == -1 );
            CHECK( sub(0) == 4 );

        }

        TEST_CASE( "Array sameSize() (total count) vs sameSizes() (per-dimension shape)", "[util]" ) {

            Array<int> a(4,5);      // 20 elements, shape (4,5)
            Array<int> b(5,4);      // 20 elements, shape (5,4) -- same count, different shape
            Array<int> c(2,10);     // 20 elements, shape (2,10)
            Array<int> d(4,6);      // 24 elements -- different count entirely

            CHECK(  a.sameSize(b)  );
            CHECK(  a.sameSize(c)  );
            CHECK(  !a.sameSize(d)  );

            CHECK(  !a.sameSizes(b)  );
            CHECK(  !a.sameSizes(c)  );
            CHECK(  !a.sameSizes(d)  );

            Array<int> aCopy(4,5);
            CHECK(  a.sameSizes(aCopy)  );
            CHECK(  a.sameSize(aCopy)  );

        }

        TEST_CASE( "Array at() returns correct values; ptr() throws on invalid offset", "[util]" ) {

            Array<int> array(4,5);
            int cnt(0);
            for( auto& v: array ) v = ++cnt;

            for( int i=0; i<4; ++i ) {
                for( int j=0; j<5; ++j ) {
                    CHECK( array.at(i,j) == array(i,j) );
                }
            }
            CHECK( array.at() == array.at(0,0) );    // 0-arg at() returns the very first element

            CHECK_THROWS_AS( array.ptr( std::vector<int64_t>{100,100} ), std::out_of_range );

        }

        TEST_CASE( "Array trim() and in-place copy(Array&) overloads", "[util]" ) {

            Array<int> array(4,5);
            int cnt(0);
            for( auto& v: array ) v = ++cnt;

            Array<int> sub( array, 1,2,1,3 );        // 2x3 sub-view, non-dense
            CHECK( !sub.dense() );

            // copy(Array<T>& out): in-place deep copy into a separate, owned target
            Array<int> target;
            sub.copy( target );
            CHECK( target.dense() );
            CHECK( target.sameSizes(sub) );
            CHECK( target == sub );
            CHECK( target.get() != array.get() );

            target(0,0) = -1;                        // must not affect the original view/array
            CHECK( sub(0,0) != -1 );

            // cross-type copy(Array<U>& out)
            Array<double> targetD;
            sub.copy( targetD );
            CHECK( targetD.sameSizes(sub) );
            for( size_t i=0; i<sub.dimSize(0); ++i ) {
                for( size_t j=0; j<sub.dimSize(1); ++j ) {
                    CHECK( targetD(i,j) == Approx( (double)sub(i,j) ) );
                }
            }

            // trim(): collapse a sub-view in-place into a dense, standalone copy of itself
            Array<int> trimmed( array, 1,2,1,3 );
            CHECK( !trimmed.dense() );
            trimmed.trim();
            CHECK( trimmed.dense() );
            CHECK( trimmed.sameSizes(sub) );
            CHECK( trimmed == sub );
            CHECK( trimmed.get() != array.get() );

        }

        TEST_CASE( "Array cross-type construction, sub-construction, and assignment", "[util]" ) {

            Array<int> src(3,4);
            int cnt(0);
            for( auto& v: src ) v = ++cnt;

            Array<double> dst( src );                // converting constructor
            CHECK( dst.sameSizes(src) );
            for( size_t i=0; i<3; ++i ) {
                for( size_t j=0; j<4; ++j ) {
                    CHECK( dst(i,j) == Approx( (double)src(i,j) ) );
                }
            }

            Array<double> dst2;
            dst2 = src;                               // cross-type operator=
            CHECK( dst2 == dst );

            Array<double> subD( src, 1,2,1,3 );        // cross-type sub-array constructor
            Array<int> subI( src, 1,2,1,3 );
            CHECK( subD.sameSizes(subI) );
            for( size_t i=0; i<subI.dimSize(0); ++i ) {
                for( size_t j=0; j<subI.dimSize(1); ++j ) {
                    CHECK( subD(i,j) == Approx( (double)subI(i,j) ) );
                }
            }

        }

        TEST_CASE( "Array-Array arithmetic operators", "[util]" ) {

            Array<int> a(3,4);
            Array<int> b(3,4);
            int cnt(0);
            for( auto& v: a ) v = ++cnt;             // 1..12
            cnt = 0;
            for( auto& v: b ) v = (++cnt)*2;         // 2,4,6...24

            Array<int> sum  = a + b;
            Array<int> diff = a - b;
            Array<int> prod = a * b;
            Array<int> quot = b / a;

            for( size_t i=0; i<3; ++i ) {
                for( size_t j=0; j<4; ++j ) {
                    CHECK( sum(i,j)  == a(i,j)+b(i,j) );
                    CHECK( diff(i,j) == a(i,j)-b(i,j) );
                    CHECK( prod(i,j) == a(i,j)*b(i,j) );
                    CHECK( quot(i,j) == b(i,j)/a(i,j) );
                }
            }

            Array<int> c = a.copy();
            c += b;
            CHECK( c == sum );
            c = a.copy(); c -= b; CHECK( c == diff );
            c = a.copy(); c *= b; CHECK( c == prod );
            Array<int> c2 = b.copy(); c2 /= a; CHECK( c2 == quot );

            Array<int> wrongShape(5,5);
            CHECK_THROWS_AS( a += wrongShape, std::invalid_argument );
            CHECK_THROWS_AS( a -= wrongShape, std::invalid_argument );
            CHECK_THROWS_AS( a *= wrongShape, std::invalid_argument );
            CHECK_THROWS_AS( a /= wrongShape, std::invalid_argument );

        }

        TEST_CASE( "Array iterator neighbourhood() and indexFromPos()", "[util]" ) {

            Array<int> array(10,10);
            int cnt(0);
            for( auto& v: array ) v = ++cnt;

            auto it = array.pos(5,5);
            Array<int> nb = it.neighbourhood(2,2);
            CHECK( nb.dimSize(0) == 5 );             // indices 3..7
            CHECK( nb.dimSize(1) == 5 );
            CHECK( nb(0,0) == array(3,3) );
            CHECK( nb(4,4) == array(7,7) );
            CHECK( nb.get() == array.get() );        // it's a view, not a copy

            auto idx = array.indexFromPos( it.pos() );
            CHECK( idx.size() == 2 );
            CHECK( idx[0] == 5 );
            CHECK( idx[1] == 5 );

            // near the boundary, the neighbourhood must clip rather than go out of bounds.
            auto itEdge = array.pos(0,0);
            Array<int> nbEdge = itEdge.neighbourhood(2,2);
            CHECK( nbEdge.dimSize(0) == 3 );         // can't go to -2, clamped to 0..2
            CHECK( nbEdge.dimSize(1) == 3 );
            CHECK( nbEdge(0,0) == array(0,0) );

        }

        TEST_CASE( "Array permuteDimensions() edge cases", "[util]" ) {

            Array<int> array(2,3,4);
            Array<int> before = array.copy();

            // too few dims given (<2): silently a no-op
            array.permuteDimensions( 0 );
            CHECK( array.dimensions() == before.dimensions() );

            // too many dims given (>nDims_): silently a no-op
            array.permuteDimensions( 0, 1, 2, 0 );
            CHECK( array.dimensions() == before.dimensions() );

            array.permuteDimensions( 0, 99 );
            CHECK( array.dimensions() == before.dimensions() );

            CHECK_THROWS_AS( array.permuteDimensions( 99, 0 ), std::out_of_range );

            Array<int> a2 = array.copy();
            size_t d0 = a2.dimSize(0), d1 = a2.dimSize(1);
            a2.permuteDimensions(0,1);
            CHECK( a2.dimSize(0) == d1 );
            CHECK( a2.dimSize(1) == d0 );

        }

        TEST_CASE( "Array fill() with a predicate gates which elements change", "[util]" ) {

            Array<int> array(10);
            int cnt(0);
            for( auto& v: array ) v = (++cnt) - 5;        // -4,-3,-2,-1,0,1,2,3,4,5

            Array<int> expected = array.copy();
            for( auto& v: expected ) if( v <= 0 ) v = 99;  // default predicate: v<=0

            array.fill( 99 );                              // use the default predicate
            CHECK( array == expected );

            // explicit predicate + a fill-function (per-element transform) variant
            Array<int> array2(10);
            cnt = 0;
            for( auto& v: array2 ) v = ++cnt;              // 1..10

            auto doubler = []( const int& v ) { return v*2; };
            auto isOdd = []( const int& v ) { return (v%2) != 0; };
            array2.fill( static_cast<int(*)(const int&)>(doubler), isOdd );

            for( int i=0; i<10; ++i ) {
                int orig = i+1;
                int expect = (orig%2) ? orig*2 : orig;     // only odd originals got doubled
                CHECK( array2(i) == expect );
            }

        }

        TEST_CASE( "Array safeDivide() divides elementwise with an epsilon cutoff", "[util]" ) {

            Array<complex_t> num(4);
            Array<double> den(4);
            num(0) = complex_t(10,0); den(0) = 2.0;         // normal division
            num(1) = complex_t(10,0); den(1) = 0.0;         // below cutoff -> zero
            num(2) = complex_t(6,8);  den(2) = 5.0;          // normal (|6+8i|/5 doesn't matter, just a/b)
            num(3) = complex_t(1,0);  den(3) = 1E-12;        // below default eps(0)? no. 1E-12>0, so NOT cut unless eps given

            Array<complex_t> result = num.copy();
            result.safeDivide( den );                        // default eps=0: only |b|>0 passes
            CHECK( result(0) == complex_t(5,0) );
            CHECK( result(1) == complex_t(0,0) );             // |0|>0 is false -> zeroed
            CHECK( result(2) == complex_t(6,8)/complex_t(5,0) );
            CHECK( std::abs(result(3) - complex_t(1,0)/complex_t(1E-12,0)) < 1E-3 );

            Array<complex_t> result2 = num.copy();
            result2.safeDivide( den, 1E-6 );                  // now the tiny denominator is also cut
            CHECK( result2(3) == complex_t(0,0) );

        }

        TEST_CASE( "Array weighted add()/addTo()/subtract()/mult()", "[util]" ) {

            Array<int> a(5), b(5);
            int cnt(0);
            for( auto& v: a ) v = ++cnt;         // 1..5
            cnt = 0;
            for( auto& v: b ) v = ++cnt;         // 1..5

            // add(Array rhs, Array weight): this[i] += rhs[i]*weight[i]
            Array<int> w(5);
            cnt = 0;
            for( auto& v: w ) v = (++cnt)%2;     // 1,0,1,0,1
            Array<int> target = a.copy();
            target.add( b, w );
            for( int i=0; i<5; ++i ) CHECK( target(i) == a(i) + b(i)*w(i) );

            // add(Array rhs, scalar weight): this[i] += rhs[i]*weight
            target = a.copy();
            target.add( b, 3 );
            for( int i=0; i<5; ++i ) CHECK( target(i) == a(i) + b(i)*3 );

            // addTo(Array& rhs, weight): rhs[i] += this[i]*weight  (roles reversed vs add())
            Array<int> dest = b.copy();
            a.addTo( dest, 2 );
            for( int i=0; i<5; ++i ) CHECK( dest(i) == b(i) + a(i)*2 );

            // subtract(Array rhs, weight): this[i] -= rhs[i]*weight
            target = a.copy();
            target.subtract( b, 2 );
            for( int i=0; i<5; ++i ) CHECK( target(i) == a(i) - b(i)*2 );

            // mult(Array rhs, weight): this[i] *= rhs[i]*weight
            target = a.copy();
            target.mult( b, 2 );
            for( int i=0; i<5; ++i ) CHECK( target(i) == a(i) * (b(i)*2) );

        }

        TEST_CASE( "Array copyFrom()/copyTo() across dense and sub-array views", "[util]" ) {

            Array<int> array(4,5);
            int cnt(0);
            for( auto& v: array ) v = ++cnt;

            // dense copyTo/copyFrom, with a type conversion (int -> double)
            std::vector<double> buf( array.nElements() );
            array.copyTo<double>( buf.data() );
            for( size_t i=0; i<buf.size(); ++i ) CHECK( buf[i] == Approx( (double)(i+1) ) );

            Array<int> target( 4, 5 );
            target.copyFrom<double>( buf.data() );
            CHECK( target == array );

            // non-dense (sub-array view) copyTo/copyFrom
            Array<int> sub( array, 1,2,1,3 );
            size_t subN = sub.nElements();
            CHECK( !sub.dense() );

            std::vector<int> subBuf( subN );
            sub.copyTo<int>( subBuf.data() );
            size_t idx(0);
            for( const int& v: sub ) CHECK( subBuf[idx++] == v );

            std::vector<int> newVals( subN );
            for( size_t i=0; i<subN; ++i ) newVals[i] = (int)(1000+i);
            sub.copyFrom<int>( newVals.data() );
            idx = 0;
            for( const int& v: sub ) CHECK( v == newVals[idx++] );

        }

        TEST_CASE( "Array cloneData() returns an independent owned copy", "[util]" ) {

            Array<int> array(3,4);
            int cnt(0);
            for( auto& v: array ) v = ++cnt;

            int* clone = array.cloneData();
            for( size_t i=0; i<array.nElements(); ++i ) CHECK( clone[i] == (int)(i+1) );

            clone[0] = -1;
            CHECK( array(0,0) != -1 );        // independent memory, not shared

            delete[] clone;

            // also works for a non-dense sub-array view
            Array<int> sub( array, 1,2,1,2 );
            int* subClone = sub.cloneData();
            size_t idx(0);
            for( const int& v: sub ) CHECK( subClone[idx++] == v );
            delete[] subClone;

        }

        TEST_CASE( "Array use_count()/empty()/valid()/mid()/getData()", "[util]" ) {

            Array<int> defaultArray;
            CHECK( defaultArray.empty() );
            CHECK( !defaultArray.valid() );          // never allocated -> null datablock

            Array<int> array(10);
            CHECK( !array.empty() );
            CHECK( array.valid() );
            CHECK( array.use_count() == 1 );

            {
                Array<int> shared = array;           // shallow copy, shares the datablock
                CHECK( array.use_count() == 2 );
                CHECK( shared.use_count() == 2 );
                CHECK( shared.getData().get() == array.get() );
            }
            CHECK( array.use_count() == 1 );          // back down after the copy went out of scope

            int cnt(0);
            for( auto& v: array ) v = ++cnt;          // 1..10
            CHECK( *array.mid() == array(5) );        // midpoint of the flat index range

        }

        TEST_CASE( "Array setBlock() copies a raw buffer in at a given index", "[util]" ) {

            Array<int> array(3,4);
            int cnt(0);
            for( auto& v: array ) v = ++cnt;
            Array<int> before = array.copy();

            int patch[2] = { -1, -2 };
            array.setBlock( patch, 2, 1, 1 );          // write 2 elements starting at (1,1)

            CHECK( array(1,1) == -1 );
            CHECK( array(1,2) == -2 );
            // everything else is unchanged
            for( size_t i=0; i<3; ++i ) {
                for( size_t j=0; j<4; ++j ) {
                    if( (i==1 && j==1) || (i==1 && j==2) ) continue;
                    CHECK( array(i,j) == before(i,j) );
                }
            }

        }

        TEST_CASE( "Array createMmap()/openMmap() round-trip through a backing file", "[util]" ) {

            bfs::path tmpFile = bfs::temp_directory_path() / bfs::unique_path( "redux_test_array_mmap_%%%%" );

            try {
                {
                    Array<double> mapped;
                    mapped.createMmap( tmpFile.string(), 3, 4 );
                    CHECK( mapped.nDimensions() == 2 );
                    CHECK( mapped.dimSize(0) == 3 );
                    CHECK( mapped.dimSize(1) == 4 );
                    int cnt(0);
                    for( auto& v: mapped ) v = (++cnt)*1.5;
                }
                // mapped went out of scope (munmap'd), but the backing file must persist
                // with the written values.
                CHECK( bfs::exists( tmpFile ) );

                Array<double> reopened;
                reopened.openMmap( tmpFile.string(), 3, 4 );
                CHECK( reopened.nDimensions() == 2 );
                CHECK( reopened.dimSize(0) == 3 );
                CHECK( reopened.dimSize(1) == 4 );
                int cnt(0);
                for( auto& v: reopened ) CHECK( v == Approx((++cnt)*1.5) );

            } catch( ... ) {
                bfs::remove( tmpFile );
                throw;
            }
            bfs::remove( tmpFile );

        }

    }

}

