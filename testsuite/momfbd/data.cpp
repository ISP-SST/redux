
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/logging/logger.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/cache.hpp"

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

using namespace redux::momfbd;
using namespace redux::util;

namespace bfs = boost::filesystem;

using namespace std;

namespace testsuite {

    namespace momfbd {
        
        namespace {
            /*MomfbdJob mjob;     // unit-wide job for testing.
            
            
            int logInit(void) {
                redux::logging::Logger& lg = mjob.getLogger();
                lg.setLevel( 5 );
                lg.setContext( "testsuite" );
                lg.addStream( cout, 0, 1 );
                return 1;
            }

            
            int initOnce(void) {
                static int init_done = logInit();
                return init_done;
            }*/
        }


        void dataTest( void ) {
            
            return;
            MomfbdJob mjob;

            //int dummy RDX_UNUSED = initOnce();
            
            {   // default-constructed job
                auto buf = sharedArray<char>( mjob.size() );
                char* ptr = buf.get();
                uint64_t count = mjob.pack( ptr );
                BOOST_CHECK_EQUAL( count, mjob.size() );
                string tmpS = string( ptr );
                redux::Job::JobPtr job = redux::Job::newJob( tmpS );
                BOOST_CHECK( job );
                count = job->unpack( ptr, false );
                BOOST_CHECK_EQUAL( count, mjob.size() );
                BOOST_CHECK_EQUAL( count, job->size() );
                //             job->info.id = ++jobCounter;
                //             job->info.step.store( Job::JSTEP_PREPROCESS );
                //             job->info.name = "job_" + to_string( job->info.id );
                //             job->info.submitTime = boost::posix_time::second_clock::universal_time();
                //             ids.push_back( jobCounter );
                //             ids[0]++;
                //             jobs.push_back( job );

            }

            {   // no image data
                PatchData pd(mjob),pd2(mjob);
                pd.id = 123;
                pd.index = Point16(1,2);
                auto buf = sharedArray<char>( pd.size() );
                char* ptr = buf.get();
                uint64_t count = pd.pack( ptr );
                BOOST_CHECK_EQUAL( count, pd.size() );
                count = pd2.unpack( ptr, false );
                BOOST_CHECK_EQUAL( count, pd.size() );
                BOOST_CHECK_EQUAL( count, pd2.size() );
                BOOST_CHECK( pd == pd2 );

                // with images
                /*pd.images.resize(10,100,120);
                float cnt(0.3);
                for(auto& it: pd.images) it = (cnt += 1);
                buf = sharedArray<char>( pd.size() );
                ptr = buf.get();
                count = pd.pack( ptr );
                BOOST_CHECK_EQUAL( count, pd.size() );
                count = pd2.unpack( ptr, false );
                BOOST_CHECK_EQUAL( count, pd.size() );
                BOOST_CHECK_EQUAL( count, pd2.size() );
                BOOST_CHECK( pd == pd2 );*/

            }
            /*




                    Array<T> tmp;
                    count = tmp.unpack( ptr,false );
                    BOOST_CHECK_EQUAL( count, array.size() );

                    BOOST_CHECK( tmp == array );*/

        }


        void patchCacheTest( void ) {

            bfs::path tmpDir = bfs::temp_directory_path() / bfs::unique_path( "redux_test_patchcache_%%%%" );
            bfs::create_directories( tmpDir );
            Cache::get().setPath( tmpDir.string() );

            try {
                MomfbdJob mjob;
                PatchData pd(mjob);
                pd.id = 42;
                pd.index = Point16(3,7);
                pd.setPath( "cachetest" );

                string expectedPath = (tmpDir / "cachetest" / ("patch_"+(string)pd.index)).string();
                BOOST_CHECK_EQUAL( pd.getFullPath(), expectedPath );

                auto buf = sharedArray<char>( pd.size() );
                uint64_t packedSize = pd.pack( buf.get() );
                BOOST_CHECK_EQUAL( packedSize, pd.size() );

                BOOST_CHECK( pd.cacheStore(true) );
                BOOST_CHECK( bfs::exists( pd.getFullPath() ) );

                BOOST_CHECK( pd.cacheLoad() );

                auto buf2 = sharedArray<char>( pd.size() );
                uint64_t packedSize2 = pd.pack( buf2.get() );
                BOOST_CHECK_EQUAL( packedSize2, packedSize );
                BOOST_CHECK( std::equal( buf.get(), buf.get()+packedSize, buf2.get() ) );

            } catch( ... ) {
                bfs::remove_all( tmpDir );
                throw;
            }
            bfs::remove_all( tmpDir );

        }


        using namespace boost::unit_test;
        void add_data_tests( test_suite* ts ) {

            ts->add( BOOST_TEST_CASE_NAME( &dataTest, "Data" ) );
            ts->add( BOOST_TEST_CASE_NAME( &patchCacheTest, "PatchCache" ) );

        }

    }

}
