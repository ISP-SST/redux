#include "redux/momfbd/momfbdjob.hpp"
#include "redux/logging/logger.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/cache.hpp"

#include "catch.hpp"

#include <boost/filesystem.hpp>

using namespace redux::momfbd;
using namespace redux::util;

namespace bfs = boost::filesystem;

using namespace std;

namespace testsuite {

    namespace momfbd {

        namespace {

            void dataTest() {

                MomfbdJob mjob;

                {   // default-constructed job
                    auto buf = sharedArray<char>( mjob.size() );
                    char* ptr = buf.get();
                    uint64_t count = mjob.pack( ptr );
                    CHECK( (count) == (mjob.size()) );
                    string tmpS = string( ptr );
                    redux::Job::JobPtr job = redux::Job::newJob( tmpS );
                    CHECK(  job  );
                    count = job->unpack( ptr, false );
                    CHECK( (count) == (mjob.size()) );
                    CHECK( (count) == (job->size()) );
                }

                {   // no image data
                    PatchData pd(mjob),pd2(mjob);
                    pd.id = 123;
                    pd.index = Point16(1,2);
                    auto buf = sharedArray<char>( pd.size() );
                    char* ptr = buf.get();
                    uint64_t count = pd.pack( ptr );
                    CHECK( (count) == (pd.size()) );
                    count = pd2.unpack( ptr, false );
                    CHECK( (count) == (pd.size()) );
                    CHECK( (count) == (pd2.size()) );
                    CHECK(  (pd == pd2)  );
                }

            }

        }

        TEST_CASE( "Data", "[momfbd]" ) {
            dataTest();
        }

        TEST_CASE( "PatchCache", "[momfbd]" ) {

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
                CHECK( (pd.getFullPath()) == (expectedPath) );

                auto buf = sharedArray<char>( pd.size() );
                uint64_t packedSize = pd.pack( buf.get() );
                CHECK( (packedSize) == (pd.size()) );

                CHECK(  pd.cacheStore(true)  );
                CHECK(  bfs::exists( pd.getFullPath() )  );

                CHECK(  pd.cacheLoad()  );

                auto buf2 = sharedArray<char>( pd.size() );
                uint64_t packedSize2 = pd.pack( buf2.get() );
                CHECK( (packedSize2) == (packedSize) );
                CHECK(  std::equal( buf.get(), buf.get()+packedSize, buf2.get() )  );

            } catch( ... ) {
                bfs::remove_all( tmpDir );
                throw;
            }
            bfs::remove_all( tmpDir );

        }

    }

}
