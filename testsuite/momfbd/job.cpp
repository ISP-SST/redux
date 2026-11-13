#include "redux/momfbd/momfbdjob.hpp"
#include "redux/logging/logger.hpp"
#include "redux/file/fileio.hpp"

#include "catch.hpp"

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>


using namespace redux::logging;
using namespace redux::momfbd;
using namespace std;
namespace bpt = boost::property_tree;


namespace testsuite {

    namespace momfbd {

        // Job::parsePropertyTree/startLog resolve OUTPUT_DIR/LOGFILE into actual
        // filesystem paths. Let's test that path concatenation works as we intend it to.
        TEST_CASE( "Job output/log path resolution", "[job]" ) {

            bpt::ptree tree;
            bpo::variables_map vm;        // left empty unless a section needs a cmd-line override

            SECTION( "no OUTPUT_DIR and no --output-dir falls back to cwd" ) {
                MomfbdJob job;
                job.Job::parsePropertyTree( vm, tree, job.getLogger() );
                CHECK( job.info.outputDir == bfs::current_path().string() );
            }

            SECTION( "relative OUTPUT_DIR is resolved against cwd" ) {
                MomfbdJob job;
                tree.put( "OUTPUT_DIR", "some/relative/dir" );
                job.Job::parsePropertyTree( vm, tree, job.getLogger() );
                CHECK( job.info.outputDir == ( bfs::current_path() / "some/relative/dir" ).string() );
            }

            SECTION( "absolute OUTPUT_DIR is used as-is" ) {
                MomfbdJob job;
                tree.put( "OUTPUT_DIR", "/some/absolute/dir" );
                job.Job::parsePropertyTree( vm, tree, job.getLogger() );
                CHECK( job.info.outputDir == "/some/absolute/dir" );
            }

            // The LOGFILE cases below exercise startLog() as well, which actually
            // opens the resolved path, so use a scratch directory instead of cwd.
            bfs::path tmpDir = bfs::temp_directory_path() / bfs::unique_path( "redux_test_job-%%%%-%%%%" );
            bfs::create_directories( tmpDir );

            SECTION( "relative LOGFILE is joined with OUTPUT_DIR" ) {
                MomfbdJob job;
                tree.put( "OUTPUT_DIR", tmpDir.string() );
                tree.put( "LOGFILE", "scan.log" );
                job.Job::parsePropertyTree( vm, tree, job.getLogger() );
                job.startLog( true );
                CHECK( bfs::exists( tmpDir / "scan.log" ) );
            }

            SECTION( "absolute LOGFILE is used as-is, ignoring OUTPUT_DIR" ) {
                MomfbdJob job;
                bfs::path otherDir = tmpDir / "unrelated_output_dir";
                bfs::path absLog = tmpDir / "scan_abs.log";
                tree.put( "OUTPUT_DIR", otherDir.string() );
                tree.put( "LOGFILE", absLog.string() );
                job.Job::parsePropertyTree( vm, tree, job.getLogger() );
                // parsePropertyTree must preserve the leading '/',
                // otherwise startLog() below would wrongly treat it as relative and
                // nest it under OUTPUT_DIR instead.
                CHECK( job.info.logFile == absLog.string() );
                job.startLog( true );
                CHECK( bfs::exists( absLog ) );
                CHECK_FALSE( bfs::exists( otherDir ) );
            }

            bfs::remove_all( tmpDir );

        }

        // OUTPUT_FILE (per-object "img_file") is, by design, always forced relative,
        // then joined under OUTPUT_DIR.
        TEST_CASE( "OUTPUT_FILE path concatenation under OUTPUT_DIR", "[job]" ) {

            using redux::file::stripLeadingSlashes;
            using redux::file::isRelative;

            SECTION( "plain relative OUTPUT_FILE nests under OUTPUT_DIR" ) {
                string outputFileName = stripLeadingSlashes( "sub/dir/name" );
                bfs::path fn( outputFileName + ".fits" );
                if( isRelative( fn ) ) fn = bfs::path( "/abs/out" ) / fn;
                CHECK( fn.string() == "/abs/out/sub/dir/name.fits" );
            }

            SECTION( "leading-slash OUTPUT_FILE is stripped, then still nests (never escapes OUTPUT_DIR)" ) {
                string outputFileName = stripLeadingSlashes( "/looks/absolute/name" );
                bfs::path fn( outputFileName + ".fits" );
                if( isRelative( fn ) ) fn = bfs::path( "/abs/out" ) / fn;
                CHECK( fn.string() == "/abs/out/looks/absolute/name.fits" );
            }

            SECTION( "different OUTPUT_DIR values each combine correctly" ) {
                string outputFileName = stripLeadingSlashes( "name" );
                bfs::path fn( outputFileName + ".momfbd" );
                CHECK( isRelative( fn ) );
                CHECK( ( bfs::path( "/first/out" ) / fn ).string() == "/first/out/name.momfbd" );
                CHECK( ( bfs::path( "/second/out" ) / fn ).string() == "/second/out/name.momfbd" );
            }

        }

    }

}
