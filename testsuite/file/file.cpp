#include "redux/file/fileio.hpp"

#include "catch.hpp"

using namespace redux::file;
using namespace std;

#ifndef RDX_TESTDATA_DIR
#error RDX_TESTDATA_DIR not set
#endif


namespace testsuite {
    namespace file {

        TEST_CASE( "Paths/Filenames", "[file]" ) {

            CHECK( isRelative("abc") );
            CHECK( !isRelative("/abc") );

            string tmp;
            bfs::path test_dir(RDX_TESTDATA_DIR);

            string my_home = getenv("HOME");
            string username = getenv("USER");
            if( !my_home.empty() ) {
                if( *my_home.rbegin() != '/' ) my_home.push_back('/');    // getHome() returns with slash, so append if needed.
                tmp = getHome();
                CHECK( tmp == my_home );
                if( !username.empty() ) {
                    CHECK( getHome(username) == my_home );
                }
            } else my_home = getHome();

            CHECK( getHome("root") == "/root/" );
            CHECK( getHome("unknown_user_zxy") == "" );

            CHECK( bfs::exists(test_dir) );
            CHECK( bfs::is_directory(test_dir) );

            // TODO: test with resolution of known symlink
            CHECK( weaklyCanonical(test_dir / "non/existing/../path/") == (test_dir/"non/path") );
            CHECK( weaklyCanonical( "~root/non/existing/../path/" ) == "/root/non/path" );

            CHECK( cleanPath( "non/existing/../path/" ) == "non/path" );
            CHECK( cleanPath( my_home+"non/existing/../path/" ) == my_home+"non/path" );
            CHECK( cleanPath( "~/non/existing/../path/" ) == my_home+"non/path" );
            CHECK( cleanPath( "~root/non/existing/../path/" ) == "/root/non/path" ); // subdirs should not be resolve due to permissions

        }

    }
}
