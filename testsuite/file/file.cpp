#include "redux/file/fileio.hpp"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace redux::file;
using namespace std;

#ifndef RDX_TESTDATA_DIR
#error RDX_TESTDATA_DIR not set
#endif


namespace testsuite {
    namespace file {
        
        void path_test( void ) {
            
            BOOST_CHECK( isRelative("abc") );
            BOOST_CHECK( !isRelative("/abc") );
            
            string tmp;
            bfs::path test_dir(RDX_TESTDATA_DIR);
            
            string my_home = getenv("HOME");
            string username = getenv("USER");
            if( !my_home.empty() ) {
                if( *my_home.rbegin() != '/' ) my_home.push_back('/');    // getHome() returns with slash, so append if needed.
                tmp = getHome();
                BOOST_TEST( tmp == my_home );
                if( !username.empty() ) {
                    BOOST_TEST( getHome(username) == my_home );
                }
            } else my_home = getHome();
            
            BOOST_TEST( getHome("root") == "/root/" );
            BOOST_TEST( getHome("unknown_user_zxy") == "" );
            
            BOOST_CHECK( bfs::exists(test_dir) );
            BOOST_CHECK( bfs::is_directory(test_dir) );
            
            //BOOST_TEST( weaklyCanonical(test_dir / "non/existing/../path/") == (test_dir/"non/path") );   // TODO test with resolution of known symlink
            BOOST_TEST( weaklyCanonical( "~root/non/existing/../path/" ) == "/root/non/path" );

            BOOST_TEST( cleanPath( "non/existing/../path/" ) == "non/path" );
            BOOST_TEST( cleanPath( my_home+"non/existing/../path/" ) == my_home+"non/path" );
            BOOST_TEST( cleanPath( "~/non/existing/../path/" ) == my_home+"non/path" );
            BOOST_TEST( cleanPath( "~root/non/existing/../path/" ) == "/root/non/path" ); // subdirs should not be resolve due to permissions
            
        }
        
        void add_ana_tests( test_suite* ts );       // defined in ana.cpp
        void add_fits_tests( test_suite* ts );      // defined in fits.cpp

        void add_tests( test_suite* ts ) {
            
            ts->add( BOOST_TEST_CASE_NAME( &path_test, "Paths/Filenames" ) );

            add_ana_tests( ts );
            add_fits_tests( ts );

        }

    }
}
