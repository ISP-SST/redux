
#include <redux/file/dirwatch.hpp>
#include <redux/util/stringutil.hpp>

#include <iostream>

#include <boost/filesystem.hpp>

using namespace redux::file;
using namespace redux::util;
using namespace std;

void dw_callback( const boost::system::error_code &ec, const DirEvent &ev ) {
    
    std::cout << "DirEvent:  file: " << ev.dirname << "/" << ev.filename << "  mask = " << bitString(ev.mask) << std::endl;
    
}

#define TEST_DIR1 "_DummyDir_" 

int main() {
    

    boost::filesystem::create_directory(TEST_DIR1);
    DirWatch dw( TEST_DIR1, IN_ALL_EVENTS); //|IN_CREATE|IN_DELETE|IN_MOVED_FROM|IN_MOVED_TO );
    //DirWatch dw( TEST_DIR1, IN_CLOSE_WRITE ); //|IN_CREATE|IN_DELETE|IN_MOVED_FROM|IN_MOVED_TO );

    dw.async_monitor( dw_callback, true );
    
/*
    int cnt(0);
    boost::system::error_code ec;
    while(cnt++ < 100) {
        DirEvent ev = dw.monitor(ec);
        cout << cnt << "  OneShot: " << flush;
        dw_callback(ec,ev);
    }
*/
    
    
    std::cin.get();
    dw.stop();

    boost::filesystem::remove_all(TEST_DIR1);
    
}
