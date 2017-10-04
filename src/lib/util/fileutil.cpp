#include "redux/util/fileutil.hpp"

#include <exception>
#include <iostream>

#include <boost/filesystem.hpp>
//#include <boost/filesystem/operations.hpp>

using namespace std;
namespace bfs = boost::filesystem;

size_t redux::util::getDirSize( string dir ) {
    size_t dirSize(0);
    bfs::path dirPath(dir);
    if( bfs::exists(dirPath) ) {
        bfs::directory_iterator end_itr;
        for ( bfs::directory_iterator dirIte(dir); dirIte != end_itr; ++dirIte ) {
            bfs::path filePath( dirIte->path() );
            try {
                if ( bfs::is_directory(dirIte->status()) ) {
                    dirSize += getDirSize( filePath.string() );
                } else {
                    dirSize += bfs::file_size( filePath );
                } 
            } catch( exception& ) {
                // if size-check failed, just ignore this item and continue
            }
        }
    }
    return dirSize;
}
