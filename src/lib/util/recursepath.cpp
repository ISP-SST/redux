#include "redux/util/recursepath.hpp"

#include <algorithm>
#include <iostream>

using namespace redux::util;
using namespace std;


RecursePath::RecursePath( bfs::path& p, functype f, int sublevels ) : nSubLevels( sublevels ), callBack(f) {
    this->operator()( p );
}


RecursePath::RecursePath( const RecursePath& rhs, int sublevels ) : nSubLevels( sublevels ) {
    
    callBack = rhs.callBack;

}


RecursePath::RecursePath( const RecursePath& rhs ) {
    
    nSubLevels = rhs.nSubLevels-1;
    callBack = rhs.callBack;

}


void RecursePath::operator()( bfs::directory_entry& p ) const {
    this->operator()( p.path() );
}


void RecursePath::operator()( const bfs::path& p ) const {

    if( !callBack( p ) || ( nSubLevels < 0 ) ) {
        return;
    }

    try {
        bfs::directory_iterator it( p );
        bfs::directory_iterator end;
        for_each( it, end, RecursePath( *this, nSubLevels - 1 ) );
    }
    catch( bfs::filesystem_error& fex ) {
        cerr << "Error in RecursePath at path \"" << p.string() << "\"" << endl;
        cerr << fex.what() << endl;
    }

}
