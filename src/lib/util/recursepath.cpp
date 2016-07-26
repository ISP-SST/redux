#include "redux/util/recursepath.hpp"

#include <iostream>

using namespace redux::util;
using namespace std;
namespace fs = boost::filesystem;


RecursePath::RecursePath( boost::filesystem::path& p, functype f, int sublevels ) : nSubLevels( sublevels ), callBack(f) {
    this->operator()( p );
}


RecursePath::RecursePath( const RecursePath& rhs, int sublevels ) : nSubLevels( sublevels ) {
    
    callBack = rhs.callBack;

}


RecursePath::RecursePath( const RecursePath& rhs ) {
    
    nSubLevels = rhs.nSubLevels-1;
    callBack = rhs.callBack;

}


void RecursePath::operator()( fs::directory_entry& p ) const {
    this->operator()( p.path() );
}


void RecursePath::operator()( const fs::path& p ) const {

    if( !callBack( p ) || ( nSubLevels < 0 ) ) {
        return;
    }

    try {
        fs::directory_iterator it( p );
        fs::directory_iterator end;
        for_each( it, end, RecursePath( *this, nSubLevels - 1 ) );
    }
    catch( fs::filesystem_error& fex ) {
        cerr << "Error in RecursePath at path \"" << p.string() << "\"" << endl;
        cerr << fex.what() << endl;
    }

}
