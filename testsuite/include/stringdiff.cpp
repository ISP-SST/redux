#include "stringdiff.hpp"

#include <algorithm>
#include <sstream>

using namespace std;


string testsuite::file::describeStringMismatch( const string& result, const string& expected ) {

    if ( result == expected ) {
        return "";
    }

    ostringstream msg;

    size_t firstMismatch(0);
    size_t resSz = result.size();
    size_t expSz = expected.size();
    if( resSz < expSz ) {
        msg << "result is shorter than expected (" << resSz << " < " << expSz <<  ")\n";
        auto mp = std::mismatch( result.begin(), result.end(), expected.begin() );
        firstMismatch = mp.first - result.begin();
    } else {
        if( resSz > expSz) {
            msg << "result is longer than expected(" << resSz << " > " << expSz <<  ")\n";
        }
        auto mp = std::mismatch( expected.begin(), expected.end(), result.begin() );
        firstMismatch = mp.first - expected.begin();
    }

    msg << "result/expected differs at position " << firstMismatch << ":";
    if( firstMismatch < resSz ) {
        msg << "  result[" << firstMismatch << "] = '" << result[firstMismatch] << "' (" << (int)result[firstMismatch] << ")";
    }
    if( firstMismatch < expSz ) {
        msg << ",  expected[" << firstMismatch << "] = '" << expected[firstMismatch] << "' (" << (int)expected[firstMismatch] << ")";
    }

    return msg.str();

}
