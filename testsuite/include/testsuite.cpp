#include "testsuite.hpp"


using namespace std;


btt::predicate_result testsuite::file::compare_strings( const string& result, const string& expected ) {

    if ( result == expected ) {
        return true;
    }

    // strings are different, generate some useful message for the test-log
    btt::predicate_result res( false );
    
    size_t firstMismatch(0);
    size_t resSz = result.size();
    size_t expSz = expected.size();
    if( resSz < expSz ) {
        res.message() << "result is shorter than expected (" << resSz << " < " << expSz <<  ")\n";
        auto mp = std::mismatch( result.begin(), result.end(), expected.begin() );
        firstMismatch = mp.first - result.begin();
    } else {
        if( resSz > expSz) {
            res.message() << "result is longer than expected(" << resSz << " > " << expSz <<  ")\n";
        }
        auto mp = std::mismatch( expected.begin(), expected.end(), result.begin() );
        firstMismatch = mp.first - expected.begin();
    }

    res.message() << "result/expected differs at position " << firstMismatch << ":";
    if( firstMismatch < resSz ) {
        res.message() << "  result[" << firstMismatch << "] = '" << result[firstMismatch] << "' (" << (int)result[firstMismatch] << ")";
    }
    if( firstMismatch < expSz ) {
        res.message() << ",  expected[" << firstMismatch << "] = '" << expected[firstMismatch] << "' (" << (int)expected[firstMismatch] << ")";
    }
    
    //res.message() << "\n      result=\"" << result << "\"\n    expected=\"" << expected << "\"\n";
    
    return res;
    
}
