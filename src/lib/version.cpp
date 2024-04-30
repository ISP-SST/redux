#include "redux/version.hpp"

#include "redux/revision.hpp"


using namespace redux;
using namespace std;

uint64_t redux::getVersionNumber ( void ) {
    return (uint64_t(reduxVersionMajor)<<32) + (reduxVersionMinor<<24) + (reduxVersionPatch<<16) + reduxVersionCommit;
}

string redux::getVersionString ( void ) {
    return to_string(reduxVersionMajor) + "." + to_string(reduxVersionMinor) + "." + to_string(reduxVersionPatch) + "-" + to_string(reduxVersionCommit);
}

string redux::getVersionString ( uint64_t v ) {
    return to_string((v>>32)&0xFF) + "." + to_string((v>>24)&0xFF) + "." + to_string((v>>16)&0xFF) + "-" + to_string(v&0xFFFF);
}

string redux::getLongVersionString ( bool includeMessage ) {
    
    string ret = to_string(reduxVersionMajor) + "." + to_string(reduxVersionMinor)
         + "." + to_string(reduxVersionPatch)+ "-" + to_string(reduxVersionCommit);
#ifdef DEBUG_
    ret += " dbg";
#endif
    if( includeMessage ) ret += " (" + string(reduxCommitMessage) + ")";
    return ret;
}
