#include "redux/version.hpp"

#include "redux/revision.hpp"

using namespace redux;
using namespace std;

int redux::getVersionNumber ( void ) {
    return (reduxVersionMajor<<16) + (reduxVersionMinor<<8) + reduxVersionPatch;
}

string redux::getVersionString ( void ) {
    return to_string(reduxVersionMajor) + "." + to_string(reduxVersionMinor) + "." + to_string(reduxVersionPatch);
}

string redux::getVersionString ( int v ) {
    return to_string((v>>16)&0xFF) + "." + to_string((v>>8)&0xFF) + "." + to_string(v&0xFF);
}

string redux::getLongVersionString ( void ) {
    
    string ret = to_string(reduxVersionMajor) + "." + to_string(reduxVersionMinor) + "." + to_string(reduxVersionPatch);
    ret += " (commit: " + to_string(reduxVersionCommit) + " - " + string(reduxCommitMessage) + ")";
    return ret;
}
