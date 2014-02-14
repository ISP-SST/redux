#include "redux/version.hpp"

#include "redux/revision.hpp"

#include <sstream>

using namespace redux;
using namespace std;

int redux::getVersionNumber ( void ) {
    return (reduxVersionMajor<<16) + (reduxVersionMinor<<8) + reduxVersionPatch;
}

std::string redux::getVersionString ( void ) {
    return to_string(reduxVersionMajor) + "." + to_string(reduxVersionMinor) + "." + to_string(reduxVersionPatch);
}

std::string redux::getVersionString ( int v ) {
    return to_string((v>>16)&0xFF) + "." + to_string((v>>8)&0xFF) + "." + to_string(v&0xFF);
}

std::string redux::getLongVersionString ( void ) {
    std::stringstream ss;
    ss << reduxVersionMajor << "." << reduxVersionMinor << "." << reduxVersionPatch;
    ss << " (commit: " << reduxVersionCommit << " - " << reduxCommitMessage << ")";
    return ss.str();
}
