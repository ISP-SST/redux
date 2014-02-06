#include "redux/version.hpp"

#include "redux/revision.hpp"

#include <sstream>

using namespace redux;

int redux::getVersionNumber ( void ) {
    return reduxVersionMajor*1000000 + reduxVersionMinor*1000 + reduxVersionPatch;
}

std::string redux::getVersionString ( void ) {
    std::stringstream ss;
    ss << reduxVersionMajor << "." << reduxVersionMinor << "." << reduxVersionPatch;
    return ss.str();
}

std::string redux::getLongVersionString ( void ) {
    std::stringstream ss;
    ss << reduxVersionMajor << "." << reduxVersionMinor << "." << reduxVersionPatch;
    ss << " (commit: " << reduxVersionCommit << " - " << reduxCommitMessage << ")";
    return ss.str();
}
