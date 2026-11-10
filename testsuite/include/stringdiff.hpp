#ifndef INCLUDE_STRINGDIFF_HPP
#define INCLUDE_STRINGDIFF_HPP

#include <string>

#include "catch.hpp"

namespace testsuite {

    namespace file {

        std::string describeStringMismatch( const std::string& result, const std::string& expected );

    }

}

#define CHECK_STRINGS_EQUAL(result, expected) \
    do { \
        INFO( testsuite::file::describeStringMismatch( (result), (expected) ) ); \
        CHECK( (result) == (expected) ); \
    } while(0)

#endif // INCLUDE_STRINGDIFF_HPP
