#ifndef INCLUDE_TESTSUITE_HPP
#define INCLUDE_TESTSUITE_HPP

#include <string>

#include <boost/test/unit_test.hpp>

namespace btt = boost::test_tools;

namespace testsuite {
    
    namespace file {

        btt::predicate_result compare_strings( const std::string& result, const std::string& expected );
        
    }

}

#endif // INCLUDE_TESTSUITE_HPP
