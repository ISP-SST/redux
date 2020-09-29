#include "redux/util/datautil.hpp"

#include "redux/util/stringutil.hpp"

#include <iostream>

using namespace redux::util;
using namespace std;


namespace redux {
    namespace util {
        
        template <typename T>       // defualt for unsigned types is to just return the value
        T Abs( const T& v ) {
            if( !std::numeric_limits<T>::is_signed ) return v;
            return std::abs( static_cast<long long int>(v) );
        }
        template <>
        float Abs( const float& v ) {
            return fabs( v );
        }
        template <>
        double Abs( const double& v ) {
            return fabs( v );
        }
        template <>
        long double Abs( const long double& v ) {
            return fabs( v );
        }
    }
}
template int8_t redux::util::Abs<int8_t>( const int8_t& );
template uint8_t redux::util::Abs<uint8_t>( const uint8_t& );
template int16_t redux::util::Abs<int16_t>( const int16_t& );
template uint16_t redux::util::Abs<uint16_t>( const uint16_t& );
template int32_t redux::util::Abs<int32_t>( const int32_t& );
template uint32_t redux::util::Abs<uint32_t>( const uint32_t& );
template int64_t redux::util::Abs<int64_t>( const int64_t& );
template uint64_t redux::util::Abs<uint64_t>( const uint64_t& );


