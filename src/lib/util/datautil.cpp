#include "redux/util/datautil.hpp"

#include "redux/momfbd/modes.hpp"
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


uint64_t redux::util::size( const momfbd::ModeList& data ) {
    static uint16_t mSize = sizeof(uint16_t) + sizeof(momfbd::ModeBase);    // modeID & modeType
    uint64_t sz = sizeof(uint64_t) + sizeof(momfbd::ModeBase);              // nModes & defaultType
    sz += data.size() * mSize;
    return sz;
}


uint64_t redux::util::pack( char* ptr, const momfbd::ModeList& in ) {
    uint64_t count = in.size();
    memcpy( ptr, reinterpret_cast<const char*>(&count), sizeof(uint64_t) );
    uint64_t totalSize = sizeof(uint64_t);
    *reinterpret_cast<momfbd::ModeBase*>(ptr+totalSize) = in.defaultType;
    totalSize += sizeof(momfbd::ModeBase);
    memcpy( ptr, reinterpret_cast<const char*>(&count), sizeof(uint64_t) );
    for( auto &m: in ) {
        *reinterpret_cast<momfbd::ModeBase*>(ptr+totalSize) = m.type;
        totalSize += sizeof(momfbd::ModeBase);
        *reinterpret_cast<uint16_t*>(ptr+totalSize) = m.mode;
        totalSize += sizeof(uint16_t);
    }
    return totalSize;
}
    
uint64_t redux::util::unpack( const char* ptr, redux::momfbd::ModeList& out, bool swap_endian ) {
    uint64_t count(0);
    memcpy( reinterpret_cast<char*>(&count), ptr, sizeof(uint64_t) );
    if( swap_endian ) swapEndian(count);
    out.resize(count);
    uint64_t totalSize = sizeof(uint64_t);
    out.defaultType = *reinterpret_cast<const momfbd::ModeBase*>(ptr+totalSize);
    totalSize += sizeof(momfbd::ModeBase);
    for( auto& m: out ) {
            m.type = *reinterpret_cast<const momfbd::ModeBase*>(ptr+totalSize);
            totalSize += sizeof(momfbd::ModeBase);
            m.mode = *reinterpret_cast<const uint16_t*>(ptr+totalSize);
            totalSize += sizeof(uint16_t);
    }
    return totalSize;
}
    
