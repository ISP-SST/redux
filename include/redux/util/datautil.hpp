#ifndef REDUX_UTIL_DATAUTIL_HPP
#define REDUX_UTIL_DATAUTIL_HPP

#include <cstring>
#include <string>

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */

        /*!  @file  datautil.hpp
         *   @brief     Collection of functions for memory/variable manipulation
         *   @details   Functions convenient for manipulation of variables, endianswapping, bit-flipping etc.
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         */

        template <typename T>
        inline char* pack(char* ptr, const T& data) {
            *reinterpret_cast<T*>(ptr) = data;
            return ptr+sizeof(T);
        }
        
        template <>
        inline char* pack<std::string>(char* ptr, const std::string& data) {
            strcpy(ptr,data.c_str());
            return ptr + data.length() + 1;
        }
        
        template <typename T>
        inline const char* unpack(const char* ptr, T& data) {
            data = *reinterpret_cast<const T*>(ptr);
            return ptr+sizeof(T);
        }
        
        template <>
        inline const char* unpack<std::string>(const char* ptr, std::string& data) {
            data = std::string(ptr);
            return ptr + data.length() + 1;
        }
        
        /*! @} */


    }  // namespace util

}  // namespace sst



#endif // REDUX_UTIL_DATAUTIL_HPP
