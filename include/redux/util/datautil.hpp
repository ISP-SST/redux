#ifndef REDUX_UTIL_DATAUTIL_HPP
#define REDUX_UTIL_DATAUTIL_HPP

#include "redux/util/endian.hpp"

#include <cstring>
#include <string>
#include <vector>

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
        inline char* pack(char* ptr, const std::vector<T>& data) {
            size_t sz = data.size()*sizeof(T);
            memcpy(ptr,reinterpret_cast<char*>(data.front()),sz);
            return ptr+sz;
        }
        
        template <typename T>
        inline char* pack(char* ptr, const T* data, size_t count) {
            size_t sz = count*sizeof(T);
            memcpy(ptr, reinterpret_cast<const char*>(data), sz);
            return ptr+sz;
         }
        
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
        inline const char* unpack(const char* ptr, std::vector<T>& data, bool swap_endian=false) {
            size_t sz = data.size()*sizeof(T);
            memcpy(reinterpret_cast<char*>(data.front()),ptr,sz);
            if(swap_endian) {
                swapEndian(&data[0],data.size());
            }
            return ptr+sz;
        }
        
        template <typename T>
        inline const char* unpack(const char* ptr, T* data, size_t count, bool swap_endian=false) {
            size_t sz = count*sizeof(T);
            memcpy(reinterpret_cast<char*>(data),ptr,sz);
            if(swap_endian) {
                swapEndian(data,count);
            }
            return ptr+sz;
        }
        
        template <typename T>
        inline const char* unpack(const char* ptr, T& data, bool swap_endian=false) {
            data = *reinterpret_cast<const T*>(ptr);
            if(swap_endian) {
                swapEndian(data);
            }
            return ptr+sizeof(T);
        }
        
        template <>
        inline const char* unpack<std::string>(const char* ptr, std::string& data, bool) {
            data = std::string(ptr);
            return ptr + data.length() + 1;
        }
        
        /*! @} */


    }  // namespace util

}  // namespace sst



#endif // REDUX_UTIL_DATAUTIL_HPP
