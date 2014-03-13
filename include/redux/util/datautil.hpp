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

        /*!  Helper for swapping elements between 2 arrays.
         *   @details This is just an iterator which calls std::swap
         */
        template <class T>
        void swap(T* a, T* b, size_t n = 1) {
            if(a == b) {
                return;
            }
            while(n--) {
                std::swap(*a++, *b++);       // uses move semantics since c++11
            }
        }


        template <typename T>
        inline char* pack(char* ptr, const std::vector<T>& in) {
            size_t sz = in.size()*sizeof(T);
            *reinterpret_cast<size_t*>(ptr) = sz;
            ptr+=sizeof(size_t);
            memcpy(ptr,reinterpret_cast<const char*>(in.data()),sz);
            return ptr+sz;
        }
        
        template <typename T>
        inline char* pack(char* ptr, const T* in, size_t count) {
            size_t sz = count*sizeof(T);
            memcpy(ptr, reinterpret_cast<const char*>(in), sz);
            return ptr+sz;
         }
        
        template <typename T>
        inline char* pack(char* ptr, const T& in) {
            *reinterpret_cast<T*>(ptr) = in;
            return ptr+sizeof(T);
        }
        
        template <>
        inline char* pack<std::string>(char* ptr, const std::string& in) {
            strcpy(ptr,in.c_str());
            return ptr + in.length() + 1;
        }
        
        template <typename T>
        inline const char* unpack(const char* ptr, std::vector<T>& out, bool swap_endian=false) {
            size_t sz = *reinterpret_cast<const size_t*>(ptr);
            ptr+=sizeof(size_t);
            if(swap_endian) swapEndian(sz);
            out.resize(sz/sizeof(T));
            memcpy(reinterpret_cast<char*>(out.data()),ptr,sz);
            if(swap_endian) {
                swapEndian(&out[0],out.size());
            }
            return ptr+sz;
        }
        
        template <typename T>
        inline const char* unpack(const char* ptr, T* out, size_t count, bool swap_endian=false) {
            size_t sz = count*sizeof(T);
            memcpy(reinterpret_cast<char*>(out),ptr,sz);
            if(swap_endian) {
                swapEndian(out,count);
            }
            return ptr+sz;
        }
        
        template <typename T>
        inline const char* unpack(const char* ptr, T& out, bool swap_endian=false) {
            out = *reinterpret_cast<const T*>(ptr);
            if(swap_endian) {
                swapEndian(out);
            }
            return ptr+sizeof(T);
        }
        
        template <>
        inline const char* unpack<std::string>(const char* ptr, std::string& out, bool) {
            out = std::string(ptr);
            return ptr + out.length() + 1;
        }
        
        /*! @} */


    }  // namespace util

}  // namespace redux



#endif // REDUX_UTIL_DATAUTIL_HPP
