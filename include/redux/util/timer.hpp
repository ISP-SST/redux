#ifndef REDUX_UTIL_TIMER_HPP
#define REDUX_UTIL_TIMER_HPP

#include <cstdint>

#include <functional>
#include <chrono>
#include <future>
#include <cstdio>
#include <iostream>

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */

        /*!  @file      timer.hpp
         *   @brief     Simple wrapper for delayed function calls
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2014
         */

        class Timer {
            
            template <class F, typename T>
            void call(F& f, T& to) { f(); if(--m_Count) post(f,to);  }
            template <class F, typename T>
            void post(F& f, T& to) {
                 if(m_Async) {
                    std::thread([this,&f,&to]() {
                        std::this_thread::sleep_for(to);
                        call(f,to);
                    }).detach();
                } else {
                    std::this_thread::sleep_for(to);
                    call(f,to);
                }
               
            }
        public:
            template <class F, typename T>
            Timer(F f, T to, int64_t count=1, bool async=false ) : m_Count(count), m_Async(async) { if(m_Count) post(f,to); }
            template <class F, typename T>
            void post(F f, T to, bool async, int64_t count=1) {
                m_Count = count;
                m_Async = async;
                 if(m_Async) {
                    std::thread([this,&f,&to]() {
                        std::this_thread::sleep_for(to);
                        call(f,to);
                    }).detach();
                }
                else {
                    std::this_thread::sleep_for(to);
                    call(f,to);
                }
               
            }
        private:
            int64_t m_Count;
            bool m_Async;
            std::chrono::duration<size_t,std::micro> m_Timeout;

        };

        /*! @} */

    }

}

#endif // REDUX_UTIL_TIMER_HPP
