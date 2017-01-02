#ifndef REDUX_UTIL_PROGRESSWATCH_HPP
#define REDUX_UTIL_PROGRESSWATCH_HPP

#include <atomic>
#include <functional>
#include <string>

#include <boost/thread.hpp>
#include <boost/thread/condition_variable.hpp>

namespace redux {

    namespace util {
        
        class ProgressWatch {
            
        public:
            ProgressWatch(int target=0, int start=0);

            void clear(void);
            void reset(void);
            void test(void);
            inline bool verify(void) { boost::lock_guard<boost::mutex> lock(mtx); return completed_; }
            void wait(void);
            
            void set(int target, int start=0);
            inline void setStart(int s) { boost::lock_guard<boost::mutex> lock(mtx); start_ = s; }
            inline void setTarget(int t) { boost::lock_guard<boost::mutex> lock(mtx); target_ = t; }
            
            void step(int step=1);
            void stepTarget(int step=1);
            inline void increase(int count=1) { step(count); };
            inline void decrease(int count=1) { step(-count); };
            inline void increaseTarget(int count=1) { stepTarget(count); }
            inline void decreaseTarget(int count=1) { stepTarget(-count); }
            
            inline void setHandler( std::function<void(void)> cb ) { boost::lock_guard<boost::mutex> lock(mtx); onCompletion = cb; }
            inline void setTicker( std::function<void(void)> cb ) { boost::lock_guard<boost::mutex> lock(mtx); onChange = cb; }
            
            // prefix
            inline ProgressWatch& operator++() { step(); return *this; }
            inline ProgressWatch& operator--() { step(-1); return *this; }
            // postfix
            //inline ProgressWatch operator++(int) { ProgressWatch tmp = *this; increase(); return tmp; }
            //inline ProgressWatch operator--(int) { ProgressWatch tmp = *this; decrease(); return tmp; }
            
            float progress(void);
            std::string dump(void);
            std::string progressString(void);
            
        private:
            int start_;
            int counter_;
            int target_;
            bool completed_;
            boost::mutex mtx;
            boost::condition_variable cv;
            std::function<void(void)> onChange;
            std::function<void(void)> onCompletion;
            
        };
    

    }  // end namespace util

}  // end namespace redux


#endif // REDUX_UTIL_PROGRESSWATCH_HPP
