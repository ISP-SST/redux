#ifndef REDUX_UTIL_PROGRESSWATCH_HPP
#define REDUX_UTIL_PROGRESSWATCH_HPP

#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <string>


namespace redux {

    namespace util {
        
        class ProgressWatch {
            
        public:
            ProgressWatch(int target=0, int start=0);

            void clear(void);
            void reset(void);
            void test(void);
            inline bool verify(void) const { return (counter_ == target_); }
            void wait(void);
            
            void set(int target, int start=0);
            inline void setTarget(int t) { target_ = t; }
            inline void increaseTarget(int t=1) { target_ += t; }
            inline void decreaseTarget(int t=1) { target_ -= t; }
            inline void setStart(int s) { start_ = s; }
            inline void tick(void) const { if(onChange) onChange(); }
            
            inline void setHandler( std::function<void(void)> cb ) { onCompletion = cb; }
            inline void setTicker( std::function<void(void)> cb ) { onChange = cb; }
            
            ProgressWatch& operator--();
            ProgressWatch& operator++();
            
            float progress(void);
            std::string progressString(void);
            
        private:
            int start_;
            std::atomic<int> counter_;
            std::atomic<int> target_;
            std::mutex mtx;
            std::condition_variable cv;
            std::function<void(void)> onChange;
            std::function<void(void)> onCompletion;
            
        };
    

    }  // end namespace util

}  // end namespace redux


#endif // REDUX_UTIL_PROGRESSWATCH_HPP
