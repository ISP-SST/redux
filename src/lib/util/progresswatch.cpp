#include "redux/util/progresswatch.hpp"

#include "redux/util/stringutil.hpp"


#include <boost/format.hpp>

using namespace redux::util;
using namespace std;


ProgressWatch::ProgressWatch(int target, int start) : start_(start), counter_(start), target_(target) {
    
}


void ProgressWatch::clear(void) {
    
    onChange = nullptr;
    onCompletion = nullptr;
    counter_ = start_;

}


void ProgressWatch::reset(void) {
    
    int previous = counter_;
    counter_ = start_;
    if( onChange && (counter_ != previous)) onChange();

}


void ProgressWatch::test(void) {
    
    int tmp = target_;
    if( counter_.compare_exchange_strong( tmp, tmp ) ) {
        cv.notify_all();
        if( onCompletion ) {
            onCompletion();
        }
    }

}


void ProgressWatch::wait(void) {
    unique_lock<mutex> lck(mtx);
    while ( counter_ != target_ ) {
        cv.wait(lck);
    }
}


void ProgressWatch::set(int end, int start) {
    
    target_ = end;
    start_ = start;
    reset();
    test();
    
}


ProgressWatch& ProgressWatch::operator--() {

    --counter_;
    if( onChange ) onChange();
    test();
    return *this;
    
}


ProgressWatch& ProgressWatch::operator++() {

    ++counter_;
    if( onChange ) onChange();
    test();
    return *this;

}


float ProgressWatch::progress(void) {

    int cnt = counter_.load();
    int tmp = target_.load();
    return (cnt - start_)*100.0/(tmp - start_);
   
}


string ProgressWatch::progressString(void) {

    return boost::str( boost::format("%03.1f%%") % progress() );
    
}