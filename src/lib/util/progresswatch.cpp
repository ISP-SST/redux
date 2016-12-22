#include "redux/util/progresswatch.hpp"

#include "redux/util/stringutil.hpp"


#include <boost/format.hpp>

using namespace redux::util;
using namespace std;


ProgressWatch::ProgressWatch(int target, int start) : start_(start), counter_(start), target_(target), completed_(false) {
    
}


void ProgressWatch::clear(void) {
    
//    cout << __LINE__ << "  " << dump() << endl;
    unique_lock<mutex> lock(mtx);
    onChange = nullptr;
    onCompletion = nullptr;
    counter_ = start_ = target_ = 0;
    completed_ = false;
}


void ProgressWatch::reset(void) {
    
    int previous = counter_;
    counter_ = start_;
    if( onChange && (counter_ != previous)) onChange();

}


void ProgressWatch::test(void) {
    

    //cout << __LINE__ << "  " << dump() << endl;
    unique_lock<mutex> lock(mtx);
    bool notify = completed_ = ((start_ != target_) && (counter_ == target_));
    if(counter_ > target_) {
        lock.unlock();
        cout << __LINE__ << "  " << dump() << "  OVERFLOW !!!" << endl;
    }
    auto completer = onCompletion;
    if( notify ) {
        onCompletion = nullptr;
        counter_ = start_;
    }
    lock.unlock();

    if( notify ) {
        cv.notify_all();
        if( completer ) completer();
    }

}


void ProgressWatch::wait(void) {
    

//    cout << __LINE__ << "  " << dump() << endl;
    unique_lock<mutex> lock(mtx);
//    cout  << "ProgressWatch::wait()  " << __LINE__ << " completed = " << completed_ << endl;
    while( !completed_ ) {
//    cout  << "ProgressWatch::wait()  " << __LINE__ << " completed = " << completed_ << endl;
        cv.wait(lock);
    }
//    cout  << "ProgressWatch::wait()  DONE" << __LINE__ << endl;

}


// void ProgressWatch::wait_for(void) {
//     unique_lock<mutex> lock(mtx);
//     while ( counter_ != target_ ) {
//         cv.wait(lock);
//     }
// 
//   std::unique_lock<std::mutex> lock(mtx);
//   while (cv.wait_for(lock,std::chrono::seconds(1))==std::cv_status::timeout) {
//     std::cout << '.';
//   }
// 
// }


void ProgressWatch::set(int end, int start) {
    
    unique_lock<mutex> lock(mtx);
    target_ = end;
    counter_ = start_ = start;
    completed_ = false;
    lock.unlock();
//cout << __LINE__ << "  " << dump() << "  SET " << end << "," << start << endl;
    
}


void ProgressWatch::step(int step) {

    unique_lock<mutex> lock(mtx);
    counter_ += step;
    bool notify = completed_ = ((start_ != target_) && (counter_ == target_));
    if(counter_ > target_) {
        lock.unlock();
        cout << __LINE__ << "  " << dump() << "  OVERFLOW !!!" << endl;
    }
    auto ticker = onChange;
    auto completer = onCompletion;
    if( notify ) {
        onCompletion = nullptr;
        counter_ = start_;
    }
    lock.unlock();
    //cout << __LINE__ << "  " << dump() << endl;

    if( ticker ) ticker();
    if( notify ) {
        cv.notify_all();
        if( completer ) {
//    cout << __LINE__ << "  " << dump() << "   Calling completer." << endl;
            completer();
        } //else cout << __LINE__ << "  " << dump() << "  NO CB !!" << endl;
    }

}


void ProgressWatch::stepTarget(int step) {

    unique_lock<mutex> lock(mtx);
    target_ += step;
    lock.unlock();
//cout << __LINE__ << "  " << dump() << "  stepTarget " << step << endl;

}


float ProgressWatch::progress(void) {

    unique_lock<mutex> lock(mtx);
    return (counter_ - start_)*100.0/(target_ - start_);
   
}


string ProgressWatch::dump(void) {
    unique_lock<mutex> lock(mtx);
    return hexString(this) + " " + to_string(start_) + "/" + to_string(counter_) + "/" + to_string(target_);
}


string ProgressWatch::progressString(void) {

    //return to_string(counter_) +"/" + to_string(target_);
    return boost::str( boost::format("%03.1f%%") % progress() );
    
}