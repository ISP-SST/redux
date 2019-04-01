#include "redux/util/semaphore.hpp"

using namespace redux::util;
using namespace std;


Semaphore::Semaphore( unsigned int cnt ) : counter( cnt ), init(cnt) {
    
}


void Semaphore::get() {
    unique_lock<mutex> lock(mtx);
    cond.wait( lock, [&]()->bool{ return counter>0; } );
    --counter;
}


void Semaphore::release(void) {
    lock_guard<mutex> lock(mtx);
    ++counter;
    cond.notify_one();
}


unsigned int Semaphore::count( void ) {
    lock_guard<mutex> lock(mtx);
    return counter;
}


unsigned int Semaphore::getInit( void ) {
    lock_guard<mutex> lock(mtx);
    return init;
}


void Semaphore::decrease( unsigned int val ) {
    lock_guard<mutex> lock(mtx);
    if( val < init ) {
        init -= val;
        if( val < counter ) counter -= val;
        else counter = 0;
    } else {
        init = 0;
        counter = 0;
    }
}


void Semaphore::increase( unsigned int val ) {
    lock_guard<mutex> lock(mtx);
    init += val;
    counter += val;
}


void Semaphore::set( unsigned int val ) {
    lock_guard<mutex> lock(mtx);
    int64_t diff = static_cast<int64_t>(val)-static_cast<int64_t>(init);
    if( diff ) {
        init += diff;
        counter += diff;
    }
}


void Semaphore::reset( void ) {
    lock_guard<mutex> lock(mtx);
    counter = init;
}


string Semaphore::getStatus(void) const {
    unsigned int i(0),c(0);
    {
        lock_guard<mutex> lock(mtx);
        i = init;
        c = counter;
    }
    return to_string(c) + "/" + to_string(i);
}


Semaphore::Scope::Scope( Semaphore& s ) : sem(s), active(true) {
    sem.get();
}


Semaphore::Scope::~Scope() {
    release();
}


void Semaphore::Scope::release( void ) {
    if( active ) {
        sem.release();
    }
    active = false;
}
