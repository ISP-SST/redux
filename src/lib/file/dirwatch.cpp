#ifdef DISABLED_FOR_NOW_BREAKS_OLDER_BOOST
#include "redux/file/dirwatch.hpp"

//#include <errno.h>
#include <iostream>

using namespace redux::file;
using namespace std;


DirWatch::DirWatch() : fd_( init() ), stream_descriptor_( service_, fd_ ), work_(),
        running_(false) {

    start();
    
}


DirWatch::DirWatch( const std::string &dirname, uint32_t mask ) : fd_( init() ), stream_descriptor_( service_, fd_ ), work_(),
        running_(false) {
            
    add_directory( dirname, mask );
    start();
    
}


void DirWatch::start( void ) {

    work_.reset( new boost::asio::io_service::work(service_) );
    begin_read();
    threads_.clear();
    for( int i = 0; i < 1; ++i ) {  // TODO multithreaded (use strand for reading/queueing)
            threads_.push_back( thread( boost::bind( &boost::asio::io_service::run, &service_) ) );
    }

    running_ = true;
}


void DirWatch::stop( void ) {

    work_.reset();
    service_.stop();
    for( auto & t : threads_ ) {
        t.join();
    }

    unique_lock<mutex> lock(event_mutex_);
    running_ = false;
    events_cond_.notify_all();
    
}


void DirWatch::add_directory( const std::string &dirname, uint32_t mask ) {
    
    int wd = inotify_add_watch(fd_, dirname.c_str(), mask);
    
    if (wd == -1) {
        boost::system::system_error err( boost::system::error_code(errno, boost::system::get_system_category()),
                                         "DirWatch::add_directory: inotify_add_watch failed.");
        boost::throw_exception(err);
    }

    unique_lock<mutex> lock(watch_descriptors_mutex_);
    watchList_.insert(watchList_t::value_type(wd, dirname));
    
    std::cout << "DirWatch::add_directory(): " << dirname << "   fd_ = " << fd_ << "   wd = " << wd << "    nWatched = " << watchList_.size() << std::endl;
    
}

void DirWatch::remove_directory( const std::string &dirname ) {
    
    unique_lock<mutex> lock(watch_descriptors_mutex_);
    watchList_t::right_map::iterator it = watchList_.right.find(dirname);
    if (it != watchList_.right.end()) {
        int ret = inotify_rm_watch(fd_, it->second);
        watchList_.right.erase(it);
        if( ret == -1 ) {
            boost::system::system_error err( boost::system::error_code(errno, boost::system::get_system_category()),
                                             "DirWatch::remove_directory: inotify_rm_watch failed.");
            boost::throw_exception(err);
        }
    }
    
    std::cout << "DirWatch::remove_directory(): " << dirname << "    nWatched = " << watchList_.size() << std::endl;
}


DirEvent DirWatch::popfront_event( boost::system::error_code &ec ) {
    
    unique_lock<mutex> lock(event_mutex_);
    while ( running_ && events_.empty() ) {
        events_cond_.wait(lock);
    }
    
    DirEvent ev;
    if ( !events_.empty() ) {
        ec = boost::system::error_code();
        ev = events_.front();
        events_.pop_front();
    } else {
        ec = boost::asio::error::operation_aborted;
    }
    return ev;
    
}


void DirWatch::pushback_event( DirEvent ev ) {
    unique_lock<mutex> lock(event_mutex_);
    if( running_ ) {
        events_.push_back(ev);
        events_cond_.notify_all();
    }
}


int DirWatch::init(void) {
    
    int fd = inotify_init();
    if (fd == -1) {
        boost::system::system_error e(boost::system::error_code(errno, boost::system::get_system_category()), "boost::asio::DirWatch::init_fd: init_inotify failed");
        boost::throw_exception(e);
    }
    return fd;
    
}


void DirWatch::begin_read(void) {

    stream_descriptor_.async_read_some( boost::asio::buffer(read_buffer_),
        boost::bind(&DirWatch::end_read, this,
        boost::asio::placeholders::error, boost::asio::placeholders::bytes_transferred));
    
}

#include "redux/util/stringutil.hpp"
using namespace redux::util;
void DirWatch::end_read( const boost::system::error_code &ec, std::size_t bytes_transferred ) {

    if (!ec) {
        pending_read_buffer_ += std::string(read_buffer_.data(), bytes_transferred);
        while (pending_read_buffer_.size() > sizeof(inotify_event)) {
            const inotify_event *iev = reinterpret_cast<const inotify_event*>(pending_read_buffer_.data());
            pushback_event( DirEvent(get_dirname(iev->wd), iev->name, iev->mask) );
            pending_read_buffer_.erase(0, sizeof(inotify_event) + iev->len);
        }
        begin_read();
    } else if ( ec != boost::asio::error::operation_aborted ) {
        boost::system::system_error e(ec);
        boost::throw_exception(e);
    }
}


std::string DirWatch::get_dirname(int wd) {
    unique_lock<mutex> lock(watch_descriptors_mutex_);
    watchList_t::left_map::iterator it = watchList_.left.find(wd);
    return it != watchList_.left.end() ? it->second : "";
}

#endif
