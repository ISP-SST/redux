#ifdef DISABLED_FOR_NOW_BREAKS_OLDER_BOOST

#ifndef REDUX_FILE_DIRWATCH_HPP
#define REDUX_FILE_DIRWATCH_HPP

#include <atomic>
#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <sys/inotify.h>
#include <thread>


#include <boost/asio.hpp>
#include <boost/array.hpp>
#include <boost/bimap.hpp>

#include <iostream>

namespace redux {

    namespace file {
        
        struct DirEvent {
            DirEvent() : mask(0) { }
            DirEvent( const std::string &d, const std::string &f, uint32_t m ) : dirname(d), filename(f), mask(m) { }
            std::string dirname;
            std::string filename;
            uint32_t mask;
        };

        class DirWatch {
            
        public:
            DirWatch();
            DirWatch( const std::string &dirname, uint32_t mask=IN_CREATE|IN_DELETE|IN_MOVED_FROM|IN_MOVED_TO );
            
            void start(void);
            void stop(void);
            
            void add_directory( const std::string &dirname, uint32_t mask=IN_CREATE|IN_DELETE|IN_MOVED_FROM|IN_MOVED_TO );
            void remove_directory( const std::string &dirname );

            DirEvent popfront_event(boost::system::error_code &ec);
            void pushback_event(DirEvent ev);
            
            template <typename CallBack>
            void async_monitor( CallBack cb, bool repeat=false );
            DirEvent monitor( boost::system::error_code &ec ) { return popfront_event( ec ); };

        private:
            int init();
            void begin_read();
            void end_read(const boost::system::error_code &ec, std::size_t bytes_transferred);
            std::string get_dirname(int wd);

            int fd_;
            boost::asio::io_service service_;
            boost::asio::posix::stream_descriptor stream_descriptor_;
            std::unique_ptr<boost::asio::io_service::work> work_;
            boost::array<char, 4096> read_buffer_;
            std::string pending_read_buffer_;
            std::mutex watch_descriptors_mutex_;
            typedef boost::bimap<int, std::string> watchList_t;
            watchList_t watchList_;
            std::mutex event_mutex_;
            std::condition_variable events_cond_;
            std::atomic<bool> running_;
            std::deque<DirEvent> events_;
            std::vector<std::thread> threads_;
            
            template<typename CallBack> friend class EventHandler;

        };

        template <typename CallBack>
        class EventHandler {
        public:
            EventHandler( DirWatch &dw, CallBack handler, bool repeat ) : dw_(dw), cb_(handler), repeat_(repeat) { }

            void operator()() const {

                boost::system::error_code ec;
                DirEvent ev = dw_.popfront_event(ec);
                cb_(ec, ev);
                if( repeat_ && dw_.running_ ) {
                    (*this)();
                }
            }

        private:
            DirWatch& dw_;
            CallBack cb_;
            bool repeat_;
        };
        
        template <typename CallBack>
        void DirWatch::async_monitor( CallBack cb, bool repeat ) {
            std::thread( EventHandler<CallBack>(*this, cb, repeat) ).detach();
        }

    }
    
}

#endif      // REDUX_FILE_DIRWATCH_HPP

#endif