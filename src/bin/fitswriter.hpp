#ifndef REDUX_FITSWRITER_HPP
#define REDUX_FITSWRITER_HPP


#include "redux/file/filefits.hpp"

//#include "frame.hpp"

#include <atomic>
#include <condition_variable>
// #include <cstdint>
#include <list>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>
namespace bpx = boost::posix_time;

namespace redux {
    
    namespace file {

    class FitsWriter {

    public:
        FitsWriter( const std::string& fn, int nThreads=1, bool compress=false );
        ~FitsWriter();

        void reset( int, int, bool );
        void makeHdr( void );
        
        template <typename T>
        void save( const redux::util::Array<T>& data, std::shared_ptr<Fits> hdr );
        
        void wait( void );
        
        void save_meta( const std::vector<std::string>& m ) { extra_meta = m; };

        void updateCard( std::string key, std::string card ) { if( hdr ) Fits::updateCard( hdr->primaryHDU.cards, key, card ); };
        void updateCard( std::string card ) { if( hdr ) Fits::updateCard( hdr->primaryHDU.cards, card ); };
        
        static void clear_gmeta(void) { std::lock_guard<std::mutex> lock(globalMtx); globalMeta.clear(); };
        static void update_gmeta( std::string card ) { std::lock_guard<std::mutex> lock(globalMtx); Fits::updateCard( globalMeta, card ); };
        template <typename T>
        static void update_gmeta( std::string key, T value, std::string comment="" ) {
            std::lock_guard<std::mutex> lock(globalMtx); Fits::updateCard( globalMeta, key, Fits::makeCard(key, value, comment) );
        }

        static std::string get_saves( void );
        static void clear_saves( void ) { std::lock_guard<std::mutex> lock(globalMtx); saves.clear(); };

    private:
        
        static std::shared_ptr<uint8_t> get_buf( size_t N );
        static void return_buf( std::shared_ptr<uint8_t>& );
        static void clear_bufs( void );
        
        static void add_save( boost::posix_time::ptime from, boost::posix_time::ptime to );
        
        void open_file( const std::string& filename, size_t sz=0 );
        void close_file( void );
        void write_exptime_table( void );
       
        void thread_run(void);
        void get_range( bpx::ptime&, bpx::ptime& );

        void push( const uint8_t*, const bpx::ptime& );
        std::shared_ptr<uint8_t> pop();

        bool running;
        bool do_fsync;
        bool do_compress;
        bool do_write;
        int index;
        int nframes;
        int nthreads;
        std::atomic<int> activeThreads;
        int bytes_per_pixel;
        int npixels;
        int fd;
        int height, width;

        int pcount, maxRowSize;
        size_t hdrEnd, dataStart;
        std::vector<std::string> extra_meta;
        std::shared_ptr<Fits> hdr;
        std::list<std::shared_ptr<uint8_t>> frames;
        std::list<boost::posix_time::ptime> times;
        std::vector<std::thread> threads;
        std::vector<int> offsets;
        std::mutex queueMtx;//, writeMtx;
        std::condition_variable cond;
        
        std::string filename;
        size_t frame_count;
        
        boost::posix_time::ptime global_first, global_last;

        static int activeCount, totalCount;
        static std::list<std::shared_ptr<uint8_t>> buffers;
        static std::vector<std::string> globalMeta;
        static std::mutex globalMtx, bufMtx;
        static std::set<std::pair<boost::posix_time::ptime,boost::posix_time::ptime>> saves;
        
    };

    }   // file

}   // redux


#endif      // REDUX_FITSWRITER_HPP
