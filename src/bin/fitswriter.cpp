#include "fitswriter.hpp"

#include "redux/util/ricecompress.hpp"

//#include "version.hpp"

#include <algorithm>
#include <atomic>
#include <functional>
#include <iostream>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <endian.h>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

using namespace redux;
using namespace redux::file;
using namespace redux::util;
using namespace std;

namespace bfs = boost::filesystem;


// Initialize static members
int FitsWriter::activeCount(0);
int FitsWriter::totalCount(0);
list<shared_ptr<uint8_t>> FitsWriter::buffers;
vector<string> FitsWriter::globalMeta;
mutex FitsWriter::globalMtx;
mutex FitsWriter::bufMtx;
set<pair<boost::posix_time::ptime,boost::posix_time::ptime>> FitsWriter::saves;

namespace {

    mutex writeMtx;

    // Pad file to a specified boundary
    void pad_to( int fd, size_t boundary, char fill=0 ) {
        size_t pos = lseek( fd, 0, SEEK_CUR) % boundary;
        if( pos ) {
            size_t pad = boundary - pos;
            unique_ptr<char[]> tmp( new char[pad] );
            memset( tmp.get(), 0, pad );
            std::fill_n( tmp.get(), pad, fill );
            size_t n = write( fd, tmp.get(), pad );
            if( n != pad ) {
                throw runtime_error("FitsWriter: Failed to write padding: " + string(strerror(errno)) );
            }
        }
    }

    void write_meta( int fd, const vector<string>& meta ) {
        
        if( meta.empty() ) return;
        
        size_t metaSize = ((meta.size()-1)/36 + 1)*2880;
        std::unique_ptr<char[]> hdu( new char[ metaSize] );
        char* hPtr = hdu.get();
        memset( hPtr, ' ', metaSize );
        for( auto& c: meta ) {
            std::copy( c.begin(), c.end(), hPtr );
            hPtr += 80;
        }
        size_t n = write( fd, hdu.get(), metaSize );
        if( n != metaSize ) {	// write offset+size for rows in binary table
            throw runtime_error("FitsWriter: write_meta failed: " + string(strerror(errno)) );
        }
    }

}


FitsWriter::FitsWriter( const string& fn, int nT, bool compress )
    : running(true), do_fsync(false), do_compress(compress), index(0), nthreads(nT), activeThreads(0),
      npixels(0), fd(-1), height(0), width(0), filename(fn), global_first(bpx::not_a_date_time), global_last(bpx::not_a_date_time) {

    bytes_per_pixel = 2;
    frame_count = 0;

    lock_guard<mutex> lock( globalMtx );
    ++totalCount;
}


FitsWriter::~FitsWriter() {

    if( do_write ) {
        if( !global_first.is_not_a_date_time() && !global_last.is_not_a_date_time() ) {
            add_save( global_first, global_last );
        }
    }
    
    lock_guard<mutex> lock( globalMtx );
    --totalCount;
    
}


void FitsWriter::makeHdr( void ) {
    
    if( !hdr ) {
        hdr.reset( new Fits() );
    }
    
    vector<string>& cards = hdr->primaryHDU.cards;
    cards.clear();
//     Fits::insertCard( cards, Fits::makeCard( "SIMPLE", true ) );
//     Fits::insertCard( cards, Fits::makeCard( "BITPIX", 16) );
//     if( do_compress ) {
//         Fits::updateCard( cards, Fits::makeCard( "NAXIS", 0 ) );
//     } else {
//         Fits::updateCard( cards, Fits::makeCard( "NAXIS", (times.size()>1)?3:2 ) );
//         Fits::updateCard( cards, Fits::makeCard( "NAXIS1", width ) );
//         Fits::updateCard( cards, Fits::makeCard( "NAXIS2", height ) );
//         if( times.size() > 1 ) {
//             Fits::updateCard( cards, Fits::makeCard( "NAXIS3", times.size() ) );
//         }
//     }
//     Fits::updateCard( cards, Fits::makeCard<string>( "EXTNAME", "Main" ) );
//     Fits::updateCard( cards, Fits::makeCard<string>( "TAB_HDUS", "TABULATIONS;DATE-BEG" ) );
//     
//     bpx::ptime first(bpx::not_a_date_time), last(bpx::not_a_date_time);
//     get_range( first, last );
// 
//     if( !first.is_not_a_date_time() && !last.is_not_a_date_time() ) {
//         string timestamp = to_iso_extended_string( first );
//         size_t pos = timestamp.find_last_of('.');
//         string timestamp_s = timestamp;
//         if( pos != string::npos ) timestamp_s = timestamp.substr(0,pos);
//         Fits::updateCard( cards, Fits::makeCard<string>( "DATE", timestamp ) );
//         Fits::updateCard( cards, Fits::makeCard<string>( "DATE-OBS", timestamp_s ) );
//         bpx::time_duration elapsed = (last - first);
//         first += elapsed/2;
//         Fits::updateCard( cards, Fits::makeCard( "DATE-AVG", to_iso_extended_string(first), "Average time of observations" ) );
//         Fits::updateCard( cards, Fits::makeCard( "DATE-END", to_iso_extended_string(last), "End time of observations" ) );
//     }
//     //Fits::insertCard( cards, Fits::makeCard<string>( "WFWFSVER", getVersionString() ) );
//     cout << "nGlobalMeta = " << globalMeta.size() << endl;
//     cout << "nExtraMeta = " << extra_meta.size() << endl;
    cards.insert( cards.end(), globalMeta.begin(), globalMeta.end() );      // copy global meta-data.
    cards.insert( cards.end(), extra_meta.begin(), extra_meta.end() );      // copy extra meta-data.

    Fits::removeCards( cards, "END" );                   // just in case it is not the last card, or if there are multiple.
    Fits::insertCard( cards, Fits::makeCard( "END" ) );

}

template <typename T>
void FitsWriter::save( const Array<T>& data, shared_ptr<Fits> hdr ) {

    size_t nD = data.nDimensions();
    if( nD < 2 || nD > 3 ) {
        throw runtime_error("FitsWriter::save only supports 2-3 dimensional data");
    }
    
    do_write = !filename.empty();
    do_compress = (do_compress && do_write && (bytes_per_pixel == 2) );
    nthreads = do_compress?nthreads:1;
    
    if( nD ) {
        nframes = data.dimSize(0);
        height = data.dimSize(1);
        width = data.dimSize(2);
    } else {
        nframes = 1;
        height = data.dimSize(0);
        width = data.dimSize(1);
    }
    
    npixels = height*width;
    
    unique_lock<mutex> wlock(writeMtx);
    index = 0;
    offsets.assign( 2*nframes, 0 );
    frames.clear();
    times.clear();
    
    vector<bpx::ptime> hdrTimes;
    if( hdr ) {
        hdrTimes = hdr->getStartTimes();
    }
    hdrTimes.resize( nframes, bpx::not_a_date_time );
    running = true;
    for( int i=0; i<nthreads; ++i ) {
        threads.push_back( thread(bind( &FitsWriter::thread_run, this )) );
    }
    //size_t blockSize = npixels*bytes_per_pixel;
    std::thread qthread([&](){
        for( int i(0); i<nframes; ++i ) {
            push( reinterpret_cast<const uint8_t*>(data.get()+i*npixels), hdrTimes[i] );
        }
    });

    if( do_write ) {
        
        open_file( filename, /*nframes*fq.frameSize+*/50*2880 );

        //bfs::path fnPath(filename);
        //Fits::updateCard( extra_meta, Fits::makeCard<string>( "FILENAME", fnPath.filename().string() ) );

        makeHdr();

        hdrEnd = ((hdr->primaryHDU.cards.size()-1)/36 + 1)*2880;		// end of primary header, and possibly start of compressed header
        dataStart = hdrEnd;
        //cout << " dataStart = " << dataStart << " nCards = " << hdr->primaryHDU.cards.size() << endl;
        if ( do_compress ) {
            dataStart += 2880+nframes*2*sizeof(int);
        }
        lseek( fd, dataStart, SEEK_SET );	// set file-pointer to where the frames should be saved.
    }
    wlock.unlock();

    qthread.join();
    wait();

    if( do_write ) {
        wlock.lock();
        close_file();
    }
    
}
template void FitsWriter::save( const Array<int16_t>&, shared_ptr<Fits> );


void FitsWriter::wait( void ) {

    try {
        {
            unique_lock<mutex> lock(queueMtx);
            while( !frames.empty() ) {
                cond.wait_for(lock,std::chrono::milliseconds(1));
            }
            running = false;
        }
        while( activeThreads ) {
            cond.notify_all();
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }

        for( auto& t: threads ) {
            if( t.joinable() ) {
                t.join();
            }
        }
        threads.clear();
    } catch ( const std::system_error& e ) {
        cerr << "FitsWriter::wait: exception caught: " << e.what() << endl;
        throw;
    }
}



void FitsWriter::thread_run(void) {

    activeThreads++;
    shared_ptr<uint8_t> f;
    shared_ptr<uint8_t> cData;
    size_t blockSize = npixels*bytes_per_pixel;
    if( do_compress ) {
        cData = get_buf( blockSize );		// compressed storage of the same size as a raw frame
    }

    while( running ) {
        if( (f = pop()) ) {
            try {
                unique_lock<mutex> lock(queueMtx);
                int ind = index;
                index += 2;
                lock.unlock();
                if( do_write ) {
                    if( do_compress ) {
                        int frameOffset(0);
                        int16_t* tmpPtr = reinterpret_cast<int16_t*>(f.get());
                        // TODO: set blocksize to width (i.e. process one row at a time) when/if the v. 4.0 FITS standard allows it.
                        int cSize = rice_comp16( tmpPtr, npixels, cData.get(), blockSize, 32 );
                        if( cSize > 0 ) {
                            lock_guard<mutex> wlock(writeMtx);
                            frameOffset = lseek( fd, 0, SEEK_CUR );
                            int count = write( fd, cData.get(), cSize );
                            if( count != cSize ) {
                                fprintf( stderr, "FitsWriter: write failed for frame #%d, count=%d  cSize=%d: %s\n", (ind/2), count, cSize, strerror(errno) );
                            }
                        } else {
                            fprintf( stderr, "FitsWriter: compression failed for frame #%d.\n", (ind/2) );
                            cSize = 0;
                        }
                        lock.lock();
                        offsets[ind] = cSize;
                        offsets[ind+1] = frameOffset;
                        lock.unlock();
                    } else {
                        if( bytes_per_pixel == 2 ) {
                            int16_t* ptr = reinterpret_cast<int16_t*>( f.get() );
                            for( int i=0; i<npixels; ++i ) ptr[i] = htobe16(ptr[i]);
                        } else if( bytes_per_pixel == 4 ) {
                            int32_t* ptr = reinterpret_cast<int32_t*>( f.get() );
                            for( int i=0; i<npixels; ++i ) ptr[i] = htobe32(ptr[i]);
                        }
                        {
                            lock_guard<mutex> wlock(writeMtx);
                            size_t count = write( fd, f.get(), blockSize );
                            if( count != blockSize ) {
                                fprintf( stderr, "FitsWriter: write failed for frame #%d, count=%zu  frameSize=%zu: %s\n", (ind/2), count, blockSize, strerror(errno) );
                            }
                        }
                    }
                }
            } catch ( exception& e ) {
                cerr << "Exception caught: " << e.what() << endl;
            } catch (... ) {
                cerr << "Unspecified exception caught." << endl;
            }
            return_buf( f );
        } else {
            unique_lock<mutex> lock(queueMtx);
            cond.wait(lock);
        }
    }

    if( cData ) {
        return_buf( cData );
    }
    
    activeThreads--;

}


void FitsWriter::get_range( bpx::ptime& first, bpx::ptime& last ) {

    first = last = bpx::not_a_date_time;
    lock_guard<mutex> lock(queueMtx);
    if( times.empty() ) return;
    
    first = last = times.front();
    for( auto& t: times ) {
        if( first.is_not_a_date_time() || (t < first) ) {
            first = t;
        }
        if( last.is_not_a_date_time() || (t > last) ) {
            last = t;
        }
    }
    
    if( global_first.is_not_a_date_time() || (first < global_first) ) {
        global_first = first;
    }
    if( global_last.is_not_a_date_time() || (last > global_last) ) {
        global_last = last;
    }

}


void FitsWriter::push( const uint8_t* data, const bpx::ptime& ts ) {
    if( data ) {
        shared_ptr<uint8_t> buf = get_buf( npixels*bytes_per_pixel );
        memcpy( buf.get(), data, npixels*bytes_per_pixel );
        {
            lock_guard<mutex> lock(queueMtx);
            frames.push_back( buf );
            times.push_back( ts );
        }
        cond.notify_all();
    }
}


shared_ptr<uint8_t> FitsWriter::pop() {
    lock_guard<mutex> lock(queueMtx);
    if( frames.empty() ) return nullptr;
    shared_ptr<uint8_t> ret = frames.front();
    frames.pop_front();
    return ret;
}


void FitsWriter::open_file( const string& filename, size_t sz ) {
    
    bfs::path fnPath(filename);
    bfs::path dirPath = fnPath.parent_path();
    maybeCreateDir( dirPath );

    fd = open( filename.c_str(), O_WRONLY|O_CREAT, 0664 ); // |O_EXCL
    if( fd < 0 ) {
        //connection->write(format("ERROR :Could not open output file '%s': %s", filename.c_str(), strerror(errno)));
        return;
    }
    
    // Allocate enough space for the whole file
    if( sz ) {
        posix_fallocate( fd, 0, sz );
    }
    
}


void FitsWriter::close_file( void ) {

    makeHdr();      // call makeHdr again to store the correct sizes.
    
    // Pad to a multiple of 2880 bytes, because FITS requires that
    pad_to( fd, 2880 );
    write_exptime_table();
    
    // Pad to a multiple of 2880 bytes, because FITS requires that
    pad_to( fd, 2880, ' ' );
    
    size_t pos = lseek( fd, 0, SEEK_CUR );		// save EOF position
    lseek( fd, 0, SEEK_SET );			        // move back to primary HDU
    
    write_meta( fd, hdr->primaryHDU.cards );

    if( do_compress ) {
        
        maxRowSize = 0;
        pcount = 0;
        
        int dataStart = INT32_MAX;
        for( size_t i=0; i<offsets.size(); i+=2 ) {
            if( offsets[i] > maxRowSize ) maxRowSize = offsets[i];
            if( offsets[i+1] && (offsets[i+1] < dataStart) ) dataStart = offsets[i+1];
            pcount += offsets[i];
        }
        
        for( size_t i=0; i<offsets.size(); i+=2 ) {
            offsets[i] = htobe32(offsets[i]);
            if( offsets[i+1] ) offsets[i+1] = htobe32(offsets[i+1]-dataStart);
        }
        
        
        const vector<string> comp_meta = {
            Fits::makeCard( "XTENSION", "BINTABLE", "binary table extension" ),
            Fits::makeCard( "BITPIX", 8, "8-bit data" ),
            Fits::makeCard( "NAXIS", 2, "2-dimensional binary table" ),
            Fits::makeCard( "NAXIS1", 8, "width of table in bytes" ),
            Fits::makeCard( "NAXIS2", nframes, "number of rows in table" ),
            Fits::makeCard( "PCOUNT", pcount, "size of special data area" ),
            Fits::makeCard( "GCOUNT", 1, "one data group (required)" ),
            Fits::makeCard( "TFIELDS", 1, "number of fields in each row" ),
            Fits::makeCard( "TTYPE1", "COMPRESSED_DATA", "label for field #1" ),
            Fits::makeCard( "TFORM1", "1PB("+to_string(maxRowSize)+")", "data format of field: variable length array" ),
            // start of mandatory keywords for tiled image compression (10.1.1 in v. 4.0 of the FITS standard)
            Fits::makeCard( "ZIMAGE", true, "extension contains compressed image" ),
            Fits::makeCard( "ZTILE1", width, "size of tiles" ),
            Fits::makeCard( "ZTILE2", height, "size of tiles" ),
            Fits::makeCard( "ZTILE3", 1, "size of tiles" ),
            Fits::makeCard( "ZCMPTYPE", "RICE_1", "compression algorithm used" ),
            Fits::makeCard( "ZNAME1", "BLOCKSIZE" ),
            Fits::makeCard( "ZVAL1", 32 ),
            Fits::makeCard( "ZNAME2", "BYTEPIX" ),
            Fits::makeCard( "ZVAL2", 2 ),
            Fits::makeCard( "ZSIMPLE", true ),
            Fits::makeCard( "ZBITPIX", 16, "bitpix of original image" ),
            Fits::makeCard( "ZNAXIS", 3, "naxis of original image" ),
            Fits::makeCard( "ZNAXIS1", width, "naxis1 of original image" ),
            Fits::makeCard( "ZNAXIS2", height, "naxis2 of original image" ),
            Fits::makeCard( "ZNAXIS3", nframes, "naxis3 of original image" ),
            Fits::makeCard( "END" )
        };

        lseek( fd, hdrEnd, SEEK_SET );		// move to end of primary header
        write_meta( fd, comp_meta );

        size_t sz = offsets.size()*sizeof(int);
        size_t n = write( fd, offsets.data(), sz );
        if( n != sz ) {	// write offset+size for rows in binary table
            cerr << "FitsWriter: Failed to write row information: " << strerror(errno) << endl;
        }
    }

    if( ftruncate( fd, pos ) ) {		// truncate file
        cerr << "FitsWriter: Failed to truncate file at position " << pos << ": " << strerror(errno) << endl;
    }
    
    if( do_fsync ) {	// Do an fsync if requested, and close the file
        fsync(fd);
    }
    
    close(fd);
    fd = -1;
    
}


void FitsWriter::write_exptime_table( void ) {

    vector<string> tmeta;
    tmeta.reserve(16);
    Fits::insertCard( tmeta, Fits::makeCard<string>( "XTENSION", "TABLE   ") );
    Fits::insertCard( tmeta, Fits::makeCard( "BITPIX", 8 ) );
    Fits::insertCard( tmeta, Fits::makeCard( "NAXIS", 2 ) );
    Fits::insertCard( tmeta, Fits::makeCard( "NAXIS1", 26 ) );
    Fits::insertCard( tmeta, Fits::makeCard( "NAXIS2", (int)times.size() ) );
    Fits::insertCard( tmeta, Fits::makeCard( "PCOUNT", 0) );
    Fits::insertCard( tmeta, Fits::makeCard( "GCOUNT", 1) );
    Fits::insertCard( tmeta, Fits::makeCard( "TFIELDS", 1) );
    Fits::insertCard( tmeta, Fits::makeCard( "TTYPE1", "DATE-BEG") );
    Fits::insertCard( tmeta, Fits::makeCard( "TBCOL1", 1) );
    Fits::insertCard( tmeta, Fits::makeCard( "TFORM1", "A26") );
    Fits::insertCard( tmeta, Fits::makeCard( "TUNIT1", "time") );
    Fits::insertCard( tmeta, Fits::makeCard( "EXTNAME", "TABULATIONS") );
    Fits::insertCard( tmeta, Fits::makeCard( "SOLARNET", 0.5) );
    Fits::insertCard( tmeta, Fits::makeCard( "OBS_HDU", 1) );
    Fits::insertCard( tmeta, Fits::makeCard( "END" ) );
    
    write_meta( fd, tmeta );
    pad_to( fd, 2880, ' ' );

    for( auto& t: times ) {
        string ts = to_iso_extended_string( t );
        ts.resize( 26, ' ' );
        if( write( fd, ts.data(), 26 ) != 26 ) {
            cerr << "write_exptime_table: Failed to write timestamp: " << strerror(errno) << endl;
        }
    }
    
    pad_to( fd, 2880, ' ' );    // After ASCII-tables, the data-secion should be padded with spaces to the next 2880 boundary.
    
}


std::shared_ptr<uint8_t> FitsWriter::get_buf( size_t N ) {
    shared_ptr<uint8_t> buf;
    {
        lock_guard<mutex> lock(bufMtx);
        if( !buffers.empty() ) {
            buf = buffers.front();
            buffers.pop_front();
        }
    }
    if(!buf) {
        buf.reset( new uint8_t[N], [](uint8_t*& p) { delete[] p; p=nullptr; });
    }
    return buf;
}


void FitsWriter::return_buf( std::shared_ptr<uint8_t>& buf ) {
    if( buf ) {
        lock_guard<mutex> block(bufMtx);
        buffers.push_back(buf);
        buf.reset();
    }
}


void FitsWriter::clear_bufs( void ) {
    lock_guard<mutex> block(bufMtx);
    buffers.clear();
}


string FitsWriter::get_saves( void ) {
    
    static const bpx::ptime epoch_time( boost::gregorian::date(1970,1,1) ); 
    string ret;
    
    lock_guard<mutex> lock( globalMtx );
    for( auto& s: saves ) {
        bpx::time_duration tmp = s.first - epoch_time;
        size_t nSecs = tmp.total_seconds();
        size_t nMicros = s.first.time_of_day().total_microseconds() - s.first.time_of_day().total_seconds()*1000000L;
        string line = boost::str( boost::format("%ld.%06ld") % nSecs % nMicros );
        tmp = s.second - epoch_time;
        nSecs = tmp.total_seconds();
        nMicros = s.second.time_of_day().total_microseconds() - s.second.time_of_day().total_seconds()*1000000L;
        line += boost::str( boost::format("    %ld.%06ld") % nSecs % nMicros );
        ret += line + "\n";
    }

    return ret;
}


void FitsWriter::add_save( boost::posix_time::ptime from, boost::posix_time::ptime to ) {
    
    lock_guard<mutex> lock( globalMtx );
    saves.emplace( make_pair(from, to) );
    
}
