#include "redux/file/fileio.hpp"

#include "redux/file/fileana.hpp"
#include "redux/file/filefits.hpp"
#include "redux/file/filemomfbd.hpp"
#include "redux/util/stringutil.hpp"

#include <iostream>
#include <future>
#include <map>
#include <mutex>

using namespace redux::file;
using namespace redux::util;
using namespace std;

namespace {

    mutex fileMutex;




    template <typename T>
    map<string, shared_ptr<T>>& getFileCache( void ) {
        static map<string, shared_ptr<T>> cache;
        return cache;
    }

    template <typename T>
    shared_ptr<T> getFile( const string& fn ) {
        static map<string, shared_ptr<T>> cache;
        {
            unique_lock<mutex> lock( fileMutex );
            auto found = cache.find( fn );
            if( found != cache.end() ) {
                return found->second;
            }
        }


    }

}


Format redux::file::readFmt( const string& filename ) {

    ifstream strm( filename, ifstream::binary );
    if( strm ) {
        uint32_t magic;
        strm.read( reinterpret_cast<char*>( &magic ), sizeof(uint32_t) );
        if( strm.good() && (strm.gcount()==sizeof(uint32_t)) ) {
            switch( magic ) {
                case Ana::MAGIC_ANA: ;
                case Ana::MAGIC_ANAR: return FMT_ANA;
                case FileMomfbd::MAGIC_MOMFBD8: ;	// Fall through
                case FileMomfbd::MAGIC_MOMFBD10: ;	// Fall through
                case FileMomfbd::MAGIC_MOMFBD11: return FMT_MOMFBD;
#ifdef REDUX_WITH_FITS
                case Fits::MAGIC_FITS: return FMT_FITS;
#endif
                //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
                default: std::ios_base::failure("readFmt needs to be implemented for this file-type: \"" +filename+"\""); 
            }
        } else throw std::ios_base::failure("Failed to read file: "+filename);
    } else throw std::ios_base::failure("Failed to open file: "+filename);

    return FMT_NONE;

}


Format redux::file::guessFmt( const string& filename ) {

    size_t pos = filename.find_last_of('.');
    if( pos != string::npos && pos < filename.length() ) {
        string ext = filename.substr(pos+1);
        std::transform(ext.begin(), ext.end(),ext.begin(), ::toupper);
        if( ext == "F0" || ext == "FZ" ) {
            return FMT_ANA;
        } else if( ext == "FITS" ) {
            return FMT_FITS;
        } else if( ext == "MOMFBD" ) {
            return FMT_MOMFBD;
        }
    }
    return FMT_NONE;

}





template <typename T>
void redux::file::getOrRead( const string& fn, shared_ptr<T>& data ) {

    //static auto& cache = getFileCache<T>();

    cout << "redux::file::getOrRead(" << fn << ")     "  << endl;

   // future<shared_ptr<T>> async( getFile<T>, fn, 100 );

    /*
        {
            unique_lock<mutex> lock( fileMutex );
            auto found = cache.find( fn );
            if( found != cache.end() ) {
                data = found->second;
                return;
            }
        }*/




}

template <typename T>
void redux::file::getOrRead2( const string& fn, shared_ptr<redux::image::Image<T>>& im ) {
//void redux::file::getOrRead(const string& fn, redux::image::Image<T>::Ptr& im) {
    cout << "redux::file::getOrRead2(" << fn << ")" << endl;
}

template void redux::file::getOrRead( const string&, typename redux::image::Image<int16_t>::Ptr& );
// template void redux::file::getOrRead<int32_t>(const string&, typename redux::image::Image<int32_t>::Ptr&);
// template void redux::file::getOrRead<float>(const string&, typename redux::image::Image<float>::Ptr&);
// template void redux::file::getOrRead<double>(const string&, typename redux::image::Image<double>::Ptr&);

template void redux::file::getOrRead2( const string&, shared_ptr<redux::image::Image<int16_t>>& );
template void redux::file::getOrRead2( const string&, shared_ptr<redux::image::Image<int32_t>>& );
template void redux::file::getOrRead2( const string&, shared_ptr<redux::image::Image<float>>& );
template void redux::file::getOrRead2( const string&, shared_ptr<redux::image::Image<double>>& );


shared_ptr<redux::file::FileMeta> redux::file::getMeta(const string& fn, bool size_only) {

    Format fmt = readFmt(fn);
    shared_ptr<redux::file::FileMeta> meta;
    
    switch(fmt) {
        case FMT_ANA: {
            meta.reset( new Ana(fn) );
            break;
        }
#ifdef REDUX_WITH_FITS
        case FMT_FITS: {
            meta.reset( new Fits(fn) );
            break;
        }
#endif
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: cout << "file::getMeta(arr) needs to be implemented for this file-type: " << fmt << "   \"" << fn << "\""  << endl;
    }
    
    return move(meta);
    
}


void redux::file::readFile( const string& filename, char* data, shared_ptr<FileMeta>& meta ) {

    Format fmt = readFmt(filename);
    switch(fmt) {
        case FMT_ANA: {
            if( !meta ) {
                meta.reset( new Ana() );
            }
            shared_ptr<Ana> hdr = static_pointer_cast<Ana>(meta);
            if( hdr ) {
                Ana::read( filename, data, hdr );
            } else cout << "file::readFile(string,char*,meta) failed to cast meta-pointer into Ana type." << endl;
            break;
        }
#ifdef REDUX_WITH_FITS
        case FMT_FITS: {
            shared_ptr<Fits> hdr;
            if( !meta ) {
                hdr.reset( new Fits() );
                meta = static_pointer_cast<FileMeta>(hdr);
            } else hdr = static_pointer_cast<Fits>(meta);
            if( hdr ) {
                hdr->read( filename );
                Fits::read( hdr, data );
            } else cout << "file::readFile(string,char*,meta) failed to cast meta-pointer into Fits type." << endl;
            break;
        }
#endif
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: cout << "file::readFile(string,char*,meta) needs to be implemented for this file-type: " << fmt << "   \"" << filename << "\"" << endl;
    }


}

template <typename T>
void redux::file::readFile( const string& filename, redux::util::Array<T>& data ) {
    Format fmt = readFmt(filename);
    switch(fmt) {
        case FMT_ANA: {
            shared_ptr<Ana> hdr(new Ana());
            Ana::read(filename,data,hdr);
            break;
        }
#ifdef REDUX_WITH_FITS
        case FMT_FITS: {
            shared_ptr<Fits> hdr(new Fits());
            Fits::read(filename,data,hdr);
            break;
        }
#endif
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: cout << "file::read(arr) needs to be implemented for this file-type: " << fmt << "   \"" << filename << "\""  << endl;
    }


}
template void redux::file::readFile( const string& filename, redux::util::Array<uint8_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<int16_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<int32_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<int64_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<float>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<double>& data );


template <typename T>
void redux::file::readFile( const string& filename, redux::image::Image<T>& image, bool metaOnly ) {
    Format fmt = readFmt(filename);
    switch(fmt) {
        case FMT_ANA: Ana::read(filename, image, metaOnly); break;
#ifdef REDUX_WITH_FITS
        case FMT_FITS: Fits::read(filename, image, metaOnly); break;
#endif
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: cout << "file::read(img) needs to be implemented for this file-type: " << fmt << "   \"" << filename << "\""  << endl;
    }


}
template void redux::file::readFile( const string& filename, redux::image::Image<uint8_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<int16_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<int32_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<int64_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<float>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<double>& data, bool );


template <typename T>
void redux::file::writeFile( const string& filename, redux::util::Array<T>& data ) {
    Format fmt = guessFmt(filename);
    switch(fmt) {
        case FMT_ANA: Ana::write(filename, data); break;
#ifdef REDUX_WITH_FITS
        case FMT_FITS: Fits::write(filename, data); break;
#endif
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: cout << "file::write(arr) needs to be implemented for this file-type: " << fmt << "   \"" << filename << "\""  << endl;
    }


}
template void redux::file::writeFile( const string& filename, redux::util::Array<uint8_t>& data );
template void redux::file::writeFile( const string& filename, redux::util::Array<int16_t>& data );
template void redux::file::writeFile( const string& filename, redux::util::Array<int32_t>& data );
template void redux::file::writeFile( const string& filename, redux::util::Array<int64_t>& data );
template void redux::file::writeFile( const string& filename, redux::util::Array<float>& data );
template void redux::file::writeFile( const string& filename, redux::util::Array<double>& data );


template <typename T>
void redux::file::writeFile( const string& filename, redux::image::Image<T>& image ) {
    Format fmt = guessFmt(filename);
    switch(fmt) {
        case FMT_ANA: Ana::write(filename, image); break;
#ifdef REDUX_WITH_FITS
        case FMT_FITS: Fits::write(filename, image); break;
#endif
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: cout << "file::write(img) needs to be implemented for this file-type: " << fmt << "   \"" << filename << "\""  << endl;
    }


}
template void redux::file::writeFile( const string& filename, redux::image::Image<uint8_t>& data );
template void redux::file::writeFile( const string& filename, redux::image::Image<int16_t>& data );
template void redux::file::writeFile( const string& filename, redux::image::Image<int32_t>& data );
template void redux::file::writeFile( const string& filename, redux::image::Image<int64_t>& data );
template void redux::file::writeFile( const string& filename, redux::image::Image<float>& data );
template void redux::file::writeFile( const string& filename, redux::image::Image<double>& data );


/*void redux::file::loadFiles( const vector<string>& filenames, char* out, size_t frameSize, uint8_t nThreads,
                             double* averages, double* times, string progressMsg ) {
    
    size_t nImages = filenames.size();
    
    atomic<size_t> imgIndex(0);
    
    vector<thread> threads;
    for( uint8_t i=0; i<nThreads; ++i ) {
        threads.push_back( std::thread(
            [&](){
                size_t myIndex;
                shared_ptr<FileMeta> myMeta;
                while( (myIndex=imgIndex.fetch_add(1)) < nImages ) {
                    char* myPtr = out + myIndex*frameSize;
                    try {
                        readFile( filenames[myIndex], myPtr, myMeta );
                        if( times ) times[ myIndex ] = myMeta->getAverageTime().time_of_day().total_nanoseconds()*1E-9;
                    } catch (const exception& e ) {
                        if( !progressMsg.empty() ) printProgress( "\nloadFiles: " + string(e.what()) + " at file #" + to_string(myIndex) + "\n", -1);
                        memset( myPtr, 0, frameSize );  // zero image and continue.
                    }
                    if( !progressMsg.empty() ) printProgress( progressMsg, (myIndex*100.0/(nImages-1)));
                }
            }));
    }
    for (auto& th : threads) th.join();

}*/


void redux::file::loadFiles( const vector<string>& filenames, char* out, size_t frameSize, uint8_t nThreads,
                             function<void(char*,size_t,shared_ptr<FileMeta>&)> postLoad) {
    
    size_t nImages = filenames.size();
    
    atomic<size_t> imgIndex(0);

    vector<thread> threads;
    for( uint8_t i=0; i<nThreads; ++i ) {
        threads.push_back( std::thread(
            [&](){
                size_t myIndex;
                shared_ptr<FileMeta> myMeta;
                while( (myIndex=imgIndex.fetch_add(1)) < nImages ) {
                    char* myPtr = out + myIndex*frameSize;
                    try {
                        readFile( filenames[myIndex], myPtr, myMeta );
                        postLoad( myPtr, myIndex, myMeta );
                    } catch (const exception& e ) {
                        printProgress( "\nloadFiles: " + string(e.what()) + " at file #" + to_string(myIndex) + "\n", -1);  // print via mutex
                        memset( myPtr, 0, frameSize );  // zero image and continue.
                    }
                }
            }));
    }
    for (auto& th : threads) th.join();
}


void redux::file::sumFiles( const std::vector<std::string>& filenames, double* out, size_t frameSize, uint8_t nThreads,
                            preSumCallback preSum ) {
    
/*    
    size_t nImages = filenames.size();
    mutex mtx;

    atomic<size_t> imgIndex(0);
    atomic<size_t> threadIndex(0);
    auto sumFunc = [=,&mtx]( double* a ) {
        std::unique_lock<mutex> lock(mtx);
        for( size_t b=0; b<nPixels; ++b ) summedData[b] += a[b];
    };
    //proc.append( std::bind( sumFunc, std::placeholders::_1 ) );
    
    std::vector<std::thread> threads;
    for( uint8_t i=0; i<nThreads; ++i ) {
        threads.push_back( std::thread(
            [&](){
                size_t myImgIndex;
                size_t myThreadIndex = threadIndex.fetch_add(1);
                double* mySumPtr = sumPtr+myThreadIndex*nPixels;
                double* myTmpPtr = tmpPtr+myThreadIndex*nPixels;
                while( (myImgIndex=imgIndex.fetch_add(1)) < nImages ) {
                    initFuncs[myImgIndex](myTmpPtr,myThreadIndex);
                    for( size_t i=0; i<nPixels; ++i ) mySumPtr[i] += myTmpPtr[i];
                    if( kw.verbose > 1 ) printProgress( statusString, (myImgIndex*100.0/(nImages-1)));
                }
                sumFunc(mySumPtr);
            }));
    }

    for (auto& th : threads) th.join();
*/
}
