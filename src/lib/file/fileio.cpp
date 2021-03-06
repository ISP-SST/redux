#include "redux/file/fileio.hpp"

#include "redux/file/fileana.hpp"
#include "redux/file/filefits.hpp"
#include "redux/file/filemomfbd.hpp"
#include "redux/util/stringutil.hpp"

#include <iostream>
#include <future>
#include <list>
#include <map>
#include <mutex>
#include <pwd.h>
#include <unistd.h>

using namespace redux::file;
using namespace redux::util;
using namespace std;

ErrorHandling redux::file::errorHandling = EH_PRINT;

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

    vector<string> txtToCards( string in ) {
        vector<string> out;
        while( in.length() >= 80 ) {
            out.push_back( in.substr(0,80) );
            in.erase( 0, 80 );
        }
        if( !in.empty() ) {
            out.push_back( in );
            // TBD:: throw or warn about not a multiple of 80?
        }
        return out;
    }
    
}


Format redux::file::readFmt( const string& filename ) {

    ifstream strm( filename, ifstream::binary );
    if( strm ) {
        uint32_t magic;
        strm.read( reinterpret_cast<char*>( &magic ), sizeof(uint32_t) );
        if( strm.good() && (strm.gcount()==sizeof(uint32_t)) ) {
            switch( magic ) {
                case Ana::MAGIC_ANA: ;	             // Fall through
                case Ana::MAGIC_ANAR: return FMT_ANA;
                case FileMomfbd::MAGIC_MOMFBD8: ;	 // Fall through
                case FileMomfbd::MAGIC_MOMFBD10: ;	 // Fall through
                case FileMomfbd::MAGIC_MOMFBD11: return FMT_MOMFBD;
#ifdef RDX_WITH_FITS
                case Fits::MAGIC_FITS: return FMT_FITS;
#endif
                //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
                default: throw std::ios_base::failure("readFmt needs to be implemented for this file-type: \"" +filename+"\""); 
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


string redux::file::filetype( const string& file ) {
    Format fmt = FMT_NONE;
    try {
        fmt = redux::file::readFmt( file );
    } catch( const std::ios_base::failure& ) {
        // silently ignore IO/read errors, and guess frmat from extension below.
    }
    if( fmt == FMT_NONE ) fmt = redux::file::guessFmt( file );
    switch(fmt) {
        case FMT_ANA: return "ANA";
        case FMT_FITS: return "FITS";
        case FMT_NCDF: return "NCDF";
        case FMT_MOMFBD: return "MOMFBD";
        default: ;
    }
    return "";
}


vector<string> redux::file::filetypes( const vector<string>& files ) {
    vector<string> filetypes;
    size_t nF = files.size();
    if( nF == 0 ) return filetypes;
    filetypes.resize( files.size(), "" );
    for( size_t i=0; i<nF; ++i ) {
        filetypes[i] = filetype( files[i] );
    }
    return filetypes;
}


string redux::file::getMetaText( string file, FileMeta::Ptr& meta, bool raw, bool all ) {
    string txt;
    try {
        if( !meta ) {
            bfs::path fn( file );
            if( bfs::exists(fn) && bfs::is_regular_file(fn) ) {
                meta = getMeta( file );
            }
        }
        if( meta ) {
            meta->getAverageTime();
            meta->getEndTime();
            meta->getStartTime();
            vector<string> hdrTexts = meta->getText( raw );
            if( !all && (hdrTexts.size() > 1) ) hdrTexts.resize(1);   // just return primary info
            for( auto& t: hdrTexts ) txt += t;
        }
    } catch( ... ) {}
    return txt;
}


vector<string> redux::file::getMetaTexts( const vector<string>& files, bool raw, bool all ) {
    vector<string> texts;
    size_t nF = files.size();
    if( nF == 0 ) return texts;
    texts.resize( files.size(), "" );
    for( size_t i(0); i<nF; ++i ) {
        texts[i] = getMetaText( files[i], raw, all );
    }
    return texts;
}


vector<string> redux::file::getMetaTextAsCards( const string& file, FileMeta::Ptr& meta, bool raw, bool all ) {
    return txtToCards( getMetaText( file, meta, raw, all ) );
}


vector<vector<string>> redux::file::getMetaTextsAsCards( const vector<string>& files, bool raw, bool all ) {
    vector<vector<string>> ret(files.size());
    for( size_t i(0); i<files.size(); ++i ) {
        ret[i] = txtToCards( getMetaText( files[i], raw, all ) );
    }
    return ret;
}


template <typename T>
void redux::file::getOrRead( const string& fn, shared_ptr<T>& data ) {

    //static auto& cache = getFileCache<T>();

    cout << "getOrRead(" << fn << ")     "  << endl;

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
    cout << "getOrRead2(" << fn << ")" << endl;
}

template void redux::file::getOrRead( const string&, typename redux::image::Image<int16_t>::Ptr& );
// template void redux::file::getOrRead<int32_t>(const string&, typename redux::image::Image<int32_t>::Ptr&);
// template void redux::file::getOrRead<float>(const string&, typename redux::image::Image<float>::Ptr&);
// template void redux::file::getOrRead<double>(const string&, typename redux::image::Image<double>::Ptr&);

template void redux::file::getOrRead2( const string&, shared_ptr<redux::image::Image<int16_t>>& );
template void redux::file::getOrRead2( const string&, shared_ptr<redux::image::Image<int32_t>>& );
template void redux::file::getOrRead2( const string&, shared_ptr<redux::image::Image<float>>& );
template void redux::file::getOrRead2( const string&, shared_ptr<redux::image::Image<double>>& );


FileMeta::Ptr redux::file::getMeta(const string& fn, bool size_only) {

    Format fmt = readFmt(fn);
    FileMeta::Ptr meta;
    
    switch(fmt) {
        case FMT_ANA: {
            meta.reset( new Ana(fn) );
            break;
        }
#ifdef RDX_WITH_FITS
        case FMT_FITS: {
            meta.reset( new Fits(fn) );
            break;
        }
#endif
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: {
            string msg = "getMeta(arr) needs to be implemented for this file-type: " + to_string(fmt)
                       + "   \"" + fn + "\"";
            throw runtime_error(msg);
        }
    }
    
    return meta;
    
}


void redux::file::readFile( const string& filename, char* data, FileMeta::Ptr& meta ) {

    try {
        Format fmt = readFmt(filename);
        switch(fmt) {
            case FMT_ANA: {
                if( !meta ) {
                    meta.reset( new Ana() );
                }
                shared_ptr<Ana> hdr = static_pointer_cast<Ana>(meta);
                if( hdr ) {
                    Ana::read( filename, data, hdr );
                } else {
                    string msg = "readFile(string,char*,FileMeta::Ptr&) failed to cast meta-pointer into Ana type.";
                    throw runtime_error(msg);
                }
                break;
            }
#ifdef RDX_WITH_FITS
            case FMT_FITS: {
                if( !meta ) {
                    meta.reset( new Fits() );
                }
                shared_ptr<Fits> hdr = static_pointer_cast<Fits>(meta);
                if( hdr ) {
                    hdr->read( filename );
                    Fits::read( hdr, data );
                } else {
                    string msg = "readFile(string,char*,FileMeta::Ptr&) failed to cast meta-pointer into Fits type.";
                    throw runtime_error(msg);
                }
                break;
            }
#endif
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
            default: {
                string msg = "readFile(string,char*,FileMeta::Ptr&) needs to be implemented for this file-type: " + to_string(fmt)
                           + "   \"" + filename + "\"";
                throw runtime_error(msg);
            }
        }
    } catch ( std::exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "FileIO Error: " << e.what() << endl;
        } else throw;
    }


}

template <typename T>
void redux::file::readFile( const string& filename, redux::util::Array<T>& data ) {
    
    try {
        Format fmt = readFmt(filename);
        switch(fmt) {
            case FMT_ANA: {
                shared_ptr<Ana> hdr(new Ana());
                Ana::read(filename,data,hdr);
                break;
            }
#ifdef RDX_WITH_FITS
            case FMT_FITS: {
                shared_ptr<Fits> hdr(new Fits());
                Fits::read(filename,data,hdr);
                break;
            }
#endif
            //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
            default: {
                string msg = "readFile(str,Array<T>) needs to be implemented for this file-type: " + to_string(fmt)
                           + "   \"" + filename + "\"";
                throw runtime_error(msg);
            }
        }
    } catch ( std::exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "FileIO Error: " << e.what() << endl;
        } else throw;
    }

}
template void redux::file::readFile( const string& filename, redux::util::Array<int8_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<uint8_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<int16_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<uint16_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<int32_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<uint32_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<int64_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<uint64_t>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<float>& data );
template void redux::file::readFile( const string& filename, redux::util::Array<double>& data );


template <typename T>
void redux::file::readFile( const string& filename, redux::image::Image<T>& image, bool metaOnly ) {
    
    try {
        Format fmt = readFmt(filename);
        switch(fmt) {
            case FMT_ANA: Ana::read(filename, image, metaOnly); break;
#ifdef RDX_WITH_FITS
            case FMT_FITS: Fits::read(filename, image, metaOnly); break;
#endif
            //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
            default: {
                string msg = "readFile(string,Image<T>,bool) needs to be implemented for this file-type: " + to_string(fmt)
                           + "   \"" + filename + "\"";
                throw runtime_error(msg);
            }
        }
    } catch ( std::exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "FileIO Error: " << e.what() << endl;
        } else throw;
    }

}
template void redux::file::readFile( const string& filename, redux::image::Image<int8_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<uint8_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<int16_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<uint16_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<int32_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<uint32_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<int64_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<uint64_t>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<float>& data, bool );
template void redux::file::readFile( const string& filename, redux::image::Image<double>& data, bool );


template <typename T>
void redux::file::writeFile( const string& filename, const redux::util::Array<T>& data ) {
    
    try {
        Format fmt = guessFmt(filename);
        switch(fmt) {
            case FMT_ANA: Ana::write(filename, data); break;
#ifdef RDX_WITH_FITS
            case FMT_FITS: Fits::write(filename, data); break;
#endif
            //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
            default: {
                string msg = "writeFile(string,Array<T>) needs to be implemented for this file-type: " + to_string(fmt)
                           + "   \"" + filename + "\"";
                throw runtime_error(msg);
            }
        }
    } catch ( std::exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "FileIO Error: " << e.what() << endl;
        } else throw;
    }


}
template void redux::file::writeFile( const string& filename, const redux::util::Array<int8_t>& data );
template void redux::file::writeFile( const string& filename, const redux::util::Array<uint8_t>& data );
template void redux::file::writeFile( const string& filename, const redux::util::Array<int16_t>& data );
template void redux::file::writeFile( const string& filename, const redux::util::Array<uint16_t>& data );
template void redux::file::writeFile( const string& filename, const redux::util::Array<int32_t>& data );
template void redux::file::writeFile( const string& filename, const redux::util::Array<uint32_t>& data );
template void redux::file::writeFile( const string& filename, const redux::util::Array<int64_t>& data );
template void redux::file::writeFile( const string& filename, const redux::util::Array<uint64_t>& data );
template void redux::file::writeFile( const string& filename, const redux::util::Array<float>& data );
template void redux::file::writeFile( const string& filename, const redux::util::Array<double>& data );


template <typename T>
void redux::file::writeFile( const string& filename, const redux::image::Image<T>& image ) {
    
    try {
        Format fmt = guessFmt(filename);
        switch(fmt) {
            case FMT_ANA: Ana::write(filename, image); break;
#ifdef RDX_WITH_FITS
            case FMT_FITS: Fits::write(filename, image); break;
#endif
            //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
            default: {
                string msg = "writeFile(string,Image<T>) needs to be implemented for this file-type: " + to_string(fmt)
                           + "   \"" + filename + "\"";
                throw runtime_error(msg);
            }
        }
    } catch ( std::exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "FileIO Error: " << e.what() << endl;
        } else throw;
    }


}
template void redux::file::writeFile( const string& filename, const redux::image::Image<int8_t>& data );
template void redux::file::writeFile( const string& filename, const redux::image::Image<uint8_t>& data );
template void redux::file::writeFile( const string& filename, const redux::image::Image<int16_t>& data );
template void redux::file::writeFile( const string& filename, const redux::image::Image<uint16_t>& data );
template void redux::file::writeFile( const string& filename, const redux::image::Image<int32_t>& data );
template void redux::file::writeFile( const string& filename, const redux::image::Image<uint32_t>& data );
template void redux::file::writeFile( const string& filename, const redux::image::Image<int64_t>& data );
template void redux::file::writeFile( const string& filename, const redux::image::Image<uint64_t>& data );
template void redux::file::writeFile( const string& filename, const redux::image::Image<float>& data );
template void redux::file::writeFile( const string& filename, const redux::image::Image<double>& data );


/*void redux::file::loadFiles( const vector<string>& filenames, char* out, size_t frameSize, uint8_t nThreads,
                             double* averages, double* times, string progressMsg ) {
    
    size_t nImages = filenames.size();
    
    atomic<size_t> imgIndex(0);
    
    vector<thread> threads;
    for( uint8_t i=0; i<nThreads; ++i ) {
        threads.push_back( std::thread(
            [&](){
                size_t myIndex;
                FileMeta::Ptr myMeta;
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
                             function<void(char*,size_t,FileMeta::Ptr&)> postLoad) {
    
    size_t nImages = filenames.size();
    
    atomic<size_t> imgIndex(0);

    mutex mtx;
    list<string> msgs;
    vector<thread> threads;
    for( uint8_t i=0; i<nThreads; ++i ) {
        threads.push_back( std::thread(
            [&](){
                size_t myIndex;
                FileMeta::Ptr myMeta;
                while( (myIndex=imgIndex.fetch_add(1)) < nImages ) {
                    char* myPtr = out + myIndex*frameSize;
                    string fn;
                    {
                        lock_guard<mutex> lock(mtx);
                        fn = filenames[myIndex];
                    }
                    try {
                        readFile( fn, myPtr, myMeta );
                        postLoad( myPtr, myIndex, myMeta );
                    } catch (const exception& e ) {
                        string msg = "file #" + to_string(myIndex) + "\"" + fn + "\": " + e.what();  
                        memset( myPtr, 0, frameSize );  // zero image and continue.
                        lock_guard<mutex> lock(mtx);
                        msgs.push_back( msg );
                    }
                }
            }));
    }
    for (auto& th : threads) th.join();
    if( msgs.size() ) {
        string msg = "FileIO Error: loadFiles(vector<string>,char*,size_t,uint8_t,callback_t)";
        for( auto& m: msgs ) msg += "\n\t  " + m;
        if( errorHandling == EH_PRINT ) {
            cerr << msg << endl;
        } else throw runtime_error( msg );
    }
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


bool redux::file::isRelative( const bfs::path &p ) {
    string s = p.string();
    return (!s.empty() && s[0] != '/');
}


string expandTilde( string p ) {

    if( !p.empty() && p[0] == '~' ) {
        string home_path;
        string subpath = p.substr( 1 );                 // drop the '~'
        if( subpath.empty() || subpath[0] == '/' ) {    // resolve own home-directory
            home_path = getHome();
        } else {                                        // resolve home-directory of named user
            string user = subpath;
            size_t pos = subpath.find_first_of('/');
            if( pos != string::npos ) {
                user = subpath.substr( 0, pos );
                subpath = subpath.substr( pos );
            }
            home_path = getHome(user);
        }
        if( !home_path.empty() ) {                      // if we found something, return it, else return input
            if( *home_path.rbegin() != '/' ) home_path.push_back('/');
            return home_path + subpath;
        }
    }
    return p;
}


string redux::file::cleanPath( string in, string base ) {

    if( in.empty() ) return in;
    
    in = expandTilde( in );       // returns unmodified if no '~'

    if( (in.length() > 1) && (in.substr(0,2) == "./") ) { // local path
        in.replace( 0, 1, bfs::current_path().string() );
    }
    
    if( in[0] != '/' ) {                                       // relative path, resolve and apply "base"
        base = expandTilde(base);
        if( !base.empty() ) {
            bfs::path base_path(base);
            if( base[0] != '/' ) base_path = bfs::current_path() / bfs::path(base);     // prepend current dir if base is relative
            in.insert( 0, base_path.string() + "/");
        }
    }
    

    bfs::path in_path(in);
    string in_filename = in_path.filename().c_str();
    if( in_filename.empty() || in_filename == "." ) {
        in_path = in_path.parent_path();
    }
#if BOOST_VERSION < 106000    // lexically_normal was introduced in 1.60
    return in_path.string();
#else
    return in_path.lexically_normal().string();
#endif
}


string redux::file::getHome( const string& username ) {
    
    string tmp;
    struct passwd pwent;
    struct passwd *pwentp;
    char buf[1024];
    if( username.empty() ) {    // resolve own home-directory
        tmp = getenv("HOME");
        if( tmp.empty() ) { // No HOME in env, try looking with pwuid
            if( !getpwuid_r( geteuid(), &pwent, buf, sizeof buf, &pwentp ) ) {
                tmp = pwent.pw_dir;
            }
        }
    } else {                // resolve home-directory of named user
        if( !getpwnam_r( username.c_str(), &pwent, buf, sizeof buf, &pwentp ) ) {
            if( username == pwent.pw_name ) {   // NOTE: non-existant user can return weird results if this is not verified
                tmp = pwent.pw_dir;
            }
        }
    }

    if( !tmp.empty() ) {    // append slash
        if( *tmp.rbegin() != '/' ) tmp.push_back('/');
    }

    return tmp;
}


bfs::path redux::file::weaklyCanonical( bfs::path p ) {
    p = expandTilde( p.string() );
    try {
#if BOOST_VERSION > 106000
        p = bfs::weakly_canonical( p );
#elif BOOST_VERSION > 104800
        bfs::path tmp = p;
        while( !tmp.empty() && !exists( tmp ) ) {
            tmp = tmp.parent_path();
        }
        if( !tmp.empty() ) {
            p = bfs::relative( p, tmp );
            p = bfs::canonical( tmp ) / p;
        }
#endif
    } catch( ... ) {};
    string p_filename = p.filename().c_str();
    if( p_filename.empty() || p_filename == "." ) {
        p = p.parent_path();
    }
#if BOOST_VERSION < 106000
    return p;
#else
    return p.lexically_normal();
#endif
}


