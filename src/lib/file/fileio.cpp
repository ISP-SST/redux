#include "redux/file/fileio.hpp"

#include "redux/file/fileana.hpp"

#include <iostream>
#include <future>
#include <mutex>

using namespace redux::file;
using namespace std;

namespace {

    mutex fileMutex;




    template <typename T>
    map<string, std::shared_ptr<T>>& getFileCache( void ) {
        static map<string, std::shared_ptr<T>> cache;
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


Format redux::file::readFmt( const std::string& filename ) {

    ifstream strm( filename, ifstream::binary );
    if( strm ) {
        uint32_t magic;
        strm.read( reinterpret_cast<char*>( &magic ), sizeof(uint32_t) );
        if( strm.good() ) {
            switch (magic) {
                case Ana::MAGIC_ANA: ;
                case Ana::MAGIC_ANAR: return FMT_ANA;
                //case Fits::MAGIC_FITS: return FMT_FITS;
                //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
                default: cout << "readFmt needs to be implemented for this file-type: \"" << filename << "\""  << endl; return FMT_NONE; 
            }
        }
    }
    return FMT_NONE;

}


Format redux::file::guessFmt( const std::string& filename ) {

    size_t pos = filename.find_last_of('.');
    if( pos != string::npos && pos < filename.length() ) {
        string ext = filename.substr(pos+1);
        if( ext == "f0" || ext == "fz" ) {
            return FMT_ANA;
        }
    }
    return FMT_NONE;

}





template <typename T>
void redux::file::getOrRead( const std::string& fn, std::shared_ptr<T>& data ) {

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
void redux::file::getOrRead2( const std::string& fn, std::shared_ptr<redux::image::Image<T>>& im ) {
//void redux::file::getOrRead(const std::string& fn, redux::image::Image<T>::Ptr& im) {
    cout << "redux::file::getOrRead2(" << fn << ")" << endl;
}

template void redux::file::getOrRead( const std::string&, typename redux::image::Image<int16_t>::Ptr& );
// template void redux::file::getOrRead<int32_t>(const std::string&, typename redux::image::Image<int32_t>::Ptr&);
// template void redux::file::getOrRead<float>(const std::string&, typename redux::image::Image<float>::Ptr&);
// template void redux::file::getOrRead<double>(const std::string&, typename redux::image::Image<double>::Ptr&);

template void redux::file::getOrRead2( const std::string&, std::shared_ptr<redux::image::Image<int16_t>>& );
template void redux::file::getOrRead2( const std::string&, std::shared_ptr<redux::image::Image<int32_t>>& );
template void redux::file::getOrRead2( const std::string&, std::shared_ptr<redux::image::Image<float>>& );
template void redux::file::getOrRead2( const std::string&, std::shared_ptr<redux::image::Image<double>>& );



template <typename T>
void redux::file::readFile( const std::string& filename, redux::util::Array<T>& data ) {
    Format fmt = readFmt(filename);
    switch(fmt) {
        case FMT_ANA: {
            shared_ptr<Ana> hdr(new Ana());
            Ana::read(filename,data,hdr);
            break;
        }
        //case Fits::MAGIC_FITS: return FMT_FITS;
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: cout << "file::read(arr) needs to be implemented for this file-type: " << fmt << "   \"" << filename << "\""  << endl;
    }


}
template void redux::file::readFile( const std::string& filename, redux::util::Array<uint8_t>& data );
template void redux::file::readFile( const std::string& filename, redux::util::Array<int16_t>& data );
template void redux::file::readFile( const std::string& filename, redux::util::Array<int32_t>& data );
template void redux::file::readFile( const std::string& filename, redux::util::Array<int64_t>& data );
template void redux::file::readFile( const std::string& filename, redux::util::Array<float>& data );
template void redux::file::readFile( const std::string& filename, redux::util::Array<double>& data );


template <typename T>
void redux::file::readFile( const std::string& filename, redux::image::Image<T>& image ) {
    Format fmt = readFmt(filename);
    switch(fmt) {
        case FMT_ANA: Ana::read(filename, image); break;
        //case Fits::MAGIC_FITS: return FMT_FITS;
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: cout << "file::read(img) needs to be implemented for this file-type: " << fmt << "   \"" << filename << "\""  << endl;
    }


}
template void redux::file::readFile( const std::string& filename, redux::image::Image<uint8_t>& data );
template void redux::file::readFile( const std::string& filename, redux::image::Image<int16_t>& data );
template void redux::file::readFile( const std::string& filename, redux::image::Image<int32_t>& data );
template void redux::file::readFile( const std::string& filename, redux::image::Image<int64_t>& data );
template void redux::file::readFile( const std::string& filename, redux::image::Image<float>& data );
template void redux::file::readFile( const std::string& filename, redux::image::Image<double>& data );


template <typename T>
void redux::file::writeFile( const std::string& filename, redux::util::Array<T>& data ) {
    Format fmt = guessFmt(filename);
    switch(fmt) {
        case FMT_ANA: Ana::write(filename,data); break;
        //case Fits::MAGIC_FITS: return FMT_FITS;
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: cout << "file::write(arr) needs to be implemented for this file-type: " << fmt << "   \"" << filename << "\""  << endl;
    }


}
template void redux::file::writeFile( const std::string& filename, redux::util::Array<uint8_t>& data );
template void redux::file::writeFile( const std::string& filename, redux::util::Array<int16_t>& data );
template void redux::file::writeFile( const std::string& filename, redux::util::Array<int32_t>& data );
template void redux::file::writeFile( const std::string& filename, redux::util::Array<int64_t>& data );
template void redux::file::writeFile( const std::string& filename, redux::util::Array<float>& data );
template void redux::file::writeFile( const std::string& filename, redux::util::Array<double>& data );


template <typename T>
void redux::file::writeFile( const std::string& filename, redux::image::Image<T>& image ) {
    Format fmt = guessFmt(filename);
    switch(fmt) {
        case FMT_ANA: Ana::write(filename, image); break;
        //case Fits::MAGIC_FITS: return FMT_FITS;
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: cout << "file::write(img) needs to be implemented for this file-type: " << fmt << "   \"" << filename << "\""  << endl;
    }


}
template void redux::file::writeFile( const std::string& filename, redux::image::Image<uint8_t>& data );
template void redux::file::writeFile( const std::string& filename, redux::image::Image<int16_t>& data );
template void redux::file::writeFile( const std::string& filename, redux::image::Image<int32_t>& data );
template void redux::file::writeFile( const std::string& filename, redux::image::Image<int64_t>& data );
template void redux::file::writeFile( const std::string& filename, redux::image::Image<float>& data );
template void redux::file::writeFile( const std::string& filename, redux::image::Image<double>& data );

