#ifndef REDUX_FILE_FILEIO_HPP
#define REDUX_FILE_FILEIO_HPP

#include "redux/image/image.hpp"
#include "redux/util/array.hpp"

#include <string>
#include <memory>
#include <thread>
#include <vector>

#include <boost/filesystem.hpp>
namespace bfs = boost::filesystem;

namespace redux {

    namespace file {
        
        /*! @defgroup file FileIO
         *  @{
         */

        enum Format : uint8_t { FMT_NONE = 0,
                                FMT_ANA,
                                FMT_FITS,
                                FMT_NCDF,
                                FMT_MOMFBD,
                                FMT_PLAIN
                              };


        Format readFmt( const std::string& );
        Format guessFmt( const std::string& );

        std::string filetype( const std::string& );
        std::vector<std::string> filetypes( const std::vector<std::string>& );

        std::string getMetaText( std::string, FileMeta::Ptr&, bool raw=false, bool all=false );
        std::vector<std::string> getMetaTextAsCards( const std::string&, FileMeta::Ptr&, bool raw=false, bool all=false );
        std::vector<std::string> getMetaTexts( const std::vector<std::string>&, bool raw=false, bool all=false );
        std::vector<std::vector<std::string>> getMetaTextsAsCards( const std::vector<std::string>&, bool raw=false, bool all=false );
        inline std::string getMetaText( std::string file, bool raw=false, bool all=false ) {
            std::string ret;
            try {
                FileMeta::Ptr meta;
                ret = getMetaText( file, meta, raw, all );
            } catch( ... ) {}
            return ret;
        }
        inline std::vector<std::string> getMetaTextAsCards( const std::string& file, bool raw=false, bool all=false ) {
            std::vector<std::string> ret;
            try {
                FileMeta::Ptr meta;
                ret = getMetaTextAsCards( file, meta, raw, all );
            } catch( ... ) {}
            return ret;
        }
        
        template <typename T>
        void getOrRead( const std::string& fn, std::shared_ptr<T>& data );

        template <typename T>
        void getOrRead2( const std::string& fn, std::shared_ptr<redux::image::Image<T>>& im );
        
        FileMeta::Ptr getMeta(const std::string& fn, bool size_only=false);

        void readFile( const std::string& fn, char* data, FileMeta::Ptr& meta );
        template <typename T>
        void readFile( const std::string& fn, redux::util::Array<T>& data );
        template <typename T>
        void readFile( const std::string& fn, redux::image::Image<T>& data, bool metaOnly=false );

        template <typename T>
        void writeFile( const std::string& fn, const redux::util::Array<T>& data );
        template <typename T>
        void writeFile( const std::string& fn, const redux::image::Image<T>& data );

        typedef std::function<void(char*,size_t,FileMeta::Ptr&)> postLoadCallback;
        typedef std::function<void(char*,size_t,FileMeta::Ptr&)> preSumCallback;
        
//         void loadFiles( const std::vector<std::string>& fn, char* out, size_t frameSize, uint8_t nThreads=std::thread::hardware_concurrency(),
//                        double* averages=nullptr, double* times=nullptr, std::string progressMsg="" );
        void loadFiles( const std::vector<std::string>& fn, char* out, size_t frameSize, uint8_t nThreads=std::thread::hardware_concurrency(),
                       postLoadCallback postLoad = postLoadCallback() );
        void sumFiles( const std::vector<std::string>& fn, double* out, size_t frameSize, uint8_t nThreads=std::thread::hardware_concurrency(),
                       preSumCallback preSum = preSumCallback() );
        
        enum ErrorHandling { EH_PRINT=1, EH_THROW };
        extern ErrorHandling errorHandling;         //<! Specify if routines should throw or print errors. Default is EH_PRINT;
        inline void setErrorHandling( ErrorHandling eh ) { errorHandling = eh; }

        bool isRelative( const bfs::path& );
        std::string cleanPath( std::string path, std::string base = "" );
        std::string getHome( const std::string& username="" );
        bfs::path weaklyCanonical( bfs::path );
        
        
        /*! @} */

    }

}

#endif // REDUX_FILE_FILEIO_HPP
