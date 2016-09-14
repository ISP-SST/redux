#ifndef REDUX_FILE_FILEIO_HPP
#define REDUX_FILE_FILEIO_HPP

#include "redux/image/image.hpp"
#include "redux/util/array.hpp"

#include <string>
#include <memory>
#include <thread>
#include <vector>

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


        Format readFmt(const std::string& );
        Format guessFmt(const std::string& );

        template <typename T>
        void getOrRead( const std::string& fn, std::shared_ptr<T>& data );

        template <typename T>
        void getOrRead2( const std::string& fn, std::shared_ptr<redux::image::Image<T>>& im );
        
        std::shared_ptr<redux::file::FileMeta> getMeta(const std::string& fn, bool size_only=false);

        void readFile( const std::string& fn, char* data, std::shared_ptr<redux::file::FileMeta>& meta );
        template <typename T>
        void readFile( const std::string& fn, redux::util::Array<T>& data );
        template <typename T>
        void readFile( const std::string& fn, redux::image::Image<T>& data, bool metaOnly=false );

        template <typename T>
        void writeFile( const std::string& fn, redux::util::Array<T>& data );
        template <typename T>
        void writeFile( const std::string& fn, redux::image::Image<T>& data );

        typedef std::function<void(char*,size_t,std::shared_ptr<redux::file::FileMeta>&)> postLoadCallback;
        typedef std::function<void(char*,size_t,std::shared_ptr<redux::file::FileMeta>&)> preSumCallback;
        
//         void loadFiles( const std::vector<std::string>& fn, char* out, size_t frameSize, uint8_t nThreads=std::thread::hardware_concurrency(),
//                        double* averages=nullptr, double* times=nullptr, std::string progressMsg="" );
        void loadFiles( const std::vector<std::string>& fn, char* out, size_t frameSize, uint8_t nThreads=std::thread::hardware_concurrency(),
                       postLoadCallback postLoad = postLoadCallback() );
        void sumFiles( const std::vector<std::string>& fn, double* out, size_t frameSize, uint8_t nThreads=std::thread::hardware_concurrency(),
                       preSumCallback preSum = preSumCallback() );




//std::shared_ptr<Image>
        /*! @} */

    }

}

#endif // REDUX_FILE_FILEIO_HPP
