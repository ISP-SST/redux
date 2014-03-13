#ifndef REDUX_FILE_FILEIO_HPP
#define REDUX_FILE_FILEIO_HPP

#include "redux/image/image.hpp"
#include "redux/util/array.hpp"

#include <string>
#include <memory>

namespace redux {

    namespace file {

        /*! @defgroup file FileIO
         *  @{
         */

        enum Format : uint8_t { FMT_NONE = 0,
                                FMT_ANA,
                                FMT_FITS,
                                FMT_NCDF,
                                FMT_PLAIN
                              };


        Format readFmt(const std::string& );
        Format guessFmt(const std::string& );

        template <typename T>
        void getOrRead( const std::string& fn, std::shared_ptr<T>& data );

        template <typename T>
        void getOrRead2( const std::string& fn, std::shared_ptr<redux::image::Image<T>>& im );

        template <typename T>
        void readFile( const std::string& fn, redux::util::Array<T>& data );
        template <typename T>
        void readFile( const std::string& fn, redux::image::Image<T>& data );

        template <typename T>
        void writeFile( const std::string& fn, redux::util::Array<T>& data );
        template <typename T>
        void writeFile( const std::string& fn, redux::image::Image<T>& data );





//std::shared_ptr<Image>
        /*! @} */

    }

}

#endif // REDUX_FILE_FILEIO_HPP
