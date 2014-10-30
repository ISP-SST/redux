#ifndef REDUX_FILE_FILEINFO_HPP
#define REDUX_FILE_FILEINFO_HPP

#include <iostream>
#include <exception>

namespace redux {

    namespace file {

        /*! @ingroup file FileIO
         *  @{
         */

        struct FileInfo {

            FileInfo ( void ) {}
            virtual ~FileInfo ( void ) = default;

        };



        template <typename T, typename U>
        inline size_t readOrThrow ( T& strm, U* out, size_t nElements, const std::string& msg = "" ) {
            size_t nBytes = nElements*sizeof ( U );
            strm.read ( reinterpret_cast<char*> ( out ), nBytes );
            if ( !strm.good() ) {
                throw std::ios_base::failure ( "Read failed: " + msg );
            }

            return nBytes;
        }



        template <typename T, typename U>
        inline size_t writeOrThrow ( T& strm, const U* data, size_t nElements, const std::string& msg = "" ) {
            size_t nBytes = nElements*sizeof ( U );
            strm.write ( reinterpret_cast<const char*> ( data ), nBytes );
            if ( !strm.good() ) {
                throw std::ios_base::failure ( "Write failed: " + msg );
            }

            return nBytes;
        }




        /*! @} */

    }

}

#endif // REDUX_FILE_FILEINFO_HPP
