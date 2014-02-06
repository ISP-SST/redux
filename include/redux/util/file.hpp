#ifndef REDUX_UTIL_FILE_HPP
#define REDUX_UTIL_FILE_HPP

#include "redux/exception.hpp"

#include <stdio.h>
#include <string>
#include <string.h>
#include <iostream>

namespace redux {

    namespace util {

        /*!  @ingroup util
         *  @{
         */

        class FileException : public redux::RecoverableException {
        public:
            FileException ( void ) : RecoverableException ( "redux::util::FileException" ) {}
            FileException ( const std::string &message ) : RecoverableException ( message ) {}

            virtual ~FileException ( void ) throw () {}
        };

        /*!  @class     File
         *   @brief     Wrapper providing RAII for std::FILE pointers
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         */
        class File {

            std::string m_FileName;
            FILE* m_File;

        public:
            File() : m_FileName(""), m_File(nullptr) {}
            File( const std::string& name, const std::string& mode = std::string( "r" ) ) : m_FileName(name), m_File(nullptr) {
                open(name,mode);
            }
            ~File() { if( m_File ) fclose( m_File ); }

            void open( const std::string& name, const std::string& mode = std::string( "r" ) ) {
                m_File = fopen( name.c_str(), mode.c_str() );
                if( !m_File ) {
                    throw FileException( std::string( "Failed to open file: " ) + name );
                }
            }
            
            void rewind(void) { ::rewind(m_File); }

            inline void seek( int64_t offset, int whence = SEEK_SET ) {
                int res = fseek( m_File, offset, whence );
                if( res < 0 ) {
                    throw FileException( std::string( "File:seekError: " ) + strerror( errno ) );
                }
            }

            inline int64_t tell( void ) {
                int64_t res = ftell( m_File );
                if( res < 0 ) {
                    throw FileException( std::string( "File:tellError: " ) + strerror( errno ) );
                }
                return res;
            }

            template <typename T>
            inline size_t readOrThrow( T* out, size_t nElements, const std::string& msg = "" ) {
                size_t count = fread( reinterpret_cast<char*>( out ), sizeof( T ), nElements, m_File );
                if( count == nElements ) {
                    return count * sizeof( T );
                }
                throw FileException( "File:readError: " + msg );
            }

            template <typename T>
            inline size_t writeOrThrow( T* data, size_t nElements, const std::string& msg = "" ) {
                size_t count = fwrite( reinterpret_cast<char*>( data ), sizeof( T ), nElements, m_File );
                if( count == nElements ) {
                    return count * sizeof( T );
                }
                throw FileException( "File:writeError: " + msg );
            }

            operator FILE*() const { return m_File; }

        };

        /*! @} */


    }  // namespace util

}  // namespace redux



#endif // REDUX_UTIL_FILE_HPP
