#ifndef REDUX_FILE_FILEANA_HPP
#define REDUX_FILE_FILEANA_HPP

#include "redux/file/fileinfo.hpp"
#include "redux/image/image.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arrayutil.hpp"

#include <fstream>


namespace redux {

    namespace file {

        /*! @ingroup file
         *  @{
         */

        /*! Container for reading/writing ANA files.
        *
        */
        struct Ana : public redux::file::FileInfo {

            enum Magic { MAGIC_ANA = 0x5555aaaa, MAGIC_ANAR = 0xaaaa5555 };
            enum TypeIndex { ANA_BYTE = 0/*uint8*/, ANA_WORD/*int16*/, ANA_LONG/*int32*/, ANA_FLOAT, ANA_DOUBLE, ANA_LONGLONG/*int64_t unsupported*/, ANA_UNDEF=255 };
            static const uint8_t typeSizes[]; // = { 1, 2, 4, 4, 8, 8 };

            typedef std::shared_ptr<Ana> Ptr;
            
            Ana( void );
            Ana( const std::string& );

            void read( std::ifstream& );
            void read( const std::string& );

            void write( std::ofstream& );

            std::string getText( void ) {
                return m_Header.txt + m_ExtendedHeader;
            }

            struct raw_header {                    // first block for ana files
                uint32_t synch_pattern;
                uint8_t subf;
                uint8_t source;
                uint8_t nhb;
                uint8_t datyp;
                uint8_t ndim;
                uint8_t free1;
                uint8_t cbytes[4];
                uint8_t free[178];
                uint32_t dim[16];
                char txt[256];
            } m_Header;

            struct compressed_header {
                uint32_t tsize, nblocks, bsize;
                uint8_t slice_size, type;
            } m_CompressedHeader;

            std::string m_ExtendedHeader;
            size_t hdrSize;

            static void readCompressed( std::ifstream& file, char* data, size_t nElements, const Ana* hdr );
            static void readUncompressed( std::ifstream& file, char* data, size_t nElements, const Ana* hdr );
            static int compressData( std::shared_ptr<uint8_t>& out, const char* data, int nElements, const std::shared_ptr<Ana>& hdr, int slice );

            static void read( const std::string& filename, char* data, std::shared_ptr<redux::file::Ana> hdr=0 );
            static void write( const std::string& filename, const char* data, std::shared_ptr<redux::file::Ana> hdr, bool compress = false, int slice=5 );


            template <typename T>
            static void read( const std::string& filename, redux::util::Array<T>& data, std::shared_ptr<redux::file::Ana>& hdr=0 );
            
            template <typename T>
            static void read( const std::string& filename, redux::image::Image<T>& data );
            
            template <typename T>
            static void write( const std::string& filename, const redux::util::Array<T>& data, std::shared_ptr<redux::file::Ana> hdr=0, int sliceSize=0 );
            
            template <typename T>
            static void write( const std::string& filename, const redux::image::Image<T>& image, int sliceSize=0 );
            

        };

        /*! @} */

    } // end namespace file

} // end namespace redux


#endif // REDUX_FILE_FILEANA_HPP
