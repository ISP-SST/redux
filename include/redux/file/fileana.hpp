#ifndef REDUX_FILE_FILEANA_HPP
#define REDUX_FILE_FILEANA_HPP

#include "redux/file/fileinfo.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arrayutil.hpp"

#include <fstream>


namespace redux {

    namespace file {

        struct Ana : public redux::util::FileInfo {

            enum Magic { MAGIC_ANA = 0x5555aaaa, MAGIC_ANAR = 0xaaaa5555 };
            enum TypeIndex { ANA_BYTE = 0/*uint8*/, ANA_WORD/*int16*/, ANA_INT/*int32*/, ANA_FLOAT, ANA_DOUBLE, ANA_LONG };
            static const uint8_t typeSizes[]; // = { 1, 2, 4, 4, 8, 8 };

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

            static void readCompressed( std::ifstream& file, char* data, size_t nElements, const std::shared_ptr<Ana>& hdr );
            static void readUncompressed( std::ifstream& file, char* data, size_t nElements, const std::shared_ptr<Ana>& hdr );
            static int compressData( std::shared_ptr<uint8_t>& out, const char* data, int nElements, const std::shared_ptr<Ana>& hdr, int slice );

            static void read( std::ifstream& file, char* data, std::shared_ptr<redux::file::Ana> hdr = 0 );

            static void write( std::ofstream& file, const char* data, std::shared_ptr<redux::file::Ana> hdr, bool compress = false, int slice = 5 );


            template <typename T>
            static void read( std::ifstream& file, redux::util::Array<T>& data = redux::util::Array<T>(), std::shared_ptr<redux::file::Ana> hdr = 0 ) {

                if( !hdr.get() ) {
                    hdr.reset( new Ana() );
                    hdr->read( file );
                }
                file.seekg( hdr->hdrSize );
                if( !file.good() ) {
                    throw std::ios_base::failure( "Seek operation failed." );
                }

                // f0 stores the dimensions with the fast index first, so swap them before allocating the array
                int nDims = hdr->m_Header.ndim;
                size_t nElements = 1;
                bool forceResize = false;
                std::vector<size_t> dimSizes( nDims, 0 );
                for( int i( 0 ); i < nDims; ++i ) {
                    dimSizes[i] = hdr->m_Header.dim[nDims - i - 1];
                    if( dimSizes[i] != data.size( i ) ) {
                        forceResize = true;
                    }
                    nElements *= hdr->m_Header.dim[nDims - i - 1];
                }

                if( forceResize ) {
                    data.reset( dimSizes );
                }

                auto tmp = std::shared_ptr<char>( new char[nElements * typeSizes[hdr->m_Header.datyp]], []( char * p ) { delete[] p; } );
                if( hdr->m_Header.subf & 1 ) {
                    readCompressed( file, tmp.get(), nElements, hdr );
                }
                else {
                    readUncompressed( file, tmp.get(), nElements, hdr );
                }

                switch( hdr->m_Header.datyp ) {
                    case( ANA_BYTE ):   data.template rawCopy<char>( tmp.get() ); break;
                    case( ANA_WORD ):   data.template rawCopy<uint16_t>( tmp.get() ); break;
                    case( ANA_INT ):    data.template rawCopy<uint32_t>( tmp.get() ); break;
                    case( ANA_FLOAT ):  data.template rawCopy<float>( tmp.get() ); break;
                    case( ANA_DOUBLE ): data.template rawCopy<double>( tmp.get() ); break;
                    default: ;
                }


            }

        };


    } // end namespace file

} // end namespace redux


#endif // REDUX_FILE_FILEANA_HPP
