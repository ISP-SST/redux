#ifndef REDUX_FILE_FILEANA_HPP
#define REDUX_FILE_FILEANA_HPP

#include "redux/file/anainfo.hpp"
#include "redux/file/exceptions.hpp"
#include "redux/util/file.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arrayutil.hpp"

#include <iostream>


namespace redux {

    namespace file {

        namespace ana {

            enum TypeIndex { ANA_BYTE = 0/*uint8*/, ANA_WORD/*int16*/, ANA_INT/*int32*/, ANA_FLOAT, ANA_DOUBLE, ANA_LONG };
            static std::tuple<uint8_t, int16_t, int32_t, float, double, int64_t> typeMap;
            static const uint8_t typeSizes[] = { 1, 2, 4, 4, 8, 8 };
            
            void readCompressed( redux::util::File& file, char* data, size_t nElements, const std::shared_ptr<AnaInfo>& hdr );
            void readUncompressed( redux::util::File& file, char* data, size_t nElements, const std::shared_ptr<AnaInfo>& hdr );
            int compressData( std::unique_ptr<uint8_t>& out, const char* data, int nElements, const std::shared_ptr<AnaInfo>& hdr, int slice );
            
        }

        void readAna( redux::util::File& file, char* data, std::shared_ptr<redux::file::AnaInfo> hdr = std::shared_ptr<AnaInfo>( 0 ) );

        void writeAna( redux::util::File& file, const char* data, std::shared_ptr<redux::file::AnaInfo> hdr, bool compress = false, int slice = 5 );


        template <typename T>
        void readAna( redux::util::File& file, redux::util::Array<T>& data = redux::util::Array<T>(), std::shared_ptr<redux::file::AnaInfo> hdr = std::shared_ptr<AnaInfo>( 0 ) ) {

            if( !hdr.get() ) {
                hdr.reset( new AnaInfo() );
                hdr->read( file );
            }

            if( fseek( file, hdr->hdrSize, SEEK_SET ) ) {
                throw DataIOException( "Seek operation failed." );
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

            auto tmp = redux::util::sharedArray<char>( nElements * ana::typeSizes[hdr->m_Header.datyp] );
            if( hdr->m_Header.subf & 1 ) {
                ana::readCompressed( file, tmp.get(), nElements, hdr );
            } else {
                ana::readUncompressed( file, tmp.get(), nElements, hdr );
            }


            switch( hdr->m_Header.datyp ) { // convert to fp_t
                case( ana::ANA_BYTE ):   data.rawCopy( tmp.get() ); break;
                case( ana::ANA_WORD ):   data.rawCopy( reinterpret_cast<uint16_t*>(tmp.get())); break;
                case( ana::ANA_INT ):    data.rawCopy( reinterpret_cast<uint32_t*>(tmp.get())); break;
                case( ana::ANA_FLOAT ):  data.rawCopy( reinterpret_cast<float*>(tmp.get())); break;
                case( ana::ANA_DOUBLE ): data.rawCopy( reinterpret_cast<double*>(tmp.get())); break;
                default: ;
            }


        }

    } // end namespace file

} // end namespace redux


#endif // REDUX_FILE_FILEANA_HPP
