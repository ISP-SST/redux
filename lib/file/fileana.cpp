#include "redux/file/fileana.hpp"

#include "redux/util/endian.hpp"

#include "redux/file/anacompress.hpp"
#include "redux/file/anadecompress.hpp"

using namespace redux::file;
using namespace redux::util;
using namespace std;

#if REDUX_BYTE_ORDER == REDUX_LITTLE_ENDIAN
static const int system_is_big_endian = 0;
#elif REDUX_BYTE_ORDER == REDUX_BIG_ENDIAN
static const int system_is_big_endian = 1;
#else
#error REDUX_BYTE_ORDER not set
#endif


void redux::file::ana::readCompressed( redux::util::File& file, char* data, size_t nElements, const shared_ptr<AnaInfo>& hdr ) {

    size_t compressedSize = hdr->m_CompressedHeader.tsize - 14;

    unique_ptr<uint8_t> tmp( new uint8_t[compressedSize] );

    size_t count = fread( tmp.get(), 1, compressedSize, file );
    if( count != compressedSize ) {
        throw DataIOException( "Failed to read the specified number of bytes." );
    }

    // fix a possible problem with chdr.nblocks
    if( hdr->m_CompressedHeader.bsize * hdr->m_CompressedHeader.nblocks > nElements ) {
        hdr->m_CompressedHeader.nblocks = static_cast<uint32_t>( nElements / hdr->m_CompressedHeader.bsize );
    }

    // consistency check
    if( hdr->m_CompressedHeader.type % 2 == hdr->m_Header.datyp ) {
        throw DataIOException( "readCompressedAna: type-mismatch between primary and compressed headers." );
    }

    int slice = hdr->m_CompressedHeader.slice_size;
    int blockSize = hdr->m_CompressedHeader.bsize;
    int nBlocks = hdr->m_CompressedHeader.nblocks;

    switch( hdr->m_CompressedHeader.type ) {
        case( 0 ):
            anadecrunch( tmp.get(), reinterpret_cast<int16_t*>( data ), slice, blockSize, nBlocks, !system_is_big_endian );
            break;
        case( 1 ):
            anadecrunch8( tmp.get(), reinterpret_cast<int8_t*>( data ), slice, blockSize, nBlocks, !system_is_big_endian );
            break;
        case( 2 ):
            anadecrunchrun( tmp.get(), reinterpret_cast<int16_t*>( data ), slice, blockSize, nBlocks, !system_is_big_endian );
            break;
        case( 3 ):
            anadecrunchrun8( tmp.get(), reinterpret_cast<int8_t*>( data ), slice, blockSize, nBlocks, !system_is_big_endian );
            break;
        case( 4 ):
            anadecrunch32( tmp.get(), reinterpret_cast<int32_t*>( data ), slice, blockSize, nBlocks, !system_is_big_endian );
            break;
        default:
            throw DataIOException( "readCompressedAna: unrecognized type of compressed data" );
    }

}


void redux::file::ana::readUncompressed( redux::util::File& file, char* data, size_t nElements, const shared_ptr<AnaInfo>& hdr ) {

    size_t nBytes = nElements * typeSizes[hdr->m_Header.datyp];

    size_t count = fread( data, 1, nBytes, file );
    if( count != nBytes ) {
        throw DataIOException( "Failed to read the specified number of bytes." );
    }

    if( hdr->m_Header.subf & 128 ) { // saved on big endian system.
        switch( typeSizes[hdr->m_Header.datyp] ) {
            case 2: swapEndian( reinterpret_cast<int16_t*>( data ), nElements ); break;
            case 4: swapEndian( reinterpret_cast<int32_t*>( data ), nElements ); break;
            case 8: swapEndian( reinterpret_cast<int64_t*>( data ), nElements ); break;
            default: ;
        }
    }

}


int redux::file::ana::compressData( unique_ptr<uint8_t>& out, const char* data, int nElements, const shared_ptr<AnaInfo>& hdr, int slice ) {

    int limit = nElements * typeSizes[hdr->m_Header.datyp];
    limit += ( limit >> 1 );
    int nx = hdr->m_Header.dim[0];
    int ny = nElements / nx;
    int runlengthflag( 0 );  // runlength unused/untested
    int res;
    uint8_t* cdata = new uint8_t[ limit ]; // allocate double size since compression can fail and generate larger data.
    out.reset( cdata );

    switch( hdr->m_Header.datyp ) {
        case( 0 ) : {
            if( runlengthflag ) res = anacrunchrun8( cdata, reinterpret_cast<const uint8_t*>( data ), slice, nx, ny, limit, system_is_big_endian );
            else res = anacrunch8( cdata, reinterpret_cast<const uint8_t*>( data ), slice, nx, ny, limit, system_is_big_endian );
            break;
        }
        case( 1 ) : {
            if( runlengthflag ) res = anacrunchrun( cdata, reinterpret_cast<const int16_t*>( data ), slice, nx, ny, limit, system_is_big_endian );
            else res = anacrunch( cdata, reinterpret_cast<const int16_t*>( data ), slice, nx, ny, limit, system_is_big_endian );
            break;
        }
        case( 2 ) : {
            if( runlengthflag ) throw DataIOException( "compressData: runlength not supported 32-bit types." );
            else res = anacrunch32( cdata, reinterpret_cast<const int32_t*>( data ), slice, nx, ny, limit, system_is_big_endian );
            break;
        }
        default: throw DataIOException( "compressData: Unsupported data type." );
    }

    return res;
}


void redux::file::readAna( redux::util::File& file, char* data, std::shared_ptr<redux::file::AnaInfo> hdr ) {

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
    unique_ptr<int> dim( new int[nDims] );
    for( int i( 0 ); i < nDims; ++i ) {
        dim.get()[i] = hdr->m_Header.dim[nDims - i - 1];
        nElements *= hdr->m_Header.dim[nDims - i - 1];
    }

    if( hdr->m_Header.subf & 1 ) ana::readCompressed( file, data, nElements, hdr );
    else ana::readUncompressed( file, data, nElements, hdr );

}


void redux::file::writeAna( redux::util::File& file, const char* data,
                            std::shared_ptr<redux::file::AnaInfo> hdr,
                            bool compress, int slice ) {

    if( !hdr.get() ) {
        throw DataIOException( "writeAna: The header object is invalid, cannot write file." );
    }

    if( hdr->m_Header.datyp > 5 || hdr->m_Header.datyp < 0 ) {
        throw DataIOException( "writeAna: hdr.datyp is not valid." );
    }

    if( hdr->m_Header.ndim > 16 ) {
        throw DataIOException( "writeAna: the ANA/f0 format does not support more dimensions than 16." );
    }

    size_t nElements = 1;
    for( uint8_t i = 0; i < hdr->m_Header.ndim; ++i ) {
        if( hdr->m_Header.dim[i] < 1 ) {
            throw DataIOException( "writeAna: dimSize < 1" );
        }
        nElements *= hdr->m_Header.dim[i];
    }

    size_t textSize = hdr->m_ExtendedHeader.length();
    size_t totalSize = nElements * ana::typeSizes[hdr->m_Header.datyp];
    size_t compressedSize = totalSize;

    hdr->m_Header.synch_pattern = MAGIC_ANA;
    hdr->m_Header.subf = 0;
    hdr->m_Header.nhb = static_cast<uint8_t>( 1 + std::max( textSize + 256, size_t( 0 ) ) / 512 );
    memset( hdr->m_Header.cbytes, 0, 4 );

    unique_ptr<uint8_t> cData;
    if( compress ) {
        if( hdr->m_Header.ndim == 2 ) {
            compressedSize = ana::compressData( cData, data, nElements, hdr, slice );
            if( compressedSize < totalSize ) {
                int tmp = static_cast<int>( compressedSize );
                hdr->m_Header.subf |= 1;
                memcpy( hdr->m_Header.cbytes, &tmp, 4 );
                memcpy( &hdr->m_CompressedHeader, cData.get(), 14 );
                totalSize = compressedSize;
            }
            else {          // compressed data larger -> store uncompressed.
                cData.reset();
                memset( &hdr->m_CompressedHeader, 0, 14 );
            }
        }
        else {
            throw DataIOException( "writeAna: compression only supported for 2D data." );
        }
    }

    rewind( file );
    hdr->write( file );

    if( cData ) fwrite( cData.get() + 14, 1, totalSize, file ); // compressed header is already written, so skip 14 bytes.
    else  fwrite( data, 1, totalSize, file );

}

