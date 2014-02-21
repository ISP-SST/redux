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

const uint8_t Ana::typeSizes[] = { 1, 2, 4, 4, 8, 8 };

Ana::Ana( void ) : hdrSize( 0 ) {
    memset( &m_Header, 0, sizeof( raw_header ) );
    memset( &m_CompressedHeader, 0, sizeof( compressed_header ) );
}

Ana::Ana( const std::string& filename ) : hdrSize( 0 ) {
    memset( &m_Header, 0, sizeof( raw_header ) );
    memset( &m_CompressedHeader, 0, sizeof( compressed_header ) );
    read( filename );
}

void Ana::read( ifstream& file ) {

    hdrSize = 0;
    file.seekg(0);

    size_t rawSize = sizeof( struct raw_header );  // = 512 by construction
    memset( &m_Header, 0, rawSize );

    file.read( reinterpret_cast<char*>( &m_Header ), rawSize );
    if( !file.good() ) {
        throw ios_base::failure( "Failed to read ANA header" );
    }

    hdrSize = file.tellg();
    bool hasMagic = ( m_Header.synch_pattern == MAGIC_ANA );
    bool hasReversedMagic = ( m_Header.synch_pattern == MAGIC_ANAR );

    if( !( hasMagic || hasReversedMagic ) ) {
        throw logic_error( "Failed to read ANA header: initial 4 bytes does not match ANA" );
    }

    if( m_Header.nhb > 1 ) {                    // read in long header

        if( m_Header.nhb > 16 ) {
            throw logic_error( "Warning: AnaHeader::read() - extended header is longer than 16 blocks!" );
        }
        else {
            size_t extraTextSize = ( m_Header.nhb - 1 ) * rawSize;
            unique_ptr<char[]> tmpText( new char [ extraTextSize + 257 ] );
            char* txtPtr = tmpText.get();
            memset( txtPtr+256, 0, extraTextSize + 1 );
            memcpy( txtPtr, m_Header.txt, 256 );
            file.read( txtPtr+256, extraTextSize );
            if( !file.good() ) {
                if( file.eof() ) {
                    throw ios_base::failure( "Failed to read ANA header: <EOF> reached" );
                }
            }
            hdrSize = file.tellg();

            char* firstNull = txtPtr;
            while( *firstNull ) {
                firstNull++;
            }

            char* ptr = firstNull;
            while( ++ptr < txtPtr + extraTextSize + 256 ) {
                if( *ptr ) {
                    swap( *ptr, *firstNull++ );
                }
            }

            m_ExtendedHeader = txtPtr;

        }
    }

    // m_Header.dim is always stored as little-endian, so swap endianess if system is big-endian.
    if( system_is_big_endian ) {
        swapEndian( &( m_Header.dim ), m_Header.ndim );
    }

    // compressed data
    if( m_Header.subf & 1 ) {
        size_t chdrSize = 14; //  <-- struct is usually padded to to next 4-byte boundary (16).
        file.read( reinterpret_cast<char*>( &m_CompressedHeader ), chdrSize );
        if( !file.good() ) {
            if( file.eof() ) {
                throw ios_base::failure( "Failed to read ANA header: <EOF> reached" );
            }
        }
        hdrSize += chdrSize;

        if( system_is_big_endian ) { // for big endian platforms
            swapEndian( &( m_CompressedHeader.tsize ) );
            swapEndian( &( m_CompressedHeader.nblocks ) );
            swapEndian( &( m_CompressedHeader.bsize ) );
        }
    }

}

void Ana::read( const std::string& filename ) {
    ifstream file( filename, ifstream::binary );
    read( file );
}

void Ana::write( ofstream& file ) {

    size_t hdrSize = sizeof( struct raw_header );   // = 512 by construction

    // file header entries are in little endian order, so swap for big-endian systems.
    if( system_is_big_endian ) {
        m_Header.subf |= 128;   // mark file as written on big-endian machine.
        swapEndian( &( m_Header.synch_pattern ) );
        swapEndian( &( m_Header.cbytes ) );
        swapEndian( &( m_Header.dim ), m_Header.ndim );
    }

    file.write( reinterpret_cast<char*>( &m_Header ), 256 );
    if( !file.good() ) {
        throw ios_base::failure( "Failed to write ANA header (m_Header)" );
    }

    if( system_is_big_endian ) { // ...and then swap back
        swapEndian( &( m_Header.synch_pattern ) );
        swapEndian( &( m_Header.cbytes ) );
        swapEndian( &( m_Header.dim ), m_Header.ndim );
    }

    // extended header
    if( ( m_Header.nhb > 1 ) && ( m_ExtendedHeader.length() > 0 ) ) {
        size_t txtSize = ( m_Header.nhb - 1 ) * hdrSize + 256;
        unique_ptr<char[]> extHdr( new char[txtSize] );
        memset( extHdr.get(), 0, txtSize );
        memcpy( extHdr.get(), m_ExtendedHeader.c_str(), m_ExtendedHeader.length() );
        file.write( extHdr.get(), txtSize );
        if( !file.good() ) {
            throw ios_base::failure( "Failed to write ANA header (extHdr)" );
        }
    }
    else {
        file.write( m_Header.txt, 256 );
        if( !file.good() ) {
            throw ios_base::failure( "Failed to write ANA header (txt)" );
        }
    }

    // compressed data
    if( m_Header.subf & 1 ) {
        if( system_is_big_endian ) {
            swapEndian( &( m_CompressedHeader.tsize ) );
            swapEndian( &( m_CompressedHeader.nblocks ) );
            swapEndian( &( m_CompressedHeader.bsize ) );
        }

        file.write( reinterpret_cast<char*>( &m_CompressedHeader ), 14 );
        if( !file.good() ) {
            throw ios_base::failure( "Failed to write ANA header (m_CompressedHeader)" );
        }

        if( system_is_big_endian ) {
            swapEndian( &( m_CompressedHeader.tsize ) );
            swapEndian( &( m_CompressedHeader.nblocks ) );
            swapEndian( &( m_CompressedHeader.bsize ) );
        }
    }


}


void redux::file::Ana::readCompressed( ifstream& file, char* data, size_t nElements, const shared_ptr<Ana>& hdr ) {

    size_t compressedSize = hdr->m_CompressedHeader.tsize - 14;

    shared_ptr<uint8_t> tmp( new uint8_t[compressedSize], []( uint8_t* p ) { delete[] p; } );

    file.read(reinterpret_cast<char*>(tmp.get()),compressedSize);
    if( !file.good() ) {
        throw ios_base::failure( "Failed to read the specified number of bytes." );
    }

    // fix a possible problem with chdr.nblocks
    if( hdr->m_CompressedHeader.bsize * hdr->m_CompressedHeader.nblocks > nElements ) {
        hdr->m_CompressedHeader.nblocks = static_cast<uint32_t>( nElements / hdr->m_CompressedHeader.bsize );
    }

    // consistency check
    if( hdr->m_CompressedHeader.type % 2 == hdr->m_Header.datyp ) {
        throw logic_error( "readCompressedAna: type-mismatch between primary and compressed headers." );
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
            throw invalid_argument( "readCompressedAna: unrecognized type of compressed data" );
    }

}


void redux::file::Ana::readUncompressed( ifstream& file, char* data, size_t nElements, const shared_ptr<Ana>& hdr ) {

    size_t nBytes = nElements * typeSizes[hdr->m_Header.datyp];

    file.read(data,nBytes);
    if( !file.good() ) {
        throw ios_base::failure( "Failed to read the specified number of bytes." );
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


int redux::file::Ana::compressData( shared_ptr<uint8_t>& out, const char* data, int nElements, const shared_ptr<Ana>& hdr, int slice ) {

    int limit = nElements * typeSizes[hdr->m_Header.datyp];
    limit += ( limit >> 1 );
    int nx = hdr->m_Header.dim[0];
    int ny = nElements / nx;
    int runlengthflag( 0 );  // runlength unused/untested
    int res;
    uint8_t* cdata = new uint8_t[ limit ]; // allocates +50% size since compression can fail and generate larger data.
    out.reset( cdata, []( uint8_t* p ) { delete[] p; } );

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
            if( runlengthflag ) throw invalid_argument( "compressData: runlength not supported for 32-bit types." );
            else res = anacrunch32( cdata, reinterpret_cast<const int32_t*>( data ), slice, nx, ny, limit, system_is_big_endian );
            break;
        }
        default: throw invalid_argument( "compressData: Unsupported data type." );
    }

    return res;
}


void redux::file::Ana::read( ifstream& file, char* data, std::shared_ptr<redux::file::Ana> hdr ) {

    if( !hdr.get() ) {
        hdr.reset( new Ana() );
        hdr->read( file );
    }

    file.seekg(hdr->hdrSize);
    if( !file.good() ) {
        throw ios_base::failure( "Seek operation failed." );
    }

    // f0 stores the dimensions with the fast index first, so swap them before allocating the array
    int nDims = hdr->m_Header.ndim;
    size_t nElements = 1;
    unique_ptr<int[]> dim( new int[nDims] );
    for( int i( 0 ); i < nDims; ++i ) {
        dim.get()[i] = hdr->m_Header.dim[nDims - i - 1];
        nElements *= hdr->m_Header.dim[nDims - i - 1];
    }

    if( hdr->m_Header.subf & 1 ) readCompressed( file, data, nElements, hdr );
    else readUncompressed( file, data, nElements, hdr );

}


void redux::file::Ana::write( ofstream& file, const char* data,
                            std::shared_ptr<redux::file::Ana> hdr,
                            bool compress, int slice ) {

    if( !hdr.get() ) {
        throw invalid_argument( "writeAna: The header object is invalid, cannot write file." );
    }

    if( hdr->m_Header.datyp > 5 || hdr->m_Header.datyp < 0 ) {
        throw invalid_argument( "writeAna: hdr.datyp is not valid." );
    }

    if( hdr->m_Header.ndim > 16 ) {
        throw invalid_argument( "writeAna: the ANA/f0 format does not support more dimensions than 16." );
    }

    size_t nElements = 1;
    for( uint8_t i = 0; i < hdr->m_Header.ndim; ++i ) {
        if( hdr->m_Header.dim[i] < 1 ) {
            throw logic_error( "writeAna: dimSize < 1" );
        }
        nElements *= hdr->m_Header.dim[i];
    }

    size_t textSize = hdr->m_ExtendedHeader.length();
    size_t totalSize = nElements * typeSizes[hdr->m_Header.datyp];
    size_t compressedSize = totalSize;


    hdr->m_Header.synch_pattern = MAGIC_ANA;
    hdr->m_Header.subf = 0;
    hdr->m_Header.nhb = static_cast<uint8_t>( 1 + std::max( textSize + 256, size_t( 0 ) ) / 512 );
    memset( hdr->m_Header.cbytes, 0, 4 );

    shared_ptr<uint8_t> cData;
    if( compress ) {
        if( hdr->m_Header.ndim == 2 ) {
            compressedSize = compressData( cData, data, nElements, hdr, slice );
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
            throw invalid_argument( "writeAna: compression only supported for 2D data." );
        }
    }

    file.seekp(0);
    hdr->write( file );

    if( cData ) file.write(reinterpret_cast<char*>(cData.get()) + 14, totalSize); // compressed header is already written, so skip 14 bytes.
    else  file.write(data, totalSize);
    
    if(!file.good())
        throw ios_base::failure("writeAna: write failed.");
    
}

