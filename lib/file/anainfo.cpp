#include "redux/file/anainfo.hpp"

#include "redux/util/endian.hpp"
#include "redux/file/exceptions.hpp"

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


AnaInfo::AnaInfo( void ) : hdrSize( 0 ) {
    memset( &m_Header, 0, sizeof( raw_header ) );
    memset( &m_CompressedHeader, 0, sizeof( compressed_header ) );
}

AnaInfo::AnaInfo( const std::string& filename ) : hdrSize( 0 ) {
    memset( &m_Header, 0, sizeof( raw_header ) );
    memset( &m_CompressedHeader, 0, sizeof( compressed_header ) );
    read( filename );
}

void AnaInfo::read( redux::util::File& file ) {

    hdrSize = 0;
    rewind( file );

    size_t rawSize = sizeof( struct raw_header );  // = 512 by construction
    memset( &m_Header, 0, rawSize );

    size_t count = fread( reinterpret_cast<char*>( &m_Header ), 1, rawSize, file );
    if( count < rawSize ) {
        throw DataIOException( "Failed to read ANA header" );
    }

    hdrSize += count;
    bool hasMagic = ( m_Header.synch_pattern == MAGIC_ANA );
    bool hasReversedMagic = ( m_Header.synch_pattern == MAGIC_ANAR );

    if( !( hasMagic || hasReversedMagic ) ) {
        throw DataIOException( "Failed to read ANA header: initial 4 bytes does not match ANA" );
    }

    if( m_Header.nhb > 1 ) {                    // read in long header

        if( m_Header.nhb > 16 ) {
            throw DataIOException( "Warning: AnaHeader::read() - extended header is longer than 16 blocks!" );
        }
        else {
            size_t extraTextSize = ( m_Header.nhb - 1 ) * rawSize;
            unique_ptr<char> tmpText( new char [ extraTextSize + 257 ] ); // add space for primary text block and a null-character.
            char* txtPtr = tmpText.get();
            memset( txtPtr+256, 0, extraTextSize + 1 );
            memcpy( txtPtr, m_Header.txt, 256 );
            count = fread( txtPtr+256, 1, extraTextSize, file );
            if( count < extraTextSize ) {
                if( feof( file ) ) {
                    throw DataIOException( "Failed to read ANA header: <EOF> reached" );
                }
            }
            hdrSize += count;

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
        size_t chdrSize = 14; // sizeof( struct compressed_header );   <-- struct is usually padded to to next 4-byte boundary (16).
        count = fread( reinterpret_cast<char*>( &m_CompressedHeader ), 1, chdrSize, file );
        if( count < chdrSize ) {
            if( feof( file ) ) {
                throw DataIOException( "Failed to read ANA header: <EOF> reached" );
            }
        }
        hdrSize += count;

        if( system_is_big_endian ) { // for big endian platforms
            swapEndian( &( m_CompressedHeader.tsize ) );
            swapEndian( &( m_CompressedHeader.nblocks ) );
            swapEndian( &( m_CompressedHeader.bsize ) );
        }
    }

}

void AnaInfo::read( const std::string& filename ) {
    File file( filename );
    read( file );
}

void AnaInfo::write( redux::util::File& file ) {

    size_t hdrSize = sizeof( struct raw_header );   // = 512 by construction

    // file header entries are in little endian order, so swap for big-endian systems.
    if( system_is_big_endian ) {
        m_Header.subf |= 128;   // mark file as written on big-endian machine.
        swapEndian( &( m_Header.synch_pattern ) );
        swapEndian( &( m_Header.cbytes ) );
        swapEndian( &( m_Header.dim ), m_Header.ndim );
    }

    size_t count = fwrite( reinterpret_cast<char*>( &m_Header ), 1, 256, file );
    if( count != 256 ) {
        throw DataIOException( "Failed to write ANA header (m_Header)" );
    }

    if( system_is_big_endian ) { // ...and then swap back
        swapEndian( &( m_Header.synch_pattern ) );
        swapEndian( &( m_Header.cbytes ) );
        swapEndian( &( m_Header.dim ), m_Header.ndim );
    }

    // extended header
    if( ( m_Header.nhb > 1 ) && ( m_ExtendedHeader.length() > 0 ) ) {
        size_t txtSize = ( m_Header.nhb - 1 ) * hdrSize + 256;
        unique_ptr<char> extHdr( new char[txtSize] );
        memset( extHdr.get(), 0, txtSize );
        memcpy( extHdr.get(), m_ExtendedHeader.c_str(), m_ExtendedHeader.length() );
        count = fwrite( extHdr.get(), 1, txtSize, file );
        if( count != txtSize ) {
            throw DataIOException( "Failed to write ANA header (extHdr)" );
        }
    }
    else {
        count = fwrite( m_Header.txt, 1, 256, file );
        if( count != 256 ) {
            throw DataIOException( "Failed to write ANA header (txt)" );
        }
    }

    // compressed data
    if( m_Header.subf & 1 ) {
        if( system_is_big_endian ) {
            swapEndian( &( m_CompressedHeader.tsize ) );
            swapEndian( &( m_CompressedHeader.nblocks ) );
            swapEndian( &( m_CompressedHeader.bsize ) );
        }

        count = fwrite( reinterpret_cast<char*>( &m_CompressedHeader ), 1, 14, file );
        if( count != 14 ) {
            throw DataIOException( "Failed to write ANA header (m_CompressedHeader)" );
        }

        if( system_is_big_endian ) {
            swapEndian( &( m_CompressedHeader.tsize ) );
            swapEndian( &( m_CompressedHeader.nblocks ) );
            swapEndian( &( m_CompressedHeader.bsize ) );
        }
    }


}


std::shared_ptr<AnaInfo> redux::file::readAnaInfo( const std::string& filename ) {

    std::shared_ptr<AnaInfo> hdr( new AnaInfo( filename ) );
    return hdr;

}

