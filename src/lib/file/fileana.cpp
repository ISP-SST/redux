#include "redux/file/fileana.hpp"

#include "redux/util/endian.hpp"
#include "redux/types.hpp"

#include "redux/file/anacompress.hpp"
#include "redux/file/anadecompress.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

using namespace redux::file;
using namespace redux::util;
using namespace redux;
using namespace std;

#if RDX_BYTE_ORDER == RDX_LITTLE_ENDIAN
static const int system_is_big_endian = 0;
#elif RDX_BYTE_ORDER == RDX_BIG_ENDIAN
static const int system_is_big_endian = 1;
#else
#error RDX_BYTE_ORDER not set
#endif

const uint8_t Ana::typeSizes[] = { 1, 2, 4, 4, 8, 8, 0, 0, 16 };

namespace {
    template <typename T> Ana::TypeIndex getDatyp( void ) { return Ana::ANA_UNDEF; }
    template <> Ana::TypeIndex getDatyp<uint8_t>( void ) { return Ana::ANA_BYTE; }
    template <> Ana::TypeIndex getDatyp<int16_t>( void ) { return Ana::ANA_WORD; }
    template <> Ana::TypeIndex getDatyp<int32_t>( void ) { return Ana::ANA_LONG; }
    template <> Ana::TypeIndex getDatyp<float  >( void ) { return Ana::ANA_FLOAT; }
    template <> Ana::TypeIndex getDatyp<double >( void ) { return Ana::ANA_DOUBLE; }
    template <> Ana::TypeIndex getDatyp<int64_t>( void ) { return Ana::ANA_LONGLONG; }
    template <> Ana::TypeIndex getDatyp<complex_t>( void ) { return Ana::ANA_COMPLEX; }
}

Ana::Ana( void ) : hdrSize( 0 ) {
    memset( &m_Header, 0, sizeof( raw_header ) );
    memset( &m_CompressedHeader, 0, sizeof( compressed_header ) );
    m_Header.datyp = ANA_UNDEF;    // just to make sure no default type is inferred.
}

Ana::Ana( const string& filename ) : hdrSize( 0 ) {
    memset( &m_Header, 0, sizeof( raw_header ) );
    memset( &m_CompressedHeader, 0, sizeof( compressed_header ) );
    read( filename );
}


Ana::~Ana( void ) {

}

void Ana::read( ifstream& file ) {

    hdrSize = 0;
    file.seekg( 0 );

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

    if( m_Header.nhb > 16 ) {
        throw logic_error( "Warning: Ana::read() - extended header is longer than 16 blocks!" );
    }
    else {
        size_t extraTextSize = ( m_Header.nhb - 1 ) * rawSize;
        unique_ptr<char[]> tmpText( new char [ extraTextSize + 257 ] );
        char* txtPtr = tmpText.get();
        memset( txtPtr + 256, 0, extraTextSize + 1 );
        memcpy( txtPtr, m_Header.txt, 256 );
        file.read( txtPtr + 256, extraTextSize );
        if( !file.good() ) {
            if( file.eof() ) {
                throw ios_base::failure( "Failed to read ANA header: <EOF> reached" );
            }
        }
        hdrSize = file.tellg();

        // Hack to strip away null-characters and copies of struct-info from the string
        // needed to correctly read broken headers written with old software.
        // The camera software at the sst used to pad the text string with copies of the
        // primary block, and the momfbd code produced a zero-padded text header.
        txtPtr += 256;  // skip part from first block
        for( int i=0; i<m_Header.nhb-1; ++i, txtPtr+=512) {
            if( !memcmp(txtPtr,reinterpret_cast<char*>( &m_Header ),256) ) {
                memset(txtPtr,0,256);
            }
        }
        txtPtr = tmpText.get(); // reset to beginning
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
        memset( m_Header.txt, 0, 256 );     // clear the .txt field, we store the whole text in m_ExtendedHeader

    }

    // m_Header.dim is always stored as little-endian, so swap endianess if system is big-endian.
    if( system_is_big_endian ) {
        swapEndian( &( m_Header.dim ), m_Header.ndim );
    }

    // compressed data
    if( m_Header.subf & 1 ) {
        size_t chdrSize = 14;
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

void Ana::read( const string& filename ) {
    ifstream file( filename, ifstream::binary );
    read( file );
}

void Ana::write( ofstream& file ) {

    size_t hdrSize = sizeof( struct raw_header );   // = 512 by construction
    boost::trim( m_ExtendedHeader );
    size_t textSize = m_ExtendedHeader.length();
    m_Header.nhb = 1;
    if( textSize > 255 ) {
        m_Header.nhb = static_cast<uint8_t>( 1+(textSize+256)/512 );
        if( m_Header.nhb > 16 ) {
            cout << "Warning: Ana::write() - this = " << hexString(this) << endl;
            cout << "Warning: Ana::write() - textSize = " << textSize << "   nhb = " << (int)m_Header.nhb << endl;
            cout << "Warning: Ana::write() - header text is too long, it will be truncated!! (max = 7935 characters)" << endl;
            cout << "Warning: Ana::write() - text = \"" << m_ExtendedHeader << "\n" << endl;
            m_Header.nhb = 16;
            textSize = 7935;
        }
    }
    
    // file header entries are in little endian order, so swap for big-endian systems.
    if( system_is_big_endian ) {
        m_Header.subf |= 128;   // mark file as written on big-endian machine.
        swapEndian( &( m_Header.synch_pattern ) );
        swapEndian( &( m_Header.cbytes ) );
        swapEndian( &( m_Header.dim ), m_Header.ndim );
    }

    file.write( reinterpret_cast<char*>( &m_Header ), 256 );
    if( !file.good() ) {
        throw ios_base::failure( "Ana::write(): Failed to write m_Header" );
    }

    if( system_is_big_endian ) { // ...and then swap back
        swapEndian( &( m_Header.synch_pattern ) );
        swapEndian( &( m_Header.cbytes ) );
        swapEndian( &( m_Header.dim ), m_Header.ndim );
    }

    if( textSize > 0 ) { // we have header text
        size_t blockSize = ( m_Header.nhb - 1 ) * hdrSize + 256;
        unique_ptr<char[]> extHdr( new char[blockSize] );
        memset( extHdr.get(), 0, blockSize );
        memcpy( extHdr.get(), m_ExtendedHeader.c_str(), min(blockSize,textSize) );
        file.write( extHdr.get(), blockSize );
        if( !file.good() ) {
            throw ios_base::failure( "Ana::write(): Failed to write additional header blocks (text)." );
        }
    }
    else {  // otherwise just write 256 zeroes
        file.write( m_Header.txt, 256 );
        if( !file.good() ) {
            throw ios_base::failure( "Ana::write(): Failed to write text block." );
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
            throw ios_base::failure( "Ana::write(): Failed to write m_CompressedHeader" );
        }

        if( system_is_big_endian ) {
            swapEndian( &( m_CompressedHeader.tsize ) );
            swapEndian( &( m_CompressedHeader.nblocks ) );
            swapEndian( &( m_CompressedHeader.bsize ) );
        }
    }


}


vector<string> redux::file::Ana::makeFitsHeader( const std::string& hdr ) {
    // TODO implement this
    return vector<string>(1,hdr);
}


// TBD: These "parsers" should probably be in the sst namespace with other telescope/hardware specific functions...
size_t redux::file::Ana::getNumberOfFrames(void) {
    
    vector<string> hdrTexts = getText();
    string hdr;
    if( !hdrTexts.empty() ) hdr = hdrTexts.front();
    boost::regex re ("(\\d+)[ .]+SUM[= ]+");
    boost::smatch match;
    if (boost::regex_search (hdr, match, re)) {
        size_t nFrames = boost::lexical_cast<size_t> (match[1]);
        if (nFrames > 1) {
            return nFrames;
        }
    }
    return 1;
    
}

bpx::ptime redux::file::Ana::getStartTime(void) {
    
    vector<string> hdrTexts = getText();
    string hdr;
    if( !hdrTexts.empty() ) hdr = hdrTexts.front();
    bpx::ptime startT;
    
    boost::regex re("Ts=([0-9.\\-: ]+)");
    boost::smatch match;
    if (boost::regex_search (hdr, match, re)) {
        startT = bpx::time_from_string(match[1]);
    }

    return startT;
    
}

bpx::ptime redux::file::Ana::getEndTime(void) {
    
    vector<string> hdrTexts = getText();
    string hdr;
    if( !hdrTexts.empty() ) hdr = hdrTexts.front();
    bpx::ptime endT;
    
    boost::regex re("Te=([0-9.\\-: ]+)");
    boost::smatch match;
    if (boost::regex_search (hdr, match, re)) {
        endT = bpx::time_from_string(match[1]);
    }

    return endT;
    
}

bpx::ptime redux::file::Ana::getAverageTime(void) {

    bpx::ptime startT = getStartTime();
    bpx::time_duration td = (getEndTime()-startT)/2;
    return (startT+td);
    
}

bpx::time_duration redux::file::Ana::getExposureTime(void) {
    
    return (getEndTime()-getStartTime());

}

vector<bpx::ptime> redux::file::Ana::getStartTimes(void) {
    size_t nFrames = getNumberOfFrames();
    if( nFrames == 0 ) return vector<bpx::ptime>();
    bpx::ptime startT = getStartTime();
    vector<bpx::ptime> ret( nFrames, startT );
    if( nFrames == 1 ) return ret;
    bpx::ptime endT = getEndTime();
    bpx::time_duration expTime = (endT-startT);
    for( size_t i=1; i<nFrames; ++i ) {
        ret[i] += expTime;
        expTime += expTime;
    }
    return ret;
}

size_t redux::file::Ana::dataSize(void) { 

    return nElements()*typeSizes[ m_Header.datyp ];
    
}


size_t redux::file::Ana::dimSize(size_t i) {
    
    if( i >= m_Header.ndim ) return 0;
    
    return m_Header.dim[ m_Header.ndim-i-1 ];
    
}


uint8_t redux::file::Ana::elementSize(void) { 
    
    return typeSizes[ m_Header.datyp ];
    
}


size_t redux::file::Ana::nElements(void) {
    
    uint8_t nDims = m_Header.ndim;
    if( nDims == 0 ) return 0;
    
    size_t ret = 1;
    for( uint8_t i=0; i < nDims; ++i ) {
        ret *= m_Header.dim[i];
    }
    
    return ret;
    
}


// IDL type-ID = ANA type-ID + 1
int redux::file::Ana::getIDLType(void) {
    
    return m_Header.datyp + 1;
    
}


void redux::file::Ana::readCompressed( ifstream& file, char* data, size_t nElements, const Ana* hdr ) {

    size_t compressedSize = hdr->m_CompressedHeader.tsize - 14;

    shared_ptr<uint8_t> tmp( new uint8_t[compressedSize+1], []( uint8_t * p ) { delete[] p; } );        // bug in anadecompress makes it go out-of-bounds by 1

    file.read( reinterpret_cast<char*>( tmp.get() ), compressedSize );
    if( !file.good() ) {
        throw ios_base::failure( "Ana::readCompressed(): Failed to read the specified number of bytes." );
    }

    uint32_t nBlocks = hdr->m_CompressedHeader.nblocks;
    uint32_t blockSize = hdr->m_CompressedHeader.bsize;
    // fix a possible problem with chdr.nblocks
    if( blockSize * nBlocks > nElements ) {
        nBlocks = static_cast<uint32_t>( nElements / blockSize );
    }

    // consistency check
    if( hdr->m_CompressedHeader.type % 2 == hdr->m_Header.datyp ) {
        throw logic_error( "Ana::readCompressed(): type-mismatch between primary and compressed headers." );
    }

    uint8_t slice = hdr->m_CompressedHeader.slice_size;

    int ret(0);
    switch( hdr->m_CompressedHeader.type ) {
        case( 0 ):
            ret = anadecrunch( tmp.get(), reinterpret_cast<int16_t*>( data ), slice, blockSize, nBlocks, !system_is_big_endian );
            break;
        case( 1 ):
            ret = anadecrunch8( tmp.get(), reinterpret_cast<int8_t*>( data ), slice, blockSize, nBlocks, !system_is_big_endian );
            break;
        case( 2 ):
            ret = anadecrunchrun( tmp.get(), reinterpret_cast<int16_t*>( data ), slice, blockSize, nBlocks, !system_is_big_endian );
            break;
        case( 3 ):
            ret = anadecrunchrun8( tmp.get(), reinterpret_cast<int8_t*>( data ), slice, blockSize, nBlocks, !system_is_big_endian );
            break;
        case( 4 ):
            ret = anadecrunch32( tmp.get(), reinterpret_cast<int32_t*>( data ), slice, blockSize, nBlocks, !system_is_big_endian );
            break;
        default:
            throw invalid_argument( "Ana::readCompressed(): unrecognized type of compressed data" );
    }
    
    if( ret < 0 ) {
        throw runtime_error( "Ana::readCompressed(): failed to decrunch data." );
    }

}


void redux::file::Ana::readUncompressed( ifstream& file, char* data, size_t nElements, const Ana* hdr ) {

    size_t nBytes = nElements * typeSizes[hdr->m_Header.datyp];

    file.read( data, nBytes );
    if( !file.good() ) {
        throw ios_base::failure( "Ana::readUncompressed(): Failed to read the specified number of bytes." );
    }

    if( hdr->m_Header.subf & 128 ) { // saved on big endian system.
        switch( typeSizes[hdr->m_Header.datyp] ) {
            case  2: swapEndian( reinterpret_cast<int16_t*>( data ), nElements ); break;
            case  4: swapEndian( reinterpret_cast<int32_t*>( data ), nElements ); break;
            case  8: swapEndian( reinterpret_cast<int64_t*>( data ), nElements ); break;
            case 16: swapEndian( reinterpret_cast<int64_t*>( data ), 2*nElements ); break;
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
    out.reset( cdata, []( uint8_t * p ) { delete[] p; } );

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
            if( runlengthflag ) throw invalid_argument( "Ana::compressData: runlength not supported for 32-bit types." );
            else res = anacrunch32( cdata, reinterpret_cast<const int32_t*>( data ), slice, nx, ny, limit, system_is_big_endian );
            break;
        }
        default: throw invalid_argument( "Ana::compressData: Unsupported data type." );
    }

    return res;
}


void redux::file::Ana::read( const string& filename, char* data, shared_ptr<redux::file::Ana>& hdr ) {

    ifstream file( filename, ifstream::binary );
    if( !file.good() ) {
        throw ios_base::failure( "Ana::read(): Failed to open file: " + filename );
    }

    if( !hdr.get() ) {
        hdr.reset( new Ana() );
    }
    hdr->read( file );

    file.seekg( hdr->hdrSize );
    if( !file.good() ) {
        throw ios_base::failure( "Ana::read(): Seek operation failed." );
    }

    size_t nElements = 1;
    for( uint8_t i(0); i < hdr->m_Header.ndim; ++i ) {
        nElements *= hdr->m_Header.dim[i];
    }

    if( hdr->m_Header.subf & 1 ) readCompressed( file, data, nElements, hdr.get() );
    else readUncompressed( file, data, nElements, hdr.get() );

}


void redux::file::Ana::write( const string& filename, const char* data, const shared_ptr<redux::file::Ana> hdr, bool compress, int slice ) {

    ofstream file( filename, ofstream::binary );
    if( !file.good() ) {
        throw ios_base::failure( "Failed to open file for writing: " + filename );
    }

    if( !hdr.get() ) {
        throw invalid_argument( "Ana::write(): The header object is invalid, cannot write file." );
    }

    if( hdr->m_Header.datyp > 5 || hdr->m_Header.datyp < 0 ) {
        throw invalid_argument( "Ana::write(): hdr.datyp is not valid." );
    }

    if( hdr->m_Header.ndim > 16 ) {
        throw invalid_argument( "Ana::write(): the ANA/f0 format does not support more dimensions than 16." );
    }

    size_t nElements = 1;
    for( uint8_t i = 0; i < hdr->m_Header.ndim; ++i ) {
        if( hdr->m_Header.dim[i] < 1 ) {
            throw logic_error( "Ana::write(): dimSize < 1" );
        }
        nElements *= hdr->m_Header.dim[i];
    }

    size_t totalSize = nElements * typeSizes[hdr->m_Header.datyp];

    hdr->m_Header.synch_pattern = MAGIC_ANA;
    hdr->m_Header.subf = 0;
    memset( hdr->m_Header.cbytes, 0, 4 );

    shared_ptr<uint8_t> cData;
    if( compress ) {
        if( hdr->m_Header.ndim == 2 ) {
            size_t compressedSize = compressData( cData, data, nElements, hdr, slice );
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
            throw invalid_argument( "Ana::write(): compression only supported for 2D data." );
        }
    }

    file.seekp( 0 );
    hdr->write( file );

    if( cData ) file.write( reinterpret_cast<char*>( cData.get() ) + 14, totalSize ); // compressed header is already written, so skip 14 bytes.
    else  file.write( data, totalSize );

    if( !file.good() )
        throw ios_base::failure( "Ana::write(): write failed: " + filename  );

}


template <typename T>
void redux::file::Ana::read( const string& filename, redux::util::Array<T>& data, shared_ptr<redux::file::Ana>& hdr ) {

    ifstream file( filename, ifstream::binary );
    if( !file.good() ) {
        throw ios_base::failure( "Ana::read() Failed to open file: " + filename );
    }

    if( !hdr.get() ) {
        hdr.reset( new Ana() );
    }
    hdr->read( file );

    file.seekg( hdr->hdrSize );
    if( !file.good() ) {
        throw ios_base::failure( "Seek operation failed." );
    }

    // f0 stores the dimensions with the fast index first, so swap them before allocating the array
    int nDims = hdr->m_Header.ndim;
    int nArrayDims = data.nDimensions();
    size_t nElements = 1;
    bool forceResize = ( nArrayDims < nDims );
    vector<size_t> dimSizes( nDims, 0 );

    for( int i( 0 ); i < nDims; ++i ) {
        dimSizes[i] = hdr->dimSize(i);
        if( !forceResize && ( dimSizes[i] != data.dimSize( nArrayDims - nDims + i ) ) ) {
            forceResize = true;
        }
        nElements *= dimSizes[i];
    }

    if( forceResize ) {
        data.resize( dimSizes );
    }

    auto tmp = shared_ptr<char>( new char[nElements * typeSizes[hdr->m_Header.datyp]], []( char * p ) { delete[] p; } );
    if( hdr->m_Header.subf & 1 ) {
        readCompressed( file, tmp.get(), nElements, hdr.get() );
    }
    else {
        readUncompressed( file, tmp.get(), nElements, hdr.get() );
    }

    switch( hdr->m_Header.datyp ) {
        case( ANA_BYTE ):   data.template copyFrom<char>( tmp.get() ); break;
        case( ANA_WORD ):   data.template copyFrom<int16_t>( tmp.get() ); break;
        case( ANA_LONG ):   data.template copyFrom<int32_t>( tmp.get() ); break;
        case( ANA_FLOAT ):  data.template copyFrom<float>( tmp.get() ); break;
        case( ANA_DOUBLE ): data.template copyFrom<double>( tmp.get() ); break;
        default: throw ios_base::failure( "Ana::read is not implemented for type: "+to_string(hdr->m_Header.datyp) );
    }

}
template void redux::file::Ana::read( const string& filename, redux::util::Array<int8_t>& data, shared_ptr<redux::file::Ana>& hdr );
template void redux::file::Ana::read( const string& filename, redux::util::Array<uint8_t>& data, shared_ptr<redux::file::Ana>& hdr );
template void redux::file::Ana::read( const string& filename, redux::util::Array<int16_t>& data, shared_ptr<redux::file::Ana>& hdr );
template void redux::file::Ana::read( const string& filename, redux::util::Array<uint16_t>& data, shared_ptr<redux::file::Ana>& hdr );
template void redux::file::Ana::read( const string& filename, redux::util::Array<int32_t>& data, shared_ptr<redux::file::Ana>& hdr );
template void redux::file::Ana::read( const string& filename, redux::util::Array<uint32_t>& data, shared_ptr<redux::file::Ana>& hdr );
template void redux::file::Ana::read( const string& filename, redux::util::Array<int64_t>& data, shared_ptr<redux::file::Ana>& hdr );
template void redux::file::Ana::read( const string& filename, redux::util::Array<uint64_t>& data, shared_ptr<redux::file::Ana>& hdr );
template void redux::file::Ana::read( const string& filename, redux::util::Array<float  >& data, shared_ptr<redux::file::Ana>& hdr );
template void redux::file::Ana::read( const string& filename, redux::util::Array<double >& data, shared_ptr<redux::file::Ana>& hdr );
template void redux::file::Ana::read( const string& filename, redux::util::Array<complex_t >& data, shared_ptr<redux::file::Ana>& hdr );


template <typename T>
void redux::file::Ana::read( const string& filename, redux::image::Image<T>& image, bool metaOnly ) {
    shared_ptr<Ana> hdr = static_pointer_cast<Ana>( image.meta );
    if( !hdr ) {
        hdr.reset( new Ana() );
        image.meta = hdr;
    }
    if( metaOnly ) {
        hdr->read( filename );
    } else {
        read( filename, image, hdr );
    }
}
template void redux::file::Ana::read( const string & filename, redux::image::Image<int8_t>& image, bool );
template void redux::file::Ana::read( const string & filename, redux::image::Image<uint8_t>& image, bool );
template void redux::file::Ana::read( const string & filename, redux::image::Image<int16_t>& image, bool );
template void redux::file::Ana::read( const string & filename, redux::image::Image<uint16_t>& image, bool );
template void redux::file::Ana::read( const string & filename, redux::image::Image<int32_t>& image, bool );
template void redux::file::Ana::read( const string & filename, redux::image::Image<uint32_t>& image, bool );
template void redux::file::Ana::read( const string & filename, redux::image::Image<int64_t>& image, bool );
template void redux::file::Ana::read( const string & filename, redux::image::Image<uint64_t>& image, bool );
template void redux::file::Ana::read( const string & filename, redux::image::Image<float  >& image, bool );
template void redux::file::Ana::read( const string & filename, redux::image::Image<double >& image, bool );
template void redux::file::Ana::read( const string & filename, redux::image::Image<complex_t >& image, bool );


template <typename T>
void redux::file::Ana::write( const string & filename, const redux::util::Array<T>& data, shared_ptr<redux::file::Ana> hdr, int sliceSize ) {

    if( !hdr.get() ) {
        hdr.reset( new Ana() );
    }
    
    hdr->m_Header.datyp = getDatyp<T>();
    if( hdr->m_Header.datyp == ANA_UNDEF ) {
        throw ios_base::failure( "Ana::write is not implemented for type: "+to_string(hdr->m_Header.datyp) );
    }
    
    auto tmpDims = data.dimensions(true);    // TBD: should we always discard dimensions of size 1, or leave it to the user ??
    int nDims = tmpDims.size();

    if( nDims == 0 ) {      // Don't write empty files, ignore silently.
        return;
    }
        
    if( nDims > 16 ) {
        throw invalid_argument( "Ana::write(): the ANA/f0 format does not support more dimensions than 16." );
    }
    
    reverse( tmpDims.begin(), tmpDims.end() );     // ANA store fast dimension first
    hdr->m_Header.ndim = 0;
    size_t nElements = 1;
    for( auto &dim: tmpDims ) {
        hdr->m_Header.dim[hdr->m_Header.ndim++] = dim;
        nElements *= dim;
    }
    
    ofstream file( filename, ifstream::binary );
    if( !file.good() ) {
        throw ios_base::failure( "Failed to open file: " + filename );
    }


    size_t textSize = hdr->m_ExtendedHeader.length();
    size_t totalSize = nElements * typeSizes[hdr->m_Header.datyp];

    hdr->m_Header.synch_pattern = MAGIC_ANA;
    hdr->m_Header.subf = 0;
    hdr->m_Header.nhb = static_cast<uint8_t>( 1 + max( textSize + 256, size_t( 0 ) ) / 512 );
    memset( hdr->m_Header.cbytes, 0, 4 );

    shared_ptr<uint8_t> cData;
    if( sliceSize > 0 ) {
        if( hdr->m_Header.ndim == 2 ) {
            size_t compressedSize = totalSize;
            if( data.dense() ) {
                compressedSize = compressData( cData, reinterpret_cast<const char*>( data.get() ), nElements, hdr, sliceSize );
            } else {
                compressedSize = compressData( cData, reinterpret_cast<const char*>( data.copy().get() ), nElements, hdr, sliceSize );
            }
            if( compressedSize < totalSize ) {
                int tmp = static_cast<int>( compressedSize );
                hdr->m_Header.subf |= 1;
                memcpy( hdr->m_Header.cbytes, &tmp, 4 );
                memcpy( &hdr->m_CompressedHeader, cData.get(), 14 );
                totalSize = compressedSize;
            }
            else {          // compressed data larger than original -> store uncompressed.
                cData.reset();
                memset( &hdr->m_CompressedHeader, 0, 14 );
            }
        }
        else {
            throw invalid_argument( "Ana::write(): compression only implemented for 2D data." );
        }
    }

    hdr->write( file );

    if( cData ) {
        file.write( reinterpret_cast<char*>( cData.get() ) + 14, totalSize ); // compressed header is already written, so skip 14 bytes.
    }
    else {
        if( data.dense() ) {
            file.write( reinterpret_cast<const char*>( data.get() ), totalSize );
        }
        else {
            file.write( reinterpret_cast<const char*>( data.copy().get() ), totalSize );
        }
    }

    if( !file.good() )
        throw ios_base::failure( "Ana::write(): write failed." );


}
template void redux::file::Ana::write( const string&, const redux::util::Array<int8_t>&, shared_ptr<redux::file::Ana>, int );
template void redux::file::Ana::write( const string&, const redux::util::Array<uint8_t>&, shared_ptr<redux::file::Ana>, int );
template void redux::file::Ana::write( const string&, const redux::util::Array<int16_t>&, shared_ptr<redux::file::Ana>, int );
template void redux::file::Ana::write( const string&, const redux::util::Array<uint16_t>&, shared_ptr<redux::file::Ana>, int );
template void redux::file::Ana::write( const string&, const redux::util::Array<int32_t>&, shared_ptr<redux::file::Ana>, int );
template void redux::file::Ana::write( const string&, const redux::util::Array<uint32_t>&, shared_ptr<redux::file::Ana>, int );
template void redux::file::Ana::write( const string&, const redux::util::Array<int64_t>&, shared_ptr<redux::file::Ana>, int );
template void redux::file::Ana::write( const string&, const redux::util::Array<uint64_t>&, shared_ptr<redux::file::Ana>, int );
template void redux::file::Ana::write( const string&, const redux::util::Array<float  >&, shared_ptr<redux::file::Ana>, int );
template void redux::file::Ana::write( const string&, const redux::util::Array<double >&, shared_ptr<redux::file::Ana>, int );
template void redux::file::Ana::write( const string&, const redux::util::Array<complex_t >&, shared_ptr<redux::file::Ana>, int );


template <typename T>
void redux::file::Ana::write( const string & filename, const redux::image::Image<T>& image, int sliceSize ) {
    write( filename, image, static_pointer_cast<redux::file::Ana>( image.meta ), sliceSize );
}
template void redux::file::Ana::write( const string&, const redux::image::Image<int8_t>&, int );
template void redux::file::Ana::write( const string&, const redux::image::Image<uint8_t>&, int );
template void redux::file::Ana::write( const string&, const redux::image::Image<int16_t>&, int );
template void redux::file::Ana::write( const string&, const redux::image::Image<uint16_t>&, int );
template void redux::file::Ana::write( const string&, const redux::image::Image<int32_t>&, int );
template void redux::file::Ana::write( const string&, const redux::image::Image<uint32_t>&, int );
template void redux::file::Ana::write( const string&, const redux::image::Image<int64_t>&, int );
template void redux::file::Ana::write( const string&, const redux::image::Image<uint64_t>&, int );
template void redux::file::Ana::write( const string&, const redux::image::Image<float  >&, int );
template void redux::file::Ana::write( const string&, const redux::image::Image<double >&, int );
template void redux::file::Ana::write( const string&, const redux::image::Image<complex_t >&, int );


template <typename T>
void redux::file::Ana::write( const string & filename, const T* data, size_t n ) {
    
    if( n == 0 ) {
        return;
    }
    
    ofstream file( filename, ifstream::binary );
    if( !file.good() ) {
        throw ios_base::failure( "Failed to open file: " + filename );
    }

    Ana tmpHdr;

    tmpHdr.m_Header.dim[0] = n;
    tmpHdr.m_Header.ndim = 1;
    tmpHdr.m_Header.datyp = getDatyp<T>();
    tmpHdr.m_Header.synch_pattern = MAGIC_ANA;
    tmpHdr.m_Header.subf = 0;
    tmpHdr.m_Header.nhb = 1;

    if( tmpHdr.m_Header.datyp == ANA_UNDEF ) {
        throw ios_base::failure( "Ana::write is not implemented for type: "+to_string(tmpHdr.m_Header.datyp) );
    }
    
    tmpHdr.write( file );

    file.write( reinterpret_cast<const char*>( data ), tmpHdr.m_Header.dim[0]*typeSizes[tmpHdr.m_Header.datyp] );

    if( !file.good() )
        throw ios_base::failure( "Ana::write(): write failed." );


}
template void redux::file::Ana::write( const string&, const int8_t*, size_t n );
template void redux::file::Ana::write( const string&, const uint8_t*, size_t n );
template void redux::file::Ana::write( const string&, const int16_t*, size_t n );
template void redux::file::Ana::write( const string&, const uint16_t*, size_t n );
template void redux::file::Ana::write( const string&, const int32_t*, size_t n );
template void redux::file::Ana::write( const string&, const uint32_t*, size_t n );
template void redux::file::Ana::write( const string&, const int64_t*, size_t n );
template void redux::file::Ana::write( const string&, const uint64_t*, size_t n );
template void redux::file::Ana::write( const string&, const float*, size_t n );
template void redux::file::Ana::write( const string&, const double*, size_t n );
template void redux::file::Ana::write( const string&, const complex_t*, size_t n );




