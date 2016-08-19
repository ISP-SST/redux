#include "redux/file/filefits.hpp"

#ifdef REDUX_WITH_FITS

#include "redux/util/arraystats.hpp"
#include "redux/util/endian.hpp"
#include "redux/types.hpp"

#include <fstream>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

using namespace redux::file;
using namespace redux::util;
using namespace redux;
using namespace std;

#if REDUX_BYTE_ORDER == REDUX_LITTLE_ENDIAN
static const int system_is_big_endian = 0;
#elif REDUX_BYTE_ORDER == REDUX_BIG_ENDIAN
static const int system_is_big_endian = 1;
#else
#error REDUX_BYTE_ORDER not set
#endif

const uint8_t Fits::typeSizes[] = { 1, 2, 4, 4, 8, 8, 0, 0, 16 };

namespace {
    template <typename T> Fits::TypeIndex getDatyp( void ) { return Fits::FITS_NOTYPE; }
    template <> Fits::TypeIndex getDatyp<uint8_t>( void ) { return Fits::FITS_BYTE; }
    template <> Fits::TypeIndex getDatyp<int16_t>( void ) { return Fits::FITS_WORD; }
    template <> Fits::TypeIndex getDatyp<int32_t>( void ) { return Fits::FITS_INT; }
    template <> Fits::TypeIndex getDatyp<float  >( void ) { return Fits::FITS_FLOAT; }
    template <> Fits::TypeIndex getDatyp<double >( void ) { return Fits::FITS_DOUBLE; }
    template <> Fits::TypeIndex getDatyp<int64_t>( void ) { return Fits::FITS_LONG; }
    template <> Fits::TypeIndex getDatyp<complex_t>( void ) { return Fits::FITS_COMPLEX; }
    
    size_t getElementSize( int dataType ) {
        switch( dataType ) {
            case( 11 ):                
            case( 12 ):                
                return 1;
            case( 20 ):                
            case( 21 ):                
                return 2;
            case( 30 ):                
            case( 31 ):                
            case( 42 ):                
                return 4;
            case( 40 ):                
            case( 41 ):                
            case( 82 ):                
                return 8;
            default: ;
        }
        return 0;
    }
    
    
    int getDataType( fitsfile* ff, int bitpix ) {
        
        int status(0);
        long bzero(0);
        float bscale(0);
        if( fits_read_key_flt( ff, "BSCALE", &bscale, nullptr, &status ) ) {
            if( status != KEY_NO_EXIST ) {
                fits_report_error( stderr, status );
                return 0;
            } else status = 0;
        }
    
        if ( bscale == 1.0 ) {
            if( fits_read_key_lng( ff, "BZERO", &bzero, nullptr, &status ) ) {
                if( status != KEY_NO_EXIST ) {
                    fits_report_error( stderr, status );
                    return 0;
                } else status = 0;
            }
        }
        
        switch( bitpix ) {
            case( 8 ): {
                if( bzero == -(1L<<7) ) return TSBYTE;
                return TBYTE;
            }
            case( 16 ): {
                if( bzero == (1L<<15) ) return TUSHORT;
                return TSHORT;
            }
            case( 32 ): {
                if( bzero == (1L<<31)) return TUINT;
                return TINT;
            }
            case( 64 ): {
                if( bzero == (1L<<63) ) return TULONG;
                return TLONG;
            }
            case( -32 ): return TFLOAT;
            case( -64 ): return TDOUBLE;
            default: return 0;
        }
    }
    
}


Fits::Fits( void ) : fitsPtr_(nullptr), status_(0)  {
    
}


Fits::Fits( const std::string& filename ) : fitsPtr_(nullptr), status_(0) {
    
    read( filename );
    
}


Fits::~Fits() {
    
    close();
    
}


void Fits::close(void) {

    if( fitsPtr_ ) {
        if( fits_close_file( fitsPtr_, &status_ ) ) {
            fits_report_error(stderr, status_);
            return;
        }
        fitsPtr_ = nullptr;
    }
}


void Fits::read( const std::string& filename ) {
    
    int nHDU, hduType;
    char data[ FLEN_VALUE ];
    char card[FLEN_CARD];
    memset( data, 0, FLEN_VALUE );
    memset( card, 0, FLEN_CARD );

    if( fitsPtr_ ) {
        close();
    }

    status_ = 0;

    if( fits_open_file( &fitsPtr_, filename.c_str(), READONLY, &status_ ) ) {
       fits_report_error(stderr, status_);
       return;
    }
    
    if( fits_get_num_hdus( fitsPtr_, &nHDU, &status_ ) ) {
       fits_report_error(stderr, status_);
       return;
    }
    
    if ( nHDU > 1 ) {
        cout << "Fits::read() only 1 HDU supported." << endl;
        return;
    }

    if( fits_get_hdu_type( fitsPtr_, &hduType, &status_ ) ) {
        fits_report_error(stderr, status_);
        return;
    }
    
    if ( hduType != IMAGE_HDU ) {
        cout << "Fits::read() only primary/image HDU implemented. " << hduType << endl;
        return;
    }
    if( fits_read_key( fitsPtr_, TINT, "BITPIX", &primaryHDU.bitpix, NULL, &status_ ) ) {
        fits_report_error(stderr, status_);
        return;
    }
    
    //cout << "bitpix: " << primaryHDU.bitpix << endl;
    
    primaryHDU.dataType = getDataType( fitsPtr_, primaryHDU.bitpix );
    primaryHDU.elementSize = getElementSize( primaryHDU.dataType );
    //cout << "dataType: " << primaryHDU.dataType << endl;
    //cout << "elementSize: " << primaryHDU.elementSize << endl;
    
    if( fits_read_key( fitsPtr_, TINT, "NAXIS", &primaryHDU.nDims, NULL, &status_ ) ) {
        fits_report_error(stderr, status_);
        return;
    }
    //cout << "naxis: " << primaryHDU.nDims << endl;
    
    primaryHDU.nElements = 0;
    if( primaryHDU.nDims > 0 ) {
        primaryHDU.nElements = 1;
        primaryHDU.dims.resize(primaryHDU.nDims);
        for(uint16_t i=0; i<primaryHDU.dims.size(); ++i) {
            string key = "NAXIS" + to_string(i+1);

            if( fits_read_key( fitsPtr_, TINT, key.c_str(), &(primaryHDU.dims[i]), NULL, &status_ ) ) {
               fits_report_error(stderr, status_);
               return;
            }
            primaryHDU.nElements *= primaryHDU.dims[i];
        }
        //cout << printArray(primaryHDU.dims,"dims") << endl;
    }
    
    //int nkeys, keypos, hdutype;
    //for (int ii = 1; !(fits_movabs_hdu(fitsPtr_, ii, &hdutype, &status_) ); ii++) 
    //{
        /* get no. of keywords */
        //if (fits_get_hdrpos(fitsPtr_, &nkeys, &keypos, &status_) ) {
            //fits_report_error(stderr, status_);
        //}
        //printf("Header listing for HDU #%d:\n", ii);
        //for (int jj = 1; jj <= nkeys; jj++)  {
        //    if ( fits_read_record(fitsPtr_, jj, card, &status_) )
        //         fits_report_error(stderr, status_);

            //printf("%s\n", card); /* print the keyword card */
        //}
        //printf("END\n\n");  /* terminate listing with END */
    //}

//     if (status_ == END_OF_FILE)   /* status values are defined in fitsioc.h */
//         status_ = 0;              /* got the expected EOF error; reset = 0  */
//     else
//        fits_report_error(stderr, status_);    /* got an unexpected error                */
    
    
}


void Fits::write( ofstream& file ) {

}

// TBD: These "parsers" should probably be in the sst namespace with other telescope/hardware specific functions...
size_t redux::file::Fits::getNumberOfFrames(void) {
    
    std::string hdr = getText();
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

bpx::ptime redux::file::Fits::getStartTime(void) {
    
    std::string hdr = getText();
    bpx::ptime startT;
    
    // TODO

    return startT;
    
}

bpx::ptime redux::file::Fits::getEndTime(void) {
    
    std::string hdr = getText();
    bpx::ptime endT;
    
    // TODO

    return endT;
    
}

bpx::ptime redux::file::Fits::getAverageTime(void) {

    bpx::ptime startT = getStartTime();
    bpx::time_duration td = (getEndTime()-startT)/2;
    return (startT+td);
    
}

bpx::time_duration redux::file::Fits::getExposureTime(void) {
    
    return (getEndTime()-getStartTime());

}


size_t redux::file::Fits::dataSize(void) { 

    return nElements()*elementSize();
    
}


size_t redux::file::Fits::dimSize(size_t i) {
    
    if( static_cast<int>(i) >= primaryHDU.nDims ) return 0;
    
    return primaryHDU.dims[ primaryHDU.nDims-i-1 ];
    
}


uint8_t redux::file::Fits::elementSize(void) { 
    
    return std::abs(primaryHDU.bitpix/8);
    
}


size_t redux::file::Fits::nElements(void) {

    if( primaryHDU.nDims == 0 ) return 0;
    
    size_t ret = 1;
    for( auto& d: primaryHDU.dims ) {
        ret *= d;
    }
    return ret;
    
}

#define IDL_TYP_UNDEF       0
#define IDL_TYP_BYTE            1
#define IDL_TYP_INT             2
#define IDL_TYP_LONG            3
#define IDL_TYP_FLOAT           4
#define IDL_TYP_DOUBLE          5
#define IDL_TYP_COMPLEX         6
#define IDL_TYP_STRING          7
#define IDL_TYP_STRUCT          8
#define IDL_TYP_DCOMPLEX        9
#define IDL_TYP_PTR     10
#define IDL_TYP_OBJREF      11
#define IDL_TYP_UINT        12
#define IDL_TYP_ULONG       13
#define IDL_TYP_LONG64      14
#define IDL_TYP_ULONG64     15

// IDL type-ID = FITS type-ID + 1
int redux::file::Fits::getIDLType(void) {
    switch( primaryHDU.dataType ) {
        case( TSBYTE ): ;
        case( TBYTE ):          return IDL_TYP_BYTE;
        case( TSHORT ):         return IDL_TYP_INT;
        case( TUSHORT ):        return IDL_TYP_UINT;
        case( TINT ):           return IDL_TYP_LONG;
        case( TUINT ):          return IDL_TYP_ULONG;
        case( TLONG ):          return IDL_TYP_LONG64;
        case( TULONG ):         return IDL_TYP_ULONG64;
        case( TFLOAT ):         return IDL_TYP_FLOAT;
        case( TDOUBLE ):        return IDL_TYP_DOUBLE;
        case( TCOMPLEX ):       return IDL_TYP_COMPLEX;
        case( TDBLCOMPLEX ):    return IDL_TYP_DCOMPLEX;
        default: return 0;
    }
}


double redux::file::Fits::getMinMaxMean( const char* data, double* Min, double* Max ){
    
    ArrayStats stats;
    size_t nEl = nElements();
    switch( primaryHDU.dataType ) {
        case( FITS_BYTE ):   stats.getMinMaxMean( data, nEl ); break;
        case( FITS_WORD ):   stats.getMinMaxMean( reinterpret_cast<const int16_t*>(data), nEl ); break;
        case( FITS_INT ):    stats.getMinMaxMean( reinterpret_cast<const int32_t*>(data), nEl ); break;
        case( FITS_FLOAT ):  stats.getMinMaxMean( reinterpret_cast<const float*>(data), nEl ); break;
        case( FITS_DOUBLE ): stats.getMinMaxMean( reinterpret_cast<const double*>(data), nEl ); break;
        default: ;
    }
    if( Min ) *Min = stats.min;
    if( Max ) *Max = stats.max;
    
    return stats.mean;
    
}


void redux::file::Fits::read( std::shared_ptr<redux::file::Fits>& hdr, char* data ) {

    if( !hdr.get() ) {
        return;
    }
    int anynull(0), status(0);
    int ret(0);
    
    bool pg_data(false);
    char card[81];
    memset(card,0,81);
    if( fits_read_keyword( hdr->fitsPtr_, "SOLARNET", card, nullptr, &status ) ) {
        if( status != KEY_NO_EXIST ) {
            fits_report_error(stderr, status);
            return;
        }
        status = 0;
        if( fits_read_key( hdr->fitsPtr_,  TSTRING, "INSTRUME", card, nullptr, &status ) ) {
            if( status != KEY_NO_EXIST ) {
                fits_report_error(stderr, status);
                return;
            }
            status = 0;
        } else {
            string strCard(card);
            if( strCard.find("Chromis-") != string::npos ) {
                pg_data = true;
            }
        }
    }

    switch( hdr->primaryHDU.dataType ) {
        case( TBYTE ): ret = fits_read_img( hdr->fitsPtr_, TBYTE, 1, hdr->nElements(), 0,
            reinterpret_cast<uint8_t*>( data ), &anynull, &status); break;
        case( TSBYTE ): ret = fits_read_img( hdr->fitsPtr_, TSBYTE, 1, hdr->nElements(), 0,
            reinterpret_cast<int8_t*>( data ), &anynull, &status); break;
        case( TSHORT ):
            ret = fits_read_img( hdr->fitsPtr_, TSHORT, 1, hdr->nElements(), 0, reinterpret_cast<int16_t*>( data ), &anynull, &status);
            if( pg_data ) {     // data is actually uint16_t, stored in little-endian form and with the 4 lowest bits zeroed.
                hdr->primaryHDU.dataType = TUSHORT;
                uint16_t* ptr = reinterpret_cast<uint16_t*>( data );
                size_t nEl = hdr->nElements();
                swapEndian( ptr, nEl );
                for( size_t i=0; i<nEl; ++i ) {
                    ptr[i] >>= 4;
                }
            }
            break;
        case( TUSHORT ): ret = fits_read_img( hdr->fitsPtr_, TUSHORT, 1, hdr->nElements(), 0,
            reinterpret_cast<uint16_t*>( data ), &anynull, &status); break;
        case( TINT ): ret = fits_read_img( hdr->fitsPtr_, TINT, 1, hdr->nElements(), 0,
            reinterpret_cast<int32_t*>( data ), &anynull, &status); break;
        case( TUINT ): ret = fits_read_img( hdr->fitsPtr_, TUINT, 1, hdr->nElements(), 0,
            reinterpret_cast<uint32_t*>( data ), &anynull, &status); break;
        case( TLONG ): ret = fits_read_img( hdr->fitsPtr_, TLONG, 1, hdr->nElements(), 0,
            reinterpret_cast<int64_t*>( data ), &anynull, &status); break;
        case( TULONG ): ret = fits_read_img( hdr->fitsPtr_, TULONG, 1, hdr->nElements(), 0,
            reinterpret_cast<uint64_t*>( data ), &anynull, &status); break;
        default: ;
    }
    
    if ( ret ) {
        fits_report_error( stderr, status );
    }

}


void redux::file::Fits::write( const std::string& filename, const char* data, const std::shared_ptr<redux::file::Fits> hdr, bool compress, int slice ) {

    // TODO
    
}


template <typename T>
void redux::file::Fits::read( const string& filename, redux::util::Array<T>& data, std::shared_ptr<redux::file::Fits>& hdr ) {

    if( !hdr.get() ) {
    //cout << "Fits::read: making new hdr" << endl;
        hdr.reset( new Fits() );
    }
    hdr->read( filename );

    // fits stores the dimensions with the fast index first, so swap them before allocating the array
    int nDims = hdr->primaryHDU.nDims;
    //cout << "Fits::read: nDims = " << nDims << endl;
    int nArrayDims = data.nDimensions();
    size_t nElements = 1;
    bool forceResize = ( nArrayDims < nDims );
    std::vector<size_t> dimSizes( nDims, 0 );
    //cout << "Fits::read: filename = " << filename << endl;
    //cout << "Fits::read: " << __LINE__ << printArray(hdr->primaryHDU.dims,"  fdims") << endl;
    //cout << "Fits::read: " << __LINE__ << printArray(data.dimensions(),"  ddims") << endl;
    for( int i( 0 ); i < nDims; ++i ) {
        dimSizes[i] = hdr->primaryHDU.dims[nDims - i - 1];
        if( !forceResize && ( dimSizes[i] != data.dimSize( nArrayDims - nDims + i ) ) ) {
            forceResize = true;
        }
        nElements *= dimSizes[i];
    }
    //cout << "Fits::read: " << __LINE__ << "  forceResize="  << forceResize << endl;
    //cout << "Fits::read: " << __LINE__ << "  nElements="  << nElements << endl;

    if( forceResize ) {
        data.resize( dimSizes );
        //cout << "Fits::read: " << __LINE__ << printArray(data.dimensions(),"  nddims") << endl;
    }

    size_t dataSize = hdr->primaryHDU.nElements * hdr->primaryHDU.elementSize;
    if( dataSize ) {
        //cout << "Fits::read: " << __LINE__ << "  dataSize = "  << dataSize << endl;
        auto tmp = std::shared_ptr<char>( new char[dataSize], []( char * p ) { delete[] p; } );
        read( hdr, tmp.get() );
        switch( hdr->primaryHDU.dataType ) {
            case( TBYTE ):   data.template copyFrom<uint8_t>( tmp.get() ); break;
            case( TSBYTE ):  data.template copyFrom<int8_t>( tmp.get() ); break;
            case( TSHORT ):  data.template copyFrom<int16_t>( tmp.get() ); break;
            case( TUSHORT ): data.template copyFrom<uint16_t>( tmp.get() ); break;
            case( TINT ):    data.template copyFrom<int32_t>( tmp.get() ); break;
            case( TUINT ):   data.template copyFrom<uint32_t>( tmp.get() ); break;
            case( TFLOAT ):  data.template copyFrom<float>( tmp.get() ); break;
            case( TDOUBLE ): data.template copyFrom<double>( tmp.get() ); break;
            default: cerr << "Fits::read: " << __LINE__ << "  unsupported data type: "  << hdr->primaryHDU.dataType << endl;
        }
        //cout << "Fits::read: " << __LINE__ << "  nElements="  << nElements << endl;
    }
}
template void redux::file::Fits::read( const string& filename, redux::util::Array<uint8_t>& data, std::shared_ptr<redux::file::Fits>& hdr );
template void redux::file::Fits::read( const string& filename, redux::util::Array<int16_t>& data, std::shared_ptr<redux::file::Fits>& hdr );
template void redux::file::Fits::read( const string& filename, redux::util::Array<int32_t>& data, std::shared_ptr<redux::file::Fits>& hdr );
template void redux::file::Fits::read( const string& filename, redux::util::Array<int64_t>& data, std::shared_ptr<redux::file::Fits>& hdr );
template void redux::file::Fits::read( const string& filename, redux::util::Array<float  >& data, std::shared_ptr<redux::file::Fits>& hdr );
template void redux::file::Fits::read( const string& filename, redux::util::Array<double >& data, std::shared_ptr<redux::file::Fits>& hdr );
template void redux::file::Fits::read( const string& filename, redux::util::Array<complex_t >& data, std::shared_ptr<redux::file::Fits>& hdr );


template <typename T>
void redux::file::Fits::read( const string& filename, redux::image::Image<T>& image, bool metaOnly ) {
    std::shared_ptr<Fits> hdr = static_pointer_cast<Fits>( image.meta );
    //cout << "Fits::read: " << __LINE__ << "  hdr ="  << hexString(hdr.get()) << "  mO = " << metaOnly << endl;
    if( !hdr ) {
        hdr.reset( new Fits() );
        //cout << "Fits::read: " << __LINE__ << "  new FitsHdr = " << hexString(hdr.get()) << endl;
        image.meta = hdr;
    }
    if( metaOnly ) {
        hdr->read( filename );
    } else {
        read( filename, image, hdr );
    }
    //cout << "Fits::read: " << __LINE__ << "  hdr ="  << hexString(hdr.get()) << "  meta = " << hexString(image.meta.get()) << endl;
}
template void redux::file::Fits::read( const string & filename, redux::image::Image<uint8_t>& image, bool );
template void redux::file::Fits::read( const string & filename, redux::image::Image<int16_t>& image, bool );
template void redux::file::Fits::read( const string & filename, redux::image::Image<int32_t>& image, bool );
template void redux::file::Fits::read( const string & filename, redux::image::Image<int64_t>& image, bool );
template void redux::file::Fits::read( const string & filename, redux::image::Image<float  >& image, bool );
template void redux::file::Fits::read( const string & filename, redux::image::Image<double >& image, bool );
template void redux::file::Fits::read( const string & filename, redux::image::Image<complex_t >& image, bool );


template <typename T>
void redux::file::Fits::write( const string & filename, const redux::util::Array<T>& data, std::shared_ptr<redux::file::Fits> hdr, int sliceSize ) {

    // TODO
    
//    if( !hdr.get() ) {
//        hdr.reset( new Fits() );
//    }
    

}
template void redux::file::Fits::write( const string&, const redux::util::Array<uint8_t>&, std::shared_ptr<redux::file::Fits>, int );
template void redux::file::Fits::write( const string&, const redux::util::Array<int16_t>&, std::shared_ptr<redux::file::Fits>, int );
template void redux::file::Fits::write( const string&, const redux::util::Array<int32_t>&, std::shared_ptr<redux::file::Fits>, int );
template void redux::file::Fits::write( const string&, const redux::util::Array<int64_t>&, std::shared_ptr<redux::file::Fits>, int );
template void redux::file::Fits::write( const string&, const redux::util::Array<float  >&, std::shared_ptr<redux::file::Fits>, int );
template void redux::file::Fits::write( const string&, const redux::util::Array<double >&, std::shared_ptr<redux::file::Fits>, int );
template void redux::file::Fits::write( const string&, const redux::util::Array<complex_t >&, std::shared_ptr<redux::file::Fits>, int );


template <typename T>
void redux::file::Fits::write( const string & filename, const redux::image::Image<T>& image, int sliceSize ) {
    write( filename, image, static_pointer_cast<redux::file::Fits>( image.meta ), sliceSize );
}
template void redux::file::Fits::write( const string&, const redux::image::Image<uint8_t>&, int );
template void redux::file::Fits::write( const string&, const redux::image::Image<int16_t>&, int );
template void redux::file::Fits::write( const string&, const redux::image::Image<int32_t>&, int );
template void redux::file::Fits::write( const string&, const redux::image::Image<int64_t>&, int );
template void redux::file::Fits::write( const string&, const redux::image::Image<float  >&, int );
template void redux::file::Fits::write( const string&, const redux::image::Image<double >&, int );
template void redux::file::Fits::write( const string&, const redux::image::Image<complex_t >&, int );


template <typename T>
void redux::file::Fits::write( const string & filename, const T* data, size_t n ) {
    
    // TODO

}
template void redux::file::Fits::write( const string&, const uint8_t*, size_t n );
template void redux::file::Fits::write( const string&, const int16_t*, size_t n );
template void redux::file::Fits::write( const string&, const int32_t*, size_t n );
template void redux::file::Fits::write( const string&, const int64_t*, size_t n );
template void redux::file::Fits::write( const string&, const float*, size_t n );
template void redux::file::Fits::write( const string&, const double*, size_t n );
template void redux::file::Fits::write( const string&, const complex_t*, size_t n );


#endif  // REDUX_WITH_FITS


