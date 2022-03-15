#include "redux/file/filefits.hpp"

#ifdef RDX_WITH_FITS

#include "redux/util/arraystats.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/ricecompress.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/types.hpp"

#include <fstream>
#include <functional>
#include <mutex>
#include <numeric>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/lexical_cast.hpp>

namespace bfs = boost::filesystem;
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

const uint8_t Fits::typeSizes[] = { 1, 2, 4, 4, 8, 8, 0, 0, 16 };

namespace {
    
    // TODO add direct read/write for compilation without cfitsio
    
    string generate_avc( void ) {
        string ret;
        for( char i(32); i<127; ++i ) ret.push_back(i);
        return ret;
    }
    const string allowed_value_chars = generate_avc();
    const string allowed_key_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ01923456789-_";
    
    const char end_card[] = "END                                                                             ";
    
    const bool use_mutex = (fits_is_reentrant() == 0);
    
    std::mutex fmtx;
    void lockFits( void ) { if(use_mutex) fmtx.lock(); }
    void unlockFits( void ) { if(use_mutex) fmtx.unlock(); }
    struct lockGuard {
        lockGuard(){ lockFits(); }
        ~lockGuard(){ unlockFits(); }
    };
    
    
    template <class T, typename... Args>
    T lockedCall( std::function<T(Args...)> f, Args...  args ) {
        lockGuard();
        return f( args... );
    }

    /* not used at the moment, so commented out to avoid compiler warnings about it
    template <typename T> Fits::TypeIndex getDatyp( void )  { return Fits::FITS_NOTYPE; }
    template <> Fits::TypeIndex getDatyp<uint8_t>( void ) { return Fits::FITS_BYTE; }
    template <> Fits::TypeIndex getDatyp<int16_t>( void ) { return Fits::FITS_WORD; }
    template <> Fits::TypeIndex getDatyp<int32_t>( void ) { return Fits::FITS_INT; }
    template <> Fits::TypeIndex getDatyp<float  >( void ) { return Fits::FITS_FLOAT; }
    template <> Fits::TypeIndex getDatyp<double >( void ) { return Fits::FITS_DOUBLE; }
    template <> Fits::TypeIndex getDatyp<int64_t>( void ) { return Fits::FITS_LONG; }
    template <> Fits::TypeIndex getDatyp<complex_t>( void )  { return Fits::FITS_COMPLEX; }
    */
    
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
    
    string statusString( const string& prefix, int status ) {
       char err_text[80];
       static std::function<void(int,char*)> fn = fits_get_errstatus;
       lockedCall( fn, status, err_text );
       return "Failed to read fits: " + prefix + "  cfitsio reports("+to_string(status)+"): " + string(err_text);
    }
    
    void printStatusError( const string& prefix, int status ) {
       cerr << statusString( prefix, status ) << endl;
    }

    void throwStatusError( const string& prefix, int status ) {
       throw ios_base::failure( statusString( prefix, status ) );
    }
    
    
    template <typename T> int getBitpix(void) { return 0; }
    template<> int getBitpix<uint8_t>(void) { return 8; }
    template<> int getBitpix<int16_t>(void) { return 16; }
    template<> int getBitpix<uint16_t>(void) { return 16; }
    template<> int getBitpix<int32_t>(void) { return 32; }
    template<> int getBitpix<uint32_t>(void) { return 32; }
    template<> int getBitpix<int64_t>(void) { return 64; }
    template<> int getBitpix<uint64_t>(void) { return 64; }
    template<> int getBitpix<float>(void)   { return -32; }
    template<> int getBitpix<double>(void)  { return -64; }
    
    
    template <typename T> long getBzero(void) { return 0; }
    template<> long getBzero<int8_t>(void) { return -(1L<<7); }
    template<> long getBzero<uint16_t>(void) { return 1L<<15; }
    template<> long getBzero<uint32_t>(void) { return 1L<<31; }
    template<> long getBzero<uint64_t>(void) { return 1L<<63; }
    
    
    long getBzero( fitsfile* ff=nullptr ) {
        long bzero(0);
        if( ff ) {
            char card[ FLEN_VALUE ];
            memset( card, 0, FLEN_VALUE );
            int status(0);
            if( fits_read_keyword( ff, "BZERO", card, nullptr, &status ) ) {
                if( status != KEY_NO_EXIST ) {
                    throwStatusError( "Fits::getDataType() BZERO", status );
                }
                if( strlen(card) ) {
                    bzero = stringTo<uint64_t>(card);
                }
            }
        }
        return bzero;
    }
    
    
    int getDataType( int bitpix, long bzero=0 ) {
                
        switch( bitpix ) {
            case( BYTE_IMG ): {
                if( bzero == -(1L<<7) ) return TSBYTE;
                return TBYTE;
            }
            case( SHORT_IMG ): {
                if( bzero == (1L<<15) ) return TUSHORT;
                return TSHORT;
            }
            case( LONG_IMG ): {
                if( bzero == (1L<<31)) return TUINT;
                return TINT;
            }
            case( LONGLONG_IMG ): {
                if( bzero == (1L<<63) ) return TULONG;
                return TLONG;
            }
            case( FLOAT_IMG ): return TFLOAT;
            case( DOUBLE_IMG ): return TDOUBLE;
            default: return 0;
        }
    }
    
    
    string trimStringValue( string str ) {
        boost::trim( str );
        size_t first = str.find_first_of("'");
        if( first != string::npos ) {
            int len = max((int)str.find_last_of("'")-(int)first-1,0);
            str = str.substr(first+1,len);
        }
        boost::replace_all( str, "''",  "'" );
        return str;
    }
    
    
    string makeKey( string key ) {  // return an 8-character string which is permitted a keyword in a FITS header.
        boost::trim( key );                                                 // reamove leading/trailing spaces
        boost::to_upper( key );                                             // make uppercase
        size_t found = key.find_first_not_of( allowed_key_chars );          // look for illegal characters
        if( (found != string::npos) && (found < 8) ) {
            throw std::domain_error("makeKey: \""+key+"\" contains illegal characters.");
        }
        key.resize( 8, ' ' );                                               // pad with spaces, or truncate, to 8 characters
        return key;
    }
    
    
    struct cicomp {  // case-insensitive comparator for the maps below.
        bool operator() ( const std::string& a, const std::string& b ) const { return redux::util::nocaseLess(a,b); }
    };
    const std::set<string, cicomp> multiKeys = { "        ", "COMMENT ", "HISTORY ", "CONTINUE" };
    
}


Fits::Fits( void ) : fitsPtr_(nullptr), status_(0)  {
    
}


Fits::Fits( const string& filename ) : fitsPtr_(nullptr), status_(0) {
    
    read( filename );
    
}


Fits::~Fits() {
    
    try {
        close();
    } catch( const std::exception& e ) {
        cerr << "Exception caught in Fits destructor: " << e.what() << endl;
    }
    
}


void Fits::close(void) {

    if( fitsPtr_ ) {
        status_ = 0;
        if( fits_close_file( fitsPtr_, &status_ ) ) {
            throwStatusError( "Fits:close:  Error: ", status_ );
        }
        fitsPtr_ = nullptr;
    }
}


void Fits::open( const std::string& filename ) {
    
    try {
        if( fitsPtr_ ) {
            if( filename != fitsPtr_->Fptr->filename ) {
                close();
            } else {
                int hdutype(ANY_HDU);
                status_ = 0;
                if( fits_movabs_hdu( fitsPtr_, 1, &hdutype, &status_ ) ) {
                    throwStatusError( "Fits:reopen: "+filename, status_ );
                }
                return;
            }
        }
    } catch( const std::exception& e ) {
        cerr << "Exception while re-opening fitsfile: " << e.what() << endl;
    }
    
    status_ = 0;
    if( fits_open_file( &fitsPtr_, filename.c_str(), READONLY, &status_ ) ) {
        throwStatusError( "Fits:open: "+filename, status_ );
    }
    
}    


namespace {

    void readRawHDU( const string& fn, size_t offset, Fits::ascii_hdu& hdu ) {
        
        // fixed values for an ASCII table.
        hdu.bitpix = 8;
        hdu.nDims = 2;
        //hdu.dataType = TSBYTE;
        hdu.dataType = TBYTE;
        hdu.elementSize = 1;
        hdu.dims.resize(2);
        
        char card[FLEN_CARD];
        memset( card, 0, FLEN_CARD );
        hdu.cards.clear();
        
        ifstream ifs( fn, ifstream::binary );
        ifs.seekg( offset );
        
        size_t count(0);
        while( ifs.good() ) {
            ifs.read( card, 80 );
            if( !strncmp(card, "END", 3) ) {
                break;
            }
            if( !strncmp( card, "TFIELDS", 7 ) ) {
                hdu.nColumns = atol( card+10 );
            }
            count += 80;
            hdu.cards.push_back(card);
        }
        if( hdu.nColumns > 0 ) {
            hdu.table_info.resize( hdu.nColumns );
        }
        for( const string& c: hdu.cards ) {
            if( !c.compare(0,6,"NAXIS1") ) {
                hdu.dims[0] = atol( c.substr(10).c_str() );
            }
            if( !c.compare(0,6,"NAXIS2") ) {
                hdu.dims[1] = atol( c.substr(10).c_str() );
            }
            if( !c.compare(0,5,"TTYPE") ) {
                int id = atoi( c.substr(5,3).c_str() );
                if( id > 0 && id <= hdu.nColumns ) {
                    hdu.table_info[id-1].columnName = trimStringValue( c.substr(10) );
                }
            }
            if( !c.compare(0,5,"TBCOL") ) {
                int id = atoi( c.substr(5,3).c_str() );
                if( id > 0 && id <= hdu.nColumns ) {
                    hdu.table_info[id-1].columnStart = atoi( c.substr(10).c_str() );
                }
            }
            if( !c.compare(0,5,"TFORM") ) {
                int id = atoi( c.substr(5,3).c_str() );
                if( id > 0 && id <= hdu.nColumns ) {
                    hdu.table_info[id-1].columnFormat = trimStringValue( c.substr(10) );
                }
            }
            if( !c.compare(0,5,"TUNIT") ) {
                int id = atoi( c.substr(5,3).c_str() );
                if( id > 0 && id <= hdu.nColumns ) {
                    hdu.table_info[id-1].columnUnit = trimStringValue( c.substr(10) );
                }
            }
        }
        
        
        
        hdu.nElements = hdu.dims[0]*hdu.dims[1];
        
        offset += (count/2880+1)*2880;      // move to data-block
        hdu.data.resize( hdu.dims[1], hdu.dims[0] );
        ifs.seekg( offset );
        ifs.read( hdu.data.get(), hdu.nElements );

    }
    
    void readPrimaryHDU( fitsfile* ff, Fits::hdu& hdu, const string& fn ){
        
        int status(0);
        if( fits_read_key( ff, TINT, "BITPIX", &hdu.bitpix, NULL, &status ) ) {
            throwStatusError( "primaryHDU:BITPIX:"+fn, status );
        }
        
        long bzero = getBzero( ff );
        hdu.dataType = getDataType( hdu.bitpix, bzero );
        hdu.elementSize = getElementSize( hdu.dataType );
        
        if( fits_read_key( ff, TINT, "NAXIS", &hdu.nDims, NULL, &status ) ) {
            throwStatusError( "primaryHDU:NAXIS:"+fn, status );
        }

        if ( hdu.nDims > 999 ) {
            throw logic_error( "Fits::read() NAXIS="+to_string(hdu.nDims)+".  file:" + fn );
        }
        
        hdu.nElements = 0;
        if( hdu.nDims > 0 ) {
            hdu.nElements = 1;
            hdu.dims.resize(hdu.nDims);
            for(uint16_t i=0; i<hdu.dims.size(); ++i) {
                string key = "NAXIS" + to_string(i+1);
                if( fits_read_key( ff, TINT, key.c_str(), &(hdu.dims[i]), NULL, &status ) ) {
                    throwStatusError( "primaryHDU:"+key+":"+fn, status );
                }
                hdu.nElements *= hdu.dims[i];
            }
        }

    }
    
    void readImageHDU( fitsfile* ff, Fits::image_hdu& hdu, const string& fn ){
        
        readPrimaryHDU( ff, dynamic_cast<Fits::hdu&>(hdu), fn );

    }
    
    void readAsciiHDU( fitsfile* ff, Fits::ascii_hdu& hdu, const string& fn ){
        
        readPrimaryHDU( ff, dynamic_cast<Fits::hdu&>(hdu), fn );
        
        // fixed values for an ASCII table.
        if( hdu.bitpix != 8 ) {
            throw ios_base::failure("Fits ASCII-table has wrong BITPIX value " + to_string(hdu.bitpix) + " (should be 8)." );
        }
        if( hdu.nDims != 2 ) {
            throw ios_base::failure("Fits ASCII-table has wrong NAXIS value " + to_string(hdu.nDims) + " (should be 2)." );
        }

        hdu.dataType = TBYTE;
        hdu.elementSize = 1;
        
        char card[FLEN_CARD];
        memset( card, 0, FLEN_CARD );
        hdu.cards.clear();

        int tmpInt;
        int status(0);
        if( fits_read_key( ff, TINT, "TFIELDS", &tmpInt, NULL, &status ) ) {
            throwStatusError( "asciiHDU:TFIELDS:"+fn, status );
        }
        hdu.nColumns = tmpInt;
        hdu.table_info.clear();
        
        for( uint16_t i=1; i<=hdu.nColumns; ++i ) {
            Fits::ascii_hdu::table_info_t ti;
            string key = "TTYPE" + to_string(i);
            if( fits_read_key( ff,  TSTRING, key.c_str(), card, nullptr, &status ) ) {
                throwStatusError( "asciiHDU:"+key+":"+fn, status );
            }
            ti.columnName = card;
            key = "TBCOL" + to_string(i);
            if( fits_read_key( ff, TINT, key.c_str(), &tmpInt, NULL, &status ) ) {
                throwStatusError( "asciiHDU:"+key+":"+fn, status );
            }
            ti.columnStart = tmpInt;
            key = "TFORM" + to_string(i);
            if( fits_read_key( ff,  TSTRING, key.c_str(), card, nullptr, &status ) ) {
                throwStatusError( "asciiHDU:"+key+":"+fn, status );
            }
            ti.columnFormat = card;
            key = "TUNIT" + to_string(i);
            if( fits_read_key( ff,  TSTRING, key.c_str(), card, nullptr, &status ) ) {
                throwStatusError( "asciiHDU:"+key+":"+fn, status );
            }
            ti.columnUnit = card;
            hdu.table_info.push_back(ti);
        }

        if( hdu.nElements == 0 ) return;

        hdu.data.resize( hdu.dims[1], hdu.dims[0] );
        hdu.data.zero();
        unsigned char* dataPtr = reinterpret_cast<unsigned char*>(hdu.data.get());

        if( ffgtbb( ff, 1, 1, hdu.nElements, dataPtr, &status ) ) {
            throwStatusError( "Fits ASCII-table read: ", status );
        }

    }
    
    void readBinaryHDU( fitsfile* ff, Fits::binary_hdu& hdu, const string& fn ){
        
        readPrimaryHDU( ff, dynamic_cast<Fits::hdu&>(hdu), fn );

    }
    
    void readAllCards( fitsfile* ff, Fits::hdu& hdu, const string& fn ){
        int status(0);
        int nkeys, keypos;
        char card[FLEN_CARD];
        memset( card, 0, FLEN_CARD );
        hdu.cards.clear();
        if ( fits_get_hdrpos( ff, &nkeys, &keypos, &status ) ) {
            throwStatusError( "readAllCards:nkeys:"+fn, status );
        }
        for (int jj = 1; jj <= nkeys; jj++)  {
            if ( fits_read_record(ff, jj, card, &status) ) {
                throwStatusError( "readAllCards:"+to_string(jj)+":"+fn, status );
            }
            memset( card+strlen(card), 32, FLEN_CARD-1-strlen(card) ); // make sure we copy all 80 chars
            hdu.cards.push_back(card);
        }
    }
    
}


void Fits::read( const string& filename ) {
    
    int nHDU, hduType;
    char data[ FLEN_VALUE ];
    char card[FLEN_CARD];
    memset( data, 0, FLEN_VALUE );
    memset( card, 0, FLEN_CARD );

    open( filename );
    
    status_ = 0;
    if( fits_get_num_hdus( fitsPtr_, &nHDU, &status_ ) ) {
        throwStatusError( "read:num_hdus:"+filename, status_ );
    }
    
    if( fits_get_hdu_type( fitsPtr_, &hduType, &status_ ) ) {
        throwStatusError( "read:hdu_type:"+filename, status_ );
    }
    
    if ( hduType != IMAGE_HDU ) {
        throw logic_error( "Fits::read() primary HDU has wrong hduType.  file:" + filename );
    }

    readPrimaryHDU( fitsPtr_, primaryHDU, filename );
    readAllCards( fitsPtr_, primaryHDU, filename );
    
    int hdutype(ANY_HDU);
    for( int ii=2; ; ++ii ) {
        status_ = 0;
        fits_movabs_hdu( fitsPtr_, ii, &hdutype, &status_ );
        if( status_ == END_OF_FILE ) {
            status_ = 0;              // reached EOF, return
            break;
        } else if( status_ == BAD_BITPIX ) {  // this is most likely one of the non-conforming early chromis ascii-tables
            status_ = 0;
            auto hdu = make_shared<ascii_hdu>();
            size_t offset = fitsPtr_->Fptr->headstart[ii-1];
            if( offset > 0 ) {
                close();
                readRawHDU( filename, offset, *hdu );   //read the HDU with direct-IO, sice cfitsio can't handle non-conforming HDUs
                extHDUs.push_back( hdu );
                open( filename );
                fits_movabs_hdu( fitsPtr_, ii, &hdutype, &status_);
                status_ = 0;
            }
            break;
        } else if( status_ ){
            throwStatusError( "read:movabs_hdu:"+to_string(ii)+":"+filename, status_ );      // got an unexpected error   
        }

        switch( hdutype ) {
            case IMAGE_HDU: {
                auto hdu = make_shared<Fits::image_hdu>();
                readImageHDU( fitsPtr_, *hdu, filename );
                readAllCards( fitsPtr_, *hdu, filename );
                if( fitsPtr_->Fptr->compressimg ) {
                    primaryHDU.dHDU = ii;
                    primaryHDU.bitpix = fitsPtr_->Fptr->zbitpix;
                    primaryHDU.nDims = fitsPtr_->Fptr->zndim;
                    primaryHDU.dims.assign( fitsPtr_->Fptr->znaxis, fitsPtr_->Fptr->znaxis+fitsPtr_->Fptr->zndim );
                    primaryHDU.nElements = 1;
                    for( auto d: primaryHDU.dims ) {
                        primaryHDU.nElements *=d;
                    }
                }
                extHDUs.push_back( hdu );
                break;
            }
            case ASCII_TBL: {
                auto hdu = make_shared<Fits::ascii_hdu>();
                readAsciiHDU( fitsPtr_, *hdu, filename );
                readAllCards( fitsPtr_, *hdu, filename );
                extHDUs.push_back( hdu );
                break;
            }
            case BINARY_TBL: {
                auto hdu = make_shared<Fits::binary_hdu>();
                readBinaryHDU( fitsPtr_, *hdu, filename);
                readAllCards( fitsPtr_, *hdu, filename );
                extHDUs.push_back( hdu );
                break;
            }
            default:
                cout << __FILE__ << ":" << __LINE__ << "   Extension is not recognized..." << endl;

        }
    }

    if( status_ == END_OF_FILE ) {
        status_ = 0;              // got the expected EOF error; reset = 0 
    } else if( status_ ) {
        throwStatusError( "read:eof:"+filename, status_ );      // got an unexpected error   
    }

}


void Fits::write( ofstream& file ) {

}


vector<string> Fits::getText( bool raw ) {
    vector<string> ret;
    string hduText;
    if( !raw && primaryHDU.dHDU && extHDUs[primaryHDU.dHDU-2]) {
        for( string k: extHDUs[primaryHDU.dHDU-2]->cards ) {
            string kk = k;
            string key = k.substr(0,8);
            if( boost::iequals( key, "ZBITPIX ") ) {
                key = "BITPIX  ";
                k.replace(0,8,key);
                updateCard( primaryHDU.cards, key, k );
            } else if( boost::iequals( key, "ZNAXIS  ") ) {
                key = "NAXIS   ";
                k.replace(0,8,key);
                updateCard( primaryHDU.cards, key, k );
            } else if( boost::iequals( key.substr(0,6), "ZNAXIS") ) {
                key = "NAXIS" + k.substr(6,1) + "  ";
                int ind = atoi(k.substr(6,1).c_str());
                k.replace(0,8,key);
                if( !updateCard( primaryHDU.cards, key, k ) ) {
                    string key2 = "NAXIS";
                    if( ind>1 ) key2 += to_string(ind-1) + "  ";
                    insertCardAfter( primaryHDU.cards, k, key2);
                }
            }
        }
        for( auto& k: primaryHDU.cards ) {
            hduText += k;
        }
        hduText += end_card;
        ret.push_back(hduText);
        for( size_t i=0; i<extHDUs.size(); ++i ) {
            if( ((int)i != (primaryHDU.dHDU-2)) && extHDUs[i] ) {
                hduText.clear();
                for( auto& k: extHDUs[i]->cards ) {
                    hduText += k;
                }
                hduText += end_card;
                ret.push_back(hduText);
            }
        }
    } else {
        for( auto& k: primaryHDU.cards ) {
            hduText += k;
        }
        hduText += end_card;
        ret.push_back(hduText);
        for( size_t i=0; i<extHDUs.size(); ++i ) {
            if( extHDUs[i] ) {
                hduText.clear();
                for( auto& k: extHDUs[i]->cards ) {
                    hduText += k;
                }
                hduText += end_card;
                ret.push_back(hduText);
            }
        }
    }
    return ret;
}


namespace redux {
    namespace file {

        template <typename T>
        string Fits::makeValue ( const T& t, bool ) {
            return alignRight( to_string(t), 20 );
        }
        template <>
        string Fits::makeValue( const string& str, bool quote ) {
            string ret = boost::trim_right_copy( str );                         // reamove trailing spaces
            size_t found = ret.find_first_not_of( allowed_value_chars );        // look for illegal characters
            if( found != string::npos ) {
                throw std::domain_error("Fits::makeValue: \""+ret+"\" contains illegal characters.");
            }
            if( quote ) ret = "'" + boost::replace_all_copy( ret, "'", "''" ) + "'";
            return ret;
        }
        template <>
        string Fits::makeValue( const bpx::ptime& pt, bool quote ) {
            string ret;
            if( pt.time_of_day() == bpx::time_duration() ) {     // Time is 00:00:00 (i.e. unset), so just print the date.
                ret = to_iso_extended_string( pt.date() );
            } else {
                ret = bpx::to_iso_extended_string(pt);
            }
            if( quote ) ret = "'" + ret + "'";
            return ret;
        }
        template <>
        string Fits::makeValue( const bool& b, bool ) {
            return alignRight( (b?"T":"F"), 20 );
        }
        
        string Fits::makeValue( const char* s, bool quote ) {
            return makeValue( string(s), quote );
        }

        template <typename T>
        T Fits::getValue( const std::string& v ) {
            return redux::util::stringTo<T>( boost::trim_copy( v ) );
        }
        template <>
        bpx::ptime Fits::getValue( const std::string& v ) {
            return stringTo<bpx::ptime>( boost::replace_all_copy( v, "'",  "" ) );
        }
        template <>
        string Fits::getValue( const std::string& v ) {
            string ret = boost::trim_copy( v );
            size_t first = ret.find_first_of("'");
            if( first != string::npos ) {
                string::const_iterator vBeg = v.cbegin()+first;
                string::const_iterator vEnd = vBeg;
                while( ++vEnd != v.cend() ) {
                    if( *vEnd == '\'' ) {
                        if( ((vEnd+1) != v.cend()) && *(vEnd+1) == '\'' ) {
                            vEnd++;
                            continue;
                        }
                        break;
                    }
                }
                // if string is enclosed by single quotes, exclude them.
                if( (vBeg != vEnd) && (*vBeg == '\'') && (*vEnd == '\'') ) {
                    vBeg++;
                    vEnd--;
                }
                ret = string( vBeg, vEnd+1 );
                
            }
            boost::replace_all( ret, "''",  "'" );  // replace double occurances with single ones
            return ret;
        }

        template <typename T>
        string Fits::makeCard( const string& key, const T& v, const string& comment, bool ) {
            string ret = makeKey(key) + "= " + makeValue ( v );
            if( comment != "" ) {
                ret += " / " + makeValue( comment, false );
            }
            ret.resize( 80, ' ' );       // pad with spaces, or truncate, to 80 characters
            return ret;
        }
        template <>
        string Fits::makeCard( const string& key, const string& v, const string& comment, bool allowLong ) {
            string ret = makeKey(key);
            string comment_separator = " / ";
            if( !v.empty() ) {
                ret += "= " + makeValue( v );
                if( !allowLong && (ret.size() > 80) ) ret[79] = '\'';               // if the string is too long, we need to add the end-quote before character 80.
                if( ret.size() < 30 ) ret.resize( 30, ' ' );                        // pad with spaces
            } else {
                comment_separator = "";
            }
            if( comment != "" ) {
                ret += comment_separator + makeValue( comment, false );
            }
            if( !allowLong || ret.size()<80 ) ret.resize( 80, ' ' );       // pad with spaces, or truncate, to 80 characters
            return ret;
        }
        template <>
        string Fits::makeCard( const string& key, const bpx::ptime& date, const string& comment, bool allowLong ) {
            return makeCard( key, bpx::to_iso_extended_string(date), comment, allowLong );
        }
        template <>
        string Fits::makeCard( const string& key, const bool& v, const string& comment, bool ) {
            string ret = makeKey(key) + "= " + makeValue(v);
            if( comment != "" ) {
                ret += " / " + makeValue( comment, false );
            }
            ret.resize( 80, ' ' );       // pad with spaces, or truncate, to 80 characters
            return ret;
        }
        


    }
}
template string Fits::makeValue<int8_t>( const int8_t&, bool );
template string Fits::makeValue<uint8_t>( const uint8_t&, bool );
template string Fits::makeValue<int16_t>( const int16_t&, bool );
template string Fits::makeValue<uint16_t>( const uint16_t&, bool );
template string Fits::makeValue<int32_t>( const int32_t&, bool );
template string Fits::makeValue<uint32_t>( const uint32_t&, bool );
template string Fits::makeValue<int64_t>( const int64_t&, bool );
template string Fits::makeValue<uint64_t>( const uint64_t&, bool );
template string Fits::makeValue<float>( const float&, bool );
template string Fits::makeValue<double>( const double&, bool );

template bool Fits::getValue<bool>( const string& );
template int8_t Fits::getValue<int8_t>( const string& );
template uint8_t Fits::getValue<uint8_t>( const string& );
template int16_t Fits::getValue<int16_t>( const string& );
template uint16_t Fits::getValue<uint16_t>( const string& );
template int32_t Fits::getValue<int32_t>( const string& );
template uint32_t Fits::getValue<uint32_t>( const string& );
template int64_t Fits::getValue<int64_t>( const string& );
template uint64_t Fits::getValue<uint64_t>( const string& );
template float Fits::getValue<float>( const string& );
template double Fits::getValue<double>( const string& );

template string Fits::makeCard( const string&, const int8_t&, const string&, bool );
template string Fits::makeCard( const string&, const uint8_t&, const string&, bool );
template string Fits::makeCard( const string&, const int16_t&, const string&, bool );
template string Fits::makeCard( const string&, const uint16_t&, const string&, bool );
template string Fits::makeCard( const string&, const int32_t&, const string&, bool );
template string Fits::makeCard( const string&, const uint32_t&, const string&, bool );
template string Fits::makeCard( const string&, const int64_t&, const string&, bool );
template string Fits::makeCard( const string&, const uint64_t&, const string&, bool );
template string Fits::makeCard( const string&, const float&, const string&, bool );
template string Fits::makeCard( const string&, const double&, const string&, bool );




string Fits::makeCard( const string& key, const char* v, const string& comment, bool allowLong ) {
    return makeCard( key, string(v), comment, allowLong );
}

string Fits::makeCard( const string& key, const string& comment, bool allowLong ) {
    string ret = makeKey(key);
    if( comment != "" ) {
        ret += makeValue( comment, false );
    }
    if( !allowLong || ret.size()<80 ) ret.resize( 80, ' ' );       // pad with spaces, or truncate, to 80 characters
    return ret;
}


string Fits::getCard( const vector<string>& cards, string key ) {
    key.resize( 8, ' ' );       // pad with spaces, or truncate, to 8 characters
    for( const auto& c: cards ) {
        if( boost::iequals( c.substr( 0, key.length() ), key ) ) {
            return c;
        }
    }
    return "";
}


void Fits::addCard( vector<string>& hdr, const string& card ) {
    
    string key = makeKey( card.substr(0,8) );
    if( multiKeys.count(key) == 0 ) {   // only allow single occurrances for non-multikeys
        for( auto& k: hdr ) {
            if( key.compare( k.substr(0,8) ) == 0 ) {
                return;
            }
        }
    }
    hdr.push_back(card);
    
}


void Fits::removeCards( vector<string>& hdr, string key ) {
    key = makeKey( key );
    hdr.erase( std::remove_if( hdr.begin(), hdr.end(), [&](string a) { return key.compare( a.substr(0,8) ) == 0; }), hdr.end() );
}


void Fits::insertCard( vector<string>& hdr, const string& card, size_t location ) {
    
    if( location < hdr.size() ) {
        hdr.insert( hdr.begin()+location, card );
    } else {
        hdr.push_back(card);
    }
    
}


void Fits::insertCardAfter( vector<string>& cards, const string& card, string after ) {
    string akey = makeKey( after );
    size_t location = string::npos;
    for( size_t i=0; i<cards.size(); ++i ) {
        if( akey.compare( cards[i].substr(0,8) ) == 0 ) {
            location = i+1;
        }
    }
    insertCard( cards, card, location );
}


void Fits::insertCardBefore( vector<string>& cards, const string& card, string before ) {
    string bkey = makeKey( before );
    size_t location = string::npos;
    for( size_t i=0; i<cards.size(); ++i ) {
        if( bkey.compare( cards[i].substr(0,8) ) == 0 ) {
            location = i;
            break;
        }
    }
    insertCard( cards, card, location );
    
}


bool Fits::updateCard( vector<string>& cards, size_t location, const string& card ) {
    
    if( location < cards.size() ) {
        cards[location] = card;
        return true;
    }
    return false;
    
}


bool Fits::updateCard( vector<string>& cards, string key, const string& card ) {
    
    key = makeKey( key );
    size_t location = string::npos;
    if( multiKeys.count(key) == 0 ) {   // never overwrite keywords which might have multiple entries
        for( size_t i=0; i<cards.size(); ++i ) {
            if( key.compare( cards[i].substr(0,8) ) == 0 ) {
                location = i;
                break;
            }
        }
    }
    return updateCard( cards, location, card );

}


bool Fits::updateCard( vector<string>& cards, const string& card ) {
    
    string key = makeKey( card.substr(0,8) );
    return updateCard( cards, key, card );

}


bool Fits::emplaceCard( vector<string>& hdr, const string& key, const string& card ) {
    
    bool updated = updateCard( hdr, key, card );
    if( !updated ) insertCard( hdr, card );
    return !updated;

}


bool Fits::emplaceCard( vector<string>& hdr, const string& card ) {
    
    bool updated = updateCard( hdr, card );
    if( !updated ) insertCard( hdr, card );
    return !updated;

}


void Fits::splitCard( const std::string& card, std::string& key, std::string& value, std::string& comment ) {
    
    key = card.substr( 0, 8 );
    
    if( (card[8] == '=') || (card[9] == ' ') ) { // Value field exists
        size_t pos = card.find_first_not_of( " ", 9 );
        if( (pos != string::npos) && (card[pos] == '\'') ) {    // it's a string value, find beginning & end
            string::const_iterator vBeg = card.cbegin()+pos;
            string::const_iterator vEnd = vBeg;
            while( ++vEnd != card.cend() ) {
                if( *vEnd == '\'' ) {
                    if( ((vEnd+1) != card.cend()) && *(vEnd+1) == '\'' ) {
                        vEnd++;
                        continue;
                    }
                    break;
                }
            }
            if( vEnd >= card.cend() ) {                         // end of string not located, assume it's the entire card.

            }
            vEnd++;
            value = string( vBeg, vEnd );
            comment = string( vEnd, card.cend() );
        } else {
            value = card.substr( 10, 20 );
            comment = card.substr( 30 );
        }
    } else {    // all after the key is comment
        value.clear();
        comment = card.substr( 8 );
    }
    
    size_t pos = comment.find_first_not_of( " /" );
    if( pos != string::npos ) comment = comment.substr(pos);
    
    boost::trim( key );
    boost::trim( value );
    boost::trim( comment );
    
}


vector<string> Fits::splitLongString( const string& card ) {
    
    vector<string> ret;

    if( (card.size() <= 80) || (card[8] != '=') || (card[9] != ' ') ) { // Nothing to do, either too short or does not have a value field
    cout <<"splitLongString: " << __LINE__ << "  cS=" << card.size() << "  C8=\"" << card[8] << "\"  C9=\"" << card[9] << "\"" << endl;
        ret.push_back( card );
        ret[0].resize( 80, ' ' );
        return ret;
    }
    
    size_t pos = card.find_first_not_of( " ", 9 );
    if( (pos == string::npos) || (card[pos] != '\'') ) {                // does not meet the FITS standard, truncate and return
    cout <<"splitLongString: " << __LINE__ << "  pos=" << pos << "  CP=\"" << card[pos] << "\"" << endl;
        ret.push_back( card );
        ret[0].resize( 80, ' ' );
        return ret;
    }
    
    string::const_iterator vBeg = card.cbegin()+pos;
    string::const_iterator vEnd = vBeg;
    while( ++vEnd != card.cend() ) {
        if( *vEnd == '\'' ) {
            if( ((vEnd+1) != card.cend()) && *(vEnd+1) == '\'' ) {
                vEnd++;
                continue;
            }
            break;
        }
    }
    
    if( vEnd >= card.cend() ) {                                         // end of string not located.
    cout <<"splitLongString: " << __LINE__ << "  vEnd=" << (vEnd-card.cbegin()) << "  cEnd=" << (card.cend()-card.cbegin()) << endl;
        ret.push_back( card );
        ret[0].resize( 80, ' ' );
        return ret;
       
    }
    
    string key = card.substr( 0, 8 );
    string value( vBeg, vEnd+1 );
    string tmpCard( card.cbegin(), vEnd+1 );
    const string cs = "&'CONTINUE  '";
    while( tmpCard.size() > 80 ) {
    cout <<"splitLongString: " << __LINE__ << "  tmpCard=\"" << tmpCard.substr(0,80) << "\"" << endl;
        tmpCard.insert( 78, cs );
    cout <<"splitLongString: " << __LINE__ << "  tmpCard=\"" << tmpCard.substr(0,80) << "\"" << endl<<endl;
        ret.push_back( tmpCard.substr(0,80) );
        tmpCard.erase( 0, 80 );
    }
    
    const string cs2 = " / ";
    string comment( vEnd+1, card.cend() );
    pos = comment.find_first_not_of( " /" );
    if( pos != string::npos ) comment = comment.substr(pos);
    cout <<"splitLongString: " << __LINE__ << "  N_tc=" << tmpCard.size() << "  N_c=" << comment.size() << endl;
    while( (tmpCard.size()+comment.size()) > 0 ) {
        size_t N_tc = tmpCard.size();
        size_t N_c = comment.size();
        if( (N_tc + N_c) < 78 ) {
            if( tmpCard == "CONTINUE  ''" ) tmpCard.append(" ");    // just to prettify by aligning the slashes for multi-line comments
            ret.push_back( tmpCard + " / " + comment );
            break;
        }
    cout <<"splitLongString: " << __LINE__ << "  tmpCard=\"" << tmpCard.substr(0,80) << "\"" << endl;
        tmpCard.insert( N_tc-1, "&" );
        size_t remainder = 80-tmpCard.size();
    cout <<"splitLongString: " << __LINE__ << "  N_tc=" << N_tc << "  N_c=" << N_c << "  remainder=" << remainder << endl;
        if( remainder > 5 ) {   // else don't bother appending the comment here
            remainder -= 3;
            tmpCard.append( cs2+comment.substr( 0, remainder ) );
    cout <<"splitLongString: " << __LINE__ << "  remainder=" << remainder << "  tmpcSz=" << tmpCard.size() << endl;
    cout <<"splitLongString: " << __LINE__ << "  tmpCard=\"" << tmpCard.substr(0,80) << "\"" << endl;
            comment.erase( 0, remainder );
        }
        tmpCard.resize( 80, ' ' );
        tmpCard.append( "CONTINUE  ''" );
    cout <<"splitLongString: " << __LINE__ << "  N=" << N_tc << "  N0=" << tmpCard.size() << "  tmpCard=\"" << tmpCard.substr(0,80) << "\"" << endl<<endl;
        ret.push_back( tmpCard.substr(0,80) );
        tmpCard.erase( 0, 80 );
    }

    for( auto& r: ret ) {
        r.resize( 80, ' ' );    // pad any too short cards to 80
    }


    cout <<"Key=\"" << key << "\"  Value=\"" << value << "\"  Comment=\"" << comment << "\"" << endl;

    return ret;
}


string Fits::mergeLongString( vector<string> cards ) {
cout << "mergeLongString: " << __LINE__ << "  nC = " << cards.size() << endl;
    mergeLongStrings( cards );
cout << "mergeLongString: " << __LINE__ << "  nC = " << cards.size() << endl;
    string ret;
    if( cards.size() ) ret = cards[0];
    return ret;
}


void Fits::splitLongStrings( vector<string>& cards ) {
    
}


void Fits::mergeLongStrings( vector<string>& cards ) {
    
}


template <typename T>
T Fits::getValue( const vector<string>& hdr, string key ) {
    
    key.resize(8,' ');       // pad with spaces, or truncate, to 8 characters

    for( string k: hdr ) {
        if( boost::iequals( k.substr(0,key.length()), key) ) {
            size_t commentStart = k.find('/');
            k = trimStringValue( k.substr( 10, commentStart-10 ) );
            try {
                return boost::lexical_cast<T>(k);
            } catch( const boost::bad_lexical_cast& ) {
                // catch and ignore, return default constructed T
            }
        }
    }

    return T();
    
}
namespace redux {
    namespace file {
        template <>
        string Fits::getValue( const vector<string>& hdr, string key ) {
            
            key.resize(8,' ');       // pad with spaces, or truncate, to 8 characters
            string ret = "";
            auto it = hdr.begin();
            for( ; it<hdr.end(); ++it ) {
                if( boost::iequals( it->substr(0,8), key) ) {
                    size_t commentStart = it->find('/');
                    ret = trimStringValue( it->substr( 10, commentStart-10 ) );
                    while( *(ret.rbegin()) == '&' ) {
                        ret.resize( ret.size()-1 );
                        if( ++it<hdr.end() ) {
                            if( boost::iequals( it->substr(0,8), "CONTINUE" )) {
                                commentStart = it->find('/');
                                ret += trimStringValue( it->substr( 10, commentStart-10 ) );
                            }
                        }
                    }
                }
            }
            return ret;
        }
        
        template <>
        bpx::ptime Fits::getValue( const vector<string>& hdr, string key ) {
            using namespace bpx;
            key.resize(8,' ');
            for( string k: hdr ) {
                if( boost::iequals( k.substr(0,8), key) ) {
                    size_t commentStart = k.find('/');
                    k = trimStringValue( k.substr( 10, commentStart-10 ) );
                    try {
                        boost::replace_first(k,"T"," ");
                        return bpx::time_from_string(k);
                    } catch( const boost::bad_lexical_cast& ) {
                        // catch and ignore, return default constructed T
                    }
                    //return bpx::from_iso_string(k);
                }
            }
            return bpx::ptime();
        }
    }
}


template <typename T>
vector<T> Fits::getTableArray( string key ) {
    
    key.resize( 8, ' ' );       // pad with spaces, or truncate, to 8 characters
    
    /*for( string k: primaryHDU.cards ) {
        if( boost::iequals( k.substr(0,8), key) ) {
            size_t commentStart = k.find('/');
    cout << "filefits: " << __LINE__ << "   s=\"" << k << "\"" << endl; 
            k = k.substr( 10, commentStart-10 );
            k.erase(remove(k.begin(), k.end(), ' '), k.end());
    cout << "filefits: " << __LINE__ << "   s=\"" << k << "\"" << endl; 
            return boost::lexical_cast<T>(k);
        }
    }*/
    for( auto& ext: extHDUs ) {
        shared_ptr<ascii_hdu> ptr = dynamic_pointer_cast<ascii_hdu>(ext);
        if( ptr ) {
            for( const ascii_hdu::table_info_t& ti: ptr->table_info ) {
                if( boost::iequals( ti.columnName.substr(0,8), key) ) {
                }
            }
    //             if( boost::iequals( k.substr(0,8), key) ) {
    //                 size_t commentStart = k.find('/');
    //                 k = k.substr( 10, commentStart );
    //                 return boost::lexical_cast<T>(k);
    //             }
        }
    }

}
namespace redux {
    namespace file {
        template <>
        vector<string> Fits::getTableArray( string key ) {
            
            key.resize(8,' ');       // pad with spaces, or truncate, to 8 characters
            vector<string> ret;
            
            int columnStart(0);
            int columnWidth(0);
            for( auto& ext: extHDUs ) {
                shared_ptr<ascii_hdu> ptr = dynamic_pointer_cast<ascii_hdu>(ext);
                if( ptr ) {
                    for( const ascii_hdu::table_info_t& ti: ptr->table_info ) {
                        if( ti.columnFormat[0] == 'A' && boost::iequals(ti.columnName.substr(0,8), key) ) {
                            columnStart = ti.columnStart;
                            columnWidth = atoi( ti.columnFormat.c_str()+1 );
                        }
                    }
                    if( columnStart && columnWidth ) {
                        for( int i=0; i<ptr->dims[1]; ++i) {
                            ret.push_back( string( ptr->data.ptr(i,columnStart-1), columnWidth ) );
                        }
                        return ret;
                    }
                }
            }
            
            return ret;
            
        }
        
        template <>
        vector<bpx::ptime> Fits::getTableArray( string key ) {
            
            using bpx::ptime;
            vector<string> tmp = getTableArray<string>( key );
            vector<ptime> ret;
            for( string it: tmp ) {
                try {
                    boost::replace_first(it,"T"," ");
                    ret.push_back(bpx::time_from_string(it));
                } catch( const boost::bad_lexical_cast& ) {
                    ret.push_back( ptime() );
                }
            }
            
            return ret;

        }
    }
}


size_t Fits::getNumberOfFrames(void) {
    
    if( primaryHDU.nDims > 2 ) {
        // TODO fix for 4+ dimensions when needed
        return primaryHDU.dims[2];
    } else if( primaryHDU.nDims < 2 ) {
        return 0;
    }

    return 1;
    
}


bpx::ptime Fits::getStartTime(void) {

    using bpx::ptime;
    ptime startT = getValue<ptime>( primaryHDU.cards, "DATE-BEG" );
    if( startT.is_special() ) {
        auto times = getTableArray<ptime>("DATE-BEG");
        for( const ptime& t: times ) {
            if( startT.is_special() || startT > t ) {
                startT = t;
            }
        }
        if( !startT.is_special() ) {
            string card = makeCard("DATE-BEG", startT, "First in table." );
            insertCardAfter( primaryHDU.cards, card, "DATE" );
        }
    }

    return startT;
    
}


bpx::ptime Fits::getEndTime(void) {
    
    using bpx::ptime;
    ptime endT = getValue<ptime>( primaryHDU.cards, "DATE-END" );
    if( endT.is_special() ) {
        auto times = getTableArray<ptime>("DATE-BEG");
        for( const ptime& t: times ) {
            if( endT.is_special() || endT < t ) {
                endT = t;
            }
        }
        if( !endT.is_special() ) {
            endT += getExposureTime();
            //  "Last in table + NAXIS3*CADENCE + XPOSURE."  ???
            string card = makeCard( "DATE-END", endT, "Last in table + XPOSURE." );
            insertCardAfter( primaryHDU.cards, card, "DATE" );
        }
    }
    return endT;
    
}


bpx::ptime Fits::getAverageTime(void) {

    using namespace bpx;
    ptime avgT = getValue<ptime>( primaryHDU.cards, "DATE-AVG" );
    if( avgT.is_special() ) {
        auto times = getTableArray<ptime>("DATE-BEG");
        time_duration sum;
        int count(0);
        for( const ptime& t: times ) {
            if( avgT.is_special() ) {
                avgT = t;
            }
            sum += (t - avgT);
            count++;
        }
        if( count ) sum /= count;
        if( !avgT.is_special() ) {
            avgT += sum + getExposureTime()/2;
            //  "Last in table + NAXIS3*CADENCE + XPOSURE."  ???
            string card = makeCard( "DATE-AVG", avgT, "Average time from table." );
            insertCardAfter( primaryHDU.cards, card, "DATE" );
        }
    }
    return avgT;
    
}


// bpx::time_duration Fits::getCadence(void) {
//     
//     float cadence = getValue<float>("CADENCE");
//     if( cadence == 0.0 ) {     // also look for old non-SolarNet keyword
//         cadence = getValue<float>("INTERVAL");
//     }
//     return bpx::microseconds( cadence*1E6 );
// 
// }


bpx::time_duration Fits::getExposureTime(void) {
    
    float exposureTime = getValue<float>( primaryHDU.cards, "XPOSURE" );
    if( exposureTime == 0.0 ) {     // also look for old non-SolarNet keyword
        exposureTime = getValue<float>( primaryHDU.cards, "EXPTIME" );
    }
    return bpx::microseconds( static_cast<long>(exposureTime*1E6) );

}


vector<bpx::ptime> Fits::getStartTimes(void){
    
    return getTableArray<bpx::ptime>("DATE-BEG");
    
}

vector<size_t> Fits::getFrameNumbers(void) {
    
    string fnstr = getValue<string>( primaryHDU.cards, "FRAMENUM" );
    size_t found = fnstr.find_first_of("-,");   // check if it is a range.
    if( found != string::npos ) {
        return redux::util::stringToUInts<size_t>(fnstr);
    }
    // else treat it as an integer.
    size_t firstFrame = 0;
    try {
        firstFrame = boost::lexical_cast<size_t>(fnstr);
    } catch( const boost::bad_lexical_cast& ) {
        // catch and ignore, use 0 as default.
    }

    size_t nFrames = getNumberOfFrames();
    vector<size_t> ret(nFrames);
    std::iota( ret.begin(), ret.end(), firstFrame );
    return ret;
    
}

// vector<bpx::ptime> Fits::getEndTimes(void){
//     
//     using namespace bpx;
//     vector<ptime> times = getTableArray<ptime>("DATE-BEG");
//     time_duration exposureTime = getExposureTime()/2;
//     for( ptime& t: times ) {
//         t += exposureTime;
//     }
//     
//     return times;
//     
// }
// 
// 
// vector<bpx::ptime> Fits::getAverageTimes(void) {
//     
//     using namespace bpx;
//     vector<ptime> times = getTableArray<ptime>("DATE-BEG");
//     time_duration halfExp = getExposureTime()/2;
//     for( ptime& t: times ) {
//         t += halfExp;
//     }
//     
//     return times;
//     
// }
// 
// 
// vector<bpx::ptime> getExposureTimes(void){
//     return vector<bpx::ptime>();   // TODO Implement when we get per-frame exposure times
// }


size_t Fits::dataSize(void) { 

    return nElements()*elementSize();
    
}


size_t Fits::dimSize(size_t i) {
    
    if( static_cast<int>(i) >= primaryHDU.nDims ) return 0;
    
    return primaryHDU.dims[ primaryHDU.nDims-i-1 ];
    
}


uint8_t Fits::elementSize(void) { 
    
    return abs(primaryHDU.bitpix/8);
    
}


size_t Fits::nElements(void) {

    return primaryHDU.nElements;
    
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
int Fits::getIDLType(void) {
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


double Fits::getMinMaxMean( const char* data, double* Min, double* Max ){
    
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


void Fits::read( shared_ptr<redux::file::Fits>& hdr, char* data ) {

    if( !hdr.get() ) {
        return;
    }
    
    int anynull(0), status(0);
    int ret(0);
    
    if( hdr->primaryHDU.dHDU ) {
        if( fits_movabs_hdu( hdr->fitsPtr_, hdr->primaryHDU.dHDU, nullptr, &status) ) {
            throwStatusError( "Fits::read(hdr,data) moving to cHDU.", status );
        }

        char* dataPtr = data;
        auto fP = hdr->fitsPtr_;
        FITSfile* fptr = fP->Fptr;
        if( fptr->zbitpix == 16 ) {      // only support int16_t at the moment
            size_t imgSize = fptr->maxtilelen;
            size_t blockSize = fptr->rice_blocksize;
            std::mutex fMtx;
            LONGLONG rowIndex(1);
            struct rowInfo {
                rowInfo() : offset(0), sz(0), idx(0) {}
                rowInfo(LONGLONG o, LONGLONG s, LONGLONG i ) : offset(o), sz(s), idx(i) {}
                LONGLONG offset, sz, idx;
                bool operator<( const struct rowInfo& rhs ) const { return offset < rhs.offset; }
            };
            std::function<bool(void)> readTile = [&status,&dataPtr,&rowIndex,&fMtx,imgSize,blockSize,fP,fptr](){
                LONGLONG rSize, rOffset;
                char* dP;
                LONGLONG thisI;
                shared_ptr<uint8_t> tmp;
                {
                    std::lock_guard<mutex> fLock( fMtx );
                    dP = dataPtr;
                    dataPtr += imgSize * 2;
                    thisI = rowIndex++;
                    if( thisI > fptr->numrows ) return false;
                    if( fits_read_descriptll( fP, fptr->cn_compressed, thisI, &rSize, &rOffset, &status ) ) {
                        printStatusError( "Fits::read(hdr,data), getting row("+to_string(thisI)+") info.", status );
                    }
                    if( (rOffset+rSize) > fptr->heapsize ) {
                        cerr << "This file has a corrupt index-table for the compressed data. Tile #" << thisI << " will be skipped!" << endl; 
                        return true;
                    }

                    tmp.reset( new uint8_t[rSize], []( uint8_t*& p ) { delete[] p; } );
                    if( fits_read_col( fP, TBYTE, fptr->cn_compressed, thisI, 1, rSize, NULL, tmp.get(), NULL, &status ) ) {
                        printStatusError( "Fits::read(hdr,data), reading row("+to_string(thisI)+").", status );
                    }
                }
                rice_decomp16( tmp.get(), rSize, reinterpret_cast<int16_t*>(dP), imgSize, blockSize );
                std::lock_guard<mutex> fLock( fMtx );
                return true;
            };
            vector<thread> threads;
            int nT = std::min<int>( 4, std::thread::hardware_concurrency() );
            for( int t(0); t<nT; ++t ) {
                threads.push_back( thread([&](){ while( readTile() ); }));
            }
            for( auto &t: threads ) t.join();
        }
        return;
    }
    
    bool pg_data(false);
    char card[81];
    memset(card,0,81);
    
    if( fits_movabs_hdu( hdr->fitsPtr_, 1, nullptr,  &status) ) {
        throwStatusError( "Fits::read(hdr,data) moving to primary HDU.", status );
    }
    
    if( fits_read_keyword( hdr->fitsPtr_, "SOLARNET", card, nullptr, &status ) ) {
        if( status != KEY_NO_EXIST ) {
            throwStatusError( "Fits::read(hdr,data) SOLARNET", status );
        }
        status = 0;
        if( fits_read_key( hdr->fitsPtr_,  TSTRING, "INSTRUME", card, nullptr, &status ) ) {
            if( status != KEY_NO_EXIST ) {
                throwStatusError( "Fits::read(hdr,data) INSTRUME", status );
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
        case( TFLOAT ): ret = fits_read_img( hdr->fitsPtr_, TFLOAT, 1, hdr->nElements(), 0,
            reinterpret_cast<float*>( data ), &anynull, &status); break;
        case( TDOUBLE ): ret = fits_read_img( hdr->fitsPtr_, TDOUBLE, 1, hdr->nElements(), 0,
            reinterpret_cast<double*>( data ), &anynull, &status); break;
		default: throw logic_error("Fits::read: Data-type not supported: " + to_string(hdr->primaryHDU.dataType) );
    }
    
    if ( ret ) {
        throwStatusError( "Fits::read(hdr,data) data (type=" + to_string(hdr->primaryHDU.dataType)+")", status );
    }

}


void Fits::write( const string& filename, const char* data, const shared_ptr<redux::file::Fits> hdr, bool compress, int slice ) {

    // TODO
    
}


template <typename T>
void Fits::read( const string& filename, redux::util::Array<T>& data, shared_ptr<redux::file::Fits>& hdr ) {

    if( !hdr.get() ) {
        hdr.reset( new Fits() );
    }
    hdr->read( filename );

    int nDims = hdr->primaryHDU.nDims;
    int nArrayDims = data.nDimensions();
    bool forceResize = ( nArrayDims < nDims );
    vector<size_t> dimSizes( hdr->primaryHDU.dims.rbegin(), hdr->primaryHDU.dims.rend() );
    if( (int)dimSizes.size() != nDims ) {
        string msg = "Fits::read() dimension mismatch:  nDims=" +to_string(nDims);
        msg += printArray(dimSizes, "\ndimSizes" );
        msg += " file:" + filename;
        throw logic_error( msg );
    }
    for( int i(0); i < nDims; ++i ) {
        if( !forceResize && ( dimSizes[i] != data.dimSize( nArrayDims - nDims + i ) ) ) {
            forceResize = true;
        }
    }

    if( forceResize ) {
        data.resize( dimSizes );
    }

    size_t dataSize = hdr->primaryHDU.nElements * hdr->primaryHDU.elementSize;
    if( dataSize ) {
        auto tmp = shared_ptr<char>( new char[dataSize], []( char * p ) { delete[] p; } );
        read( hdr, tmp.get() );
        switch( hdr->primaryHDU.dataType ) {
            case( TBYTE ):   data.template copyFrom<uint8_t>( tmp.get() ); break;
            case( TSBYTE ):  data.template copyFrom<int8_t>( tmp.get() ); break;
            case( TSHORT ):  data.template copyFrom<int16_t>( tmp.get() ); break;
            case( TUSHORT ): data.template copyFrom<uint16_t>( tmp.get() ); break;
            case( TINT ):    data.template copyFrom<int32_t>( tmp.get() ); break;
            case( TUINT ):   data.template copyFrom<uint32_t>( tmp.get() ); break;
            case( TLONG ):   data.template copyFrom<int64_t>( tmp.get() ); break;
            case( TULONG ):  data.template copyFrom<uint64_t>( tmp.get() ); break;
            case( TFLOAT ):  data.template copyFrom<float>( tmp.get() ); break;
            case( TDOUBLE ): data.template copyFrom<double>( tmp.get() ); break;
            default: string msg = "Fits::read: unsupported data type: "  + to_string(hdr->primaryHDU.dataType);
                throw logic_error( msg );
        }
    }
}
template void Fits::read( const string& filename, redux::util::Array<int8_t>& data, shared_ptr<redux::file::Fits>& hdr );
template void Fits::read( const string& filename, redux::util::Array<uint8_t>& data, shared_ptr<redux::file::Fits>& hdr );
template void Fits::read( const string& filename, redux::util::Array<int16_t>& data, shared_ptr<redux::file::Fits>& hdr );
template void Fits::read( const string& filename, redux::util::Array<uint16_t>& data, shared_ptr<redux::file::Fits>& hdr );
template void Fits::read( const string& filename, redux::util::Array<int32_t>& data, shared_ptr<redux::file::Fits>& hdr );
template void Fits::read( const string& filename, redux::util::Array<uint32_t>& data, shared_ptr<redux::file::Fits>& hdr );
template void Fits::read( const string& filename, redux::util::Array<int64_t>& data, shared_ptr<redux::file::Fits>& hdr );
template void Fits::read( const string& filename, redux::util::Array<uint64_t>& data, shared_ptr<redux::file::Fits>& hdr );
template void Fits::read( const string& filename, redux::util::Array<float  >& data, shared_ptr<redux::file::Fits>& hdr );
template void Fits::read( const string& filename, redux::util::Array<double >& data, shared_ptr<redux::file::Fits>& hdr );
template void Fits::read( const string& filename, redux::util::Array<complex_t >& data, shared_ptr<redux::file::Fits>& hdr );


template <typename T>
void Fits::read( const string& filename, redux::image::Image<T>& image, bool metaOnly ) {
    shared_ptr<Fits> hdr = static_pointer_cast<Fits>( image.meta );
    if( !hdr ) {
        hdr.reset( new Fits() );
        image.meta = hdr;
    }
    if( metaOnly ) {
        hdr->read( filename );
    } else {
        read( filename, image, hdr );
    }
    hdr->close();
}
template void Fits::read( const string & filename, redux::image::Image<int8_t>& image, bool );
template void Fits::read( const string & filename, redux::image::Image<uint8_t>& image, bool );
template void Fits::read( const string & filename, redux::image::Image<int16_t>& image, bool );
template void Fits::read( const string & filename, redux::image::Image<uint16_t>& image, bool );
template void Fits::read( const string & filename, redux::image::Image<int32_t>& image, bool );
template void Fits::read( const string & filename, redux::image::Image<uint32_t>& image, bool );
template void Fits::read( const string & filename, redux::image::Image<int64_t>& image, bool );
template void Fits::read( const string & filename, redux::image::Image<uint64_t>& image, bool );
template void Fits::read( const string & filename, redux::image::Image<float  >& image, bool );
template void Fits::read( const string & filename, redux::image::Image<double >& image, bool );
template void Fits::read( const string & filename, redux::image::Image<complex_t >& image, bool );


template <typename T>
void Fits::write( const string& filename, const redux::util::Array<T>& data, shared_ptr<redux::file::Fits> hdr, int sliceSize ) {
    
    if( !hdr.get() ) {
        hdr.reset( new Fits() );
    }
    

    if( sliceSize ) {   // compression
        
        // TODO  add compression
        
    } else {
        
        // FIXME: fugly hack since cfitsio doesn't have a clear way to overwrite files.
        bfs::path fn(filename);
        if( bfs::exists( bfs::path(filename) ) ) {
            bfs::remove(fn);
        }
    
        if( !hdr->fitsPtr_ ) {
            fits_create_file( &hdr->fitsPtr_, filename.c_str() , &hdr->status_ );
            if( hdr->status_ ) {
                throwStatusError( "Fits::write() create_file", hdr->status_ );
            }
        }
        
        
        int bitpix = getBitpix<T>();
        long long bzero = getBzero<T>();
        int datatype = getDataType( bitpix, bzero );
        int nDims = data.nDimensions();
        if( nDims > 9 ) {
            cerr << "Fits::write()  Array has > 9 dimensions !!" << endl;
        }
        long naxes[9] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        copy_n( data.dimensions().rbegin(), nDims, naxes );

        fits_create_img( hdr->fitsPtr_, bitpix, nDims, naxes, &hdr->status_ );
        if( hdr->status_ ) {
            throwStatusError( "Fits::write() create_img", hdr->status_ );
        }
        
        if( bzero ) {
            string card;
            if( bzero == -128 ) {    // for signed 8-bit data
                card = Fits::makeCard( "BZERO", bzero, "Data is signed integer" );
            } else {
                card = Fits::makeCard( "BZERO", static_cast<uint64_t>(bzero), "Data is unsigned integer" );
            }
            Fits::addCard( hdr->primaryHDU.cards, card );
        }
        
        Fits::removeCards( hdr->primaryHDU.cards, "END" );    // just in case it is not the last card, or if there are multiple.
        for( auto& c: hdr->primaryHDU.cards ) {
            fits_update_card( hdr->fitsPtr_, c.substr(0,8).c_str(), c.c_str(), &hdr->status_ ); 
            if( hdr->status_ ) {
                throwStatusError( "Fits::write() card: "+c, hdr->status_ );
            }
        }
        
        fits_write_img( hdr->fitsPtr_, datatype, 1, data.nElements(), (void*)data.get(), &hdr->status_);
        if( hdr->status_ ) {
            throwStatusError( "Fits::write() write_img: ", hdr->status_ );
        }
        
        for( auto& eh: hdr->extHDUs ) {
            shared_ptr<ascii_hdu> ahdu = dynamic_pointer_cast<ascii_hdu>(eh);
            if( ahdu ) {
                int tfields = 0;
                const char* extname = ahdu->name.empty()?nullptr:ahdu->name.c_str();
                vector<const char*> ttype, tform, tunit;
                for( auto& ti: ahdu->table_info ) {
                    tfields++;
                    ttype.push_back( ti.columnName.empty()?nullptr:ti.columnName.c_str() );
                    tform.push_back( ti.columnFormat.empty()?nullptr:ti.columnFormat.c_str() );
                    tunit.push_back( ti.columnUnit.empty()?nullptr:ti.columnUnit.c_str() );
                }
                fits_create_tbl( hdr->fitsPtr_, ASCII_TBL, 0, tfields, const_cast<char**>(ttype.data()),
                                 const_cast<char**>(tform.data()), const_cast<char**>(tunit.data()), extname, &hdr->status_ );
                if( hdr->status_ ) {
                    throwStatusError( "Fits::write() create_ascii_tbl: ", hdr->status_ );
                }
                
                if( ahdu->data.nDimensions() == 2 ) {   // TODO: fix for generic case, this is hardcoded for timestamps @ SST
                    size_t nRows = ahdu->data.dimSize(0);
                    vector<char*> ptrs;
                    for( size_t i(0); i<nRows; ++i ) {
                        ptrs.push_back( ahdu->data.ptr(i,0) );
                    }
                    fits_write_col( hdr->fitsPtr_, TSTRING, 1, 1, 1, nRows, ptrs.data(), &hdr->status_ );
                    if( hdr->status_ ) {
                        throwStatusError( "Fits::write() write_ascii_tbl: ", hdr->status_ );
                    }

                }
                
            }
            shared_ptr<binary_hdu> bhdu = dynamic_pointer_cast<binary_hdu>(eh);
            if( bhdu ) {
                // TODO
            }
        }
        
    }
    
    hdr->close();

}
template void Fits::write( const string&, const redux::util::Array<int8_t>&, shared_ptr<redux::file::Fits>, int );
template void Fits::write( const string&, const redux::util::Array<uint8_t>&, shared_ptr<redux::file::Fits>, int );
template void Fits::write( const string&, const redux::util::Array<int16_t>&, shared_ptr<redux::file::Fits>, int );
template void Fits::write( const string&, const redux::util::Array<uint16_t>&, shared_ptr<redux::file::Fits>, int );
template void Fits::write( const string&, const redux::util::Array<int32_t>&, shared_ptr<redux::file::Fits>, int );
template void Fits::write( const string&, const redux::util::Array<uint32_t>&, shared_ptr<redux::file::Fits>, int );
template void Fits::write( const string&, const redux::util::Array<int64_t>&, shared_ptr<redux::file::Fits>, int );
template void Fits::write( const string&, const redux::util::Array<uint64_t>&, shared_ptr<redux::file::Fits>, int );
template void Fits::write( const string&, const redux::util::Array<float  >&, shared_ptr<redux::file::Fits>, int );
template void Fits::write( const string&, const redux::util::Array<double >&, shared_ptr<redux::file::Fits>, int );
template void Fits::write( const string&, const redux::util::Array<complex_t >&, shared_ptr<redux::file::Fits>, int );


template <typename T>
void Fits::write( const string & filename, const redux::image::Image<T>& image, int sliceSize ) {
    write( filename, image, static_pointer_cast<redux::file::Fits>( image.meta ), sliceSize );
}
template void Fits::write( const string&, const redux::image::Image<int8_t>&, int );
template void Fits::write( const string&, const redux::image::Image<uint8_t>&, int );
template void Fits::write( const string&, const redux::image::Image<int16_t>&, int );
template void Fits::write( const string&, const redux::image::Image<uint16_t>&, int );
template void Fits::write( const string&, const redux::image::Image<int32_t>&, int );
template void Fits::write( const string&, const redux::image::Image<uint32_t>&, int );
template void Fits::write( const string&, const redux::image::Image<int64_t>&, int );
template void Fits::write( const string&, const redux::image::Image<uint64_t>&, int );
template void Fits::write( const string&, const redux::image::Image<float  >&, int );
template void Fits::write( const string&, const redux::image::Image<double >&, int );
template void Fits::write( const string&, const redux::image::Image<complex_t >&, int );


template <typename T>
void Fits::write( const string & filename, const T* data, size_t n ) {
    
    // TODO add wrapper for POD

}
template void Fits::write( const string&, const int8_t*, size_t n );
template void Fits::write( const string&, const uint8_t*, size_t n );
template void Fits::write( const string&, const int16_t*, size_t n );
template void Fits::write( const string&, const uint16_t*, size_t n );
template void Fits::write( const string&, const int32_t*, size_t n );
template void Fits::write( const string&, const uint32_t*, size_t n );
template void Fits::write( const string&, const int64_t*, size_t n );
template void Fits::write( const string&, const uint64_t*, size_t n );
template void Fits::write( const string&, const float*, size_t n );
template void Fits::write( const string&, const double*, size_t n );
template void Fits::write( const string&, const complex_t*, size_t n );


#endif  // RDX_WITH_FITS


