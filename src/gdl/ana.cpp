#include "ana.hpp"

#include "redux/file/fileana.hpp"
#include "redux/file/fileio.hpp"

#include <iostream>
#include <stdexcept>

/*   Examples
    a=f0('/scratch/Data/6302_pf-scan/camXXV_ff21Sep2008.fp3.2447.2943.lc0.0347263')
    fzread,img,'/scratch/Data/6302_pf-scan/camXXV_ff21Sep2008.fp3.2447.2943.lc0.0347263',hdr
    print,fzhead('/scratch/Data/6302_pf-scan/camXXV_ff21Sep2008.fp3.2447.2943.lc0.0347263')
    fzwrite,img,'/scratch/Data/6302_pf-scan/newfile',hdr
    fzwrite,img,'/scratch/Data/6302_pf-scan/newfile',hdr,compress=10
    fcwrite,img,'/scratch/Data/6302_pf-scan/newfile',hdr
*/

using namespace redux::file;
using namespace std;

namespace {

    static const uint8_t ana_sizes[] = { 1, 2, 4, 4, 8, 8, 0, 0, 16 };
    static const string ana_type_names[] = {"uint8_t", "int16_t", "int32_t", "float", "double", "long", "undef", "undef", "complex"};

/*    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
        IDL_INT compress;
        IDL_INT help;
        IDL_INT verbose;
    } KW_RESULT;

    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { ( char* )"COMPRESS", IDL_TYP_INT, 1, IDL_KW_ZERO, 0, ( char* )IDL_KW_OFFSETOF( compress ) },
        { ( char* )"HELP",     IDL_TYP_INT, 1, IDL_KW_ZERO, 0, ( char* )IDL_KW_OFFSETOF( help ) },
        { ( char* )"VERBOSE",  IDL_TYP_INT, 1, 0,           0, ( char* )IDL_KW_OFFSETOF( verbose ) },
        { NULL }
    };

    void ana_help( void ) {

        cout << "ANA_HELP\n"
                "       Syntax:   img=f0('/path/to/file/data.f0')\n"
                "                 hdr=fzhead('/path/to/file/data.f0')\n"
                "                 fzread,img,'/path/to/file/data.f0',hdr,/KEYWORDS\n"
                "                 fzwrite,img,'/path/to/file/data.f0',hdr,/KEYWORDS\n"
                "                 fcwrite,img,'/path/to/file/data.fz',hdr,/KEYWORDS\n"
                "       Accepted Keywords:\n"
                "              HELP                Display this help.\n"
                "              COMPRESS=[0-31]     Compress data using the specified slice-size, (writing only)\n"
                "                                      default is 0 (no compression) for fzwrite and 5 for fcwrite.\n"
                "              VERBOSE={0,1,2}     Verbosity, default is 0 (only error output)." << endl;

    }
*/

}

namespace {
    const string fzhead_help[]={
        "Usage: hdr = rdx_fzhead(filename [,/KEYWORDS])",
        "   Accepted Keywords:",
        "      FITS                Return header as a list of 80-character \"cards\".",
        "      HELP                Display this info."
    };
    constexpr int fzhead_help_sz = sizeof(fzhead_help) / sizeof(string);
}

BaseGDL* redux::rdx_fzhead_gdl( EnvT* e ) {

    static int helpIx = e->KeywordIx("HELP");
    static int fitsIx = e->KeywordIx("FITS");
    
    if( e->KeywordSet( helpIx ) ) {
        e->Help( fzhead_help, fzhead_help_sz );
    }
    
    SizeT nParam RDX_UNUSED = e->NParam( 1 );
    DStringGDL* files = e->GetParAs<DStringGDL>(0);
    SizeT nEl = files->N_Elements();
    if( nEl > 1 ) {    // TODO: implement support for multiple files
        e->Throw( "redux::rdx_fzhead_gdl: Only a single file supported at the moment." );
    }

    bool fitsSet = e->KeywordSet( fitsIx );
    vector<string> texts;
    for( SizeT i(0); i<nEl; ++i ) {
        FileMeta::Ptr meta;
        if( fitsSet ) {
            texts = redux::file::getMetaTextAsCards( (*files)[i], meta, true );
        } else {
            texts.push_back(redux::file::getMetaText( (*files)[i], meta, true ));
        }
    }
    
    DStringGDL* ret = new DStringGDL( texts.size(), BaseGDL::NOZERO );
    for( SizeT i(0); i<texts.size(); ++i ) {
        (*ret)[i] = texts[i];
    }

    return ret;    

}


namespace {
    const string f0_help[]={
        "Usage: data = rdx_f0(filename [,/KEYWORDS])",
        "   Accepted Keywords:",
        "      FITS                Return header as a list of 80-character \"cards\".",
        "      HEADER              (output) Return header.",
        "      HELP                Display this info."
    };
    constexpr int f0_help_sz = sizeof(f0_help) / sizeof(string);
}

BaseGDL* redux::rdx_f0_gdl( EnvT* e ) {

    static int fitsIx = e->KeywordIx("FITS");
    static int headIx = e->KeywordIx("HEADER");
    static int helpIx = e->KeywordIx("HELP");
    
    bool helpSet = e->KeywordSet( helpIx );
    if( helpSet ) {
        e->Help( f0_help, f0_help_sz );
    }

    SizeT nParam RDX_UNUSED = e->NParam( 1 );
    DStringGDL* files = e->GetParAs<DStringGDL>(0);
    SizeT nEl = files->N_Elements();

    if( nEl > 1 ) {    // TODO: implement support for multiple files
        e->Throw( "redux::rdx_f0_gdl: Only a single file supported at the moment." );
    }
    

    bool fitsSet = e->KeywordSet( fitsIx );
    bool headSet = e->KeywordPresent( headIx );

    BaseGDL* ret = nullptr;
    vector<string> texts;
    for( SizeT i(0); i<nEl; ++i ) {
        FileMeta::Ptr meta;
        if( headSet ) {
            if( fitsSet ) {
                texts = redux::file::getMetaTextAsCards( (*files)[i], meta, true );
            } else {
                texts.push_back(redux::file::getMetaText( (*files)[i], meta, true ));
            }
        } else {
            meta = getMeta( (*files)[i] );
        }

        size_t nDims = meta->nDims();
        if( !nDims ) {
            e->Throw( "redux::rdx_f0_gdl: File contains no data: " + (*files)[i] );
            //Warning( "redux::rdx_f0_gdl: File contains no data: " + (*files)[i] );
            //ret = new DLongGDL( 0 );
            //return ret;
        }
        vector<SizeT> fileDims;
        for( size_t i=0; i<nDims; ++i ) {
            fileDims.push_back( meta->dimSize(i) );
        }
        std::reverse( fileDims.begin(), fileDims.end() );   // IDL has reversed dimensions
        dimension gdlDims( fileDims.data(), nDims );
        int dataTYpe = meta->getIDLType();
        void* data = nullptr;
        switch( dataTYpe ) {
            case( GDL_BYTE ): ret = new DByteGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_INT ): ret = new DIntGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_LONG ): ret = new DLongGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_FLOAT ): ret = new DFloatGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_DOUBLE ): ret = new DDoubleGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_COMPLEX ): ret = new DComplexGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_COMPLEXDBL ): ret = new DComplexDblGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_UINT ): ret = new DUIntGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_ULONG ): ret = new DLongGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_LONG64 ): ret = new DLong64GDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_ULONG64 ): ret = new DULong64GDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            default: e->Throw( "redux::rdx_f0_gdl: Unrecognized dataType: " + (*files)[i] );
        }
        readFile( (*files)[i], reinterpret_cast<char*>(data), meta );
    }
    
    if( headSet ) {
        DStringGDL* hdr = new DStringGDL( texts.size(), BaseGDL::NOZERO );
        for( SizeT i(0); i<texts.size(); ++i ) {
            (*hdr)[i] = texts[i];
        }
        e->SetKW( headIx, hdr );
    }

    return ret;    

}


namespace {
    const string fzread_help[]={
        "Usage: rdx_fzread, data, filename [,header] [,/KEYWORDS]",
        "   Accepted Keywords:",
        "      FITS                Return header as a list of 80-character \"cards\".",
        "      HELP                Display this info."
    };
    constexpr int fzread_help_sz = sizeof(fzread_help) / sizeof(string);
}

void redux::rdx_fzread_gdl( EnvT* e ) {

    static int fitsIx = e->KeywordIx("FITS");
    static int helpIx = e->KeywordIx("HELP");

    if( e->KeywordSet( helpIx ) ) {
        e->Help( f0_help, f0_help_sz );
    }

    SizeT nParam = e->NParam( 2 );
    
    BaseGDL*& ret = e->GetPar(0);
    
    DString filename;
    e->AssureStringScalarPar( 1, filename );

    bool fitsSet = e->KeywordSet( fitsIx );
    bool getHdr(false);
    
    if( nParam > 2 ) {
        getHdr = true;
    }

    vector<string> texts;
    FileMeta::Ptr meta;
    if( getHdr ) {
        if( fitsSet ) {
            texts = redux::file::getMetaTextAsCards( filename, meta, true );
        } else {
            texts.push_back(redux::file::getMetaText( filename, meta, true ));
        }
    } else {
        meta = getMeta( filename );
    }

    size_t nDims = meta->nDims();
    if( !nDims ) {
        e->Throw( "redux::rdx_fzread_gdl: File contains no data: " + filename );
    }
    vector<SizeT> fileDims;
    for( size_t i=0; i<nDims; ++i ) {
        fileDims.push_back( meta->dimSize(i) );
    }
    std::reverse( fileDims.begin(), fileDims.end() );   // IDL has reversed dimensions
    dimension gdlDims( fileDims.data(), nDims );
    int dataTYpe = meta->getIDLType();
    void* data = nullptr;
    switch( dataTYpe ) {
        case( GDL_BYTE ): ret = new DByteGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
        case( GDL_INT ): ret = new DIntGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
        case( GDL_LONG ): ret = new DLongGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
        case( GDL_FLOAT ): ret = new DFloatGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
        case( GDL_DOUBLE ): ret = new DDoubleGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
        case( GDL_COMPLEX ): ret = new DComplexGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
        case( GDL_COMPLEXDBL ): ret = new DComplexDblGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
        case( GDL_UINT ): ret = new DUIntGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
        case( GDL_ULONG ): ret = new DLongGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
        case( GDL_LONG64 ): ret = new DLong64GDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
        case( GDL_ULONG64 ): ret = new DULong64GDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
        default: e->Throw( "redux::rdx_fzread_gdl: Unrecognized dataType: " + filename );
    }
    readFile( filename, reinterpret_cast<char*>(data), meta );
    
    if( getHdr ) {
        DStringGDL* hdr = new DStringGDL( texts.size(), BaseGDL::NOZERO );
        for( SizeT i(0); i<texts.size(); ++i ) {
            (*hdr)[i] = texts[i];
        }
        e->SetPar( 2, hdr );
    }

}


namespace {
    const string fzwrite_help[]={
        "Usage: rdx_fzwrite, data, filename [,header] [,/KEYWORDS]",
        "       rdx_fcwrite, data, filename [,header] [,/KEYWORDS]",
        "   Accepted Keywords:",
        "      COMPRESS=[0-31]     Compress data using the specified slice-size, (default=[0,5] for [fzwrite,fcwrite])",
        "      HELP                Display this info."
    };
    constexpr int fzwrite_help_sz = sizeof(fzwrite_help) / sizeof(string);
}

void redux::rdx_fzwrite_gdl( EnvT* e ) {
    
    static int compIx = e->KeywordIx("COMPRESS");
    static int helpIx = e->KeywordIx("HELP");

    if( e->KeywordSet( helpIx ) ) {
        e->Help( fzwrite_help, fzwrite_help_sz );
    }

    SizeT nParam = e->NParam( 2 );
    
    BaseGDL* data = e->GetPar(0);
    
    DString filename;
    e->AssureStringScalarPar( 1, filename );

    bool hasHdr(false);
    
    DString headerText;
    if( nParam > 2 ) {
        e->AssureStringScalarPar( 2, headerText );
    }

    int verbosity(0);   // TBD: should this be kept?
    DLong slice(0);
    string procname = e->GetProName();
    if( procname == "RDX_FCWRITE" ) {
        slice = 5;
    }
    e->AssureLongScalarKWIfPresent( compIx, slice );
    slice = std::min( std::max( slice, 0 ), 31 );
    bool compress = slice > 0;

    void* dataPtr = data->DataAddr();

    try {
        Ana::Ptr header( new Ana() );
        
        const dimension& dim = data->Dim();
        SizeT nDims = dim.Rank();
        string dimString = "";
        for( SizeT i=0; i<nDims; ++i ) {
            header->m_Header.dim[i] = dim[i];
            if(i) dimString += "x";
            dimString += to_string(header->m_Header.dim[i]);
        }
        header->m_Header.ndim = nDims;
        header->m_Header.datyp = data->Type() - 1;    // ANA type-ID = IDL type-ID - 1

        if( !headerText.empty() ) {
            header->m_ExtendedHeader = headerText;
        }
        if( verbosity > 0 ) {
            cout << "Writing file: \"" << filename << "\"" << endl;
            if( verbosity > 1 ) {
                cout << "        type: " << ana_type_names[ header->m_Header.datyp ] << endl;
                cout << "       nDims: " << nDims << "   (" << dimString << ")" << endl;
                cout << "    compress: " << compress << endl;
                cout << "       slice: " << slice << endl;
                cout << "      header: " << header->m_ExtendedHeader << endl;
            }
        }

        Ana::write( filename, reinterpret_cast<char*>(dataPtr), header, compress, slice );
    } catch( exception& ex ) {
        e->Throw( "redux::rdx_fzwrite: Failed to write Ana file: " + filename + " Reason: " + ex.what() );
    }
    

}


extern "C" {

    BaseGDL* rdx_fzhead( EnvT* e ) {  return redux::rdx_fzhead_gdl(e); }
    BaseGDL* rdx_f0( EnvT* e ) {  return redux::rdx_f0_gdl(e); }

    void rdx_fzread( EnvT* e ) {  return redux::rdx_fzread_gdl(e); }
    void rdx_fzwrite( EnvT* e ) {  return redux::rdx_fzwrite_gdl(e); }
    void rdx_fcwrite( EnvT* e ) {  return redux::rdx_fzwrite_gdl(e); }
    
}
