#include "ana.hpp"

#include "redux/file/fileana.hpp"

#include <boost/filesystem.hpp>

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

namespace bfs=boost::filesystem;

namespace {

    static const uint8_t ana_sizes[] = { 1, 2, 4, 4, 8, 8, 0, 0, 16 };
    static const string ana_type_names[] = {"uint8_t", "int16_t", "int32_t", "float", "double", "long", "undef", "undef", "complex"};

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
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


}


IDL_VPTR redux::fzhead( int argc, IDL_VPTR* argv ) {

    if( argc < 1 ) {
        cout << "FZHEAD: needs 1 argument." << endl;
        ana_help();
        return IDL_GettmpInt(0);    // return a dummy
    }
    IDL_VPTR name_raw = argv[0];

    IDL_ENSURE_SIMPLE( name_raw );
    IDL_ENSURE_STRING( name_raw );

    string name = IDL_VarGetString( name_raw );
    string hdrText;
    try {
        if( !bfs::is_regular_file(name) ) throw std::runtime_error("Not a readable file: " + std::string(name) );
        Ana header( name );
        vector<string> hdrTexts = header.getText();
        if( !hdrTexts.empty() ) hdrText = hdrTexts.front();
    } catch (exception& e) {
        cout << "Failed to read ANA Header from file \"" << name << "\": " << e.what() << endl;
    }

    size_t textSize = hdrText.length();

    IDL_VPTR result = IDL_Gettmp();
    result->type = IDL_TYP_STRING;
    result->value.str.s = new char[ textSize + 1 ];
    if(textSize) {
        memcpy( ( result->value.str.s ), hdrText.c_str(), textSize );
    }
    result->value.str.s[ textSize ] = 0;
    result->value.str.slen = textSize;
    result->value.str.stype = 0; //IDL_V_DYNAMIC;    // flag as dynamic (container will be deleted when destructed)

    return result;

}


IDL_VPTR redux::f0( int argc, IDL_VPTR* argv, char* argk ) {

    KW_RESULT kw;
    kw.help = 0;
    kw.verbose = 0;
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, kw_pars, ( IDL_VPTR* )0, 255, &kw );

    if( nPlainArgs < 1 ) {
        cout << "F0: needs 1 argument." << endl;
        ana_help();
        return IDL_GettmpInt(0);    // return a dummy
    }
    
    IDL_VPTR name_raw = argv[0];
    
    if( kw.help ) {
        ana_help();
        return IDL_GettmpInt(0);    // return a dummy
    }
    
    IDL_KW_FREE;
    
    IDL_ENSURE_SIMPLE( name_raw );
    IDL_ENSURE_STRING( name_raw );

    char* name = IDL_VarGetString( name_raw );

    int verbosity = std::min( std::max( ( int )kw.verbose, 0 ), 8 );
    if( verbosity > 0 ) {
        cout << "Loading file: \"" << name << "\"" << endl;
    }

    try {
        if( !bfs::is_regular_file(name) ) throw std::runtime_error("Not a readable file: " + std::string(name) );
        
        Ana::Ptr header( new Ana(name) );

        int nDims = header->m_Header.ndim;
        if( nDims < 1 || nDims > 8 ) {        // IDL supports up to 8 dimensional data, ANA up to 16.
            throw logic_error("Dimension error.");
        }
        size_t totalSize = 1;
        IDL_ARRAY_DIM dim;
        // ANA files are stored with the fast index first, so reverse the dimensions.
        string dimString = "";
        for( int i=0; i<nDims; ++i ) {
            dim[i] = header->m_Header.dim[i];
            if(i) dimString += "x";
            dimString += to_string(dim[i]);
            totalSize *= dim[i];
        }
        totalSize *= ana_sizes[header->m_Header.datyp];
        
        char* data = new char[totalSize];

        Ana::read(name,data,header);
        if( verbosity > 1 ) {
            cout << "        type: " << ana_type_names[ header->m_Header.datyp ] << endl;
            cout << "       nDims: " << nDims << "   (" << dimString << ")" << endl;
            if( header->m_Header.subf & 1 ) { // compressed
                cout << "       slice: " << ( int )header->m_CompressedHeader.slice_size << endl;
            }
            cout << "      header: " << header->m_ExtendedHeader << endl;
        }

        // IDL type-ID = ANA type-ID + 1
        return IDL_ImportArray( nDims, dim, header->m_Header.datyp + 1, ( UCHAR* )data, redux::util::castAndDelete<char>, NULL );
    } catch (exception& e) {
        cout << "Failed to read Ana file: " << name << endl << "Reason: " << e.what() << endl;
        return IDL_GettmpInt(-1);    // return a dummy
    }
    


}


void redux::fzread( int argc, IDL_VPTR* argv, char* argk ) {

    KW_RESULT kw;
    kw.verbose = 0;
    kw.help = 0;
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, kw_pars, ( IDL_VPTR* )0, 255, &kw );

    if( nPlainArgs < 2 ) {
        cout << "FZREAD: needs (at least) 2 arguments." << endl;
        ana_help();
        return;
    }

    IDL_VPTR imageVar  = argv[0];
    IDL_VPTR fileName  = argv[1];

    if( kw.help ) {
        ana_help();
        return;
    }

    IDL_ENSURE_SIMPLE( imageVar );

    IDL_ENSURE_SIMPLE( fileName );
    IDL_ENSURE_STRING( fileName );

    int verbosity = std::min( std::max( ( int )kw.verbose, 0 ), 8 );

    char *name = IDL_VarGetString( fileName );

    if( verbosity > 0 ) {
        cout << "Loading file: \"" << name << "\"" << endl;
    }

    try {
        if( !bfs::is_regular_file(name) ) throw std::runtime_error("Not a readable file: " + std::string(name) );
        Ana::Ptr header( new Ana(name) );

        int nDims = header->m_Header.ndim;
        if( nDims < 1 || nDims > 8 ) {        // IDL supports up to 8 dimensional data, ANA up to 16.
            throw logic_error("Dimension error.");
        }
        size_t totalSize = 1;
        IDL_ARRAY_DIM dim;
        string dimString = "";
        for( int i=0; i<nDims; ++i ) {
            dim[i] = header->m_Header.dim[i];
            if(i) dimString += "x";
            dimString += to_string(dim[i]);
            totalSize *= dim[i];
        }
        totalSize *= ana_sizes[header->m_Header.datyp];

        char* data = new char[totalSize];

        Ana::read(name,data,header);

        if( verbosity > 1 ) {
            cout << "        type: " << ana_type_names[ header->m_Header.datyp ] << endl;
            cout << "       nDims: " << nDims << "   (" << dimString << ")" << endl;
            if( header->m_Header.subf & 1 ) { // compressed
                cout << "       slice: " << ( int )header->m_CompressedHeader.slice_size << endl;
            }
            cout << "      header: " << header->m_ExtendedHeader << endl;
        }

        char *var_name = IDL_VarName( imageVar );

        // IDL type-ID = ANA type-ID + 1
        argv[0] = IDL_ImportNamedArray( var_name, nDims, dim, header->m_Header.datyp + 1, ( UCHAR* )data, redux::util::castAndDelete<char>, NULL );

        vector<string> hdrTexts = header->getText();
        string hdrText;
        if( !hdrTexts.empty() ) hdrText = hdrTexts.front();
        size_t textSize = hdrText.length();
        if( textSize && argc>2 ) {
            IDL_VPTR hdrVar    = argv[2];
            IDL_ENSURE_SIMPLE( hdrVar );
            IDL_ALLTYPES h;
            h.str.stype = IDL_V_DYNAMIC;
            h.str.s = new char[ textSize + 1 ];
            memcpy( h.str.s, hdrText.c_str(), textSize );
            h.str.s[ textSize ] = 0;
            h.str.slen = textSize;
            IDL_StoreScalar(hdrVar, IDL_TYP_STRING, &h );
        }
    } catch (exception& e) {
        cout << "Failed to read Ana file: " << name << endl << "Reason: " << e.what() << endl;
        return;
    }

    IDL_KW_FREE;

}

void redux::fzwrite( int argc, IDL_VPTR* argv, char* argk ) {

    KW_RESULT kw;
    kw.help = 0;
    kw.verbose = 0;
    kw.compress = 0;
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, kw_pars, ( IDL_VPTR* )0, 255, &kw );

    if( nPlainArgs < 2 ) {
        cout << "FZWRITE: needs (at least) 2 arguments." << endl;
        ana_help();
        return;
    }
    
    IDL_VPTR imageVar  = argv[0];
    IDL_VPTR fileName  = argv[1];

    if( kw.help ) {
        ana_help();
        return;
    }

    IDL_ENSURE_SIMPLE( imageVar );
    IDL_ENSURE_ARRAY( imageVar );

    IDL_ENSURE_SIMPLE( fileName );
    IDL_ENSURE_STRING( fileName );

    char* headerText = nullptr;
    if(nPlainArgs > 2) {
        IDL_VPTR hdrVar    = argv[2];
        IDL_ENSURE_SIMPLE( hdrVar );
        IDL_ENSURE_STRING( hdrVar );
        headerText = IDL_VarGetString( hdrVar );
    }
    int slice = std::min( std::max( static_cast<int>(kw.compress), 0 ), 31 );
    bool compress = slice > 0;
    int verbosity = std::min( std::max( static_cast<int>(kw.verbose), 0 ), 8 );

    char *name = IDL_VarGetString( fileName );

    char* data( NULL );
    IDL_MEMINT ndim;
    IDL_VarGetData( imageVar, &ndim, ( char** )&data, 0 );

    try {

        Ana::Ptr header( new Ana() );

        int nDims = imageVar->value.arr->n_dim;
        string dimString = "";
        for( int i=0; i<nDims; ++i ) {
            header->m_Header.dim[i] = imageVar->value.arr->dim[i];
            if(i) dimString += "x";
            dimString += to_string(header->m_Header.dim[i]);
        }

        header->m_Header.ndim = imageVar->value.arr->n_dim;
        header->m_Header.datyp = imageVar->type - 1;    // ANA type-ID = IDL type-ID - 1

        if( headerText ) {
            header->m_ExtendedHeader = headerText;
        }
        if( verbosity > 0 ) {
            cout << "Writing file: \"" << name << "\"" << endl;
            if( verbosity > 1 ) {
                cout << "        type: " << ana_type_names[ header->m_Header.datyp ] << endl;
                cout << "       nDims: " << nDims << "   (" << dimString << ")" << endl;
                cout << "    compress: " << compress << endl;
                cout << "       slice: " << slice << endl;
                cout << "      header: " << header->m_ExtendedHeader << endl;
            }
        }

        Ana::write( name, data, header, compress, slice );
    } catch (exception& e) {
        cout << "Failed to write Ana file: " << name << endl << "Reason: " << e.what() << endl;
        return;
    }

    IDL_KW_FREE;

}

void redux::fcwrite( int argc, IDL_VPTR* argv, char* argk ) {

    KW_RESULT kw;
    kw.help = 0;
    kw.verbose = 0;
    kw.compress = 5;
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, kw_pars, ( IDL_VPTR* )0, 255, &kw );

    if( nPlainArgs < 2 ) {
        cout << "FCWRITE: needs (at least) 2 arguments." << endl;
        ana_help();
        return;
    }
    
    IDL_VPTR imageVar  = argv[0];
    IDL_VPTR fileName  = argv[1];

    if( kw.help ) {
        ana_help();
        return;
    }

    IDL_ENSURE_SIMPLE( imageVar );
    IDL_ENSURE_ARRAY( imageVar );

    IDL_ENSURE_SIMPLE( fileName );
    IDL_ENSURE_STRING( fileName );

    char* headerText = nullptr;
    if(nPlainArgs > 2) {
        IDL_VPTR hdrVar    = argv[2];
        IDL_ENSURE_SIMPLE( hdrVar );
        IDL_ENSURE_STRING( hdrVar );
        headerText = IDL_VarGetString( hdrVar );
    }
    int slice = std::min( std::max( static_cast<int>(kw.compress), 0 ), 31 );
    bool compress = slice > 0;
    int verbosity = std::min( std::max( static_cast<int>(kw.verbose), 0 ), 8 );

    char *name = IDL_VarGetString( fileName );

    char* data( NULL );
    IDL_MEMINT ndim;
    IDL_VarGetData( imageVar, &ndim, ( char** )&data, 0 );

    try {
        
        Ana::Ptr header( new Ana() );
        
        int nDims = imageVar->value.arr->n_dim;
        for( int i=0; i<nDims; ++i ) {
            header->m_Header.dim[i] = imageVar->value.arr->dim[i];
        }

        header->m_Header.ndim = imageVar->value.arr->n_dim;
        header->m_Header.datyp = imageVar->type - 1;    // ANA type-ID = IDL type-ID - 1

        if(headerText) {
            header->m_ExtendedHeader = headerText;
        }
        
        if( verbosity > 0 ) {
            cout << "Writing file: \"" << name << "\"" << endl;
            if( verbosity > 1 ) {
                cout << "    compress: " << compress << endl;
                cout << "       slice: " << slice << endl;
                cout << "       nDims: " << nDims << endl;
                cout << "      header: " << header->m_ExtendedHeader << endl;
                cout << "        type: " << ana_type_names[ header->m_Header.datyp ] << endl;
            }
        }

        Ana::write( name, data, header, compress, slice );
    } catch (exception& e) {
        cout << "Failed to write Ana file: " << name << endl << "Reason: " << e.what() << endl;
        return;
    }

    IDL_KW_FREE;

}


extern "C" {

    int IDL_Load( void ) {

        static IDL_SYSFUN_DEF2 function_addr[] = {
            { {( IDL_VPTR( * )() )redux::f0}, ( char* )"F0", 0, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { {( IDL_VPTR( * )() )redux::fzhead}, ( char* )"FZHEAD", 0, 1, 0, 0 },
        };

        static IDL_SYSFUN_DEF2 procedure_addr[] = {
            { {( IDL_SYSRTN_GENERIC )redux::fzread}, ( char* )"FZREAD", 0, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { {( IDL_SYSRTN_GENERIC )redux::fzwrite}, ( char* )"FZWRITE", 0, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { {( IDL_SYSRTN_GENERIC )redux::fcwrite}, ( char* )"FCWRITE", 0, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
        };

        /* Register our routine. The routines must be specified exactly the same as in testmodule.dlm. */
        return IDL_SysRtnAdd( function_addr, TRUE, IDL_CARRAY_ELTS( function_addr ) ) &&
               IDL_SysRtnAdd( procedure_addr, FALSE, IDL_CARRAY_ELTS( procedure_addr ) );

    }

}
