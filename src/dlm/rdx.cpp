#include "rdx.hpp"

#include <iomanip>
#include <iostream>

using namespace std;


namespace {

    const char *var_names[] = { "undef", "int8", "int16", "int32", "float32",
                                "float64", "complex", "string", "struct", "double complex",
                                "int32", "int32", "uint16", "uint32", "int64", "uint64"
                              };

    const size_t var_sizes[] = { 0, sizeof ( UCHAR ), sizeof ( IDL_INT ), sizeof ( IDL_LONG ), sizeof ( float ),
                                 sizeof ( double ), sizeof ( IDL_COMPLEX ), sizeof ( IDL_STRING ), 0, sizeof ( IDL_DCOMPLEX ),
                                 sizeof ( IDL_ULONG ), sizeof ( IDL_ULONG ), sizeof ( IDL_UINT ), sizeof ( IDL_ULONG ), sizeof ( IDL_LONG64 ), sizeof ( IDL_ULONG64 )
                               };

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT help;
        IDL_INT indent;
    } KW_RESULT;

    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { ( char* ) "HELP",    IDL_TYP_INT, 1, IDL_KW_ZERO, 0, ( char* ) IDL_KW_OFFSETOF ( help ) },
        { ( char* ) "INDENT",  IDL_TYP_INT, 1, 0,           0, ( char* ) IDL_KW_OFFSETOF ( indent ) },
        { NULL }
    };

    void print_help ( void ) {

        cout << "RDX HELP\n"
             "       Syntax:   struct_info,data,/KEYWORDS\n"
             "                 sz=struct_info(data,/KEYWORDS)\n"
             "       Accepted Keywords:\n"
             "              HELP              Display this help.\n"
             "              INDENT=[0-30]     Indentation for each level in the structure, default is 2." << endl;

    }


}


//
// print the layout of an IDL structure and return the total data-size
//
size_t redux::dumpStruct ( IDL_VPTR data, int current, int indent ) {

    size_t sz = 0;
    if ( data->type == IDL_TYP_STRUCT ) {
        if ( current < 0 ) {
            cout << "           TAG           TYPE                     OFFSET         SIZE           " << endl;
            sz += dumpStruct ( data, indent, indent );
                cout << string ( 45, ' ' ) << "Total Size:         " << to_string ( sz ) << endl;
            return sz;
        }


        IDL_StructDefPtr structDef = data->value.s.sdef;
        //uint8_t *buf = data->value.s.arr->data;
        int nTags = IDL_StructNumTags ( structDef );
        int count;

        for ( int t = 0; t < nTags; ++t ) {
            char *name = IDL_StructTagNameByIndex ( structDef, t, 0, 0 );
            IDL_VPTR v;
            IDL_MEMINT offset = IDL_StructTagInfoByIndex ( structDef, t, 0, &v );
            string type = string ( var_names[v->type] );
            count = 1;
            if ( v->flags & IDL_V_ARR ) {
                type.append ( "(" );
                for ( int d = 0; d < v->value.arr->n_dim; ++d ) {
                    if ( d ) type.append ( "," );
                    count *= v->value.arr->dim[d];
                    type.append ( to_string ( ( int ) v->value.arr->dim[d] ) );
                }
                type.append ( ")" );
            }

            cout.setf ( std::ios::left );
            if ( v->type == IDL_TYP_STRUCT ) {
                cout << std::setw(25) << (string(current, ' ') + name);
                cout << std::setw(25) << type;
                //cout.setf ( std::ios::hex );
                //cout << std::setw(20) << (void*)( buf + offset );
                //cout.setf ( std::ios::dec );
                cout << std::setw(15) << to_string ( (size_t)offset ) << endl;
                sz += count * dumpStruct ( v, current + indent, indent );
            } else {
                sz += count * var_sizes[v->type];
                cout << std::setw(25) << (string(current,' ')+name);
                cout << std::setw(25) << type;
                //cout.setf ( std::ios::hex );
                //cout << std::setw(20) << (void*)( buf + offset );
                //cout.setf ( std::ios::dec );
                cout << std::setw(15) << to_string ( (size_t)offset );
                cout << std::setw(15) << to_string ( count * var_sizes[v->type] ) << endl; // << endl;
            }

        }

    }

    return sz;
}


void redux::printStruct ( int argc, IDL_VPTR argv[], char* argk ) {

    if ( argc < 1 ) {
        print_help();
        return;
    }

    IDL_VPTR data = argv[0];

    KW_RESULT kw;
    kw.help = 0;
    kw.indent = 2;
    ( void ) IDL_KWProcessByOffset ( argc, argv, argk, kw_pars, ( IDL_VPTR* ) 0, 255, &kw );

    if ( kw.help ) {
        print_help();
        return;
    }

    int indent = std::min ( std::max ( ( int ) kw.indent, 0 ), 30 );


    dumpStruct ( data, -1, indent );
}


IDL_VPTR redux::structInfo ( int argc, IDL_VPTR argv[], char* argk ) {
    if ( argc < 1 ) {
        print_help();
        return IDL_Gettmp();
    }
    IDL_VPTR data = argv[0];

    KW_RESULT kw;
    kw.help = 0;
    kw.indent = 2;
    ( void ) IDL_KWProcessByOffset ( argc, argv, argk, kw_pars, ( IDL_VPTR* ) 0, 255, &kw );

    if ( kw.help ) {
        print_help();
        return IDL_Gettmp();
    }

    int indent = std::min ( std::max ( ( int ) kw.indent, 0 ), 30 );

    return IDL_GettmpULong ( dumpStruct ( data, -1, indent ) );
}


extern "C" {

    int IDL_Load ( void ) {

        static IDL_SYSFUN_DEF2 function_addr[] = {
            { { ( IDL_VPTR ( * ) () ) redux::structInfo}, ( char* ) "STRUCT_INFO", 0, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
        };

        static IDL_SYSFUN_DEF2 procedure_addr[] = {
            { { ( IDL_SYSRTN_GENERIC ) redux::printStruct}, ( char* ) "STRUCT_INFO", 0, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
        };

        return IDL_SysRtnAdd ( function_addr, TRUE, IDL_CARRAY_ELTS ( function_addr ) ) &&
               IDL_SysRtnAdd ( procedure_addr, FALSE, IDL_CARRAY_ELTS ( procedure_addr ) );

    }

}
