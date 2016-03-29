#include "idlutil.hpp"

#include <redux/revision.hpp>
#include <redux/version.hpp>
#include <redux/util/arrayutil.hpp>
#include <redux/util/arraystats.hpp>
#include <redux/util/datautil.hpp>
#include <redux/util/stringutil.hpp>

#include <algorithm>
#include <iostream>


using namespace redux::util;
using namespace redux;

using namespace std;

int IdlContainer::registerRoutine( IDL_SYSFUN_DEF2 def, int is_function, info_func info, void_func cleanup ) {

    auto result = get().routines.emplace( def.name, RoutineInfo(def, info) );
    if( result.second ) {
        result.first->second.is_function = is_function;
        result.first->second.cleanup = cleanup;
    }
    
    return 0;
    
}


void IdlContainer::cleanup( void ) {
    
    for( auto& r: routines ) {
        if( r.second.cleanup ) r.second.cleanup();
    }

}


void IdlContainer::exit( void ) {
    
    get().cleanup();
    get().routines.clear();
    
}


int IdlContainer::load( void ) {
    
    int ret(1);

    for( auto& r: get().routines ) {
        ret &= IDL_SysRtnAdd( &r.second.def, r.second.is_function, 1);
    }
    
    return ret;
    
}


void IdlContainer::reset( void ) {
    
   get().cleanup();
   
}


void IdlContainer::info( int argc, IDL_VPTR* argv, char* argk ) {
    
    int verb(0);
    string out;
    
    if( argc ) {
        
        IDL_VPTR arg = argv[0];
        IDL_ENSURE_SIMPLE( arg );
        IDL_ENSURE_SCALAR( arg );
        
        if ( arg->type == IDL_TYP_STRING ) {
            string fragment = IDL_VarGetString( arg );
            for( auto& r: get().routines ) {
                if( r.second.info && contains(r.first, fragment, true) ) {
                    out += r.second.info(2);
                }
            }
        } else {
            verb = IDL_LongScalar( arg );
            for( auto& r: get().routines ) {
                if( r.second.info ) {
                    out += r.second.info(verb);
                }
            }
        }
    } else {
        for( auto& r: get().routines ) {
            if( r.second.info ) {
                out += r.second.info(0);
            }
        }
    }
    
    cout << "\nRedux DLM";
    
    if( verb > 2 ) {
        cout << "\n   Version:  " << getLongVersionString() << endl;
        cout << "   Commited: " << reduxCommitTime << endl;
        cout << "   Compiled: " << reduxBuildTime << endl;
    } else {
        cout << " - " << getVersionString() << endl;
    }
    
    if( verb <= 2 ) cout << endl << out << endl;
    else cout << endl;
    
}


namespace redux {
    template <typename T> UCHAR idlType(void) { return 0; }
    template<> UCHAR idlType<int16_t>(void) { return IDL_TYP_INT; }
    template<> UCHAR idlType<int32_t>(void) { return IDL_TYP_LONG; }
    template<> UCHAR idlType<float>(void) { return IDL_TYP_FLOAT; }
    template<> UCHAR idlType<double>(void) { return IDL_TYP_DOUBLE; }
}


template <typename T>
void redux::addToRaw( const UCHAR* in, T* out, int64_t nElements, UCHAR type ) {

    if( !nElements ) return;
    
    switch( type ) {
        case( IDL_TYP_BYTE ): {
            for( int64_t i=0; i<nElements; ++i) out[i] += in[i];
            break;
        }
        case( IDL_TYP_INT ): {
            auto inPtr = reinterpret_cast<const IDL_INT*>( in );
            for( int64_t i=0; i<nElements; ++i) out[i] += inPtr[i];
            break;
        }
        case( IDL_TYP_LONG ): {
            auto inPtr = reinterpret_cast<const IDL_LONG*>( in );
            for( int64_t i=0; i<nElements; ++i) out[i] += inPtr[i];
            break;
        }
        case( IDL_TYP_FLOAT ): {
            auto inPtr = reinterpret_cast<const float*>( in );
            for( int64_t i=0; i<nElements; ++i) out[i] += inPtr[i];
            break;
        }
        case( IDL_TYP_DOUBLE ): {
            auto inPtr = reinterpret_cast<const double*>( in );
            for( int64_t i=0; i<nElements; ++i) out[i] += inPtr[i];
            break;
        }
        default: ;
    }

}
template void redux::addToRaw(const UCHAR*, uint8_t*, int64_t, UCHAR);
template void redux::addToRaw(const UCHAR*, int16_t*, int64_t, UCHAR);
template void redux::addToRaw(const UCHAR*, int32_t*, int64_t, UCHAR);
template void redux::addToRaw(const UCHAR*, float*, int64_t, UCHAR);
template void redux::addToRaw(const UCHAR*, double*, int64_t, UCHAR);


template <typename T>
void redux::copyFromRaw( const T* in, UCHAR* out, int64_t nElements, UCHAR type ) {

    if( !nElements || (type == IDL_TYP_UNDEF) ) return;

    if( idlType<T>() == type ) {     // matching type
        memcpy(out, in, nElements*sizeof(T) );
    } else {
        switch( type ) {
            case( IDL_TYP_BYTE ): {
                std::copy(in, in+nElements, out);
                break;
            }
            case( IDL_TYP_INT ): {
                std::copy(in, in+nElements, reinterpret_cast<IDL_INT*>( out ));
                break;
            }
            case( IDL_TYP_LONG ): {
                std::copy(in, in+nElements, reinterpret_cast<IDL_LONG*>( out ));
                break;
            }
            case( IDL_TYP_FLOAT ): {
                std::copy(in, in+nElements, reinterpret_cast<float*>( out ));
                break;
            }
            case( IDL_TYP_DOUBLE ): {
                std::copy(in, in+nElements, reinterpret_cast<double*>( out ));
                break;
            }
            default: memset( out, 0, nElements*IDL_TypeSizeFunc(type) );
        }
    }

}
template void redux::copyFromRaw(const uint8_t*, UCHAR*, int64_t, UCHAR);
template void redux::copyFromRaw(const int16_t*, UCHAR*, int64_t, UCHAR);
template void redux::copyFromRaw(const int32_t*, UCHAR*, int64_t, UCHAR);
template void redux::copyFromRaw(const float*, UCHAR*, int64_t, UCHAR);
template void redux::copyFromRaw(const double*, UCHAR*, int64_t, UCHAR);


template <typename T>
void redux::copyToRaw( const UCHAR* in, T* out, int64_t nElements, UCHAR type ) {

    if( !nElements ) return;
    
    if( idlType<T>() == type ) {     // matching type
        memcpy( out, in, nElements*sizeof(T) );
    } else {
        switch( type ) {
            case( IDL_TYP_BYTE ): {
                std::copy(in, in+nElements, out);
                break;
            }
            case( IDL_TYP_INT ): {
                auto inPtr = reinterpret_cast<const IDL_INT*>( in );
                std::copy(inPtr, inPtr+nElements, out);
                break;
            }
            case( IDL_TYP_LONG ): {
                auto inPtr = reinterpret_cast<const IDL_LONG*>( in );
                std::copy(inPtr, inPtr+nElements, out);
                break;
            }
            case( IDL_TYP_FLOAT ): {
                auto inPtr = reinterpret_cast<const float*>( in );
                std::copy(inPtr, inPtr+nElements, out);
                break;
            }
            case( IDL_TYP_DOUBLE ): {
                auto inPtr = reinterpret_cast<const double*>( in );
                std::copy(inPtr, inPtr+nElements, out);
                break;
            }
            default: memset(out,0,nElements*sizeof(T));
        }
    }


}
template void redux::copyToRaw(const UCHAR*, uint8_t*, int64_t, UCHAR);
template void redux::copyToRaw(const UCHAR*, int16_t*, int64_t, UCHAR);
template void redux::copyToRaw(const UCHAR*, int32_t*, int64_t, UCHAR);
template void redux::copyToRaw(const UCHAR*, float*, int64_t, UCHAR);
template void redux::copyToRaw(const UCHAR*, double*, int64_t, UCHAR);


template <typename T>
shared_ptr<T> redux::castOrCopy( IDL_VPTR array ) {
    
    shared_ptr<T> ret(nullptr);

    if ( (array->flags & IDL_V_ARR) ) {
        int64_t nElements = array->value.arr->n_elts;
        if( nElements ) {
            if( idlType<T>() == array->type ) {     // matching type, no copying
                ret.reset( reinterpret_cast<T*>( array->value.arr->data), [](T* p){} );
            } else {
                ret.reset( new T[ nElements ], [](T* p){ delete[] p; } );
                switch( array->type ) {
                    case( IDL_TYP_BYTE ): {
                        auto data = reinterpret_cast<const UCHAR*>( array->value.arr->data );
                        std::copy(data, data+nElements, ret.get());
                        break;
                    }
                    case( IDL_TYP_INT ): {
                        auto data = reinterpret_cast<const IDL_INT*>( array->value.arr->data );
                        std::copy(data, data+nElements, ret.get());
                        break;
                    }
                    case( IDL_TYP_LONG ): {
                        auto data = reinterpret_cast<const IDL_LONG*>( array->value.arr->data );
                        std::copy(data, data+nElements, ret.get());
                        break;
                    }
                    case( IDL_TYP_FLOAT ): {
                        auto data = reinterpret_cast<const float*>( array->value.arr->data );
                        std::copy(data, data+nElements, ret.get());
                        break;
                    }
                    case( IDL_TYP_DOUBLE ): {
                        auto data = reinterpret_cast<const double*>( array->value.arr->data );
                        std::copy(data, data+nElements, ret.get());
                        break;
                    }
                    default: memset(ret.get(),0,nElements*sizeof(T));
                }
            }
        }
    }
        
    return std::move(ret);

}
template shared_ptr<uint8_t> redux::castOrCopy<uint8_t>(IDL_VPTR array);
template shared_ptr<int16_t> redux::castOrCopy<int16_t>(IDL_VPTR array);
template shared_ptr<int32_t> redux::castOrCopy<int32_t>(IDL_VPTR array);
template shared_ptr<float> redux::castOrCopy<float>(IDL_VPTR array);
template shared_ptr<double> redux::castOrCopy<double>(IDL_VPTR array);


template <typename T>
void redux::copyInto( IDL_VPTR array, T* out, size_t outY, size_t outX, size_t posY, size_t posX ) {

    if ( (array->flags & IDL_V_ARR) && (array->value.arr->n_dim == 2) ) {
        size_t sizeX = array->value.arr->dim[0];
        size_t sizeY = array->value.arr->dim[1];
        if( sizeX && sizeY ) {
            if( idlType<T>() == array->type ) {
                redux::util::copyInto( reinterpret_cast<T*>( array->value.arr->data), sizeY, sizeX, out, outY, outX, posY, posX );
            } else {
                switch( array->type ) {
                    case( IDL_TYP_BYTE ): {
                        redux::util::copyInto( reinterpret_cast<const UCHAR*>( array->value.arr->data), sizeY, sizeX, out, outY, outX, posY, posX );
                        break;
                    }
                    case( IDL_TYP_INT ): {
                        redux::util::copyInto( reinterpret_cast<const IDL_INT*>( array->value.arr->data), sizeY, sizeX, out, outY, outX, posY, posX );
                        break;
                    }
                    case( IDL_TYP_LONG ): {
                        redux::util::copyInto( reinterpret_cast<const IDL_LONG*>( array->value.arr->data), sizeY, sizeX, out, outY, outX, posY, posX );
                        break;
                    }
                    case( IDL_TYP_FLOAT ): {
                        redux::util::copyInto( reinterpret_cast<const float*>( array->value.arr->data), sizeY, sizeX, out, outY, outX, posY, posX );
                        break;
                    }
                    case( IDL_TYP_DOUBLE ): {
                        redux::util::copyInto( reinterpret_cast<const double*>( array->value.arr->data), sizeY, sizeX, out, outY, outX, posY, posX );
                        break;
                    }
                    default: memset( out, 0, outX*outY*sizeof(T) );
                }
            }
        }
    }
}
//template void redux::copyInto( IDL_VPTR, uint8_t*, size_t, size_t, size_t, size_t );
//template void redux::copyInto( IDL_VPTR, int16_t*, size_t, size_t, size_t, size_t );
//template void redux::copyInto( IDL_VPTR, int32_t*, size_t, size_t, size_t, size_t );
//template void redux::copyInto( IDL_VPTR, float*, size_t, size_t, size_t, size_t );
template void redux::copyInto( IDL_VPTR, double*, size_t, size_t, size_t, size_t );


template <typename T>
void redux::copyFromIDL(const UCHAR* in, T* out, size_t nElements, UCHAR IDLtype) {

    switch( IDLtype ) {
        case( IDL_TYP_BYTE ): {
            std::copy(in, in+nElements, out);
            break;
        }
        case( IDL_TYP_INT ): {
            const IDL_INT* cin = reinterpret_cast<const IDL_INT*>( in );
            std::copy(cin, cin+nElements, out);
            break;
        }
        case( IDL_TYP_LONG ): {
            const IDL_LONG* cin = reinterpret_cast<const IDL_LONG*>( in );
            std::copy(cin, cin+nElements, out);
            break;
        }
        case( IDL_TYP_FLOAT ): {
            const float* cin = reinterpret_cast<const float*>( in );
            std::copy(cin, cin+nElements, out);
            break;
        }
        case( IDL_TYP_DOUBLE ): {
            const double* cin = reinterpret_cast<const double*>( in );
            std::copy(cin, cin+nElements, out);
            break;
        }
        default: cout << "copyFromIDL():  type = " << (int)IDLtype << " has not been implemented." << endl;
    }

}
template void redux::copyFromIDL<UCHAR>(const UCHAR*, UCHAR*, size_t, UCHAR);
template void redux::copyFromIDL<int16_t>(const UCHAR*, int16_t*, size_t, UCHAR);
template void redux::copyFromIDL<int32_t>(const UCHAR*, int32_t*, size_t, UCHAR);
template void redux::copyFromIDL<float>(const UCHAR*, float*, size_t, UCHAR);
template void redux::copyFromIDL<double>(const UCHAR*, double*, size_t, UCHAR);


template <typename T>
void redux::copyToIDL(const T* in, UCHAR* out, size_t nElements, UCHAR IDLtype) {

    switch( IDLtype ) {
        case( IDL_TYP_BYTE ): {
            std::copy(in, in+nElements, out);
            break;
        }
        case( IDL_TYP_INT ): {
            std::copy(in, in+nElements, reinterpret_cast<IDL_INT*>( out ));
            break;
        }
        case( IDL_TYP_LONG ): {
            std::copy(in, in+nElements, reinterpret_cast<IDL_LONG*>( out ));
            break;
        }
        case( IDL_TYP_FLOAT ): {
            std::copy(in, in+nElements, reinterpret_cast<float*>( out ));
            break;
        }
        case( IDL_TYP_DOUBLE ): {
            std::copy(in, in+nElements, reinterpret_cast<double*>( out ));
            break;
        }
        default:  cout << "copyToIDL():  type = " << (int)IDLtype << " has not been implemented." << endl;
    }

}
template void redux::copyToIDL<UCHAR>(const UCHAR*, UCHAR*, size_t, UCHAR);
template void redux::copyToIDL<int16_t>(const int16_t*, UCHAR*, size_t, UCHAR);
template void redux::copyToIDL<int32_t>(const int32_t*, UCHAR*, size_t, UCHAR);
template void redux::copyToIDL<float>(const float*, UCHAR*, size_t, UCHAR);
template void redux::copyToIDL<double>(const double*, UCHAR*, size_t, UCHAR);


template <typename T, typename U>
void redux::subtractFromIDL(const UCHAR* in1, const T* in2, U* out, size_t nElements, UCHAR IDLtype) {

    switch( IDLtype ) {
        case( IDL_TYP_BYTE ): {
            std::transform(in1, in1+nElements, in2, out, []( const UCHAR& a, const T&b ){ return a-b;} );
            break;
        }
        case( IDL_TYP_INT ): {
            std::transform(in1, in1+nElements, in2, reinterpret_cast<IDL_INT*>( out ),
                           []( const UCHAR& a, const T&b ){ return a-b;} );
            break;
        }
        case( IDL_TYP_LONG ): {
            std::transform(in1, in1+nElements, in2, reinterpret_cast<IDL_LONG*>( out ),
                           []( const UCHAR& a, const T&b ){ return a-b;} );
            break;
        }
        case( IDL_TYP_FLOAT ): {
            std::transform(in1, in1+nElements, in2, reinterpret_cast<float*>( out ),
                           []( const UCHAR& a, const T&b ){ return a-b;} );
            break;
        }
        case( IDL_TYP_DOUBLE ): {
            std::transform(in1, in1+nElements, in2, reinterpret_cast<double*>( out ),
                           []( const UCHAR& a, const T&b ){ return a-b;} );
            break;
        }
        default:  cout << "subtractFromIDL():  type = " << (int)IDLtype << " has not been implemented." << endl;
    }

}
template void redux::subtractFromIDL<float,float>(const UCHAR*, const float*, float*, size_t, UCHAR);
template void redux::subtractFromIDL<float,double>(const UCHAR*, const float*, double*, size_t, UCHAR);
template void redux::subtractFromIDL<double,float>(const UCHAR*, const double*, float*, size_t, UCHAR);
template void redux::subtractFromIDL<double,double>(const UCHAR*, const double*, double*, size_t, UCHAR);


double redux::getMinMaxMean( const UCHAR* data, int64_t nElements, UCHAR IDLtype, double* Min, double* Max ) {

    ArrayStats stats;
    switch( IDLtype ) {
        case( IDL_TYP_BYTE ):   stats.getMinMaxMean( data, nElements ); break;
        case( IDL_TYP_INT ):    stats.getMinMaxMean( reinterpret_cast<const IDL_INT*>(data), nElements ); break;
        case( IDL_TYP_LONG ):   stats.getMinMaxMean( reinterpret_cast<const IDL_LONG*>(data), nElements ); break;
        case( IDL_TYP_FLOAT ):  stats.getMinMaxMean( reinterpret_cast<const float*>(data), nElements ); break;
        case( IDL_TYP_DOUBLE ): stats.getMinMaxMean( reinterpret_cast<const double*>(data), nElements ); break;
        default: ;
    }
    if( Min ) *Min = stats.min;
    if( Max ) *Max = stats.max;
    
    return stats.mean;
    
}


// print the layout of an IDL structure and return the total data-size
size_t printStruct (IDL_VPTR data, int current, int indent) {

    size_t sz = 0;

    if (data->type == IDL_TYP_STRUCT) {
        if (current < 0) {
            cout << "           TAG           TYPE                     OFFSET         SIZE           " << endl;
            sz += printStruct (data, indent, indent);
            cout << string (45, ' ') << "Total Size:         " << to_string (sz) << endl;
            return sz;
        }

        IDL_StructDefPtr structDef = data->value.s.sdef;

        int nTags = IDL_StructNumTags (structDef);
        int count;

        for (int t = 0; t < nTags; ++t) {
            char *name = IDL_StructTagNameByIndex (structDef, t, 0, 0);
            IDL_VPTR v;
            IDL_MEMINT offset = IDL_StructTagInfoByIndex (structDef, t, 0, &v);
            string type = IDL_TypeNameFunc(v->type);
            count = 1;

            if (v->flags & IDL_V_ARR) {
                type.append ("(");
                for (int d = 0; d < v->value.arr->n_dim; ++d) {
                    if (d) type.append (",");

                    count *= v->value.arr->dim[d];
                    type.append (to_string ( (int) v->value.arr->dim[d]));
                }
                type.append (")");
            }

            cout.setf (std::ios::left);

            if (v->type == IDL_TYP_STRUCT) {
                cout << std::setw (25) << (string (current, ' ') + name);
                cout << std::setw (25) << type;
                cout << std::setw (15) << to_string ( (size_t) offset) << endl;
                sz += count * printStruct(v, current + indent, indent);
            } else {
                sz += count * IDL_TypeSizeFunc(v->type);
                cout << std::setw (25) << (string (current, ' ') + name);
                cout << std::setw (25) << type;
                cout << std::setw (15) << to_string ( (size_t) offset);
                cout << std::setw (15) << to_string (count * IDL_TypeSizeFunc(v->type)) << endl;   // << endl;
            }

        }

    }

    return sz;
}


typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    IDL_INT indent;
} KW_RESULT;

static IDL_KW_PAR kw_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "HELP",      IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (help) },
    { (char*) "INDENT",    IDL_TYP_INT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (indent) },
    { NULL }
};

string structinfo_help( int lvl ) {
    
    string ret = "RDX_STRUCTINFO";
    if( lvl > 0 ) {
        
        ret += ((lvl > 1)?"\n":"     ");          // newline if lvl>1
        ret += "   Syntax:   rdx_structinfo, data, /KEYWORDS\n";
    
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      INDENT              Indentation for each level in the structure (2).\n";
        }
    
    } else ret += "\n";

    return ret;
    
}


void redux::structinfo( int argc, IDL_VPTR argv[], char* argk ) {

    if( argc < 1 ) {
        cout << structinfo_help(2) << endl;
        return;
    }

    IDL_VPTR data = argv[0];

    KW_RESULT kw;
    kw.help = 0;
    kw.indent = 2;
    (void) IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if (kw.help) {
        cout << structinfo_help(2) << endl;
        return;
    }

    int indent = std::min (std::max ( (int) kw.indent, 0), 30);

    printStruct( data, -1, indent);
    
}


namespace {
    static int dummy UNUSED = 
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)structinfo, (char*)"RDX_STRUCTINFO", 0, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },  0, structinfo_help );
}

