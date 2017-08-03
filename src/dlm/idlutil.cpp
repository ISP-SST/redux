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

    string name = def.name + to_string(is_function);
    auto result = get().routines.emplace( name, RoutineInfo(def, info) );
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


#define xstringize(s) stringize(s)
#define stringize(s) #s
int IdlContainer::load( void ) {
    
    int ret(1);
    
    string dlmVer = xstringize(RDX_IDL_HDR_VER);
    string idlVer = string(IDL_SysvVersion.release.s);
    size_t nChars = min( dlmVer.length(), idlVer.length() );
    if( strncmp( idlVer.c_str(), dlmVer.c_str(), nChars ) ) {
        cout << colorString("Warning",YELLOW) << ": The redux DLM was not compiled using the same IDL-header (idl_export.h) as the current session." << endl;
        cout << "IDL Version: " << IDL_SysvVersion.release.s << endl;
        cout << "Header Version: " << dlmVer << endl;
        cout << "This might lead to undefined behaviour if there is a mismatch in the data/functions used, YOU HAVE BEEN WARNED!" << endl;
    }
    for( auto& r: get().routines ) {
        ret &= IDL_SysRtnAdd( &r.second.def, r.second.is_function, 1);
    }
    
    return ret;
    
}


void IdlContainer::reset( void ) {
    
   get().cleanup();
   
}

typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    IDL_VPTR version;
} KW_INFO;

static IDL_KW_PAR kw_info[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "HELP",      IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*)IDL_KW_OFFSETOF2(KW_INFO,help) },
    { (char*) "VERSION",   IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*)IDL_KW_OFFSETOF2(KW_INFO,version) },
    { NULL }
};

string info_help( int lvl ) {
    
    string ret = "RDX";
    if( lvl > 0 ) {
        
        ret += ((lvl > 1)?"\n":"     ");          // newline if lvl>1
        ret += "   Syntax:   rdx, 1/2, /KEYWORDS\n";
        ret += "             rdx, str, /KEYWORDS\n";
    
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      VERSION             Return the DLM version as a string.\n";
        }
    
    } else ret += "\n";

    return ret;
    
}


void IdlContainer::info( int argc, IDL_VPTR* argv, char* argk ) {
    
    KW_INFO kw;
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, kw_info, (IDL_VPTR*)0, 255, &kw );

    if( kw.help ) {
        cout << info_help(2) << endl;
        return;
    }

    if( kw.version ) {
        string ver = getLongVersionString(false);
        IDL_VPTR tmp = IDL_StrToSTRING( (char*)ver.c_str() );
        IDL_VarCopy( tmp, kw.version );
        return;
    }
    
    int verb(0);
    string out;
    
    if( nPlainArgs ) {
        
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
        cout << endl;
        cout << "   Version:  " << getLongVersionString() << endl;
        cout << "   Commited: " << reduxCommitTime << endl;
        cout << "   Compiled: " << reduxBuildTime;
        cout << " (using idl_exports " << xstringize(RDX_IDL_HDR_VER) << ")" << endl;
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
 
    string printValue( IDL_VPTR v, UCHAR* data, size_t maxLength=30 ) {
        
        if( v->type == IDL_TYP_PTR ) {
            IDL_HVID tmpHVID = *((IDL_HVID*)(data));
            IDL_HEAP_VPTR heapPtr = IDL_HeapVarHashFind(tmpHVID);
            if( heapPtr ) {
                v = &(heapPtr->var);
                if( heapPtr->var.flags & IDL_V_ARR ) {
                    data = reinterpret_cast<UCHAR*>(v->value.arr->data);
                } else data = reinterpret_cast<UCHAR*>(&(v->value));
            }
        }
        
        if( v->flags & IDL_V_ARR ) {
            string ret;
            size_t nElements = std::min<size_t>( 5, v->value.arr->n_elts );
            IDL_ALLTYPES* tmpAll = reinterpret_cast<IDL_ALLTYPES*>(data);
            UCHAR* arrayData = tmpAll->arr->data;
            switch( v->type ) {
                case( IDL_TYP_BYTE ):   ret = printArray( arrayData, nElements, "" ); break;
                case( IDL_TYP_INT ):    ret = printArray( reinterpret_cast<IDL_INT*>(arrayData), nElements, "" ); break;
                case( IDL_TYP_LONG ):   ret = printArray( reinterpret_cast<IDL_LONG*>(arrayData), nElements, "" ); break;
                case( IDL_TYP_FLOAT ):  ret = printArray( reinterpret_cast<float*>(arrayData), nElements, "" ); break;
                case( IDL_TYP_DOUBLE ): ret = printArray( reinterpret_cast<double*>(arrayData), nElements, "" ); break;
                case( IDL_TYP_STRING ): {
                    ret = "\"";
                    IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>(arrayData);
                    for( size_t i=0; i<nElements; ++i ) {
                        if( strptr[i].slen && strptr[i].s ) {
                            if( i ) ret += ", ";
                            ret += strptr[i].s;
                        }
                    }
                    ret += "\"";
                    break;
                }
                default: ret =  "";
            }
            return ret.substr( 0, maxLength );
        } else {
            IDL_ALLTYPES* tmpAll = reinterpret_cast<IDL_ALLTYPES*>(data);
            switch( v->type ) {
                case( IDL_TYP_BYTE ):   return to_string( (int)tmpAll->c );
                case( IDL_TYP_INT ):    return to_string( tmpAll->i );
                case( IDL_TYP_LONG ):   return to_string( tmpAll->l );
                case( IDL_TYP_FLOAT ):  return to_string( tmpAll->f );
                case( IDL_TYP_DOUBLE ): return to_string( tmpAll->d );
                case( IDL_TYP_STRING ): {
                    if( tmpAll->str.slen && tmpAll->str.s ) {
                        return "\"" + string( tmpAll->str.s ) + "\"";
                    } else return "\"\"";
                }
                default: return "";
            }
        }
    }
    
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


double redux::getMinMaxMean( const UCHAR* data, int64_t nElements, UCHAR IDLtype, double* Min, double* Max, bool* hasInf ) {

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
    if( hasInf ) *hasInf = stats.hasInfinity;
    
    return stats.mean;
    
}


// print the layout of an IDL structure and return the total data-size
size_t printStruct (IDL_VPTR data, int current, int indent, IDL_INT max_length) {

    size_t sz = 0;
    if( data->type == IDL_TYP_OBJREF ) {    // if it is an object-reference, fetch the heap-pointer and call again.
        IDL_HEAP_VPTR heapPtr = IDL_HeapVarHashFind(data->value.hvid);
        return printStruct( &(heapPtr->var), current, indent, max_length );
    }
    
    if (data->type == IDL_TYP_STRUCT) {
        if (current < 0) {
            cout << "           TAG           TYPE                     OFFSET         SIZE           VALUE" << endl;
            sz += printStruct (data, indent, indent, max_length );
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
            UCHAR dataType = v->type;
            string typeName = string(" ")+IDL_TypeNameFunc(dataType);
            size_t dataSize = IDL_TypeSizeFunc(dataType);
            UCHAR* dataPtr = (data->value.s.arr->data + offset);
            
            if( v->flags & IDL_V_ARR ) {
                v->value.arr->data = dataPtr;
            }

            
            if( dataType == IDL_TYP_PTR ) {
                IDL_HVID tmpHVID = *((IDL_HVID*)(dataPtr));
                IDL_HEAP_VPTR heapPtr = IDL_HeapVarHashFind(tmpHVID);
                if( !heapPtr ) continue;
                v = &(heapPtr->var);
                dataType = heapPtr->var.type;
                dataPtr = reinterpret_cast<UCHAR*>(&(heapPtr->var.value));
                typeName = string("*")+IDL_TypeNameFunc(dataType);
            }
            
            count = 1;
            if( v->flags & IDL_V_ARR ) {
                typeName.append ("(");
                for (int d = 0; d < v->value.arr->n_dim; ++d) {
                    if (d) typeName.append (",");
                    count *= v->value.arr->dim[d];
                    typeName.append (to_string ( (int) v->value.arr->dim[d]));
                }
                typeName.append (")");
                dataPtr = reinterpret_cast<UCHAR*>(&(v->value));
            }

            cout.setf (std::ios::left);

            cout << std::setw (25) << (string (current, ' ') + name);
            cout << std::setw (25) << typeName;
            cout << std::setw (15) << to_string ( (size_t) offset);
            if( dataType == IDL_TYP_STRUCT ) {
                cout << endl;
                sz += count * printStruct(v, current + indent, indent, max_length);
            } else {
                cout << std::setw (15) << to_string (count * dataSize);
                cout << std::setw (30) << printValue( v, dataPtr, max_length );
                cout << endl;
                sz += count * dataSize;
            }

        }

    }

    return sz;
}


typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    IDL_INT indent;
    IDL_INT max_length;
} KW_RESULT;

static IDL_KW_PAR kw_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "HELP",      IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (help) },
    { (char*) "INDENT",    IDL_TYP_INT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (indent) },
    { (char*) "MAX_LENGTH",IDL_TYP_INT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (max_length) },
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

    KW_RESULT kw;
    kw.help = 0;
    kw.indent = 2;
    kw.max_length = 90;
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if( nPlainArgs < 1 ) {
        cout << structinfo_help(2) << endl;
        return;
    }

    IDL_VPTR data = argv[0];

    if (kw.help) {
        cout << structinfo_help(2) << endl;
        return;
    }

    int indent = std::min (std::max ( (int) kw.indent, 0), 30);

    printStruct( data, -1, indent, kw.max_length );
    
}


size_t getStructDefSize( IDL_StructDefPtr structDef, vector<pair<IDL_MEMINT,size_t>>& strOffsets ) {
    
    size_t sz(1);    // \0
    if( structDef->id ) {
        sz += structDef->id->len;
    }
    sz += sizeof(int);    // nTags
    for( int i = 0; i < structDef->ntags; ++i ) {
        sz += 1;    // \0
        if( structDef->tags[i].id ) {
            sz += structDef->tags[i].id->len;
        }
        sz += sizeof(IDL_MEMINT);   // offset
        sz += 2*sizeof( UCHAR );    // all IDL variables a type id & flags
        IDL_VPTR v = &(structDef->tags[i].var);
        size_t count(1);
        sz += sizeof( UCHAR );      // n_dim
        if( v->flags & IDL_V_ARR ) { // it's an array, also store n_dim + dim[]
            sz += v->value.arr->n_dim * sizeof(IDL_MEMINT);
            count = v->value.arr->n_elts;
        }
        if( v->type == IDL_TYP_STRING ) {
            strOffsets.push_back( make_pair( structDef->tags[i].offset, count ) );
        } else if( v->flags & IDL_V_STRUCT ) {
            IDL_StructDefPtr tagStructDef = v->value.s.sdef;
            vector<pair<IDL_MEMINT,size_t>> myStrOffsets;
            size_t tagSize = getStructDefSize( tagStructDef, myStrOffsets );
            sz += tagSize;
            for( auto& it: myStrOffsets ) it.first += structDef->tags[i].offset;
            strOffsets.insert( strOffsets.end(), myStrOffsets.begin(), myStrOffsets.end() );
        }
    }
    return sz;
}


uint64_t packStructDef( char* ptr, IDL_StructDefPtr structDef, vector<pair<IDL_MEMINT,size_t>>& strOffsets) {
    
    uint64_t count(0);    // \0
    if( structDef->id ) {
        string tmp( structDef->id->name );
        count += pack( ptr+count, tmp );
    } else count += pack( ptr+count, UCHAR(0) );

    count += pack( ptr+count, structDef->ntags );
    for( int i=0; i < structDef->ntags; ++i ) {

        if( structDef->tags[i].id ) {
            string tmp( structDef->tags[i].id->name );
            count += pack( ptr+count, tmp );
        } else count += pack( ptr+count, UCHAR(0) );
        count += pack( ptr+count, structDef->tags[i].offset );
        IDL_VPTR v = &(structDef->tags[i].var);
        count += pack( ptr+count, v->type );
        count += pack( ptr+count, v->flags );

        size_t count2(1);

        if( v->flags & IDL_V_ARR ) { // it's an array, also store n_dim + dim[]
            count += pack( ptr+count, v->value.arr->n_dim );
            memcpy( ptr+count, v->value.arr->dim, v->value.arr->n_dim*sizeof(IDL_MEMINT) );
            count += v->value.arr->n_dim * sizeof(IDL_MEMINT);
            count2 = v->value.arr->n_elts;

        } else count += pack( ptr+count, UCHAR(0) );    // scalar
        
        if( v->type == IDL_TYP_STRING ) {
            strOffsets.push_back( make_pair( structDef->tags[i].offset, count2 ) );
        } else if( v->flags & IDL_V_STRUCT ) {
            IDL_StructDefPtr tagStructDef = v->value.s.sdef;
            vector<pair<IDL_MEMINT,size_t>> myStrOffsets;
            count += packStructDef( ptr+count, tagStructDef, myStrOffsets );
            for( auto& it: myStrOffsets ) it.first += structDef->tags[i].offset;
            strOffsets.insert( strOffsets.end(), myStrOffsets.begin(), myStrOffsets.end() );
        }
    }

    return count;
    
}


uint64_t unpackStructDef( const char* ptr, IDL_StructDefPtr& structDef ) {
    
    string structName;
    uint64_t count(0);
    count += unpack( ptr+count, structName );

    int nTags;
    count += unpack( ptr+count, nTags );
    vector<IDL_STRUCT_TAG_DEF> tags;
    
    IDL_STRUCT_TAG_DEF tmp;
    for( int i=0; i < nTags; ++i ) {
        string tagName;
        count += unpack( ptr+count, tagName );
        tmp.name = const_cast<char*> ( tagName.c_str() );
        IDL_MEMINT tagOffset;
        UCHAR tagType;
        count += unpack( ptr+count, tagOffset );
        count += unpack( ptr+count, tagType );
        count += unpack( ptr+count, tmp.flags );
        
        UCHAR nTagDims;
        IDL_MEMINT* tagDims = new IDL_MEMINT[IDL_MAX_ARRAY_DIM];
        void* tagTypePtr = reinterpret_cast<void*>(tagType);
        count += unpack( ptr+count, nTagDims );
        tmp.dims = nullptr;
        if( nTagDims ) {
            memcpy( tagDims+1, ptr+count, nTagDims*sizeof(IDL_MEMINT) );
            count += nTagDims * sizeof(IDL_MEMINT);
            tagDims[0] = nTagDims;
            tmp.dims = tagDims;
        }
        if( tagType == IDL_TYP_STRUCT ) {
            IDL_StructDefPtr tagStructDef;
            count += unpackStructDef( ptr+count, tagStructDef );
            tagTypePtr = (void*)tagStructDef;
        }
        tmp.type = tagTypePtr;
        tags.push_back(tmp);
    }
    memset( &tmp, 0, sizeof(tmp) );
    tags.push_back(tmp);
    char* namePtr = nullptr;
    if( !structName.empty() ) namePtr = const_cast<char*> ( structName.c_str() );
    structDef = IDL_MakeStruct ( namePtr, tags.data() );
    for( auto& t: tags ) {
        if( t.dims ) {
            delete[] t.dims;
            t.dims = nullptr;
        }
    }
    return count;
}


size_t redux::getVarSize( IDL_VPTR v ) {
    
    size_t sz = 2*sizeof( UCHAR ); // all IDL variables a type id & flags
    UCHAR dataType = v->type;
    
    if( dataType == IDL_TYP_OBJREF ) {    // if it is an object-reference, fetch the heap-pointer and call again.
        IDL_HEAP_VPTR heapPtr = IDL_HeapVarHashFind( v->value.hvid );
        cout << "This is an objref" << endl;
        return getVarSize( &(heapPtr->var) );
    }

    if( dataType == IDL_TYP_STRUCT ) {
        
        if( !v->value.s.arr ) {
            cout << "weird, all structs should have the arr..." << endl;
        }
        
        // data block
        sz += sizeof( IDL_MEMINT );         // arr_len
        sz += v->value.s.arr->arr_len;
        sz += sizeof( UCHAR );              // n_dim
        sz += (v->value.s.arr->n_dim+1) * sizeof(IDL_MEMINT);    // dim + elt_len

        // struct def
        IDL_StructDefPtr structDef = v->value.s.sdef;
        vector< pair<IDL_MEMINT, size_t> > strOffsets;
        size_t structDefSz = getStructDefSize( structDef, strOffsets );
        sz += structDefSz;
        size_t sSize(0);
        size_t sCount = v->value.s.arr->n_elts;
        sz += sizeof(size_t);   // nStrings
        for( size_t c=0; c<sCount; ++c ) {
            for( auto& it: strOffsets ) {
                size_t offset = c*v->value.s.arr->elt_len + it.first;
                IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>( v->value.s.arr->data+offset );
                for( size_t i=0; i<it.second; ++i ) {
                    sSize += strptr[i].slen+1;
                }
                sz += it.second*sizeof(size_t);    // offsets
            }
        }

        sz += sSize;
        return sz;
    }
    
    int count = 1;
    if( v->flags & IDL_V_STRUCT ) {
        // structs should be dealt with above
    } else if( v->flags & IDL_V_ARR ) { // it's an array, also store n_dim + dim[]
        sz += sizeof( UCHAR );      // n_dim
        sz += v->value.arr->n_dim * sizeof(IDL_MEMINT);
        if( dataType == IDL_TYP_STRING ) {
            IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>( v->value.arr->data );
            for( int i=0; i<v->value.arr->n_elts; ++i ) {
                sz += strptr[i].slen+1;
            }
            count = 0;  // don't add arraySz*elSz below
        } else {
            count = v->value.arr->n_elts;
        }
    } else if( dataType == IDL_TYP_STRING ) {
        sz += v->value.str.slen+1;
        return sz;
    }

    sz += count * IDL_TypeSizeFunc(dataType);

    return sz;
    
}


uint64_t redux::packVar( char* ptr, IDL_VPTR v ) {
    
    using redux::util::pack;
    
    uint64_t count = pack( ptr, v->type );
    count += pack( ptr+count, v->flags );
    
    if( v->type == IDL_TYP_STRUCT ) {
        
        // data block
        count += pack( ptr+count, v->value.s.arr->arr_len );
        memcpy( ptr+count, v->value.s.arr->data, v->value.s.arr->arr_len );
        count += v->value.s.arr->arr_len;
        count += pack( ptr+count, v->value.s.arr->n_dim );
        memcpy( ptr+count, v->value.s.arr->dim, v->value.s.arr->n_dim*sizeof(IDL_MEMINT) );
        count += v->value.s.arr->n_dim*sizeof(IDL_MEMINT);
        count += pack( ptr+count, v->value.s.arr->elt_len );
        
        // struct def
        IDL_StructDefPtr structDef = v->value.s.sdef;
        vector< pair<IDL_MEMINT, size_t> > strOffsets;
        count += packStructDef( ptr+count, structDef, strOffsets );
        size_t nStrings(0);
        for( auto& it: strOffsets ) nStrings += it.second;
        nStrings *= v->value.s.arr->n_elts;
        count += pack( ptr+count, nStrings );
        for( int c=0; c<v->value.s.arr->n_elts; ++c ) {
            for( auto& it: strOffsets ) {
                size_t offset = c*v->value.s.arr->elt_len + it.first;
                IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>( v->value.s.arr->data+offset );
                for( size_t i=0; i<it.second; ++i ) {
                    string tmp( strptr[i].s );
                    count += pack( ptr+count, offset );
                    count += pack( ptr+count, tmp );
                }
            }
        }
        
        return count;
    }
    
    if( v->flags & IDL_V_STRUCT ) {
        // structs should be dealt with above
    } else if( v->flags & IDL_V_ARR ) { // it's an array, also store n_dim + dim[]
        count += pack( ptr+count, v->value.arr->n_dim );
        for( UCHAR i=0; i<v->value.arr->n_dim; ++i ) {
            count += pack( ptr+count, v->value.arr->dim[i] );
        }
        if( v->type == IDL_TYP_STRING ) {
            IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>( v->value.arr->data );
            for( int i=0; i<v->value.arr->n_elts; ++i ) {
                count += pack( ptr+count, string( strptr[i].s ) );
            }
        } else {
            memcpy( ptr+count, v->value.arr->data, v->value.arr->arr_len );
            count += v->value.arr->arr_len;
        }
        return count;
    }
    
    // scalar values
    if( v->type == IDL_TYP_STRING ) {
        string tmp( v->value.str.s );
        count += pack( ptr+count, tmp );
    } else {
        memcpy( ptr+count, &(v->value), IDL_TypeSizeFunc( v->type ) );
        count += IDL_TypeSizeFunc( v->type );
    }
    return count;
}


uint64_t redux::unpackVar( const char* ptr, IDL_VPTR& v ) {
    
    using redux::util::unpack;
    
    UCHAR dataType, flags;
    uint64_t count = unpack( ptr, dataType );
    count += unpack( ptr+count, flags );

    if( dataType == IDL_TYP_STRUCT ) {
        
        // data block
        IDL_MEMINT blockSize;
        count += unpack( ptr+count, blockSize );
        unique_ptr<UCHAR[]> tmpData( new UCHAR[blockSize] );
        memcpy( tmpData.get(), ptr+count, blockSize );
        count += blockSize;
        UCHAR nDims;
        count += unpack( ptr+count, nDims );
        IDL_ARRAY_DIM dims;
        memcpy( dims, ptr+count, nDims*sizeof(IDL_MEMINT) );
        count += nDims*sizeof(IDL_MEMINT);
        IDL_MEMINT elt_len;
        count += unpack( ptr+count, elt_len );
        
        // struct def
        IDL_StructDefPtr structDef = v->value.s.sdef;
        count += unpackStructDef( ptr+count, structDef );
        
        size_t nStrings;
        count += unpack( ptr+count, nStrings );
        for( size_t i=0; i<nStrings; ++i) {
            size_t offset;
            string tmp;
            count += unpack( ptr+count, offset );
            count += unpack( ptr+count, tmp );
            IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>( tmpData.get()+offset );
            IDL_StrStore( strptr, const_cast<char*>(tmp.c_str()) );
        }
        IDL_VPTR tmp = IDL_ImportArray( nDims, dims, IDL_TYP_STRUCT, tmpData.release(), 0, structDef );
        *v = {0};
        IDL_VarCopy( tmp, v );

        return count;
    }
    
    if( flags & IDL_V_STRUCT ) {
        // structs should be dealt with above
    } else if( flags & IDL_V_ARR ) { // it's an array, also store n_dim + dim[]
        UCHAR nDim;
        IDL_ARRAY_DIM dims;
        count += unpack( ptr+count, nDim );
        for( UCHAR i=0; i<nDim; ++i ) {
            count += unpack( ptr+count, dims[i] );
        }
        if( !nDim ) return count;
        IDL_VPTR tmpV;
        IDL_MakeTempArray( dataType, nDim, dims, IDL_ARR_INI_NOP, &tmpV );
        if( dataType == IDL_TYP_STRING ) {
             IDL_STRING* resptr = reinterpret_cast<IDL_STRING*>( tmpV->value.arr->data );
             for( int i=0; i<tmpV->value.arr->n_elts; ++i ) {
                string tmp;
                count += unpack( ptr+count, tmp );
                IDL_StrStore(resptr+i,const_cast<char*>(tmp.c_str()));
            }
        } else {
            memcpy( tmpV->value.arr->data, ptr+count, tmpV->value.arr->arr_len );
            count += tmpV->value.arr->arr_len;
        }
        IDL_VarCopy( tmpV, v );
        return count;
    }
    
    // scalar values
    if( dataType == IDL_TYP_STRING ) {
        string tmp;
        count += unpack( ptr+count, tmp );
        IDL_VPTR tmpV = IDL_StrToSTRING( const_cast<char*>(tmp.c_str()) );
        IDL_VarCopy( tmpV, v );
    } else {
        v->type = dataType;
        v->flags = 0;   // should always be 0 here
        memcpy( &(v->value), ptr+count, IDL_TypeSizeFunc(dataType) );
        count += IDL_TypeSizeFunc( dataType );
    }
    return count;
}


namespace {
    static int dummy RDX_UNUSED = 
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)structinfo, (char*)"RDX_STRUCTINFO", 0, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },  0, structinfo_help );
}

