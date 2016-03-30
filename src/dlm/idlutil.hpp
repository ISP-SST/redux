#ifndef REDUX_DLM_IDLUTIL_HPP
#define REDUX_DLM_IDLUTIL_HPP

#include <functional>
#include <map>
#include <memory>
#include <stdio.h>          // has to be included before idl_export.h, otherwise FILE is undefined.
#include <string>

#include <idl_export.h>

namespace redux {
    
    class IdlContainer {
        typedef std::function<std::string(int)> info_func;
        typedef std::function<void(void)> void_func;
        struct RoutineInfo {
            RoutineInfo( IDL_SYSFUN_DEF2& d, info_func i ) : def(d), info(i), is_function(0), cleanup(nullptr) {}
            IDL_SYSFUN_DEF2 def;
            info_func info;
            int is_function;
            void_func cleanup;
        };
    public:
        static IdlContainer& get(void) { static IdlContainer ic; return ic; };
        static int registerRoutine( IDL_SYSFUN_DEF2 def, int is_function, info_func info=nullptr, void_func cleanup=nullptr );
        static void exit(void);
        static int load(void);
        static void reset(void);
        static void info( int argc, IDL_VPTR argv[], char* argk );
        
    private:
        IdlContainer() {};
        
        void cleanup(void);
        
        std::map< std::string, RoutineInfo > routines;
    };

    template <typename T> UCHAR idlType(void);
    
    template <typename T> void addToRaw( const UCHAR* in, T* out, int64_t nElements, UCHAR type );
    template <typename T> void addToRaw( const IDL_VPTR in, T* out  ) {
        if ( in->flags & IDL_V_ARR ) {
            addToRaw( in->value.arr->data, out, in->value.arr->n_elts, in->type );
        }
    }

    template <typename T, typename U, typename V, typename W>
    void applyDarkAndGain( const T* img, U* out, V* dd, W* gg, size_t nElements ) {
        if( dd && gg ) {
            for( size_t i=0; i<nElements; ++i ) out[i] = (img[i]-dd[i])*gg[i];
        } else if ( dd ) {
            for( size_t i=0; i<nElements; ++i ) out[i] = (img[i]-dd[i]);
        } else if ( gg ) {
            for( size_t i=0; i<nElements; ++i ) out[i] = img[i]*gg[i];
        }
    }
    template <typename T, typename U, typename V>
    void applyDarkAndGain( const UCHAR* in, T* out, U* dd, V* gg, int64_t nElements, UCHAR type ) {
        if( !nElements ) return;
        switch( type ) {
            case( IDL_TYP_BYTE ): {
                applyDarkAndGain( in, out, dd, gg, nElements );
                break;
            }
            case( IDL_TYP_INT ): {
                auto inPtr = reinterpret_cast<const IDL_INT*>( in );
                applyDarkAndGain( inPtr, out, dd, gg, nElements );
                break;
            }
            case( IDL_TYP_LONG ): {
                auto inPtr = reinterpret_cast<const IDL_LONG*>( in );
                applyDarkAndGain( inPtr, out, dd, gg, nElements );
                break;
            }
            case( IDL_TYP_FLOAT ): {
                auto inPtr = reinterpret_cast<const float*>( in );
                applyDarkAndGain( inPtr, out, dd, gg, nElements );
                break;
            }
            case( IDL_TYP_DOUBLE ): {
                auto inPtr = reinterpret_cast<const double*>( in );
                applyDarkAndGain( inPtr, out, dd, gg, nElements );
                break;
            }
            default: ;
        }
    }
    
    template <typename T, typename U, typename V>
    void applyDarkAndGain( const IDL_VPTR in, T* out, U* dd, V* gg ) {
        if ( in->flags & IDL_V_ARR ) {
            applyDarkAndGain( in->value.arr->data, out, dd, gg, in->value.arr->n_elts, in->type );
        }
    }


    template <typename T> void copyFromRaw( const T* in, UCHAR* out, int64_t nElements, UCHAR type );
    template <typename T> void copyFromRaw( const T* in, IDL_VPTR out ) {
        if ( in->flags & IDL_V_ARR ) {
            copyFromRaw( in, out->value.arr->data, in->value.arr->n_elts, in->type );
        }
    }
    
    template <typename T> void copyToRaw( const UCHAR* in, T* out, int64_t nElements, UCHAR type );
    template <typename T> void copyToRaw( const IDL_VPTR in, T* out  ) {
        if ( in->flags & IDL_V_ARR ) {
            copyToRaw( in->value.arr->data, out, in->value.arr->n_elts, in->type );
        }
    }

    template <typename T> std::shared_ptr<T> castOrCopy( IDL_VPTR array );
    
    template <typename T>
    void copyInto( IDL_VPTR array, T* out, size_t outY, size_t outX, size_t posY, size_t posX );
    
    template <typename T>
    void copyFromIDL(const UCHAR* in, T* out, size_t nElements, UCHAR IDLtype);
    template <typename T>
    void copyToIDL(const T* in, UCHAR* out, size_t nElements, UCHAR IDLtype);
    template <typename T, typename U>
    void subtractFromIDL(const UCHAR* in1, const T* in2, U* out, size_t nElements, UCHAR IDLtype);

    double getMinMaxMean( const UCHAR* data, int64_t nElements, UCHAR IDLtype, double* Min=nullptr, double* Max=nullptr, bool* hasInf=nullptr );
    
    void structinfo( int argc, IDL_VPTR argv[], char* argk );
    
}

#endif  // REDUX_DLM_IDLUTIL_HPP
