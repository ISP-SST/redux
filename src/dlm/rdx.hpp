#ifndef REDUX_DLM_RDX_HPP
#define REDUX_DLM_RDX_HPP

#include <stdio.h>          // has to be included before idl_export.h, otherwise FILE is undefined.
#include <idl_export.h>

namespace redux {

    size_t dumpStruct( IDL_VPTR data, int current=0, int indent=0 );
    void printStruct( int argc, IDL_VPTR argv[], char* argk );
    IDL_VPTR structInfo( int argc, IDL_VPTR argv[], char* argk );

}

#endif  // REDUX_DLM_RDX_HPP
