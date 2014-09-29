#ifndef REDUX_DLM_ANA_HPP
#define REDUX_DLM_ANA_HPP

#include <stdio.h>          // has to be included before idl_export.h, otherwise FILE is undefined.
#include <idl_export.h>

namespace redux {

    IDL_VPTR fzhead( int argc, IDL_VPTR* argv );
    IDL_VPTR f0( int argc, IDL_VPTR* argv, char* argk );
    void fzread( int argc, IDL_VPTR* argv, char* argk );
    void fzwrite( int argc, IDL_VPTR* argv, char* argk );
    void fcwrite( int argc, IDL_VPTR* argv, char* argk );

}

#endif  // REDUX_DLM_ANA_HPP
