#ifndef REDUX_DLM_GRIDMATCH_HPP
#define REDUX_DLM_GRIDMATCH_HPP

#include <stdio.h>          // has to be included before idl_export.h, otherwise FILE is undefined.
#include <idl_export.h>

namespace redux {

    IDL_VPTR gridmatch( int argc, IDL_VPTR argv[] );
    IDL_VPTR stretch( int argc, IDL_VPTR argv[] );
}

#endif  // REDUX_DLM_GRIDMATCH_HPP
