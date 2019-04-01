#ifndef REDUX_DLM_MOMFBD_HPP
#define REDUX_DLM_MOMFBD_HPP

#include <stdio.h>          // has to be included before idl_export.h, otherwise FILE is undefined.
#include <idl_export.h>

namespace redux {

    IDL_VPTR momfbd_read( int argc, IDL_VPTR* argv, char* argk );
    IDL_VPTR momfbd_header( int argc, IDL_VPTR* argv, char* argk );
    void momfbd_write( int argc, IDL_VPTR* argv, char* argk );
    IDL_VPTR momfbd_mozaic( int argc, IDL_VPTR *argv, char *argk );

}

#endif  // REDUX_DLM_MOMFBD_HPP
