#ifndef REDUX_DLM_OPENCV_HPP
#define REDUX_DLM_OPENCV_HPP

#include <stdio.h>          // has to be included before idl_export.h, otherwise FILE is undefined.
#include <idl_export.h>

namespace redux {

    IDL_VPTR img_align( int argc, IDL_VPTR* argv, char* argk );
    IDL_VPTR img_project( int argc, IDL_VPTR* argv, char* argk );

}

#endif  // REDUX_DLM_OPENCV_HPP
