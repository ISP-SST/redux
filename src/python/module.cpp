#include "redux/python/common.hpp"

using namespace redux::python;

PYBIND11_MODULE(pyredux, m) {
    m.doc() = "Python bindings for libredux.";

    pybind11::register_exception<FileNotFoundError>( m, "FileNotFoundError", PyExc_FileNotFoundError );

    init_io( m );
    init_momfbd( m );
    init_image( m );
    init_cfg( m );
    init_mosaic( m );

}
