#ifndef REDUX_PYTHON_COMMON_HPP
#define REDUX_PYTHON_COMMON_HPP

#include <pybind11/pybind11.h>
#include <boost/filesystem.hpp>
#include <stdexcept>
#include <string>

namespace redux {

    namespace python {

        // Initialize/bind the submodules.
        void init_io( pybind11::module_& m );
        void init_momfbd( pybind11::module_& m );
        void init_image( pybind11::module_& m );
        void init_cfg( pybind11::module_& m );
        void init_mosaic( pybind11::module_& m );


        struct FileNotFoundError : std::runtime_error {
            using std::runtime_error::runtime_error;
        };


        inline void require_file( const std::string& path ) {
            if( !boost::filesystem::is_regular_file( path ) ) {
                throw FileNotFoundError( "No such file: \"" + path + "\"" );
            }
        }


    }  // namespace python

}  // namespace redux

#endif  // REDUX_PYTHON_COMMON_HPP
