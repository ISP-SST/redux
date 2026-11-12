# Optional python bindings for some libredux tools/functions.
set(PYBIND11_FINDPYTHON ON)
find_package(Python COMPONENTS Interpreter Development)
find_package(pybind11)

if( Python_FOUND )
    append_libs_unique( RDX_CURRENT_LIBRARIES "${PYTHON_LIBRARY}" )
endif()
