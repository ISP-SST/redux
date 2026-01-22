#find_package( pybind11 REQUIRED )

set(PYBIND11_FINDPYTHON ON)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

append_libs_unique( RDX_CURRENT_LIBRARIES "${PYTHON_LIBRARY}" )
