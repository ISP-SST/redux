USE_REDUX(redux)
USE_EXTERNAL(fits boost pybind)


if( RDX_WITH_FFTW3 )
    add_definitions(-DRDX_WITH_FFTW3)
endif()

if( NOT Python_FOUND OR NOT pybind11_FOUND )
    message(STATUS "pyredux python module will not be built (requires Python3-dev and pybind11).")
    return()
endif()

file(GLOB PYREDUX_CPP "*.cpp")

pybind11_add_module(pyredux ${PYREDUX_CPP})

# Link against the redux library + dependencies
target_link_libraries( pyredux PRIVATE ${RDX_CURRENT_LIBRARIES} )
#target_link_options(pyredux PRIVATE "-Wl,--as-needed")
