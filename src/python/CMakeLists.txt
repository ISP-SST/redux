USE_REDUX(redux)
USE_EXTERNAL(fits boost pybind)

file(GLOB PYREDUX_CPP "*.cpp")

pybind11_add_module(pyredux ${PYREDUX_CPP})

# Link against the redux library + dependencies
target_link_libraries( pyredux PRIVATE ${RDX_CURRENT_LIBRARIES} )
#target_link_options(pyredux PRIVATE "-Wl,--as-needed")
