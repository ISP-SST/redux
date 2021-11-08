#

 
#find_package(OpenMP)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
#endif()

# The built-in find_package call seems to have some issues on older cmake versions,
# so we use a hard compiler-flag for now.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

