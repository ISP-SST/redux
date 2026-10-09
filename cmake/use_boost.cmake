if( POLICY CMP0167 )
    # CMake >= 3.30 has removed the bundled FindBoost module; without explicitly
    # selecting NEW here, find_package() tries the (now-removed) module first and
    # fails outright instead of falling back to Boost's own exported config package.
    cmake_policy( SET CMP0167 NEW )
endif()

find_package( Boost 1.66 QUIET COMPONENTS date_time filesystem program_options serialization system thread regex )

if( Boost_FOUND )
    set( Boost_VERSION "${Boost_VERSION_MAJOR}.${Boost_VERSION_MINOR}.${Boost_VERSION_PATCH}" CACHE STRING "Version" FORCE )
    append_libs_unique( RDX_CURRENT_LIBRARIES "${Boost_LIBRARIES}" )
    append_paths_unique( RDX_CURRENT_INCLUDES "${Boost_INCLUDE_DIRS}" )
    append_paths_unique( RDX_CURRENT_LIBDIRS "${Boost_LIBRARY_DIRS}" )
else()
    message(STATUS "Boost >= 1.66 (with components: date_time filesystem program_options serialization system thread regex) not found."
                    " Try your system's equivalent of \"apt-get install libboost-all-dev\","
                    " or point cmake at a non-standard install via -DBOOST_ROOT=/path/to/boost or the BOOST_ROOT environment variable.")
endif()
