USE_REDUX(redux)

add_definitions(-DBOOST_ALL_DYN_LINK)

if( REDUX_WITH_FFTW3 )
    add_definitions(-DREDUX_WITH_FFTW3)
endif()

if( Boost_VERSION VERSION_LESS "1.41" )
    message(STATUS "The redux binaries will not be built (only tested with boost >= 1.41).")
    return()
endif()


REDUX_ADD_EXECUTABLE( rdx_sub rsub.cpp )
REDUX_ADD_EXECUTABLE( rdx_stat rstat.cpp )
REDUX_ADD_EXECUTABLE( rdx_del rdel.cpp )
REDUX_ADD_EXECUTABLE( rdx_ctl rctl.cpp )

REDUX_ADD_EXECUTABLE( reduxd reduxd.cpp )

if( REDUX_SKIP_GUI )
    return()
endif()

USE_REDUX(reduxgui redux)

# Check that we found Qt4
if (NOT QT4_FOUND)
    message(WARNING "QT4 not found: can't build redux-gui!")
    return()
endif()

set(CMAKE_AUTOMOC ON)
REDUX_ADD_EXECUTABLE( rdx_gui redux-gui.cpp )

