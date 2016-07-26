USE_REDUX(redux)

add_definitions(-DBOOST_ALL_DYN_LINK)

if( REDUX_WITH_FFTW3 )
    add_definitions(-DREDUX_WITH_FFTW3)
endif()


REDUX_ADD_EXECUTABLE( rsub jsub.cpp )
REDUX_ADD_EXECUTABLE( rstat jstat.cpp )
REDUX_ADD_EXECUTABLE( rdel jdel.cpp )

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
REDUX_ADD_EXECUTABLE( redux-gui redux-gui.cpp )

