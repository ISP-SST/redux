
USE_REDUX(redux)

add_definitions(-DBOOST_ALL_DYN_LINK)

REDUX_ADD_EXECUTABLE( jsub jsub.cpp )
REDUX_ADD_EXECUTABLE( jstat jstat.cpp )
REDUX_ADD_EXECUTABLE( jdel jdel.cpp )

REDUX_ADD_EXECUTABLE( reduxd reduxd.cpp )

