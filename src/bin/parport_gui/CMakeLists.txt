USE_REDUX(reduxgui redux)

set(CMAKE_AUTOMOC ON)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

redux_add_executable(parport-gui main.cpp parportpin.cpp parportframe.cpp parportwidget.cpp )

