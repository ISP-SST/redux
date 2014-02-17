
if (NOT QT4_FOUND)
    set( QT_REQUIRED_VERSION "4.8")
    # FindQt4 only serches in the environment path "QTDIR", so we allow cmake to change it
    # to let us build against different versions easily
    if( IS_DIRECTORY "${QTDIR}" )
        set( ENV{QTDIR} ${QTDIR} )
    endif()
    # Attempt to locate QT
    find_package( Qt4 ${QT_REQUIRED_VERSION} COMPONENTS QTCORE QTMAIN QTGUI QTUITOOLS QTDESIGNER QTNETWORK REQUIRED )
endif()

if( QT4_FOUND )
    include(${QT_USE_FILE})
    add_definitions(-DWITH_QT)
    # If found, add relevant paths for current target.
    append_libs_unique( REDUX_CURRENT_LIBRARIES "${QT_LIBRARIES}" )
    list( APPEND REDUX_CURRENT_INCLUDES ${QT_INCLUDES} )
    list( APPEND REDUX_CURRENT_LIBDIRS ${QT_LIBRARY_DIR} )
endif()
