#include <QtGui/QApplication>

#include <redux/gui/reduxmw.hpp>
#include "parportwidget.hpp"

using namespace redux::gui;
using namespace gui;


int main ( int argc, char *argv[] ) {

    QApplication a ( argc, argv );

    QCoreApplication::setOrganizationName ( "redux" );
    QCoreApplication::setOrganizationDomain ( "redux" );
    QCoreApplication::setApplicationName ( "parport-gui" );
    QCoreApplication::setApplicationVersion ( "alpha" );

    ReduxMW mw;
    mw.setWindowTitle ( "parport-gui" );
    ParPortWidget* w = new ParPortWidget(&mw);
    mw.setCentralWidget ( w );
    mw.show();

    return a.exec();
}
