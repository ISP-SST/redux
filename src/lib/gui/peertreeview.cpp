#include "redux/gui/peertreeview.hpp"

using namespace redux::gui;
using namespace redux;


PeerTreeView::PeerTreeView( QWidget *parent ) : QTreeView( parent ) {
    createActions();
    setExpandsOnDoubleClick( false );
    setRootIsDecorated( false );
    setAlternatingRowColors( true );
    setSelectionMode( QAbstractItemView::ExtendedSelection );
}


void PeerTreeView::createActions( void ) {

    addAct = new QAction( tr( "&Add Host" ), this );
    //addAct->setShortcut(tr("Ctrl+A"));
    addAct->setStatusTip( tr( "Add Host" ) );
    connect( addAct, SIGNAL( triggered() ), this, SLOT( adding() ) );

    removeAct = new QAction( tr( "&Remove" ), this );
    //removeAct->setShortcut(tr("Ctrl+R"));
    removeAct->setStatusTip( tr( "Remove Host" ) );
    connect( removeAct, SIGNAL( triggered() ), this, SLOT( removing() ) );

    launchAct = new QAction( tr( "&Launch Daemon" ), this );
    //launchAct->setShortcut(tr("Ctrl+L"));
    launchAct->setStatusTip( tr( "Launch Host" ) );
    connect( launchAct, SIGNAL( triggered() ), this, SLOT( launching() ) );

    killAct = new QAction( tr( "&Kill" ), this );
    //killAct->setShortcut(tr("Ctrl+K"));
    killAct->setStatusTip( tr( "Kill Host" ) );
    connect( killAct, SIGNAL( triggered() ), this, SLOT( killing() ) );

}


void PeerTreeView::contextMenuEvent( QContextMenuEvent* event ) {

    QMenu menu( this );
    menu.addAction( addAct );
    menu.addAction( removeAct );
    menu.addAction( launchAct );
    menu.addAction( killAct );
    menu.exec( event->globalPos() );

}
