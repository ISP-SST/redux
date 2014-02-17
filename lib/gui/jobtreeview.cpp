#include "redux/gui/jobtreeview.hpp"

#include "redux/job.hpp"

using namespace redux::gui;
using namespace redux;

JobTreeView::JobTreeView( QWidget *parent ) : QTreeView( parent ) {
    createActions();
    setExpandsOnDoubleClick( true );
    setRootIsDecorated( false );
    setAlternatingRowColors( true );

    setSelectionMode( QAbstractItemView::ExtendedSelection );

}


void JobTreeView::createActions( void ) {

    //addAct = new QAction ( tr ( "&Add Job" ), this );
    //addAct->setShortcut(tr("Ctrl+A"));
    //addAct->setStatusTip ( tr ( "Add Job" ) );
    //connect ( addAct, SIGNAL ( triggered() ), this, SLOT ( adding() ) );
    ajAction* tmp;
    /*
        for (int i=1; i<=REDUX::DEBUGJOB; ++i ) {
           tmp = new ajAction ( this, i );
           tmp->setText( (REDUX::jobNames[i]).c_str() );

           connect ( tmp, SIGNAL ( ajTriggered(int) ), this, SLOT ( addJob( int ) ) );
           addJobActs.append( tmp );
        }*/

    cloneAct = new QAction( tr( "&Clone Job" ), this );
    //cloneAct->setShortcut(tr("Ctrl+C"));
    cloneAct->setStatusTip( tr( "Clone Job" ) );
    connect( cloneAct, SIGNAL( triggered() ), this, SLOT( cloning() ) );

    removeAct = new QAction( tr( "&Remove Job" ), this );
    //removeAct->setShortcut(tr("Ctrl+R"));
    removeAct->setStatusTip( tr( "Remove Job" ) );
    connect( removeAct, SIGNAL( triggered() ), this, SLOT( removing() ) );

    submitAct = new QAction( tr( "&Submit Job" ), this );
    //launchAct->setShortcut(tr("Ctrl+L"));
    submitAct->setStatusTip( tr( "Submit Job" ) );
    connect( submitAct, SIGNAL( triggered() ), this, SLOT( submit() ) );

    moveUpAct = new QAction( tr( "&Move Up" ), this );
    //killAct->setShortcut(tr("Ctrl+K"));
    moveUpAct->setStatusTip( tr( "Move Up in List" ) );
    connect( moveUpAct, SIGNAL( triggered() ), this, SLOT( moveUp() ) );

    moveDownAct = new QAction( tr( "&Move Down" ), this );
    //killAct->setShortcut(tr("Ctrl+K"));
    moveDownAct->setStatusTip( tr( "Move Down in List" ) );
    connect( moveDownAct, SIGNAL( triggered() ), this, SLOT( moveDown() ) );

}


void JobTreeView::contextMenuEvent( QContextMenuEvent* event ) {

    QMenu menu( this );
    QMenu submenu( this );

    for( int i = 0; i < addJobActs.size(); ++i ) {
        submenu.addAction( addJobActs[i] );
    }

    //cout << "BLA " << currentIndex().row() << " - " << currentIndex().column() << " - " << (int)currentIndex().isValid() << endl;
    submenu.setTitle( "Add Job" );
    //menu.addAction ( addAct );
    menu.addMenu( &submenu );
    if( selectedIndexes().size() > 0 ) menu.addAction( cloneAct );
    menu.addAction( removeAct );
    menu.addAction( submitAct );
    menu.addAction( moveUpAct );
    menu.addAction( moveDownAct );
    menu.exec( event->globalPos() );

}
