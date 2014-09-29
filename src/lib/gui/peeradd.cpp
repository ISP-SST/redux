#include "redux/gui/peeradd.hpp"

//#include "redux/network/protocol.hpp"

using namespace redux::gui;
using namespace redux;

HostAdd::HostAdd( Peer::HostInfo* host, bool mst, QWidget* par ) : QDialog( par ) {

    newHost = host;
    parent = par;

    mainLayout = new QGridLayout( this );

    hostNameLabel = new QLabel( this );
    hostNameLabel->setText( tr( "Name (IP)" ) );
    mainLayout->addWidget( hostNameLabel, 0, 0, 1, 1 );
    hostName = new QLineEdit( this );
    hostName->setText( "localhost" );
    mainLayout->addWidget( hostName, 0, 1, 1, 1 );

    hostPortLabel = new QLabel( this );
    hostPortLabel->setText( tr( "Port" ) );
    mainLayout->addWidget( hostPortLabel, 0, 2, 1, 1 );
    hostPort = new QSpinBox( this );
    hostPort->setMaximum( 65536 );
    hostPort->setValue( 30000 );
    hostPort->setToolTip( QString( "The port entered will be tried first.\n" ) +
                          QString( "If binding to that port fails, " ) +
                          QString( "it will\nattempt the following 20 " ) +
                          QString( "ports before\nfailure is reported." ) ); ;
    mainLayout->addWidget( hostPort, 0, 3, 1, 1 );

    hostThreadsLabel = new QLabel( this );
    hostThreadsLabel->setText( tr( "Threads" ) );
    mainLayout->addWidget( hostThreadsLabel, 1, 0, 1, 1 );
    hostThreads = new QSpinBox( this );
    hostThreads->setMaximum( 64 );
    hostThreads->setValue( 8 );
    hostThreads->setToolTip( QString( "The entered value will be a maximum.\n" ) +
                             QString( "If the computer has fewer cores, " ) +
                             QString( "the number of threads will be " ) +
                             QString( "reduced accordingly." ) ); ;
    mainLayout->addWidget( hostThreads, 1, 1, 1, 1 );


    buttonBox = new QDialogButtonBox( QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, this );
    connect( buttonBox, SIGNAL( accepted() ), this, SLOT( accepting() ) );
    //connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
    connect( buttonBox, SIGNAL( rejected() ), this, SLOT( reject() ) );
    mainLayout->addWidget( buttonBox, 3, 1, 1, 3 );


    if( mst ) createMasterLayout();
    else createSlaveLayout();

    createActions();

    //setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed );

}

HostAdd::~HostAdd( ) { }

void HostAdd::init( void ) {

}

void HostAdd::contextMenuEvent( QContextMenuEvent* event ) {

    QMenu menu( this );
    menu.addAction( addAct );
    menu.addAction( removeAct );
    menu.addAction( moveUpAct );
    menu.addAction( moveDownAct );
    menu.exec( event->globalPos() );

}

void HostAdd::createMasterLayout( void ) {

    setWindowTitle( "Add Master" );
    //newHost->hostType = TYPE_MASTER;

    hostTunnelLabel = new QLabel( this );
    hostTunnelLabel->setText( tr( "Use SSH Tunnel" ) );
    mainLayout->addWidget( hostTunnelLabel, 1, 2, 1, 1 );
    hostTunnelCB = new QCheckBox( this );
    hostTunnelCB->setToolTip( QString( "Check this box to use SSH tunneling\n" ) +
                              QString( "to connect to this master." ) ); ;
    mainLayout->addWidget( hostTunnelCB, 1, 3, 1, 1 );


}

void HostAdd::createSlaveLayout( void ) {

    setWindowTitle( "Add Slave" );
    //newHost->hostType = TYPE_SLAVE;

}

void HostAdd::accepting( void ) {

    //newHost->hostName = hostName->text().toStdString();
    //newHost->listeningPort = hostPort->value();
    //newHost->nThreads = hostThreads->value();

    accept();
}


void HostAdd::createActions( void ) {

    addAct = new QAction( tr( "&Add" ), this );
    addAct->setShortcut( tr( "Ctrl+A" ) );
    addAct->setStatusTip( tr( "Add Filter at end." ) );
    //connect(addAct, SIGNAL(triggered()), this, SLOT(AddFilter()));

    removeAct = new QAction( tr( "&Remove" ), this );
    removeAct->setShortcut( tr( "Ctrl+R" ) );
    removeAct->setStatusTip( tr( "Remove Filter" ) );
    //connect(removeAct, SIGNAL(triggered()), this, SLOT(open()));

    moveUpAct = new QAction( tr( "Move &Up" ), this );
    moveUpAct->setShortcut( tr( "Ctrl+U" ) );
    moveUpAct->setStatusTip( tr( "Move Filter Up" ) );
    //connect(moveUpAct, SIGNAL(triggered()), this, SLOT(save()));

    moveDownAct = new QAction( tr( "Move &Down" ), this );
    moveDownAct->setShortcut( tr( "Ctrl+D" ) );
    moveDownAct->setStatusTip( tr( "Move Filter Down" ) );
    //connect(moveDownAct, SIGNAL(triggered()), this, SLOT(print()));

}
