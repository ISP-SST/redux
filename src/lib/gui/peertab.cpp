#include "redux/gui/peertab.hpp"

using namespace redux::gui;
using namespace redux;

PeerTab::PeerTab( QWidget* par, const char* name ) : QWidget( par ), modified( false ) {

    parent = par;

    createLayout();
    createActions();
    this->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );

}



void PeerTab::init( void ) {


}

void PeerTab::loadHost( int active ) {

    activeHost = active;
    loadValues();

}

void PeerTab::contextMenuEvent( QContextMenuEvent* event ) {

    QMenu menu( this );
    menu.addAction( addAct );
    menu.addAction( removeAct );
    menu.addAction( moveUpAct );
    menu.addAction( moveDownAct );
    menu.exec( event->globalPos() );

}

void PeerTab::createLayout( void ) {

    tabLayout = new QGridLayout( this );

    globalFrame = new QFrame( this );
    globalFrame->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );
    globalFrame->setFrameShape( QFrame::Box );  //StyledPanel);
    globalLayout = new QGridLayout( globalFrame );

    globalLabel = new QLabel( globalFrame );
    globalLabel->setText( tr( "Daemon Status:" ) );
    globalStatusLabel = new QLabel( globalFrame );
    globalStatusLabel->setText( tr( "OFFLINE" ) );

    globalLogLabel = new QLabel( globalFrame );
    globalLogLabel->setText( tr( "System Log:" ) );
    globalLognameLabel = new QLabel( globalFrame );
    globalLognameLabel->setText( tr( "/path/to/sysLog.log" ) );

    globalVerbLabel = new QLabel( globalFrame );
    globalVerbLabel->setText( tr( "Verbosity" ) );
    //*globalVerbLabel,*globalReplaceLabel;
    globalVerbCombo = new QComboBox( globalFrame );
    globalVerbCombo->addItem( "0" );
    globalVerbCombo->addItem( "1" );
    globalVerbCombo->addItem( "2" );
    globalVerbCombo->addItem( "3" );
    globalVerbCombo->addItem( "4" );
    globalVerbCombo->addItem( "5" );
    globalVerbCombo->addItem( "6" );
    globalVerbCombo->addItem( "7" );
    globalVerbCombo->addItem( "8" );
    globalVerbCombo->setToolTip( QString( "0 - Silent\n" ) +
                                 QString( "1 - Fatal Errors only\n" ) +
                                 QString( "2 - All Errors\n" ) +
                                 QString( "3 - Errors & Warnings\n" ) +
                                 QString( "4 - Standard Messages\n" ) +
                                 QString( "5 - More Detail\n" ) +
                                 QString( "6 - All Messages\n" ) +
                                 QString( "7 - Debug\n" ) +
                                 QString( "8 - Debug Details" ) );

    globalLayout->addWidget( globalLabel, 0, 0, 1, 1 );
    globalLayout->addWidget( globalStatusLabel, 0, 3, 1, 2 );
    globalLayout->addWidget( globalLogLabel, 1, 0, 1, 1 );
    globalLayout->addWidget( globalLognameLabel, 1, 1, 1, 1 );
    globalLayout->addWidget( globalVerbLabel, 1, 3, 1, 1 );
    globalLayout->addWidget( globalVerbCombo, 1, 4, 1, 1 );
    globalLayout->setRowStretch( 1, 10 );
    globalLayout->setColumnStretch( 1, 10 );

    masterFrame = new QFrame( this );
    masterFrame->setFrameShape( QFrame::Box );  //StyledPanel);
    masterLayout = new QGridLayout( masterFrame );

    masterLabel = new QLabel( masterFrame );
    masterLabel->setText( tr( "Master" ) );
    masterLayout->addWidget( masterLabel, 0, 0, 1, 1 );
    //masterLogLabel = new QLabel( netWidget );
    //masterLogLabel->setText(tr("SysLog"));
    //masterLabel->setFixedHeight ( 10 );
    //slavesLabel = new QLabel( netWidget );
    //slavesLabel->setText(tr("Slave List"));
    //jobsLabel = new QLabel( netWidget );
    //jobsLabel->setText(tr("Job List"));
    //masterLayout->addWidget( masterLabel, 0, 0, 1, 1 );

    slaveFrame = new QFrame( this );
    slaveFrame->setFrameShape( QFrame::Box );  //StyledPanel);
    slaveLayout = new QGridLayout( slaveFrame );

    slaveLabel = new QLabel( slaveFrame );
    slaveLabel->setText( tr( "Slave" ) );
    //slaveLogLabel = new QLabel( slaveFrame );
    //slaveLogLabel->setText(tr("after"));

    slaveLayout->addWidget( slaveLabel, 0, 0, 1, 1 );


    //slaveFrame->hide();

    //tabLayout->addWidget( globalLabel, 0, 0, 1, 1 );
    tabLayout->addWidget( globalFrame, 0, 0, 1, 3 );
    tabLayout->addWidget( masterFrame, 1, 0, 1, 3 );
    tabLayout->addWidget( slaveFrame, 2, 0, 1, 3 );
    //tabLayout->addWidget( slaveLogLabel, 4, 3, 1, 1 );
    tabLayout->setRowStretch( 3, 10 );
    tabLayout->setColumnStretch( 2, 10 );

}


void PeerTab::createActions( void ) {

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

void PeerTab::loadValues( void ) {

    globalStatusLabel->setTextFormat( Qt::RichText );
    globalStatusLabel->setText( tr( "Testing something <font color=\"lime\">ONLINE</font> different..." ) );
    //globalStatusLabel->setForegroundRole ( QPalette::Highlight );
    globalLognameLabel->setText( ( std::string( ":" ) ).c_str() );
    //globalLognameLabel->setForegroundRole ( QPalette::Highlight );

    globalVerbCombo->setCurrentIndex( 0 );

}
