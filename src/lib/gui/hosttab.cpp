#include "redux/gui/hosttab.hpp"

#include "redux/job.hpp"
#include "redux/network/protocol.hpp"
#include "redux/util/endian.hpp"

using namespace redux::gui;
using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;

HostTab::HostTab( Host& m, TcpConnection::Ptr ptr, QWidget* par ) :
    QWidget( par ), Host(m.info), connection(ptr),
    parent( par ),
    jobModel( jobs ),
    hostModel( hosts ) {

    static Host::HostInfo myInfo;
    swap_endian = ( myInfo.littleEndian != info.littleEndian );

    createLayout();
    createActions();

    connect( this, SIGNAL( infoChanged() ), this, SLOT( updateInfo() ), Qt::AutoConnection );

    this->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );

    getInfo();

}


HostTab::~HostTab() {


}


void HostTab::updateInfo( void ) {

    laValueLabel->setText( "" );

    string tmp = "";

    tcLabel->setText( tmp.c_str() );

    struct timespec ts;
    clock_gettime( CLOCK_REALTIME, &ts );

    hostModel.refreshTree();
    hostList->reset();

    jobModel.refreshTree();
    jobList->reset();

}


void HostTab::getJobs( void ) {

    jobs.clear();

    if( !connection || !*connection ) {
        return;
    }

    uint8_t cmd = CMD_GET_JOBLIST;
    boost::asio::write(connection->socket(),boost::asio::buffer(&cmd,1));

    uint64_t blockSize;
    shared_ptr<char> buf = connection->receiveBlock( blockSize );
    
    if( !blockSize ) return;

    const char* ptr = buf.get();
    uint64_t count(0);
    try {
        while( count < blockSize ) {
            string tmpS = string( ptr+count );
            Job::JobPtr job = Job::newJob( tmpS );
            if( job ) {
                count += job->unpack( ptr+count, swap_endian );
                jobs.insert( job );
            }
            else throw invalid_argument( "Unrecognized Job tag: \"" + tmpS + "\"" );
        }
    }
    catch( const exception& e ) {
        cout << "addJobs: Exception caught while parsing block: " << e.what() << endl;
    }
    if( count != blockSize ) {
        cout << "addJobs: Parsing of datablock failed,  count = " << count << "  blockSize = " << blockSize << " bytes." << endl;
    }

}


void HostTab::getHosts( void ) {

    hosts.clear();

    if( !connection || !*connection ) {
        return;
    }

    uint8_t cmd = CMD_PSTAT;
    boost::asio::write(connection->socket(),boost::asio::buffer(&cmd,1));

    size_t blockSize;
    uint64_t count(0);
    shared_ptr<char> buf = connection->receiveBlock( blockSize );
    
    if( !blockSize ) return;

    const char* ptr = buf.get();
    try {
        Host host;
        while( count < blockSize ) {
            count += host.unpack( ptr+count, swap_endian );
            hosts.insert( Host::Ptr(new Host(host.info)) );
        }
    }
    catch( const exception& e ) {
        cout << "getPeers: Exception caught while parsing block: " << e.what() << endl;
    }
    
    if( count != blockSize ) {
        cout << "getPeers: Parsing of datablock failed,  count = " << count << "  blockSize = " << blockSize << " bytes." << endl;
    }


}


void HostTab::getInfo( void ) {

    if( !connection || !*connection ) {
        cout << "broken connection" << endl;
        return;
    }

    getJobs();
    getHosts();

    emit infoChanged();

}


void HostTab::contextMenuEvent( QContextMenuEvent* event ) {

    QMenu menu( this );
    menu.addAction( addAct );
    menu.addAction( removeAct );
    menu.addAction( moveUpAct );
    menu.addAction( moveDownAct );
    menu.exec( event->globalPos() );

}


void HostTab::createLayout( void ) {

    tabLayout = new QGridLayout( this );

    statFrame = new QFrame( this );
    statFrame->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );
    statFrame->setFrameShape( QFrame::Box );  //StyledPanel);
    statLayout = new QGridLayout( statFrame );

    hsLabel = new QLabel( statFrame );
    hsLabel->setText( tr( "Host Status:" ) );
    hsValueLabel = new QLabel( statFrame );
    laLabel = new QLabel( statFrame );
    laLabel->setText( tr( "Load Average:" ) );
    laValueLabel = new QLabel( statFrame );
    threadsLabel = new QLabel( statFrame );
    threadsLabel->setText( tr( "Threads/Cores" ) );
    tcLabel = new QLabel( statFrame );
    tcLabel->setToolTip( QString( "1/2/3/4\n1: Computation Threads\n2: IO Threads\n3: Total Threads\n4: Max number of IO/Comp. threads." ) );
    uptimeLabel = new QLabel( statFrame );
    uptimeLabel->setText( tr( "Uptime" ) );
    uptime = new QLabel( statFrame );

    globalLogLabel = new QLabel( statFrame );
    globalLogLabel->setText( tr( "SysLog:" ) );
    globalLognameLabel = new QLabel( statFrame );
    globalLognameLabel->setText( tr( "/path/to/sysLog.log" ) );

    updateIntervalLabel = new QLabel( statFrame );
    updateIntervalLabel->setText( tr( "Update Interval" ) );
    updateInterval = new QSpinBox( this );
    updateInterval->setMaximum( 255 );
    connect( updateInterval, SIGNAL( valueChanged( int ) ), this, SLOT( updateModified( int ) ) ); //, Qt::AutoConnection ); // DirectConnection );

    connect( &timer, SIGNAL( timeout() ), this, SLOT( getInfo() ) );
    timer.setSingleShot( false );
    timer.start();
    updateInterval->setValue( 5 );
    //globalLayout->addWidget ( globalVerbLabel, 1, 3, 1, 1 );
    //globalLayout->addWidget ( globalVerbCombo, 1, 4, 1, 1 );


    /*globalVerbLabel = new QLabel ( globalFrame );
    globalVerbLabel->setText ( tr ( "Verbosity" ) );
    // *globalVerbLabel,*globalReplaceLabel;
    globalVerbCombo = new QComboBox ( globalFrame );
    globalVerbCombo->addItem ( "0" );
    globalVerbCombo->addItem ( "1" );
    globalVerbCombo->addItem ( "2" );
    globalVerbCombo->addItem ( "3" );
    globalVerbCombo->addItem ( "4" );
    globalVerbCombo->addItem ( "5" );
    globalVerbCombo->addItem ( "6" );
    globalVerbCombo->addItem ( "7" );
    globalVerbCombo->addItem ( "8" );
    globalVerbCombo->setToolTip ( QString ( "0 - Silent\n" ) +
                                  QString ( "1 - Fatal Errors only\n" ) +
                                  QString ( "2 - All Errors\n" ) +
                                  QString ( "3 - Errors & Warnings\n" ) +
                                  QString ( "4 - Standard Messages\n" ) +
                                  QString ( "5 - More Detail\n" ) +
                                  QString ( "6 - All Messages\n" ) +
                                  QString ( "7 - Debug\n" ) +
                                  QString ( "8 - Debug Details" ) );
    */
    statLayout->addWidget( hsLabel, 0, 0, 1, 1 );
    statLayout->addWidget( hsValueLabel, 0, 1, 1, 1 );
    statLayout->addWidget( laLabel, 0, 3, 1, 1 );
    statLayout->addWidget( laValueLabel, 0, 4, 1, 1 );
    statLayout->addWidget( threadsLabel, 1, 3, 1, 1 );
    statLayout->addWidget( tcLabel, 1, 4, 1, 1 );
    statLayout->addWidget( uptimeLabel, 2, 3, 1, 1 );
    statLayout->addWidget( uptime, 2, 4, 1, 1 );
    statLayout->addWidget( globalLogLabel, 3, 0, 1, 1 );
    statLayout->addWidget( globalLognameLabel, 3, 1, 1, 1 );
    statLayout->addWidget( updateIntervalLabel, 3, 3, 1, 1 );
    statLayout->addWidget( updateInterval, 3, 4, 1, 1 );

    hostsFrame = new QFrame( statFrame );
    hostsFrame->setFrameShape( QFrame::Box );  //StyledPanel);
    hostsLayout = new QGridLayout( hostsFrame );

    //hostsLabel = new QLabel ( hostsFrame );
    //hostsLabel->setText ( tr ( "No Slaves" ) );
    //hostsLayout->addWidget ( hostsLabel, 0, 0, 1, 1 );

    hostList = new HostTreeView( hostsFrame );
    hostsLayout->addWidget( hostList, 1, 0, 1, 5 );
    //hostModel.setRootObject( &( sysConfig.hosts ) );
    hostList->setModel( &hostModel );
    hostList->setSelectionBehavior( QAbstractItemView::SelectRows );
    //connect( hostList, SIGNAL ( addLog() ), this, SLOT ( open() ) );

    jobList = new JobTreeView( hostsFrame );
    hostsLayout->addWidget( jobList, 2, 0, 1, 5 );
    //jobModel.setRootObject( jobs );
    jobList->setModel( &jobModel );
    jobList->setSelectionBehavior( QAbstractItemView::SelectRows );
    //connect( hostList, SIGNAL ( addLog() ), this, SLOT ( open() ) );

    jobsFrame = new QFrame( statFrame );
    jobsFrame->setFrameShape( QFrame::Box );  //StyledPanel);
    jobsLayout = new QGridLayout( jobsFrame );

    statLayout->addWidget( hostsFrame, 6, 0, 1, 5 );
    statLayout->addWidget( jobsFrame, 7, 0, 1, 5 );

    statLayout->setRowStretch( 1, 10 );
    statLayout->setColumnStretch( 1, 10 );

    tabLayout->addWidget( statFrame, 0, 0, 1, 3 );

    tabLayout->setRowStretch( 3, 10 );
    tabLayout->setColumnStretch( 2, 10 );

}


void HostTab::createActions( void ) {

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


