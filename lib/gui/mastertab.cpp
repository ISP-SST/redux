#include "redux/gui/mastertab.hpp"

#include "redux/job.hpp"
#include "redux/network/protocol.hpp"
#include "redux/util/endian.hpp"

//#include <unistd.h>

using namespace redux::gui;
using namespace redux::network;
using namespace redux::util;
using namespace redux;
using namespace std;

MasterTab::MasterTab( Peer& m, QWidget* par ) :
    QWidget( par ),
    master( m ),
    parent( par ),
    modified( false ),
    jobModel( jobs ),
    peerModel( peers ) {

    static Peer::HostInfo myInfo;
    swap_endian = ( myInfo.littleEndian != master.host.littleEndian );

//     size_t blockSize = 2 + sizeof( size_t ) + myInfo.size();
//     char* tmp = new char[ blockSize ];
//     memset( tmp, 0 , blockSize );
//
//     char* ptr = tmp;
//     *ptr++ = CMD_CONNECT;
//     // *ptr++ = myEndian;
//     *( ( size_t* )ptr ) = blockSize - 1;
//     ptr += sizeof( size_t );
//     ptr = myInfo.pack( ptr );
//     int ret2 = 0;//sock->send( tmp, blockSize );
//     //if ( ret != blockSize )
//     //cout << "reduxd::conf() " << ret2 << " - " << blockSize << endl;
//
//
//     if( ret2 != blockSize ) {
//         cout << ( string( "MasterTab::MasterTab(1): Sending " ) + std::to_string( blockSize ) + string( " bytes, but only " ) + std::to_string( ret2 ) + string( " verified." ) ) << endl;
//     }
//
//     bool differentEndian = false;
//     blockSize = 1 + sizeof( size_t );
//     tmp = new char[ blockSize ];
//     ret2 = 0; //sock->peek( tmp, blockSize );
//     if( ret2 != blockSize ) {
//         cout << ( string( "MasterTab::MasterTab(2): Expected " ) + std::to_string( blockSize ) + string( " bytes, but received " ) + std::to_string( ret2 ) + string( "." ) ) << endl;
//         delete[] tmp;
//         return;
//     }
//
//     if( *tmp != 0 ) { //myEndian ) {
//         differentEndian = true;
//         //swapEndian( ( size_t* )( tmp + 1 ) );
//     }
//     blockSize = *( ( size_t* )( tmp + 1 ) );
//     delete[] tmp;
//     tmp = new char[ blockSize ];
//     memset( tmp, 0 , blockSize );
//     ret2 = 0; //sock->recv( tmp, blockSize );
//     if( ret2 != blockSize ) {
//         cout << ( string( "MasterTab::MasterTab(3): Expected " ) + std::to_string( blockSize ) + string( " bytes, but received " ) + std::to_string( ret2 ) + string( "." ) ) << endl;
//     }
//
//     ptr = tmp + 1 + sizeof( size_t );
//     //ptr = masterInfo.unpack( ptr, differentEndian );
//     delete[] tmp;
// //     if( ptr != ( tmp + blockSize ) ) {
// //         cout << ( string( "MasterTab::MasterTab(4): Unpacking of received HostInfo failed." ) ) << endl;
// //     }
//
//     //masterInfo.IP = sock->getPeerIP();

    createLayout();
    createActions();
    //cout << "MasterTab::MasterTab() 2" << endl;
    connect( this, SIGNAL( infoChanged() ), this, SLOT( updateInfo() ), Qt::AutoConnection );

    this->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );

    cout << "MasterTab::MasterTab() 3" << endl;
    getInfo();
    cout << "MasterTab::MasterTab() 4" << endl;
}

MasterTab::~MasterTab() {
    //cout << "MasterTab::~MasterTab()" << endl;
    //myTimer.exit();
//     if( conn ) {
//         *conn << CMD_DISCONNECT;
//         //cout << "MasterTab::~MasterTab() 2" << endl;
//     }
//     //cout << "MasterTab::~MasterTab() 3" << endl;

    //if ( pool ) myTimer.exit();

}



void MasterTab::init( void ) {


}

void MasterTab::loadHost( int active ) {

    activeHost = active;
    loadValues();

}

void MasterTab::updateInfo( void ) {

    //cout << "MasterTab::updateInfo()" << endl;
    laValueLabel->setText( "" );
//     if( masterInfo.loadAvg > masterInfo.nThreads[0] ) laValueLabel->setForegroundRole( QPalette::LinkVisited );
//     else laValueLabel->setForegroundRole( QPalette::WindowText );

    //updateInterval->setValue( myTimer.getInterval() );
    string tmp = ""; // toString( ( int )masterInfo.nThreads[1] ) + string( "/" ) + toString( ( int )masterInfo.nThreads[2] ) + string( "/" )
    //+ toString( ( int )masterInfo.nThreads[3] ) + string( "/" ) + toString( ( int )masterInfo.nThreads[0] );
    tcLabel->setText( tmp.c_str() );

    struct timespec ts;
    clock_gettime( CLOCK_REALTIME, &ts );
    //gettimeofday ( &ts, NULL );
    //uptime->setText( tsToString( tsSubtract( ts, masterInfo.startTime ) ).c_str() );

    // globalLognameLabel->setText( sysConfig.sysLog.file.name().c_str() );

    //cout << "MasterTab::updateInfo(5)" << endl;
    peerModel.refreshTree();
    peerList->reset();

    //if ( !sysConfig.hosts.size() ) hostList->hide();
    //else hostList->show();
    //cout << "MasterTab::updateInfo(7)" << endl;

    jobModel.refreshTree();
    jobList->reset();


}

bool MasterTab::getJobs( void ) {

    jobs.clear();   // TODO do something smarter...

    if( !master.conn || !*master.conn ) {
        return true;
    }

    cout << "getJobs 1" << endl;
    //*master.conn << CMD_GET_JOBLIST;
    uint8_t cmd = CMD_GET_JOBLIST;
    boost::asio::write(master.conn->socket(),boost::asio::buffer(&cmd,1));

    size_t blockSize;
    shared_ptr<char> buf = master.receiveBlock( blockSize, swap_endian );
   
    if( !blockSize ) return true;

    const char* cptr = buf.get();
    char* end = buf.get() + blockSize;
    bool modified = false;
    try {
        while( cptr < end ) {
            string tmpS = string( cptr );
            Job::JobPtr job = Job::newJob( tmpS );
            if( job ) {
                cptr = job->unpack( cptr, swap_endian );
                auto ret = jobs.insert( job );
                if( !ret.second ) {
                    //jobs.erase(ret.first);
                    //jobs.insert( job );
                    cout << "getJobs:  duplicate" << endl;
                    //ret.first->stat = job->stat;
                }
                else {
                    modified = true;
                }
            }
            else throw invalid_argument( "Unrecognized Job tag: \"" + tmpS + "\"" );
        }
    }
    catch( const exception& e ) {
        cout << "addJobs: Exception caught while parsing block: " << e.what() << endl;
    }
    if( cptr != end ) {
        cout << "addJobs: Parsing of datablock failed, there was a missmatch of " << ( cptr - end ) << " bytes." << endl;
    }

    cout << "getJobs E" << endl;
    return modified;

}


bool MasterTab::getPeers( void ) {

    peers.clear();   // TODO do something smarter...

    if( !master.conn || !*master.conn ) {
        return true;
    }
    
    cout << "getPeers 1" << endl;

    //*master.conn << CMD_PSTAT;
    uint8_t cmd = CMD_PSTAT;
    boost::asio::write(master.conn->socket(),boost::asio::buffer(&cmd,1));

    size_t blockSize;
    shared_ptr<char> buf = master.receiveBlock( blockSize, swap_endian );
    
    if( !blockSize ) return true;

    const char* cptr = buf.get();
    char* end = buf.get() + blockSize;
    bool modified = false;
    try {
        Peer peer;
        while( cptr < end ) {
            cptr = peer.unpack( cptr, swap_endian );
            auto ret = peers.insert( peer );
            if( !ret.second ) { // already in the set, replace with new version
                peers.erase( ret.first );
                peers.insert( peer );
            }
            modified = true;
        }
    }
    catch( const exception& e ) {
        cout << "getPeers: Exception caught while parsing block: " << e.what() << endl;
    }
    if( cptr != end ) {
        cout << "getPeers: Parsing of datablock failed, there was a missmatch of " << ( cptr - end ) << " bytes." << endl;
        return false;
    }

    cout << "getPeers E" << endl;
    return modified;

}



void MasterTab::getInfo( void ) {

cout << "MasterTab::getInfo(1)" << endl;

    if( !master.conn || !*master.conn ) {
        cout << "broken connection" << endl;
        return;
    }

    bool modified = getJobs() | getPeers();

    if( modified ) {
        emit infoChanged();
    }

cout << "MasterTab::getInfo(2)" << endl;
}

void MasterTab::contextMenuEvent( QContextMenuEvent* event ) {

    QMenu menu( this );
    menu.addAction( addAct );
    menu.addAction( removeAct );
    menu.addAction( moveUpAct );
    menu.addAction( moveDownAct );
    menu.exec( event->globalPos() );

}

void MasterTab::createLayout( void ) {

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

    peerList = new PeerTreeView( hostsFrame );
    hostsLayout->addWidget( peerList, 1, 0, 1, 5 );
    //hostModel.setRootObject( &( sysConfig.hosts ) );
    peerList->setModel( &peerModel );
    peerList->setSelectionBehavior( QAbstractItemView::SelectRows );
    //connect( hostList, SIGNAL ( addLog() ), this, SLOT ( open() ) );

    jobList = new JobTreeView( hostsFrame );
    hostsLayout->addWidget( jobList, 2, 0, 1, 5 );
    //jobModel.setRootObject( jobs );
    jobList->setModel( &jobModel );
    jobList->setSelectionBehavior( QAbstractItemView::SelectRows );
    //connect( hostList, SIGNAL ( addLog() ), this, SLOT ( open() ) );

    //hostList->hide();

    jobsFrame = new QFrame( statFrame );
    jobsFrame->setFrameShape( QFrame::Box );  //StyledPanel);
    jobsLayout = new QGridLayout( jobsFrame );

    //jobsLabel = new QLabel ( jobsFrame );
    //jobsLabel->setText ( tr ( "No Jobs" ) );
    //jobsLayout->addWidget ( jobsLabel, 0, 0, 1, 1 );


    statLayout->addWidget( hostsFrame, 6, 0, 1, 5 );
    statLayout->addWidget( jobsFrame, 7, 0, 1, 5 );

    //globalLayout->addWidget ( globalVerbLabel, 1, 3, 1, 1 );
    //globalLayout->addWidget ( globalVerbCombo, 1, 4, 1, 1 );
    statLayout->setRowStretch( 1, 10 );
    statLayout->setColumnStretch( 1, 10 );

    //masterLogLabel = new QLabel( netWidget );
    //masterLogLabel->setText(tr("SysLog"));
    //masterLabel->setFixedHeight ( 10 );
    //slavesLabel = new QLabel( netWidget );
    //slavesLabel->setText(tr("Slave List"));
    //jobsLabel = new QLabel( netWidget );
    //jobsLabel->setText(tr("Job List"));
    //masterLayout->addWidget( masterLabel, 0, 0, 1, 1 );


    //slaveFrame->hide();

    //tabLayout->addWidget( globalLabel, 0, 0, 1, 1 );
    tabLayout->addWidget( statFrame, 0, 0, 1, 3 );

    tabLayout->setRowStretch( 3, 10 );
    tabLayout->setColumnStretch( 2, 10 );

}


void MasterTab::createActions( void ) {

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

void MasterTab::loadValues( void ) {

    laLabel->setTextFormat( Qt::RichText );
    laLabel->setText( tr( "Testing something <font color=\"lime\">ONLINE</font> different..." ) );
    //globalStatusLabel->setForegroundRole ( QPalette::Highlight );
    //globalLognameLabel->setText ( (myNet->logFile.hostName + string(":") + myNet->logFile.fileName).c_str() );
    //globalLognameLabel->setForegroundRole ( QPalette::Highlight );

    //globalVerbCombo->setCurrentIndex( myNet->logVerbosity );

}

