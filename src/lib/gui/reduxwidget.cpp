#include "redux/gui/reduxwidget.hpp"

#include "redux/network/tcpconnection.hpp"

#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace redux::gui;
using namespace redux::network;
using namespace redux;
using namespace std;

ReduxWidget::ReduxWidget( QMainWindow* par, const char* name ) : QWidget( par ), info( NULL ) {

    parent = par;

    createLayout();
    createActions();

}

ReduxWidget::~ReduxWidget( ) {

}


void ReduxWidget::setInfo( QTextEdit* d ) {
    if( info ) {
        mainLayout->removeWidget( info );
        info = NULL;
    }
    if( d ) {
        info = d;
        mainLayout->addWidget( info, 5, 0, 1, 4 );
    }
}

void ReduxWidget::createLayout() {

    mainLayout = new QGridLayout( this );

    // *************** Side Panel ***************
    sidePanel = new QWidget( this );
    sidePanelLayout = new QGridLayout( sidePanel );

    sidePanel->setMaximumWidth( 250 );
    sidePanel->setMinimumWidth( 250 );
    sidePanel->setMinimumHeight( 500 );
    sidePanel->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Expanding );

    masterFrame = new QFrame( sidePanel );
    masterFrame->setFrameShape( QFrame::StyledPanel );
    masterLayout = new QGridLayout( masterFrame );

    //masterLabel = new QLabel( netWidget );
    //masterLabel->setText(tr("Master"));
    //masterLogLabel = new QLabel( netWidget );
    //masterLogLabel->setText(tr("SysLog"));
    //masterLabel->setFixedHeight ( 10 );
    //slavesLabel = new QLabel( netWidget );
    //slavesLabel->setText(tr("Slave List"));
    //jobsLabel = new QLabel( netWidget );
    //jobsLabel->setText(tr("Job List"));
    //masterLayout->addWidget( masterLabel, 0, 0, 1, 1 );

    masterName = new QLineEdit( this );
    masterName->setText( tr( "localhost" ) );
    //masterLayout->addWidget( masterLabel, 0, 0, 1, 1 );
    masterLayout->addWidget( masterName, 0, 0, 1, 1 );

    masterPort = new QSpinBox( this );
    masterPort->setMaximum( 65535 );
    masterPort->setMinimum( 1024 );
    masterPort->setValue( 30000 );
    masterLayout->addWidget( masterPort, 0, 2, 1, 1 );

    //masterCB = new QCheckBox ( this );
    //masterCB->setText ( tr ( "SSH tunnel" ) );
    //masterLayout->addWidget ( masterCB, 0, 3, 1, 1 );

    //masterLog = new QLineEdit ( this );
    //masterLog->setText ( tr ( "./logfile.txt" ) );
    //masterLayout->addWidget( masterLogLabel, 1, 0, 1, 1 );
    //masterLayout->addWidget ( masterLog, 1, 0, 1, 1 );

    masterButton = new QPushButton( this );
    masterButton->setText( tr( "connect" ) );
    masterLayout->addWidget( masterButton, 2, 2, 1, 1 );


    sidePanelLayout->addWidget( masterFrame, 0, 0, 1, 1 );

    connect( masterButton, SIGNAL( clicked() ), this, SLOT( tryConnection() ) );


    sidePanelLayout->setRowStretch( 8, 10 );
    mainLayout->addWidget( sidePanel, 0, 0, 3, 1 );
    // *************** End of Side Panel ***************

    // *************** Tabs  ***************
    mainTabs = new QTabWidget( this );
    //masterTab = new MasterTab ( this );
    mainTabs->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );
    mainTabs->setTabsClosable( true );
    //mainTabs->addTab ( ReduxTab, "Redux" );
    mainLayout->addWidget( mainTabs, 0, 1, 2, 2 );
    connect( mainTabs, SIGNAL( tabCloseRequested( int ) ), this, SLOT( closeTab( int ) ) );
    //   plotFrame = new QFrame( this );
    //   plotLayout = new QGridLayout( plotFrame );
    // plotFrame->setFrameShape( QFrame::StyledPanel );
    // mainLayout->addWidget( plotFrame, 0, 1, 3, 3 );
    // *************** End of Tabs ***************

    //mainLayout->addWidget( mainTabs, 0, 1, 2, 2 );
    mainLayout->setRowStretch( 1, 10 );
    mainLayout->setColumnStretch( 1, 10 );

    // connect( this, SIGNAL( setsChanged() ), &logModel, SLOT( refreshTree() ) );
    // connect( logTree, SIGNAL( doubleClicked( QModelIndex ) ), this, SLOT( openLog( QModelIndex ) ) );
    //connect(this, SIGNAL( setsChanged() ), this, SLOT( refreshSetTab() ) );

    return;




}

void ReduxWidget::createActions( void ) {

    addAct = new QAction( tr( "&Add" ), this );
    addAct->setShortcut( tr( "Ctrl+A" ) );
    addAct->setStatusTip( tr( "Add Filter at end." ) );
    //connect(addAct, SIGNAL(triggered()), this, SLOT(AddFilter()));

    removeAct = new QAction( tr( "&Remove" ), this );
    removeAct->setShortcut( tr( "Ctrl+R" ) );
    removeAct->setStatusTip( tr( "Remove Filter" ) );
    //connect(removeAct, SIGNAL(triggered()), this, SLOT(open()));

}

void ReduxWidget::loadFile( const char* filename ) {

    ifstream my_file( filename );
    string buffer;

    if( !my_file ) { dumpMsg( "Can not open file: " + QString( filename ) + " !!", Qt::red ); }
    else {
        dumpMsg( "Loading file: " + QString( filename ) + " ..." );
        while( !my_file.eof() ) {
            if( getline( my_file, buffer ) ) {
                /*if ( buffer.find ( "Atmosphere:" ) == 0 )
                {
                   //atmosphere.parseLine( buffer.substr(11,buffer.length()).c_str() );
                   //refreshAtmosphere();
                }
                else if ( buffer.find ( "LambdaRegion:" ) == 0 )
                {
                   //spectrum.parseRegion( buffer.substr(13,buffer.length()).c_str() );
                }
                else if ( buffer.find ( "SpectralLine:" ) == 0 )
                {
                   //spectrum.parseLine( buffer.substr(13,buffer.length()).c_str() );
                }
                else if ( buffer.find ( "Instrument:" ) == 0 )
                {
                   //instrument.parseLine( buffer.substr(11,buffer.length()).c_str() );
                }*/
            }
        }
        my_file.close();
    }
}

void ReduxWidget::readSettings() {

    QSettings settings;
    settings.beginGroup( "Redux" );
    //QPoint pos = settings.value("position", QPoint(200, 200)).toPoint();
    //QSize size = settings.value("size", QSize(400, 400)).toSize();
    settings.endGroup();
    //resize(size);
    //move(pos);


}

void ReduxWidget::writeSettings() {

    QSettings settings;
    settings.beginGroup( "Redux" );
    //settings.setValue("position", pos());
    //settings.setValue("size", size());
    settings.endGroup();

}

void ReduxWidget::open() {

    QStringList files = QFileDialog::getOpenFileNames( this, tr( "Open File" ), "", tr( "Log Files (*_lg*)" ) );

    int sz = PATH_MAX + 1; // N.B. It is possible to construct a path longer than PATH_MAX on most systems,
    // so this is really not fool-proof...
    char* buf = new char[sz];

    QStringList::Iterator it = files.begin();
    while( it != files.end() ) {
        memset( buf, 0, sz * sizeof( char ) );
        if( ( ! it->isEmpty() ) && realpath( it->toStdString().c_str(), buf ) ) {
            dumpMsg( QString( "Opening LogFile : " ) + *it );
            /*LogFile* tmpLog = new LogFile ( buf );
            bool skip = false;
            for ( uint i=0; i<myLogs.size(); ++i)
            if ( !(*(myLogs[i]) != *tmpLog) )
               skip = true;
             cout << *tmpLog << endl;
                 if ( ! skip ) myLogs.push_back( tmpLog );
                 else delete tmpLog;
                 //myLog.load();
                 */
        }
        ++it;
    }
    delete[] buf;

    //emit setsChanged();
    //logTree->reset();

    //emit jobsChanged();
    //jobTree->reset();
    //for (int i=0; i<myJobs.size(); ++i) cout << *myJobs[i];
    //cout << myNet;
    //cout << dumpXML(true) << endl;
}

void ReduxWidget::close() {

    QStringList files = QFileDialog::getOpenFileNames( this, tr( "Open File" ), "", tr( "Log Files (*_lg*)" ) );

    //int sz = PATH_MAX + 1; // N.B. It is possible to construct a path longer than PATH_MAX on most systems,
    // so this is really not fool-proof...
//     char buf[sz];
// 
//     QStringList::Iterator it = files.begin();
//     while( it != files.end() ) {
//         memset( buf, 0, sz * sizeof( char ) );
//         /*if ( (! it->isEmpty()) && realpath ( it->toStdString().c_str(), buf ) ) {
//         dumpMsg( QString("Opening LogFile : ")+*it);
//          LogFile* tmpLog = new LogFile ( buf );
//          bool skip = false;
//          for ( uint i=0; i<myLogs.size(); ++i)
//             if ( !(*(myLogs[i]) != *tmpLog) )
//                skip = true;
//             if ( ! skip ) myLogs.push_back( tmpLog );
//             else delete tmpLog;
//             //myLog.load();
// 
//             }*/
//         ++it;
//     }

    //emit setsChanged();
    //logTree->reset();

    //emit jobsChanged();
    //jobTree->reset();
    //for (int i=0; i<myJobs.size(); ++i) cout << *myJobs[i];
    //cout << myNet;
    //cout << dumpXML(true) << endl;
}

void ReduxWidget::save( void ) {

    dumpMsg( QString( "ReduxWidget::save() : " ), QColor( rand() % 256, rand() % 256, rand() % 256 ) );
    return;
    //mainPlot = new DisplayWidget();
    //mainPlot->show();

    return;
    /* QString fileName = QFileDialog::getSaveFileName ( this, tr ( "Save File" ), "", tr ( "Config Files (*.xml)" )  );
    if ( fileName.isEmpty() )
       return;

    QFile file ( fileName );
    if ( !file.open ( QFile::WriteOnly | QFile::Text ) )
    {
       QMessageBox::warning ( this, tr ( "Application" ),
                              tr ( "Cannot write file %1:\n%2." )
                              .arg ( fileName )
                              .arg ( file.errorString() ) );
       return;
    }

    QTextStream out ( &file );
    QApplication::setOverrideCursor ( Qt::WaitCursor );
    out << QString ( dumpXML ( true ).c_str() );
    QApplication::restoreOverrideCursor();
    */
    //setCurrentFile(fileName);
    //statusBar()->showMessage ( tr ( "File saved" ), 2000 );
    return;


}


void ReduxWidget::tryConnection( void ) {

    auto conn = TcpConnection::newPtr( ioservice );
    conn->connect( masterName->text().toStdString(), to_string( masterPort->value() ) );

    if( *conn ) {
        dumpMsg( string( "Connected to " ) + conn->socket().remote_endpoint().address().to_string() + string( ":" ) + masterPort->text().toStdString(), Qt::blue );

        uint8_t cmd = CMD_CONNECT;
        Host host;
        boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));
        boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));
        if( cmd == CMD_AUTH ) {
            // implement
        }
        if( cmd == CMD_CFG ) {  // handshake requested
            *conn << myInfo;
            *conn >> host.info;
            boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));       // ok or err
            //*conn >> cmd;       // ok or err
        }
        if( cmd != CMD_OK ) {
            dumpMsg( string( "Handshake with server failed: server replied " )+to_string((int)cmd), Qt::red );
            return;
        }

        mainTabs->addTab( new HostTab( host, conn, this ), masterName->text() + QString( ":" ) + masterPort->text() );

    }

}

void ReduxWidget::closeTab( int ind ) {

    delete mainTabs->widget( ind );

}


void ReduxWidget::openLog( const QModelIndex& ind ) {


}

void ReduxWidget::addHost( void ) {

    /*HostInfo newHost;
    bool mst = false;

    if ( myNet.hosts.size() < 1 ) mst = true;

    HostAdd ha ( &newHost, mst, this );
    ha.setModal ( true );

    if ( ha.exec() )
    {
       myNet.hosts.push_back ( newHost );
       cout << "exec() returned TRUE." << endl;
    }

    emit hostsChanged();
    hostTree->reset();
    */
}

void ReduxWidget::removeHost( void ) {

    /*   QModelIndexList selection = hostTree->selectedIndexes();
    vector<int> indices;
    int j;
    for ( auto &item: selection ) {
       if ( item.column() != 0 ) continue;
       j = item.row();
       cout << "selected row = " << j << endl;
       if ( item.parent().flags() ) j++;
       indices.push_back ( j );
    }
    std::sort ( indices.begin(), indices.end() );

    for ( auto it = indices.end()-1; it >= indices.begin(); --it )
    {
       myNet.hosts.erase ( myNet.hosts.begin() + *it );
       if ( *it < activeHost ) activeHost--;
       else if ( *it == activeHost ) activeHost = -1;
    }

    emit hostsChanged();
    hostTree->reset();
    */
}

void ReduxWidget::launchHost( void ) {

    cout << "launchHost()." << endl;
    //startStatusThread();

}

void ReduxWidget::killHost( void ) {
    cout << "killHost()." << endl;
    //stopStatusThread();
}

void ReduxWidget::refreshPeerTab( void ) {
    /*   if ( activeHost < myNet.hosts.size() ) {
       mainTabs->setTabText ( mainTabs->indexOf ( hostTab ), myNet.hosts[activeHost].hostName.c_str() );
       hostTab->loadHost ( activeHost );
       mainTabs->setCurrentWidget( hostTab );
    }
    */
}

void ReduxWidget::openJob( const QModelIndex& ind ) {
    activeJob = ind.row();
    //if ( ind.parent().flags() ) activeJob++;
    /*int jobInd = mainTabs->indexOf ( jobTab );

    emit jobsChanged();
    jobTree->reset();
    */
}

void ReduxWidget::addJob( void ) { cout << "addJob()." << endl; }
void ReduxWidget::cloneJob( void ) { cout << "cloneJob()." << endl; }

void ReduxWidget::moveJobUp( void ) {

    /*   QModelIndexList selection = jobTree->selectedIndexes();
    vector<int> indices;
    int j;
    for ( auto &item: selection ) {
       if ( item.column() != 0 ) continue;
       j = item.row();
       if ( item.parent().flags() ) j++;
       indices.push_back ( j );
    }
    std::sort ( indices.begin(), indices.end() );

    for ( auto it = indices.begin(); it != indices.end(); ++it ) {
       if (*it > 0) ROYAC::swap( &myJobs[*it], &myJobs[*it - 1]);
    }

    activeJob--;
    emit jobsChanged();
    jobTree->reset();
     */
}

void ReduxWidget::moveJobDown( void ) {

    /*   QModelIndexList selection = jobTree->selectedIndexes();
    vector<int> indices;
    int j;
    for ( auto &item: selection ) {
       if ( item.column() != 0 ) continue;
       j = item.row();
       if ( item.parent().flags() ) j++;
       indices.push_back ( j );
    }

    std::sort( indices.begin(), indices.end() );
    std::reverse( indices.begin(), indices.end() );

    for ( auto it = indices.begin(); it != indices.end(); ++it ) {
       if (*it < myJobs.size()-1) ROYAC::swap( &myJobs[*it], &myJobs[*it + 1]);
    }

    activeJob++;
    emit jobsChanged();
    jobTree->reset();
    */
}

void ReduxWidget::removeJob( void ) {

    /*   QModelIndexList selection = jobTree->selectedIndexes();
    vector<int> indices;
    int j;
    for ( auto &item: selection ) {
       if ( item.column() != 0 ) continue;
       j = item.row();
       if ( item.parent().flags() ) j++;
       indices.push_back ( j );
    }
    std::sort ( indices.begin(), indices.end() );

    for ( auto it = indices.end()-1; it >= indices.begin(); --it )
    {
       cout << "removing job #" << *it << endl;
       delete myJobs[*it];
       myJobs.erase ( myJobs.begin() + *it );
       if ( *it < activeJob ) activeJob--;
       else if ( *it == activeJob ) activeJob = -1;
    }

    emit jobsChanged();
    jobTree->reset();
    */
}

void ReduxWidget::submitJob( void ) { cout << "submitJob()." << endl; }

void ReduxWidget::refreshJobTab( void ) {
    /*   if ( activeJob < myJobs.size() ) {
       mainTabs->setTabText ( mainTabs->indexOf ( jobTab ), myJobs[activeJob]->label.c_str() );
       jobTab->loadJob ( activeJob );
       mainTabs->setCurrentWidget( jobTab );
    }
    */
}

void ReduxWidget::dumpMsg( const QString m, QColor c ) {

    if( info ) {
        QColor tmpColor = info->textColor();
        info->setTextColor( c );
        info->append( m );
        info->setTextColor( tmpColor );
        //dumper->update();
        info->repaint();
    }
    else emit msg( m, c );

}

