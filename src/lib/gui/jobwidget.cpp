#include "redux/gui/jobwidget.hpp"

#include "redux/gui/jobmodel.hpp"

#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
using std::string;
using std::ifstream;
using std::getline;
using std::cout;
using std::endl;

// #include "MainWindow.h"
//
// #include <REDUX/Model.h>
//#include <ROYAC.h>
//using ROYAC::PI;

using namespace redux::gui;
using namespace redux;
using namespace std;

JobWidget::JobWidget( QMainWindow* parent, const char* name )
    : QWidget( parent ), activeHost( 0 ), jobModel(jobs), peerModel(peers), statusThread( NULL ) {
        
    myMainWindow = parent;

    createLayout();
    createActions();
    //readSettings();
    emit hostsChanged();
    peerTree->reset();

    //startStatusThread();
}

JobWidget::~JobWidget( ) {
    stopStatusThread();
    //writeSettings();
    //for (int i=0; i<myJobs.size(); ++i) delete myJobs[i];

}

void JobWidget::contextMenuEvent( QContextMenuEvent* event ) {

    QMenu menu( peerTree );
    menu.addAction( addAct );
    menu.addAction( removeAct );
    menu.exec( event->globalPos() );

}


void JobWidget::createLayout() {

    mainLayout = new QGridLayout( this );

    jobView = new QTextEdit( this );
    jobDW = new QDockWidget;
    jobView->setReadOnly( TRUE );
    jobDW->setObjectName( "jobDW" );
    jobDW->setWindowTitle( "Job-Log" );
    jobDW->setWidget( jobView );
    myMainWindow->addDockWidget( Qt::BottomDockWidgetArea, jobDW );
    //myMainWindow->viewMenu->addAction(jobDW->toggleViewAction());

    sysView = new QTextEdit( this );
    sysDW = new QDockWidget;
    sysView->setReadOnly( TRUE );
    sysDW->setObjectName( "sysDW" );
    sysDW->setWindowTitle( "Sys-Log" );
    sysDW->setWidget( sysView );
    myMainWindow->addDockWidget( Qt::BottomDockWidgetArea, sysDW );
    //myMainWindow->viewMenu->addAction(sysDW->toggleViewAction());

    //tabWidget = new QTabWidget;
    //setCentralWidget(tabWidget);

    netWidget = new QWidget( this );
    netDW = new QDockWidget;
    netDW->setObjectName( "netDW" );
    netDW->setWindowTitle( "Network" );
    netDW->setWidget( netWidget );

    netWidget->setMaximumWidth( 350 );
    netWidget->setMinimumWidth( 350 );
    netWidget->setMinimumHeight( 300 );
    netWidget->setSizePolicy( QSizePolicy::Fixed, QSizePolicy::Expanding );

    netLayout = new QGridLayout( netWidget );

    masterFrame = new QFrame( this );
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
    masterName->setText( tr( "polar" ) );
    //masterLayout->addWidget( masterLabel, 0, 0, 1, 1 );
    masterLayout->addWidget( masterName, 0, 0, 1, 1 );

    masterPort = new QSpinBox( this );
    masterPort->setMaximum( 65536 );
    masterPort->setValue( 30000 );
    masterLayout->addWidget( masterPort, 0, 2, 1, 1 );

    masterCB = new QCheckBox( this );
    masterCB->setText( tr( "SSH tunnel" ) );
    masterLayout->addWidget( masterCB, 0, 3, 1, 1 );

    masterLog = new QLineEdit( this );
    masterLog->setText( tr( "./logfile.txt" ) );
    //masterLayout->addWidget( masterLogLabel, 1, 0, 1, 1 );
    masterLayout->addWidget( masterLog, 1, 0, 1, 1 );

    masterButton = new QPushButton( this );
    masterButton->setText( tr( "connect" ) );
    masterLayout->addWidget( masterButton, 2, 3, 1, 1 );

    //hostModel.setRootObject( &myNet );
    peerTree = new PeerTreeView( this );
    peerTree->setModel( &peerModel );
    peerTree->setSelectionBehavior( QAbstractItemView::SelectRows );

    //jobModel.setRootObject( jobs );
    jobTree = new JobTreeView( this );
    jobTree->setModel( &jobModel );
    jobTree->setSelectionBehavior( QAbstractItemView::SelectRows );

    netLayout->addWidget( masterFrame, 0, 0, 3, 3 );
    //netLayout->addWidget( slavesLabel, 4, 0, 1, 1 );
    netLayout->addWidget( peerTree, 3, 0, 3, 3 );
    //netLayout->addWidget( jobsLabel, 8, 0, 1, 1 );
    netLayout->addWidget( jobTree, 6, 0, 3, 3 );

    netLayout->setRowStretch( 4, 10 );
    netLayout->setRowStretch( 7, 10 );

    connect( peerTree, SIGNAL( doubleClicked( QModelIndex ) ), this, SLOT( openHost( QModelIndex ) ) );
    connect( this, SIGNAL( hostsChanged() ), &peerModel, SLOT( refreshTree() ) );
    connect( this, SIGNAL( hostsChanged() ), this, SLOT( refreshPeerTab() ) );

    connect( jobTree, SIGNAL( doubleClicked( QModelIndex ) ), this, SLOT( openJob( QModelIndex ) ) );
    connect( this, SIGNAL( jobsChanged() ), &jobModel, SLOT( refreshTree() ) );
    connect( this, SIGNAL( jobsChanged() ), this, SLOT( refreshJobTab() ) );

    connect( peerTree, SIGNAL( addHost() ), this, SLOT( addHost() ) );
    connect( peerTree, SIGNAL( removeHost() ), this, SLOT( removeHost() ) );
    connect( peerTree, SIGNAL( launchHost() ), this, SLOT( launchHost() ) );
    connect( peerTree, SIGNAL( killHost() ), this, SLOT( killHost() ) );

    connect( jobTree, SIGNAL( addJob() ), this, SLOT( addJob() ) );
    connect( jobTree, SIGNAL( removeJob() ), this, SLOT( removeJob() ) );
    connect( jobTree, SIGNAL( cloneJob() ), this, SLOT( cloneJob() ) );
    connect( jobTree, SIGNAL( submitJob() ), this, SLOT( submitJob() ) );
    connect( jobTree, SIGNAL( moveJobUp() ), this, SLOT( moveJobUp() ) );
    connect( jobTree, SIGNAL( moveJobDown() ), this, SLOT( moveJobDown() ) );

    myMainWindow->addDockWidget( Qt::LeftDockWidgetArea, netDW );
    //((MainWindow*)myMainWindow)->viewMenu->addAction(netDW->toggleViewAction());

    mainTabs = new QTabWidget( this );
    peerTab = new PeerTab( this );
    jobTab = new JobTab( jobs, this );
    //tabsDW = new QDockWidget;
    //tabsDW->setObjectName("tabsDW");
    //tabsDW->setWindowTitle("Tabs");
    //tabsDW->setWidget( mainTabs );
    //hostTab->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );
    //tabsDW->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );
    mainTabs->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );
    mainTabs->addTab( peerTab, "Host" );
    mainTabs->addTab( jobTab, "Job" );
    //hostTab->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );  //jobTab->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );

    mainLayout->addWidget( mainTabs, 0, 0, 3, 3 );
    mainLayout->setRowStretch( 1, 10 );
    mainLayout->setColumnStretch( 1, 10 );

    //myMainWindow->addDockWidget(Qt::RightDockWidgetArea, tabsDW);
    // myMainWindow->setCentralWidget( mainTabs );
    //logview = new QTextEdit( this );
    //plot = new PlotWidget( this );

    //controlTabs->setMaximumWidth(248);
    //controlTabs->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);

    //logview->setMaximumHeight(60);
    //logview->setReadOnly( TRUE );
    //logview->setTabChangesFocus( TRUE );
    //logview->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

    //mainLayout->addWidget( mainTabs, 0, 0, 1, 3 );
    //gridLayout->addWidget( plot, 0, 1, 1, 1 );
    //mainLayout->addWidget( logview, 1, 0, 1, 3 );
    //   gridLayout->addWidget( test, 1, 0, 1, 1 );

    //gridLayout->setRowStretch( 0, 0 );
    //gridLayout->setRowStretch( 1, 0 );
    //mainLayout->setColumnMinimumWidth( 0, 248 );
    //mainLayout->setColumnMinimumWidth( 1, 300 );
    //gridLayout->setColumnStretch( 0, 0 );
    //gridLayout->setColumnStretch( 1, 0 );

    //configTab = new ConfigTab( this );
    //mainTabs->addTab(configTab, "Config");

    //specTab = new QWidget( this );
    //mainTabs->addTab(specTab, "Spectrum");

    //instrTab = new QWidget( this );
    //mainTabs->addTab(instrTab, "Instr.");



}

void JobWidget::createActions( void ) {

    addAct = new QAction( tr( "&Add" ), this );
    addAct->setShortcut( tr( "Ctrl+A" ) );
    addAct->setStatusTip( tr( "Add Filter at end." ) );
    //connect(addAct, SIGNAL(triggered()), this, SLOT(AddFilter()));

    removeAct = new QAction( tr( "&Remove" ), this );
    removeAct->setShortcut( tr( "Ctrl+R" ) );
    removeAct->setStatusTip( tr( "Remove Filter" ) );
    //connect(removeAct, SIGNAL(triggered()), this, SLOT(open()));

}


void JobWidget::loadFile( const char* filename ) {

    ifstream my_file( filename );
    string buffer;

    if( !my_file ) { printMsg( "Can not open file: " + QString( filename ) + " !!", Qt::red ); }
    else {
        printMsg( "Loading file: " + QString( filename ) + " ..." );
        while( !my_file.eof() ) {
            if( getline( my_file, buffer ) ) {
                if( buffer.find( "Atmosphere:" ) == 0 ) {
                    //atmosphere.parseLine( buffer.substr(11,buffer.length()).c_str() );
                    //refreshAtmosphere();
                }
                else if( buffer.find( "LambdaRegion:" ) == 0 ) {
                    //spectrum.parseRegion( buffer.substr(13,buffer.length()).c_str() );
                }
                else if( buffer.find( "SpectralLine:" ) == 0 ) {
                    //spectrum.parseLine( buffer.substr(13,buffer.length()).c_str() );
                }
                else if( buffer.find( "Instrument:" ) == 0 ) {
                    //instrument.parseLine( buffer.substr(11,buffer.length()).c_str() );
                }
            }
        }
        my_file.close();
    }
}

void JobWidget::readSettings() {

    QSettings settings;
    /*settings.beginGroup("JobWidget");
    QPoint pos = settings.value("position", QPoint(200, 200)).toPoint();
    QSize size = settings.value("size", QSize(400, 400)).toSize();
    settings.endGroup();
    resize(size);
    move(pos);*/

    settings.beginGroup( "MerlinJobLog" );
    bool bl = settings.value( "visible", true ).toBool();
    jobDW->setVisible( bl );
    bl = settings.value( "floating", false ).toBool();
    jobDW->setFloating( bl );
    QPoint pos = settings.value( "position", QPoint( 200, 200 ) ).toPoint();
    QSize size = settings.value( "size", QSize( 400, 400 ) ).toSize();
    Qt::DockWidgetArea area = ( Qt::DockWidgetArea ) settings.value( "dockingarea", Qt::BottomDockWidgetArea ).toInt();
    settings.endGroup();

    jobDW->resize( size );
    jobDW->move( pos );
    myMainWindow->addDockWidget( area, jobDW );

    settings.beginGroup( "MerlinSysLog" );
    bl = settings.value( "visible", true ).toBool();
    sysDW->setVisible( bl );
    bl = settings.value( "floating", false ).toBool();
    sysDW->setFloating( bl );
    pos = settings.value( "position", QPoint( 200, 200 ) ).toPoint();
    size = settings.value( "size", QSize( 400, 400 ) ).toSize();
    area = ( Qt::DockWidgetArea ) settings.value( "dockingarea", Qt::BottomDockWidgetArea ).toInt();
    settings.endGroup();

    sysDW->resize( size );
    sysDW->move( pos );
    myMainWindow->addDockWidget( area, sysDW );

}

void JobWidget::writeSettings() {

    QSettings settings;
    /*settings.beginGroup("MainWindow");
    settings.setValue("position", pos());
    settings.setValue("size", size());
    settings.endGroup();*/

    settings.beginGroup( "MerlinJobLog" );
    settings.setValue( "visible", jobDW->isVisible() );
    settings.setValue( "floating", jobDW->isFloating() );
    settings.setValue( "dockingarea", myMainWindow->dockWidgetArea( jobDW ) );
    settings.setValue( "position", jobDW->pos() );
    settings.setValue( "size", jobDW->size() );
    settings.endGroup();

    settings.beginGroup( "MerlinSysLog" );
    settings.setValue( "visible", sysDW->isVisible() );
    settings.setValue( "floating", sysDW->isFloating() );
    settings.setValue( "dockingarea", myMainWindow->dockWidgetArea( sysDW ) );
    settings.setValue( "position", sysDW->pos() );
    settings.setValue( "size", sysDW->size() );
    settings.endGroup();


}

void JobWidget::open() {

    QString fileName = QFileDialog::getOpenFileName( this, tr( "Open File" ), "", tr( "Config Files (*.xml)" ) );

//   if ( !fileName.isEmpty() && doc.LoadFile ( fileName.toStdString() ) ) {

//      int i=0;
//      TiXmlNode* node = 0;
//      TiXmlElement* element = 0;
    //cout << "Test1" << endl;
    /*if ( node = doc.FirstChild ( "NETWORK" ) ) myNet.parseXML ( node );
    else cout << "No <NETWORK> tag found in xml-file !!!" << endl;

    if ( node = doc.FirstChild ( "PARAMETERS" ) ) defaultModel.parseParameters ( node );
    if ( node = doc.FirstChild ( "LINES" ) ) defaultModel.parseLines ( node );
    //cout << defaultModel.dumpXML ( 2,true ) << endl;
    cout << "Test2" << endl;

    defaultGenetic.setParents ( NULL, &defaultModel );
    defaultMerlin.setParents ( NULL, &defaultModel );
    if ( node = doc.FirstChild ( "GENETIC" ) ) defaultGenetic.parseXML ( node );
    if ( node = doc.FirstChild ( "MERLIN" ) ) defaultMerlin.parseXML ( node );
    //cout << defaultGenetic.dumpXML ( 0, true ) << endl;
    //cout << defaultMerlin.dumpXML ( 0,true ) << endl;


    //JobConfig defaultJob ( defaultModel, defaultGenetic, defaultMerlin );
    //JobConfig* tmpJob;
    if ( ! ( node = doc.FirstChild ( "INVERSION" ) ) )
       cout << "No <INVERSION> tag found in xml-file !!!" << endl;
    else { //cout << myNet.dumpXML(2,true) << endl;
       while ( node ) {
     Job* tmpJob = new Inversion(); // defaultModel, &defaultGenetic, &defaultMerlin );// = defaultJob;
          tmpJob->parseXML ( node );
          myJobs.push_back ( tmpJob );
          //cout << "node->Value() " << ++i << "  " << node->Value() << endl;
          //cout << tmpJob.dumpXML(2,true) << endl;
          node = node->NextSibling ( "INVERSION" );
       }
    }*/

//if (node) cout << "node->ToText() " << ++i << "  " << node->Value() << endl;
    //else cout << "nothing..." << endl;
    //while ( node ) {
    //cout << "node->Value() " << ++i << "  " << node->Value() << endl;
    //node = node->NextSibling( "GENETIC" );
    //cout << "node->ToText() " << ++i << "  " << node->Value() << endl;
    //}
    //assert( node );
    //element = node->ToElement();
    //assert( element  );

    // Going to the toy store is now our second priority...
    // So set the "priority" attribute of the first item in the list.
    //node = element->FirstChildElement();  // This skips the "PDA" comment.
    //assert( node );
    //itemElement = node->ToElement();
    //assert( itemElement  );
    //itemElement->SetAttribute( "priority", 2 );

    //cfgFile.parse( myJobs );
    //cfgFile.parse ( myNet );


    //cout << myNet.dumpXML(2, true) << endl;
//      myNet.clean();
    // }

    emit hostsChanged();
    peerTree->reset();

    emit jobsChanged();
    jobTree->reset();
    //for (int i=0; i<myJobs.size(); ++i) cout << *myJobs[i];
    //cout << myNet;
    //cout << dumpXML(true) << endl;
}

void JobWidget::save( void ) {

    QString fileName = QFileDialog::getSaveFileName( this, tr( "Save File" ), "", tr( "Config Files (*.xml)" ) );
    if( fileName.isEmpty() )
        return;

    QFile file( fileName );
    if( !file.open( QFile::WriteOnly | QFile::Text ) ) {
        QMessageBox::warning( this, tr( "Application" ),
                              tr( "Cannot write file %1:\n%2." )
                              .arg( fileName )
                              .arg( file.errorString() ) );
        return;
    }

    QTextStream out( &file );
    QApplication::setOverrideCursor( Qt::WaitCursor );
    //out << QString( dumpXML( true ).c_str() );
    QApplication::restoreOverrideCursor();

    //setCurrentFile(fileName);
    //statusBar()->showMessage ( tr ( "File saved" ), 2000 );
    return;


}

void JobWidget::openHost( const QModelIndex& ind ) {

    activeHost = ind.row();
    if( ind.parent().flags() ) activeHost++;

    emit hostsChanged();
    peerTree->reset();

}

void JobWidget::addHost( void ) {

//     HostInfo* newHost = new HostInfo;
//     bool mst = false;
//
//     if( myNet.hosts.size() < 1 ) mst = true;
//
//     HostAdd ha( newHost, mst, this );
//     ha.setModal( true );
//
//     if( ha.exec() ) {
//         myNet.hosts.push_back( newHost );
//         cout << "exec() returned TRUE." << endl;
//     }

    emit hostsChanged();
    peerTree->reset();

}

void JobWidget::removeHost( void ) {

    QModelIndexList selection = peerTree->selectedIndexes();
    vector<uint> indices;
    uint j;
    for( QList<QModelIndex>::iterator i = selection.begin(); i < selection.end(); i++ ) {
        if( i->column() != 0 ) continue;
        j = i->row();
        cout << "selected row = " << j << endl;
        if( i->parent().flags() ) j++;
        indices.push_back( j );
    }
    std::sort( indices.begin(), indices.end() );

    vector<uint>::iterator it;
    for( it = indices.end() - 1; it >= indices.begin(); --it ) {
        //myNet.hosts.erase( myNet.hosts.begin() + *it );
        if( *it < activeHost ) activeHost--;
        else if( *it == activeHost ) activeHost = -1;
    }

    emit hostsChanged();
    peerTree->reset();

}

void JobWidget::launchHost( void ) {

    cout << "launchHost()." << endl;
    startStatusThread();

}

void JobWidget::killHost( void ) {
    cout << "killHost()." << endl;
    stopStatusThread();
}

void JobWidget::refreshPeerTab( void ) {
//     if( activeHost < myNet.hosts.size() ) {
//         mainTabs->setTabText( mainTabs->indexOf( hostTab ), myNet.hosts[activeHost]->hostName.c_str() );
//         hostTab->loadHost( activeHost );
//         mainTabs->setCurrentWidget( hostTab );
//     }

}

void JobWidget::openJob( const QModelIndex& ind ) {
    activeJob = ind.row();
    //if ( ind.parent().flags() ) activeJob++;
    //int jobInd = mainTabs->indexOf( jobTab );

    emit jobsChanged();
    jobTree->reset();

}

void JobWidget::addJob( void ) { cout << "addJob()." << endl; }
void JobWidget::cloneJob( void ) { cout << "cloneJob()." << endl; }

void JobWidget::moveJobUp( void ) {

    QModelIndexList selection = jobTree->selectedIndexes();
    vector<int> indices;
    int j;
    for( QList<QModelIndex>::iterator i = selection.begin(); i < selection.end(); i++ ) {
        if( i->column() != 0 ) continue;
        j = i->row();
        if( i->parent().flags() ) j++;
        indices.push_back( j );
    }
    std::sort( indices.begin(), indices.end() );

    vector<int>::iterator it;
    for( it = indices.begin(); it < indices.end(); ++it ) {
        // if( *it > 0 ) swap( &myJobs[*it], &myJobs[*it - 1] );
    }

    activeJob--;
    emit jobsChanged();
    jobTree->reset();

}

void JobWidget::moveJobDown( void ) {

    QModelIndexList selection = jobTree->selectedIndexes();
    vector<int> indices;
    int j;
    for( QList<QModelIndex>::iterator i = selection.begin(); i < selection.end(); i++ ) {
        if( i->column() != 0 ) continue;
        j = i->row();
        if( i->parent().flags() ) j++;
        indices.push_back( j );
    }

    std::sort( indices.begin(), indices.end() );
    std::reverse( indices.begin(), indices.end() );

    vector<int>::iterator it;
    for( it = indices.begin(); it < indices.end(); ++it ) {
        //if( *it < myJobs.size() - 1 ) ROYAC::swap( &myJobs[*it], &myJobs[*it + 1] );
    }

    activeJob++;
    emit jobsChanged();
    jobTree->reset();

}

void JobWidget::removeJob( void ) {

    QModelIndexList selection = jobTree->selectedIndexes();
    vector<uint> indices;
    uint j;
    for( QList<QModelIndex>::iterator i = selection.begin(); i < selection.end(); i++ ) {
        if( i->column() != 0 ) continue;
        j = i->row();
        if( i->parent().flags() ) j++;
        indices.push_back( j );
    }
    std::sort( indices.begin(), indices.end() );

    vector<uint>::iterator it;
    for( it = indices.end() - 1; it >= indices.begin(); --it ) {
        cout << "removing job #" << *it << endl;
        //jobs.erase( jobs.begin() + *it );
        if( *it < activeJob ) activeJob--;
        else if( *it == activeJob ) activeJob = -1;
    }

    emit jobsChanged();
    jobTree->reset();

}

void JobWidget::submitJob( void ) { cout << "submitJob()." << endl; }

void JobWidget::refreshJobTab( void ) {
    if( activeJob < jobs.size() ) {
       // mainTabs->setTabText( mainTabs->indexOf( jobTab ), jobs[activeJob]->info.name.c_str() );
        jobTab->loadJob( activeJob );
        mainTabs->setCurrentWidget( jobTab );
    }

}

void JobWidget::printMsg( const QString &msg, QColor c ) {
    //logview->setTextColor( c );
    //logview->append( msg );
}
