#include "redux/gui/reduxmw.hpp"

#include <QtGui>

using namespace redux::gui;
using namespace redux;


ReduxMW::ReduxMW(QWidget* cw) {

    createLayout();
    createActions();
    createMenus();
    //createToolBars();
    //createStatusBar();

    readSettings();

    setCentralWidget( cw );
}

void ReduxMW::closeEvent( QCloseEvent *event ) {

    if( maybeSave() ) {
        writeSettings();
        event->accept();
    }
    else {
        event->ignore();
    }


}

void ReduxMW::about() {

    QLabel *icon = new QLabel;
    icon->setPixmap( QPixmap( ":/images/network.png" ) );

    QLabel *text = new QLabel;
    text->setWordWrap( true );
    text->setText( "<p>The <b>REDUX</b> gui is developed and maintained"
                   " by <a href=\"mailto:hillberg@astro.su.se\">Tomas Hillberg</a></p>"
                   "<p>Please report any bugs.</p>" );

    QPushButton *quitButton = new QPushButton( "OK" );

    QHBoxLayout *topLayout = new QHBoxLayout;
    topLayout->setMargin( 10 );
    topLayout->setSpacing( 10 );
    topLayout->addWidget( icon );
    topLayout->addWidget( text );

    QHBoxLayout *bottomLayout = new QHBoxLayout;
    bottomLayout->addStretch();
    bottomLayout->addWidget( quitButton );
    bottomLayout->addStretch();

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addLayout( topLayout );
    mainLayout->addLayout( bottomLayout );

    QDialog aboutDialog( this );
    aboutDialog.setModal( true );
    aboutDialog.setWindowTitle( tr( "About the REDUX GUI" ) );
    aboutDialog.setLayout( mainLayout );

    connect( quitButton, SIGNAL( clicked() ), &aboutDialog, SLOT( close() ) );

    aboutDialog.exec();

}

void ReduxMW::createLayout() {

    //cw = new ReduxWidget( this );
    output = new QTextEdit( this );
    output->setReadOnly( TRUE );
    //output->setMaximumHeight( 80 );
    output->setMinimumHeight( 20 );
    //output->append( QString("Bla") );
    //cw->setInfo( output );

    outputDW = new QDockWidget(tr("Output"), this);
    outputDW->setWidget(output);

    //setCentralWidget( cw );

}


void ReduxMW::createActions() {
    newAct = new QAction( QIcon( ":/images/filenew.xpm" ), tr( "&New" ), this );
    newAct->setShortcut( tr( "Ctrl+N" ) );
    newAct->setStatusTip( tr( "Create a new file" ) );
    //connect(newAct, SIGNAL(triggered()), this, SLOT(newFile()));

    openAct = new QAction( QIcon( ":/images/fileopen.xpm" ), tr( "&Open..." ), this );
    openAct->setShortcut( tr( "Ctrl+O" ) );
    openAct->setStatusTip( tr( "Open an existing file" ) );
    //connect( openAct, SIGNAL( triggered() ), cw, SLOT( open() ) );

    saveAct = new QAction( QIcon( ":/images/filesave.xpm" ), tr( "&Save" ), this );
    saveAct->setShortcut( tr( "Ctrl+S" ) );
    saveAct->setStatusTip( tr( "Save the document to disk" ) );
    //connect( saveAct, SIGNAL( triggered() ), cw, SLOT( save() ) );

    saveAsAct = new QAction( tr( "Save &As..." ), this );
    saveAsAct->setStatusTip( tr( "Save the document under a new name" ) );
    //connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveAs()));

    exitAct = new QAction( tr( "E&xit" ), this );
    exitAct->setShortcut( tr( "Ctrl+Q" ) );
    exitAct->setStatusTip( tr( "Exit the application" ) );
    connect( exitAct, SIGNAL( triggered() ), this, SLOT( close() ) );

    cutAct = new QAction( QIcon( ":/images/editcut.xpm" ), tr( "Cu&t" ), this );
    cutAct->setShortcut( tr( "Ctrl+X" ) );
    cutAct->setStatusTip( tr( "Cut the current selection's contents to the "
                              "clipboard" ) );
    //connect(cutAct, SIGNAL(triggered()), textEdit, SLOT(cut()));

    copyAct = new QAction( QIcon( ":/images/editcopy.xpm" ), tr( "&Copy" ), this );
    copyAct->setShortcut( tr( "Ctrl+C" ) );
    copyAct->setStatusTip( tr( "Copy the current selection's contents to the "
                               "clipboard" ) );
    //connect(copyAct, SIGNAL(triggered()), textEdit, SLOT(copy()));

    pasteAct = new QAction( QIcon( ":/images/editpaste.xpm" ), tr( "&Paste" ), this );
    pasteAct->setShortcut( tr( "Ctrl+V" ) );
    pasteAct->setStatusTip( tr( "Paste the clipboard's contents into the current "
                                "selection" ) );
    //connect(pasteAct, SIGNAL(triggered()), textEdit, SLOT(paste()));


    aboutAct = new QAction( tr( "&About" ), this );
    aboutAct->setStatusTip( tr( "Show the application's About box" ) );
    connect( aboutAct, SIGNAL( triggered() ), this, SLOT( about() ) );

    aboutQtAct = new QAction( tr( "About &Qt" ), this );
    aboutQtAct->setStatusTip( tr( "Show the Qt library's About box" ) );
    connect( aboutQtAct, SIGNAL( triggered() ), qApp, SLOT( aboutQt() ) );

    cutAct->setEnabled( false );
    copyAct->setEnabled( false );
    //connect(textEdit, SIGNAL(copyAvailable(bool)),
    //     cutAct, SLOT(setEnabled(bool)));
    //connect(textEdit, SIGNAL(copyAvailable(bool)),
    //      copyAct, SLOT(setEnabled(bool)));
}

void ReduxMW::createMenus() {
    fileMenu = menuBar()->addMenu( tr( "&File" ) );
    fileMenu->addAction( newAct );
    fileMenu->addAction( openAct );
    fileMenu->addAction( saveAct );
    fileMenu->addAction( saveAsAct );
    fileMenu->addSeparator();
    fileMenu->addAction( exitAct );

    editMenu = menuBar()->addMenu( tr( "&Edit" ) );
    editMenu->addAction( cutAct );
    editMenu->addAction( copyAct );
    editMenu->addAction( pasteAct );

    viewMenu = menuBar()->addMenu( tr( "&View" ) );
    viewMenu->addAction(outputDW->toggleViewAction());
    
    menuBar()->addSeparator();

    helpMenu = menuBar()->addMenu( tr( "&Help" ) );
    //helpMenu->addAction( aboutAct );
    helpMenu->addAction( aboutQtAct );
}

void ReduxMW::createToolBars() {
    fileToolBar = addToolBar( tr( "File" ) );
    fileToolBar->addAction( newAct );
    fileToolBar->addAction( openAct );
    fileToolBar->addAction( saveAct );

    editToolBar = addToolBar( tr( "Edit" ) );
    editToolBar->addAction( cutAct );
    editToolBar->addAction( copyAct );
    editToolBar->addAction( pasteAct );
}

void ReduxMW::createStatusBar() {
    statusBar()->showMessage( tr( "Ready" ) );
}

void ReduxMW::readSettings() {

    QSettings settings;
    settings.beginGroup( "ReduxMW" );
    QPoint pos = settings.value( "position", QPoint( 200, 200 ) ).toPoint();
    QSize size = settings.value( "size", QSize( 400, 400 ) ).toSize();
    //bool bl = settings.value( "outputVisible", true ).toBool();
    settings.endGroup();
    resize( size );
    move( pos );
    //output->setVisible( bl );
    //showOutput->setChecked( bl );

    settings.beginGroup ( "Output" );
    bool bl = settings.value ( "visible", true ).toBool();
    outputDW->setVisible ( bl );
    bl = settings.value ( "floating", false ).toBool();
    outputDW->setFloating ( bl );
    Qt::DockWidgetArea dwArea = ( Qt::DockWidgetArea ) settings.value ( "dockingarea", Qt::BottomDockWidgetArea ).toInt();
    pos = settings.value ( "position", QPoint ( 200, 200 ) ).toPoint();
    size = settings.value ( "size", QSize ( 400, 400 ) ).toSize();
    if(~bl) {
        addDockWidget ( dwArea, outputDW );
    }
    outputDW->resize ( size );
    outputDW->move ( pos );
    settings.endGroup();


}

void ReduxMW::writeSettings() {

    QSettings settings;
    settings.beginGroup( "ReduxMW" );
    settings.setValue( "position", pos() );
    settings.setValue( "size", size() );
    settings.endGroup();

    settings.beginGroup ( "Output" );
    settings.setValue ( "visible", outputDW->isVisible() );
    bool bl = outputDW->isFloating();
    settings.setValue ( "floating", bl );
    settings.setValue ( "dockingarea", dockWidgetArea ( outputDW ) );
    settings.setValue ( "position", outputDW->pos() );
    settings.setValue ( "size", outputDW->size() );
    settings.endGroup();

}

bool ReduxMW::maybeSave() {
    /*if (cw->isModified()) {
              int ret = QMessageBox::warning(this, tr("Application"),
                          tr("The document has been modified.\n"
                          "Do you want to save your changes?"),
                          QMessageBox::Yes | QMessageBox::Default,
                          QMessageBox::No,
                          QMessageBox::Cancel | QMessageBox::Escape);
              if (ret == QMessageBox::Yes)
              return save();
              else if (ret == QMessageBox::Cancel)
              return false;
          }*/
    return true;
}

void ReduxMW::printMsg( const QString &msg, QColor c ) {
    output->setTextColor( c );
    output->append( msg );
}

