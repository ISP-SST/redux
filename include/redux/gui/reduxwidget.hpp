#ifndef REDUX_GUI_REDUXWIDGET_HPP
#define REDUX_GUI_REDUXWIDGET_HPP

#include "redux/gui/hosttab.hpp"

#include <string>

#include <QSettings>
#include <QTextEdit>
#include <QTreeView>
#include <QtGui>


namespace redux {

    namespace gui {

        class ReduxWidget : public QWidget {
            
            Q_OBJECT

        public:
            ReduxWidget( QMainWindow* par = NULL, const char* name = 0 );
            ~ReduxWidget();

            void loadFile( const char* filename );
            void setInfo( QTextEdit* d = NULL );

            void readSettings();
            void writeSettings();

        public slots:
            void open();
            void close();
            void save();
            void tryConnection();
            void closeTab( int );

            void openLog( const QModelIndex& );
            void addHost( void );
            void removeHost( void );
            void launchHost( void );
            void killHost( void );
            void refreshPeerTab( void );

            void openJob( const QModelIndex& );
            void addJob( void );
            void cloneJob( void );
            void moveJobUp( void );
            void moveJobDown( void );
            void removeJob( void );
            void submitJob( void );
            void refreshJobTab( void );

            void dumpMsg( const QString, QColor c = Qt::black );
            void dumpMsg( std::string s, QColor c = Qt::black ) { dumpMsg( QString( s.c_str() ), c ); };

        protected:

            //void contextMenuEvent ( QContextMenuEvent* event );

        signals:
            void setsChanged( void );
            void msg( const QString msg, QColor c = Qt::black );

        private:

            boost::asio::io_service ioservice;
            Host::HostInfo myInfo;
            
            uint activeHost, activeJob;

            void createLayout( void );
            void createActions( void );

            QTextEdit *info;
            QFrame *infoFrame, *plotFrame;

            QWidget *sidePanel;

            QFrame *masterFrame;

            QLabel *masterLabel, *masterLogLabel, *slavesLabel, *jobsLabel;
            QLineEdit *masterName, *masterLog;
            QSpinBox *masterPort;
            QPushButton *masterButton;
            QCheckBox *masterCB;
            QGridLayout *masterLayout, *netLayout;

            HostTab *masterTab;

            QTabWidget *mainTabs;
            QMainWindow *parent;

            QAction *addAct, *removeAct;

            QGridLayout *mainLayout, *sidePanelLayout, *plotLayout;


        };


    }  // namespace gui

}  // namespace redux

#endif // REDUX_GUI_REDUXWIDGET_HPP
