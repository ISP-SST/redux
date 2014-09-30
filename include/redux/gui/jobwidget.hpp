#ifndef REDUX_GUI_JOBWIDGET_HPP
#define REDUX_GUI_JOBWIDGET_HPP

#include "redux/gui/peermodel.hpp"
#include "redux/gui/peertab.hpp"
#include "redux/gui/peertreeview.hpp"
#include "redux/gui/peeradd.hpp"
#include "redux/gui/jobmodel.hpp"
#include "redux/gui/jobtab.hpp"
#include "redux/gui/jobtreeview.hpp"
#include "redux/gui/statusthread.hpp"

#include "redux/job.hpp"

#include <string>

//#include <QMainWindow>
#include <QSettings>
#include <QTextEdit>
#include <QDockWidget>
#include <QTreeView>
#include <QtGui>

//#include <qwidget.h>
//#include <qtabwidget.h>
//#include <qlabel.h>
//#include <QDoubleSpinBox>
//#include <QHBoxLayout>
//#include <QVBoxLayout>
//#include <qcombobox.h>

//#include <tinyxml.h>
//#include <XMLFile.h>
// #include <REDUX/Job.h>
// #include <REDUX/InvMerlin.h>
// #include <REDUX/InvGenetic.h>
//#include <NetConfig.h>
//#include <HostInfo.h>

//#include "JobAdd.h"
//#include "MainWindow.h"
//class MainWindow;
//#include "MERLIN_MilneEddington.h"
//#include "MERLIN_Model.h"
//#include "MERLIN_Profile.h"

//#include "PlotWidget.h"
//#include "GUI_Atmosphere.h"
//#include "GUI_Region.h"
//#include "GUI_AtomicLine.h"
//#include "GUI_Filter.h"
//#include "GUI_StrayLight.h"

//#include "ConfigTab.h"

namespace redux {

    namespace gui {

        class JobWidget : public QWidget {
            Q_OBJECT

        public slots:
            void open();
            void save();

            void openHost( const QModelIndex& );
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

        signals:
            void hostsChanged( void );
            void jobsChanged( void );

        public:
            JobWidget( QMainWindow* parent = 0, const char* name = 0 );
            ~JobWidget();

            void loadFile( const char* filename );

            void readSettings();
            void writeSettings();

        protected:

            void contextMenuEvent( QContextMenuEvent* event );

        private:

//      TiXmlDocument doc;

            uint activeHost, activeJob;

            void createLayout( void );
            void createActions( void );
            //std::string dumpXML( bool comments = false );
            void startStatusThread( void ) {
                if( ! statusThread ) {
                    statusThread = new StatusThread( jobs );
                    statusThread->start();
                }
            };
            void stopStatusThread( void ) {
                if( statusThread ) {
                    statusThread->stop();
                    statusThread->wait();
                    delete statusThread;
                    statusThread = NULL;
                }
            };

            QDockWidget *jobDW, *sysDW, *netDW, *tabsDW;
            QTextEdit *jobView, *sysView;
            QFrame *masterFrame;

            QLabel *masterLabel, *masterLogLabel, *slavesLabel, *jobsLabel;
            QLineEdit *masterName, *masterLog;
            QSpinBox *masterPort;
            QPushButton *masterButton;
            QCheckBox *masterCB;
            //QComboBox *masterLogVerb,*masterLogReplace,*masterLogBuffer;

            QWidget *netWidget;

            JobTreeView *jobTree;
            JobModel jobModel;
            JobTab *jobTab;
            Job::JobSet jobs;

            PeerTreeView *peerTree;
            PeerModel peerModel;
            PeerTab *peerTab;
            std::set<Peer> peers;
            
            QTabWidget *mainTabs;
            QMainWindow *myMainWindow;

            StatusThread *statusThread;

            //Model defaultModel;
            //MerlinTask defaultMerlin;
            //GeneticTask defaultGenetic;

            //ROYAC::XMLFile cfgFile;

            QAction *addAct, *removeAct;
            //QWidget *specTab,*instrTab;

            //ConfigTab *configTab;

            QGridLayout *netLayout, *mainLayout, *masterLayout;
            //QVBoxLayout *specLayout,*atmLayout,*instrLayout;


            void printMsg( const QString &msg, QColor c = Qt::black );


        };

    }   // gui

}   // redux

#endif // REDUX_GUI_JOBWIDGET_HPP
