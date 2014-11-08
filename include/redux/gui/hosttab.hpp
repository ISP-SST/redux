#ifndef REDUX_GUI_HOSTTAB_HPP
#define REDUX_GUI_HOSTTAB_HPP

#include "redux/gui/hosttreeview.hpp"
#include "redux/gui/hostmodel.hpp"
#include "redux/gui/jobtreeview.hpp"
#include "redux/gui/jobmodel.hpp"
#include "redux/network/host.hpp"

#include <mutex>

#include <QWidget>
#include <QAction>
#include <QMenu>
#include <QContextMenuEvent>
#include <QtGui>


// #include <Socket.h>
// #include <ThreadPool.h>
//
// #include <SysConfig.h>
//#include <REDUX/DaemonTasks.h>
//#include <REDUX/Jobs/Jobs.h>
/*
#include <pthread.h>

#include "PeerTreeView.h"
#include "PeerModel.h"
#include "JobTreeView.h"
#include "JobModel.h"*/

/*QT_BEGIN_NAMESPACE
class QAction;
class QActionGroup;
class QLabel;
class QMenu;
QT_END_NAMESPACE
*/

using redux::network::Host;
using redux::network::TcpConnection;

namespace redux {

    namespace gui {

        class HostTab : public QWidget, public redux::network::Host {
            
            Q_OBJECT

            TcpConnection::Ptr connection;

            QTimer timer;
            
            struct classcomp {
                bool operator() (const int& lhs, const int& rhs) const
                {return lhs<rhs;}
            };

            Job::JobSet jobs;
            Host::Set hosts;

        public:
            HostTab( Host&, TcpConnection::Ptr, QWidget* par=0 );
            virtual ~HostTab();

            void getJobs( void );
            void getHosts( void );

        signals:
            void infoChanged( void );

        private slots:
            void updateModified( int i ) { if( i ) { timer.setInterval( i * 1000 ); timer.start(); } else timer.stop(); };
            void updateInfo( void );
            void getInfo( void );

        protected:
            QWidget* parent;

            void contextMenuEvent( QContextMenuEvent* event );

        protected slots:


        private:
            
            bool swap_endian;

            void createLayout( void );
            void createActions( void );
            void loadValues( void );

            QAction *addAct, *removeAct, *moveUpAct, *moveDownAct;

            QGridLayout *tabLayout, *statLayout, *hostsLayout, *jobsLayout;
            QFrame *statFrame, *hostsFrame, *fsFrame, *jobsFrame;

            QLabel *hsLabel, *hsValueLabel, *laLabel, *laValueLabel, *globalLogLabel, *globalLognameLabel;
            QLabel *globalVerbLabel, *globalReplaceLabel, *updateIntervalLabel, *uptime, *uptimeLabel;
            QLabel *threadsLabel, *tcLabel;
            QComboBox *globalVerbCombo;
            QSpinBox *updateInterval;

            QLabel *hostsLabel, *masterLogLabel, *masterNameLabel, *slavesLabel, *jobsLabel;
            QLabel *slaveLogLabel, *slaveNameLabel;

            JobTreeView* jobList;
            JobModel jobModel;

            HostTreeView* hostList;
            HostModel hostModel;

        };

    }   // gui

}   // redux


#endif // REDUX_GUI_HOSTTAB_HPP
