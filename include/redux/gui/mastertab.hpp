#ifndef REDUX_GUI_MASTERTAB_HPP
#define REDUX_GUI_MASTERTAB_HPP

#include "redux/gui/peertreeview.hpp"
#include "redux/gui/peermodel.hpp"
#include "redux/gui/jobtreeview.hpp"
#include "redux/gui/jobmodel.hpp"
#include "redux/network/peer.hpp"

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

using redux::network::Peer;

namespace redux {

    namespace gui {

        class MasterTab : public QWidget {
            Q_OBJECT

            Peer master;

            QTimer timer;

            Job::JobSet jobs;
            std::set<Peer> peers;

        public:
            MasterTab( Peer&, QWidget* par=0 );
            virtual ~MasterTab();

            void init( void );
            void loadHost( int );

            //void setPool( ThreadPool* p ) { if( true || pool || !p ) return; pool = p; };  // pool->addTask( &myTimer ); };
            bool getJobs( void );
            bool getPeers( void );

        signals:
            void infoChanged( void );

        private slots:
            void hostEdited( void ) { modified = true; };
            void updateModified( int i ) { if( i ) { timer.setInterval( i * 1000 ); timer.start(); } else timer.stop(); };
            void updateInfo( void );
            void getInfo( void );

        protected:
            QWidget* parent;

            void contextMenuEvent( QContextMenuEvent* event );

        protected slots:


        private:
            
            bool swap_endian;
            
            bool modified;
            int activeHost;

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

            PeerTreeView* peerList;
            PeerModel peerModel;

        };

    }   // gui

}   // redux


#endif // REDUX_GUI_MASTERTAB_HPP
