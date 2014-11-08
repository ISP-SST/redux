#ifndef REDUX_GUI_PEERTAB_HPP
#define REDUX_GUI_PEERTAB_HPP

#include "redux/network/host.hpp"

#include <QWidget>
#include <QAction>
#include <QMenu>
#include <QContextMenuEvent>
#include <QtGui>


//#include <REDUX/NetConfig.h>
/*QT_BEGIN_NAMESPACE
class QAction;
class QActionGroup;
class QLabel;
class QMenu;
QT_END_NAMESPACE
*/

namespace redux {

    namespace gui {


        class PeerTab : public QWidget {
            Q_OBJECT

        public:
            PeerTab( QWidget* par = 0, const char* name = 0 );
            virtual ~PeerTab() {};

            void init( void );
            void loadHost( int );

        private slots:
            void hostEdited( void ) { modified = true; };

        protected:
            QWidget* parent;

            void contextMenuEvent( QContextMenuEvent* event );

        protected slots:


        private:

            bool modified;
            int activeHost;

            void createLayout( void );
            void createActions( void );
            void loadValues( void );

            QAction *addAct, *removeAct, *moveUpAct, *moveDownAct;

            QGridLayout *tabLayout, *globalLayout, *masterLayout, *slaveLayout;
            QFrame *globalFrame, *masterFrame, *fsFrame, *slaveFrame, *jobFrame;

            QLabel *globalLabel, *globalStatusLabel, *globalLogLabel, *globalLognameLabel;
            QLabel *globalVerbLabel, *globalReplaceLabel;
            QComboBox *globalVerbCombo;

            QLabel *masterLabel, *masterLogLabel, *masterNameLabel, *slavesLabel, *jobsLabel;
            QLabel *slaveLabel, *slaveLogLabel, *slaveNameLabel;

        };

    }   // gui

}   // redux

#endif // REDUX_GUI_PEERTAB_HPP
