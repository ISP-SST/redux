#ifndef REDUX_GUI_JOBTAB_HPP
#define REDUX_GUI_JOBTAB_HPP

#include "redux/job.hpp"

#include <QWidget>
#include <QAction>
#include <QMenu>
#include <QContextMenuEvent>
#include <QtGui>

//#include <REDUX/Job.h>
/*QT_BEGIN_NAMESPACE
class QAction;
class QActionGroup;
class QLabel;
class QMenu;
QT_END_NAMESPACE
*/

namespace redux {

    namespace gui {

        class JobTab : public QWidget {
            Q_OBJECT

        public:
            JobTab( Job::JobSet&, QWidget* par=0, const char* name=0 );
            virtual ~JobTab() {};

            void init( void );
            void loadJob( int );

        private slots:
            void JobEdited( void ) { modified = true; };

        protected:
            QWidget* parent;

            void contextMenuEvent( QContextMenuEvent* event );

        protected slots:


        private:

            Job::JobSet& jobs;

            bool modified;
            int activeJob;

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

#endif // REDUX_GUI_JOBTAB_HPP
