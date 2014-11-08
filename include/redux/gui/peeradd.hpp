#ifndef REDUX_GUI_PEERADD_HPP
#define REDUX_GUI_PEERADD_HPP

#include "redux/network/host.hpp"

#include <QDialog>
#include <QAction>
#include <QMenu>
#include <QContextMenuEvent>
#include <QtGui>

// #include <REDUX/NetConfig.h>
// #include <REDUX/HostInfo.h>
/*QT_BEGIN_NAMESPACE
class QAction;
class QActionGroup;
class QLabel;
class QMenu;
QT_END_NAMESPACE
*/
using redux::network::Host;

namespace redux {

    namespace gui {

        class HostAdd : public QDialog {
            Q_OBJECT

        public:
            HostAdd( Host::HostInfo* host, bool mst = false, QWidget* par = 0 );
            virtual ~HostAdd();

            void init( void );

        private slots:
            void accepting( void );

        protected:
            QWidget* parent;

            void contextMenuEvent( QContextMenuEvent* event );

        protected slots:


        private:

            void createMasterLayout( void );
            void createSlaveLayout( void );
            void createActions( void );

            Host::HostInfo* newHost;

            QAction *addAct, *removeAct, *moveUpAct, *moveDownAct;

            QGridLayout *mainLayout;
            QFrame *globalFrame, *masterFrame, *fsFrame, *slaveFrame, *jobFrame;

            QDialogButtonBox *buttonBox;

            QLabel *globalLabel, *globalStatusLabel, *globalLogLabel, *globalLognameLabel;
            QLabel *globalVerbLabel, *globalReplaceLabel;
            QComboBox *globalVerbCombo;

            QLineEdit *hostName;
            QSpinBox *hostPort, *hostThreads;
            QCheckBox *hostTunnelCB;
            QLabel *hostLabel, *hostLogLabel, *hostNameLabel,
                   *hostPortLabel, *hostThreadsLabel, *hostTunnelLabel;
            QLabel *slaveLabel, *slaveLogLabel, *slaveNameLabel;

        };

    }   // gui

}   // redux

#endif // REDUX_GUI_PEERADD_HPP
