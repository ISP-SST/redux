#ifndef REDUX_GUI_REDUXMW_HPP
#define REDUX_GUI_REDUXMW_HPP

//#include <iostream>

#include <QMainWindow>
#include <QCloseEvent>
//#include <QTabWidget>
//#include <QGridLayout>

//#include <ThreadPool.h>

//#include "ReduxWidget.h"

//class MerlinWidget;
class ToolBar;
class QMenu;
//class QSignalMapper;
class QAction;
class QTextEdit;
//class QTabWidget;

//#include "PolRT_Widget.h"

namespace redux {

    namespace gui {

        class ReduxMW : public QMainWindow {

            Q_OBJECT

        public:
            ReduxMW(QWidget* cw=nullptr);
            ~ReduxMW() { };


        public slots:
            void printMsg( const QString &msg, QColor c = Qt::black );
            //void actionTriggered(QAction *action);

        protected:
            void closeEvent( QCloseEvent *event );
            
        private slots:
            void about();

        private:
            
            QTextEdit *output;
            QDockWidget *outputDW;

            QMenu *fileMenu, *editMenu, *viewMenu, *helpMenu;

            QToolBar *fileToolBar, *editToolBar;

            QAction *newAct, *openAct, *saveAct, *saveAsAct, *exitAct;
            QAction *cutAct, *copyAct, *pasteAct;
            QAction *aboutAct, *aboutQtAct;

            void createLayout();
            void createActions();
            void createMenus();
            void createToolBars();
            void createStatusBar();
            void readSettings();
            void writeSettings();
            bool maybeSave();


        }; // end of ReduxMW class

    }   // gui

}   // redux

#endif // REDUX_GUI_REDUXMW_HPP
