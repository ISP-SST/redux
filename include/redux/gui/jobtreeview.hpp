#ifndef REDUX_GUI_JOBTREEVIEW_HPP
#define REDUX_GUI_JOBTREEVIEW_HPP

#include <QTreeView>
#include <QAction>
#include <QMenu>
#include <QContextMenuEvent>


namespace redux {

    namespace gui {

        class ajAction : public QAction {

            Q_OBJECT

            int type;

        public:
            ajAction( QWidget *p, int i ) : QAction( p ), type( i ) {
                connect( this, SIGNAL( triggered() ), this, SLOT( qTriggered() ) );
            };

        private slots:
            void qTriggered( void ) { emit ajTriggered( type ); };

        signals:
            void ajTriggered( int );


        };

        class JobTreeView : public QTreeView {
            Q_OBJECT

            void createActions( void );

            QAction *addAct, *cloneAct, *removeAct, *moveUpAct, *moveDownAct, *submitAct;
            QList<ajAction*> addJobActs;

        protected:
            void contextMenuEvent( QContextMenuEvent* event );

        private slots:
            void addJob( int i ) { }; //cout << "P=" << i << endl;  };
            void cloning( void ) { emit cloneJob(); };
            void removing( void ) { emit removeJob(); };
            void submit( void ) { emit submitJob(); };
            void killing( void ) { emit killJob(); };
            void moveUp( void ) { emit moveJobUp(); };
            void moveDown( void ) { emit moveJobDown(); };

        signals:
            void cloneJob( void );
            void removeJob( void );
            void submitJob( void );
            void killJob( void );
            void moveJobUp( void );
            void moveJobDown( void );

        public:
            JobTreeView( QWidget *parent = 0 );
            ~JobTreeView() { for( int i = 0; i < addJobActs.size(); ++i ) delete addJobActs[i]; };

            QModelIndexList selectedIndexes( void ) const { return QTreeView::selectedIndexes(); };
            virtual void reset( void ) { QTreeView::reset(); expandAll(); for( int i = 0; i < 5; ++i ) resizeColumnToContents( i ); };

        };

    }   // gui

}   // redux

#endif // REDUX_GUI_JOBTREEVIEW_HPP
