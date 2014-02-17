#ifndef REDUX_GUI_PEERTREEVIEW_HPP
#define REDUX_GUI_PEERTREEVIEW_HPP

#include <QTreeView>
#include <QAction>
#include <QMenu>
#include <QContextMenuEvent>

namespace redux {

    namespace gui {

        class PeerTreeView : public QTreeView {
            Q_OBJECT

            void createActions( void );

            QAction *addAct, *removeAct, *launchAct, *killAct;

        protected:
            void contextMenuEvent( QContextMenuEvent* event );

        private slots:
            void adding( void ) { emit addHost(); };
            void removing( void ) { emit removeHost(); };
            void launching( void ) { emit launchHost(); };
            void killing( void ) { emit killHost(); };

        signals:
            void addHost( void );
            void removeHost( void );
            void launchHost( void );
            void killHost( void );

        public:
            PeerTreeView( QWidget *parent = 0 );
            ~PeerTreeView() {};

            QModelIndexList selectedIndexes( void ) const { return QTreeView::selectedIndexes(); };
            virtual void reset( void ) { QTreeView::reset(); for( int i = 0; i < 5; ++i ) resizeColumnToContents( i ); }; //expandAll(); };

        };

    }   // gui

}   // redux

#endif // REDUX_GUI_PEERTREEVIEW_HPP
