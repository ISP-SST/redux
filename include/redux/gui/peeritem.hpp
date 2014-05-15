#ifndef REDUX_GUI_PEERITEM_HPP
#define REDUX_GUI_PEERITEM_HPP

#include <QList>
#include <QVariant>


namespace redux {

    namespace gui {

        class PeerItem {
            
        public:
            PeerItem( const QList<QVariant> &data, PeerItem *parent = 0, int hi = -1 );
            ~PeerItem();

            void appendChild( PeerItem *child );

            int getHostIndex( void ) { return hostIndex; };
            PeerItem *child( int row );
            int childCount() const;
            int columnCount() const;
            QVariant data( int column ) const;
            int row() const;
            PeerItem *parent();
            void clear( void ) { childItems.clear(); }  //qDeleteAll(childItems); };

        private:
            QList<PeerItem*> childItems;
            QList<QVariant> itemData;
            PeerItem *parentItem;
            int hostIndex;
        };

    }   // gui

}   // redux

#endif // REDUX_GUI_PEERITEM_HPP
