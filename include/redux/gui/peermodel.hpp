#ifndef REDUX_GUI_PEERMODEL_HPP
#define REDUX_GUI_PEERMODEL_HPP

#include "redux/gui/peeritem.hpp"
#include "redux/network/peer.hpp"

#include <vector>

#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>

using redux::network::Peer;

namespace redux {

    namespace gui {

        class PeerItem;

        class PeerModel : public QAbstractItemModel {
            Q_OBJECT

        public:
            PeerModel( const std::set<Peer>&, QObject *parent = 0 );
            ~PeerModel();

            void setData( const QString &data );

            QVariant data( const QModelIndex &index, int role ) const;
            Qt::ItemFlags flags( const QModelIndex &index ) const;
            QVariant headerData( int section, Qt::Orientation orientation,
                                 int role = Qt::DisplayRole ) const;
            QModelIndex index( int row, int column,
                               const QModelIndex &parent = QModelIndex() ) const;
            QModelIndex parent( const QModelIndex &index ) const;
            int rowCount( const QModelIndex &parent = QModelIndex() ) const;
            int columnCount( const QModelIndex &parent = QModelIndex() ) const;

        public slots:
            void refreshTree( void );

        private:
            void setupModelData( const QStringList &lines, PeerItem *parent );

            const std::set<Peer>& peers;

            PeerItem *rootItem;

        };

    }   // gui

}   // redux


#endif // REDUX_GUI_PEERMODEL_HPP
