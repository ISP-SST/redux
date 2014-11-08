#ifndef REDUX_GUI_PEERMODEL_HPP
#define REDUX_GUI_PEERMODEL_HPP

#include "redux/gui/hostitem.hpp"

#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>

using redux::network::Host;

namespace redux {

    namespace gui {

        class HostModel : public QAbstractItemModel {
            Q_OBJECT

        public:
            HostModel( const Host::Set&, QObject *parent = 0 );
            ~HostModel();

            QVariant data( const QModelIndex &index, int role ) const;
            Qt::ItemFlags flags( const QModelIndex &index ) const;
            QVariant headerData( int section, Qt::Orientation orientation, int role = Qt::DisplayRole ) const;
            QModelIndex index( int row, int column, const QModelIndex &parent = QModelIndex() ) const;
            QModelIndex parent( const QModelIndex &index ) const;
            int rowCount( const QModelIndex &parent = QModelIndex() ) const;
            int columnCount( const QModelIndex &parent = QModelIndex() ) const;

        public slots:
            void refreshTree( void );

        private:
            const Host::Set& hosts;

            HostItem::Ptr rootItem;

        };

    }   // gui

}   // redux


#endif // REDUX_GUI_PEERMODEL_HPP
