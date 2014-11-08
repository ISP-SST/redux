#include "redux/gui/hostmodel.hpp"

#include "redux/gui/hostitem.hpp"
#include "redux/util/stringutil.hpp"

using namespace redux::gui;
using namespace redux;


HostModel::HostModel( const Host::Set& p, QObject *parent ) : QAbstractItemModel( parent ), hosts( p ) {
    
    QList<QVariant> rootData;
    rootData << "Hostname" << "PID" << "T/C" << "Load" << "Uptime" << "Last Seen" << "";
    rootItem.reset( new HostItem( rootData ) );
    
}


HostModel::~HostModel() {

}


int HostModel::columnCount( const QModelIndex &parent ) const {

    if( parent.isValid() ) {
        HostItem* hiPtr = static_cast<HostItem*>( parent.internalPointer() );
        if(hiPtr) {
            return hiPtr->columnCount();
        }
    }
    if(rootItem) {

        return rootItem->columnCount();
    }
    
    return 0;
    
}


QVariant HostModel::data( const QModelIndex &index, int role ) const {
    
    if( role != Qt::DisplayRole ){
        return QVariant();
    }
      
    if( index.isValid() ) {
        HostItem *item = static_cast<HostItem*>( index.internalPointer() );
        if(item) {
            return item->data( index.column() );
        }
    }
     
    return QVariant();
    
}


Qt::ItemFlags HostModel::flags( const QModelIndex &index ) const {
    if( !index.isValid() )
        return 0;

    return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}


QVariant HostModel::headerData( int section, Qt::Orientation orientation, int role ) const {
    if( orientation == Qt::Horizontal && role == Qt::DisplayRole )
        return rootItem->data( section );

    return QVariant();
}


QModelIndex HostModel::index( int row, int column, const QModelIndex &parent ) const {
    
    if( !hasIndex( row, column, parent ) ) {
        return QModelIndex();
    }
    
    HostItem::Ptr parentItem;

    if( !parent.isValid() ) {
        parentItem = rootItem;
    } else {
        HostItem* hiPtr = static_cast<HostItem*>( parent.internalPointer() );
        if(hiPtr) {
            parentItem = hiPtr->shared();
        }
    }
    
    if(parentItem) {
        HostItem::Ptr childItem = parentItem->child( row );
        if( childItem ) {
            return createIndex( row, column, childItem.get() );
        }
    }
    
    return QModelIndex();
    
}


QModelIndex HostModel::parent( const QModelIndex &index ) const {
    
    if( index.isValid() ) {
        HostItem* hiPtr = static_cast<HostItem*>(index.internalPointer());
        HostItem::Ptr childItem = hiPtr->shared();
        if(childItem) {
            HostItem::Ptr parentItem = childItem->parent();
            if( parentItem != rootItem ) {
                return createIndex( parentItem->indexOf(childItem), 0, parentItem.get() );
            }
        }
    }
    
    return QModelIndex();
    
}


int HostModel::rowCount( const QModelIndex &parent ) const {

    if( parent.column() > 0 ) {
        return 0;
    }
    
    if( parent.isValid() ) {
        HostItem* pPtr = static_cast<HostItem*>(parent.internalPointer());
        if( pPtr ) {
            HostItem::Ptr parentItem = pPtr->shared();
            return parentItem->childCount();
        }
    }
    
    return rootItem->childCount();
    
}


void HostModel::refreshTree( void ) {
    
    static Host::Set cache;

    for( auto &it: hosts ) {
        rootItem->appendChild( it );
        cache.erase(it);
    }

    for( auto &it: cache ) { // delete items that was removed between updates
        rootItem->removeChild( it );
    }

    cache = hosts;      // remember until next update
    
}

