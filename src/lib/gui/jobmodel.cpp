#include "redux/gui/jobmodel.hpp"


using namespace redux::gui;
using namespace redux;

JobModel::JobModel( const Job::JobSet& j, QObject *parent ) : QAbstractItemModel( parent ), jobs(j) {
    QList<QVariant> rootData;
    rootData << "Job" << "ID" << "Label" << "Type" << "Priority" << "User" << "Status";
    rootItem.reset( new JobItem( rootData ) );
    //setupModelData(data.split(QString("\n")), rootItem);
}


JobModel::~JobModel() {

}


int JobModel::columnCount( const QModelIndex &parent ) const {

    if( parent.isValid() ) {
        JobItem* jiPtr = static_cast<JobItem*>( parent.internalPointer() );
        if(jiPtr) {
            return jiPtr->columnCount();
        }
    }
    if(rootItem) {

        return rootItem->columnCount();
    }
    
    return 0;
    
}


QVariant JobModel::data( const QModelIndex &index, int role ) const {

    if( role != Qt::DisplayRole ){
        return QVariant();
    }
      
    if( index.isValid() ) {
        JobItem *item = static_cast<JobItem*>( index.internalPointer() );
        if(item) {
            return item->data( index.column() );
        }
    }
     
    return QVariant();
    
}


Qt::ItemFlags JobModel::flags( const QModelIndex &index ) const {
    if( !index.isValid() )
        return 0;

    return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}


QVariant JobModel::headerData( int section, Qt::Orientation orientation,
                               int role ) const {
    if( orientation == Qt::Horizontal && role == Qt::DisplayRole )
        return rootItem->data( section );

    return QVariant();
}


QModelIndex JobModel::index( int row, int column, const QModelIndex &parent ) const {
    
    if( !hasIndex( row, column, parent ) )
        return QModelIndex();

    JobItem::Ptr parentItem;

    if( !parent.isValid() ) {
        parentItem = rootItem;
    } else {
        JobItem* jiPtr = static_cast<JobItem*>( parent.internalPointer() );
        if(jiPtr) {
            parentItem = jiPtr->shared();
        }
    }
    
    if(parentItem) {
        JobItem::Ptr childItem = parentItem->child( row );
        if( childItem ) {
            return createIndex( row, column, childItem.get() );
        }
    }
    
    return QModelIndex();

}


QModelIndex JobModel::parent( const QModelIndex &index ) const {
   
    if( index.isValid() ) {
        JobItem* jiPtr = static_cast<JobItem*>(index.internalPointer());
        JobItem::Ptr childItem = jiPtr->shared();
        if(childItem) {
            JobItem::Ptr parentItem = childItem->parent();
            if( parentItem != rootItem ) {
                return createIndex( parentItem->indexOf(childItem), 0, parentItem.get() );
            }
        }
    }
    
    return QModelIndex();
    
}


int JobModel::rowCount( const QModelIndex &parent ) const {

    if( parent.column() > 0 ) {
        return 0;
    }
    
    if( parent.isValid() ) {
        JobItem* pPtr = static_cast<JobItem*>(parent.internalPointer());
        if( pPtr ) {
            JobItem::Ptr parentItem = pPtr->shared();
            return parentItem->childCount();
        }
    }
    
    return rootItem->childCount();
    
}



void JobModel::refreshTree( void ) {
    
    static Job::JobSet cache;

    for( auto &it: jobs ) {
        rootItem->appendChild( it );
        cache.erase(it);
    }

    for( auto &it: cache ) { // delete items that was removed between updates
        rootItem->removeChild( it );
    }

    cache = jobs;      // remember until next update

}
