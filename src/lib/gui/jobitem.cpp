#include "redux/gui/jobitem.hpp"

#include <QStringList>

using namespace redux::gui;
using namespace redux;

JobItem::JobItem( const QList<QVariant> &data, JobItem *parent, int ji ) {
    parentItem = parent;
    itemData = data;
    jobIndex = ji;
}

JobItem::~JobItem() {
    qDeleteAll( childItems );
}

void JobItem::appendChild( JobItem *item ) {
    childItems.append( item );
}

JobItem *JobItem::child( int row ) {
    return childItems.value( row );
}

int JobItem::childCount() const {
    return childItems.count();
}

int JobItem::columnCount() const {
    return itemData.count();
}

QVariant JobItem::data( int column ) const {
    return itemData.value( column );
}

JobItem *JobItem::parent() {
    return parentItem;
}

int JobItem::row() const {
    if( parentItem && parentItem->childItems.size() > 0 )
        return parentItem->childItems.indexOf( const_cast<JobItem*>( this ) );

    return 0;
}
