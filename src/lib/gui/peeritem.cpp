#include "redux/gui/peeritem.hpp"

#include <QStringList>

#include <iostream>

using namespace redux::gui;
using namespace redux;

PeerItem::PeerItem( const QList<QVariant> &data, PeerItem *parent, int hi ) {
    //std::cout << "PeerItem::PeerItem()" << ROYAC::hexString(this) << std::endl;
    parentItem = parent;
    itemData = data;
    hostIndex = hi;
}

PeerItem::~PeerItem() {
    qDeleteAll( childItems );
}

void PeerItem::appendChild( PeerItem *item ) {
    childItems.append( item );
}

PeerItem *PeerItem::child( int row ) {
    return childItems.value( row );
}

int PeerItem::childCount() const {
    return childItems.count();
}

int PeerItem::columnCount() const {
    return itemData.count();
}

QVariant PeerItem::data( int column ) const {
    return itemData.value( column );
}

PeerItem *PeerItem::parent() {
    return parentItem;
}

int PeerItem::row() const {
    if( parentItem && parentItem->childItems.size() > 0 )
        return parentItem->childItems.indexOf( const_cast<PeerItem*>( this ) );

    return 0;
}
