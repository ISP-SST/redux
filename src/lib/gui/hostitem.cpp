#include "redux/gui/hostitem.hpp"

#include <boost/date_time/posix_time/time_formatters.hpp>
#include <boost/format.hpp>

using namespace redux::gui;
using namespace std;

HostItem::HostItem(const QList<QVariant> &data) : itemData(data), hostIndex(-1) {}


HostItem::HostItem(const Host::Ptr& h, Ptr parent, int hi) : host(h), parentItem(parent), hostIndex(hi) {
    reset();
}


HostItem::~HostItem() {
    
}


void HostItem::reset(Host::Ptr h) {

    if(h && h != host) host = h;

    if(!host) return;

    using namespace boost::posix_time;
    
    string tmp;
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    
    itemData.clear();
    tmp = host->info.name; // + string( " (" ) + ipString( ( *rootObject )[i]->IP ) + string( ":" );
    itemData << tmp.c_str();
    itemData << to_string(host->info.pid).c_str();
    tmp = to_string(host->status.nThreads) + string("/") + to_string(host->info.nCores);       // + string( "/" )
    itemData << tmp.c_str();
    itemData << boost::str(boost::format("%3.1f%%") % host->status.load[0]).c_str();     //to_string( host->status.loadAvg ).c_str();
    time_duration elapsed = (second_clock::universal_time() - host->info.startedAt);
    itemData << to_simple_string(elapsed).c_str();//to_iso_extended_string( host->info.startedAt ).c_str(); //tsToString( tsSubtract( ts, ( *rootObject )[i]->startTime ) ).c_str();
    //elapsed = (second_clock::local_time() - host->status.lastSeen);
    itemData << to_iso_extended_string( host->status.lastSeen ).c_str(); //tsToString( tsSubtract( ts, ( *rootObject )[i]->lastSeen ) ).c_str();

    itemData << Host::TypeNames[host->info.peerType].c_str();


}


void HostItem::appendChild(const Host::Ptr& h) {
    
    for( auto& child: children ) {
        if (*(child->host) == *h) {
            child->host->info = h->info;
            child->host->status = h->status;
            child->reset();
            return;
        }
    }
    
    children.push_back( newItem(h,shared()) );
    
}


void HostItem::removeChild( const Host::Ptr& h ) {
 
    children.remove_if( [this,&h](const HostItem::Ptr& item){ return (*(item->host) == *h); });

}


HostItem::Ptr HostItem::child(int row) {
    
    if( row < 0 || row > (int)children.size() ) {
        return nullptr;
    }
    
    auto it = children.begin();
    while(row--) it++;
    
    return *it;
    
}


int HostItem::childCount() const {
    return children.size();
}


int HostItem::columnCount() const {
    return itemData.count();
}


int HostItem::indexOf( const Ptr& child ) {
    
    int i(-1);
    for( auto& it: children ) {
        ++i;
        if (it == child) {
            return i;
        }
    }
    
    return i;
    
}


QVariant HostItem::data(int column) const {
    return itemData.value(column);
}


HostItem::Ptr HostItem::parent() {
    return parentItem;
}

