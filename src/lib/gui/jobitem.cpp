#include "redux/gui/jobitem.hpp"

#include <QStringList>

using namespace redux::gui;
using namespace redux;
using namespace std;

JobItem::JobItem(const QList<QVariant> &data) : itemData(data), jobIndex(-1) {}


JobItem::JobItem(const Job::JobPtr& j, Ptr parent, int ji) : job(j), parentItem(parent), jobIndex(ji) {
    reset();
}


JobItem::~JobItem() {

}


void JobItem::reset(Job::JobPtr j) {

    if(j && j != job) job = j;

    if(!job) return;

    using namespace boost::posix_time;
    
    string tmp;
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    
        itemData.clear();
        itemData << std::to_string( 0 ).c_str();           // rownumber+1 = jobnumber
        itemData << std::to_string( job->info.id ).c_str(); // job ID
        itemData << job->info.name.c_str();             // job-label
        itemData << job->info.typeString.c_str();        // job-type
        itemData << std::to_string( job->info.priority ).c_str(); // priority
        itemData << job->info.user.c_str();              // username
        tmp = Job::stateString( job->info.state ) + " - " + job->stepString( job->info.step );
        itemData << tmp.c_str();     // current state
/*
    itemData.clear();
    tmp = host->info.name; // + string( " (" ) + ipString( ( *rootObject )[i]->IP ) + string( ":" );
    itemData << tmp.c_str();
    itemData << to_string(host->info.pid).c_str();
    tmp = to_string(host->status.nThreads) + string("/") + to_string(host->info.nCores);       // + string( "/" )
    itemData << tmp.c_str();
    itemData << boost::str(boost::format("%3.1f%%") % host->status.loadAvg).c_str();     //to_string( host->status.loadAvg ).c_str();
    time_duration elapsed = (second_clock::local_time() - host->info.startedAt);
    itemData << to_simple_string(elapsed).c_str();//to_iso_extended_string( host->info.startedAt ).c_str(); //tsToString( tsSubtract( ts, ( *rootObject )[i]->startTime ) ).c_str();
    elapsed = (second_clock::local_time() - host->status.lastSeen);
    itemData << to_simple_string(elapsed).c_str();//to_iso_extended_string( host->lastSeen ).c_str(); //tsToString( tsSubtract( ts, ( *rootObject )[i]->lastSeen ) ).c_str();

    itemData << Host::TypeNames[host->info.peerType].c_str();

*/
}


void JobItem::appendChild( const Job::JobPtr& j ) {
    
    for( auto& child: children ) {
        if (!(*(child->job) != *j)) {
            child->job = j;
            //child->job->info = h->info;
            //child->job->status = h->status;
            child->reset();
            return;
        }
    }
    
    children.push_back( newItem(j,shared()) );
    
}


void JobItem::removeChild( const Job::JobPtr& j ) {
 
    children.remove_if( [this,&j](const JobItem::Ptr& item){ return !(*(item->job) != *j); });

}



JobItem::Ptr JobItem::child( int row ) {
    
    if( row < 0 || row > (int)children.size() ) {
        return nullptr;
    }
    
    auto it = children.begin();
    while(row--) it++;
    
    return *it;
    
}


int JobItem::childCount() const {
    return children.size();
}


int JobItem::columnCount() const {
    return itemData.count();
}


int JobItem::indexOf( const Ptr& child ) {
    
    int i(-1);
    for( auto& it: children ) {
        ++i;
        if (it == child) {
            return i;
        }
    }
    
    return i;
    
}


QVariant JobItem::data( int column ) const {
    return itemData.value( column );
}


JobItem::Ptr JobItem::parent() {
    return parentItem;
}
