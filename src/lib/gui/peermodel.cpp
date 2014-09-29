#include "redux/gui/peermodel.hpp"

#include "redux/gui/peeritem.hpp"

#include <QtGui>
/*
#include <DataTools.h>
#include <StringTools.h>*/

#include <boost/date_time/posix_time/time_formatters.hpp>
#include <boost/format.hpp>

using namespace redux::gui;
using namespace redux;
using namespace std;

PeerModel::PeerModel( const std::set<Peer>& p, QObject *parent ) : QAbstractItemModel( parent ), peers( p ) {
    QList<QVariant> rootData;
    rootData << "Hostname" << "PID" << "T/C" << "Load" << "Uptime" << "Last Seen" << "";
    rootItem = new PeerItem( rootData );
    //setupModelData(data.split(QString("\n")), rootItem);
}

PeerModel::~PeerModel() {
    delete rootItem;
}

void PeerModel::setData( const QString &data ) {
    //QList<QVariant> rootData;
    //rootData << "Name" << "Type" << "Status";
    //rootItem = new PeerItem(rootData);
    setupModelData( data.split( QString( "\n" ) ), rootItem );
}

int PeerModel::columnCount( const QModelIndex &parent ) const {
    if( parent.isValid() )
        return static_cast<PeerItem*>( parent.internalPointer() )->columnCount();
    else
        return rootItem->columnCount();
}

QVariant PeerModel::data( const QModelIndex &index, int role ) const {
    if( !index.isValid() )
        return QVariant();

    if( role != Qt::DisplayRole )
        return QVariant();

    PeerItem *item = static_cast<PeerItem*>( index.internalPointer() );

    return item->data( index.column() );
}

Qt::ItemFlags PeerModel::flags( const QModelIndex &index ) const {
    if( !index.isValid() )
        return 0;

    return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}

QVariant PeerModel::headerData( int section, Qt::Orientation orientation, int role ) const {
    if( orientation == Qt::Horizontal && role == Qt::DisplayRole )
        return rootItem->data( section );

    return QVariant();
}

QModelIndex PeerModel::index( int row, int column, const QModelIndex &parent ) const {
    if( !hasIndex( row, column, parent ) )
        return QModelIndex();

    PeerItem *parentItem;

    if( !parent.isValid() )
        parentItem = rootItem;
    else
        parentItem = static_cast<PeerItem*>( parent.internalPointer() );

    PeerItem *childItem = parentItem->child( row );
    if( childItem )
        return createIndex( row, column, childItem );
    else
        return QModelIndex();
}

QModelIndex PeerModel::parent( const QModelIndex &index ) const {
    if( !index.isValid() )
        return QModelIndex();

    PeerItem *childItem = static_cast<PeerItem*>( index.internalPointer() );
    PeerItem *parentItem = childItem->parent();

    if( parentItem == rootItem )
        return QModelIndex();

    return createIndex( parentItem->row(), 0, parentItem );
}

int PeerModel::rowCount( const QModelIndex &parent ) const {

    PeerItem *parentItem;
    if( parent.column() > 0 )
        return 0;

    if( !parent.isValid() )
        parentItem = rootItem;
    else
        parentItem = static_cast<PeerItem*>( parent.internalPointer() );

    return parentItem->childCount();
}

void PeerModel::setupModelData( const QStringList &lines, PeerItem *parent ) {

    QList<PeerItem*> parents;
    QList<int> indentations;
    parents << parent;
    indentations << 0;

    int number = 0;

    while( number < lines.count() ) {
        int position = 0;
        while( position < lines[number].length() ) {
            if( lines[number].mid( position, 1 ) != " " )
                break;
            position++;
        }

        QString lineData = lines[number].mid( position ).trimmed();

        if( !lineData.isEmpty() ) {
            // Read the column data from the rest of the line.
            QStringList columnStrings = lineData.split( "\t", QString::SkipEmptyParts );
            QList<QVariant> columnData;
            for( int column = 0; column < columnStrings.count(); ++column )
                columnData << columnStrings[column];

            if( position > indentations.last() ) {
                // The last child of the current parent is now the new parent
                // unless the current parent has no children.

                if( parents.last()->childCount() > 0 ) {
                    parents << parents.last()->child( parents.last()->childCount() - 1 );
                    indentations << position;
                }
            }
            else {
                while( position < indentations.last() && parents.count() > 0 ) {
                    parents.pop_back();
                    indentations.pop_back();
                }
            }

            // Append a new item to the current parent's list of children.
            parents.last()->appendChild( new PeerItem( columnData, parents.last() ) );
        }

        number++;
    }
}

void PeerModel::refreshTree( void ) {
    
    using namespace boost::posix_time;

    QList<PeerItem*> parents;
    QList<int> indentations;
    parents << rootItem;
    indentations << 0;

    rootItem->clear();

    if( peers.empty() ) return;

    QList<QVariant> columnData;
    /*columnData << ( ROYAC::toString ( (*rootObject)[0].hostName ) + string ( ":" ) +
    ROYAC::toString ( (int)(*rootObject)[0].listeningPort ) ).c_str();
     columnData << ROYAC::toString ( (int)(*rootObject)[0].nThreads[0] ).c_str();
    columnData << "OFFLINE";
    PeerItem *masterItem = new PeerItem ( columnData, rootItem, 0 );
    rootItem->appendChild ( masterItem );
    */
    //cout << "muuu" << rootObject->size() << endl;

    string tmp;
    struct timespec ts;
    clock_gettime( CLOCK_REALTIME, &ts );
    std::stringstream ssElapsed;
    //ssDuration << duration;

    //std::string str = ssDuration.str();
    
    //gettimeofday ( &tv, NULL );
    int i=0;
    for( auto &it: peers ) {
        ssElapsed.clear();
//         //cout << "muuu            " << i << " - " << hexString((*rootObject)[i]) << endl;;
        columnData.clear();
        tmp = it.host.name; // + string( " (" ) + ipString( ( *rootObject )[i]->IP ) + string( ":" );
        //tmp += toString( ( int )( *rootObject )[i]->listeningPort ) + string( ")" );
        //cout << "yyyyyyy         " << i << endl;
        columnData << tmp.c_str();
        columnData << to_string( it.host.pid ).c_str();
        tmp = to_string( it.stat.nThreads ) + string( "/" ) + to_string( it.host.nCores ); // + string( "/" )
           //   + to_string( ( int )( *rootObject )[i]->nThreads[3] ) + string( "/" ) + to_string( ( int )( *rootObject )[i]->nThreads[0] );
        //cout << "xxxxxxx         " << i << endl;
        columnData << tmp.c_str();
        columnData << boost::str( boost::format( "%3.1f%%" ) % it.stat.loadAvg ).c_str(); //to_string( it.stat.loadAvg ).c_str();
        //cout << "wwwwwww         " << i << endl;
        // job->info.submitTime = boost::posix_time::second_clock::local_time();
        time_duration elapsed = (second_clock::local_time()-it.host.startedAt);
        //to_simple_string(elapsed)
        //ssElapsed << elapsed;
        columnData << to_simple_string(elapsed).c_str();//to_iso_extended_string( it.host.startedAt ).c_str(); //tsToString( tsSubtract( ts, ( *rootObject )[i]->startTime ) ).c_str();
        
        //string tmp;
        //tmp.resize ( 25, 0 );
        //strftime ( & (tmp[0]), 24, "%Y-%m-%d %H:%M:%S.\0", localtime ( &msg->timestamp.tv_sec ) );
        elapsed = (second_clock::local_time()-it.lastSeen);
        columnData << to_simple_string(elapsed).c_str();//to_iso_extended_string( it.lastSeen ).c_str(); //tsToString( tsSubtract( ts, ( *rootObject )[i]->lastSeen ) ).c_str();
        //cout << "rrrrrrr         " << i << endl;
        columnData << Peer::TypeNames[it.host.peerType].c_str(); //  ( *rootObject )[i]->stateString().c_str();
        //cout << "nuuuu         " << i << endl;
        //masterItem->appendChild ( new PeerItem ( columnData, masterItem, i ) );
        rootItem->appendChild( new PeerItem( columnData, rootItem, i ) );
        //cout << "gnuuuu      " << i << endl;
    }


    /*    while (number < lines.count()) {
            int position = 0;
            while (position < lines[number].length()) {
                if (lines[number].mid(position, 1) != " ")
                    break;
                position++;
            }

            QString lineData = lines[number].mid(position).trimmed();

            if (!lineData.isEmpty()) {
                // Read the column data from the rest of the line.
                QStringList columnStrings = lineData.split("\t", QString::SkipEmptyParts);
                QList<QVariant> columnData;
                for (int column = 0; column < columnStrings.count(); ++column)
                    columnData << columnStrings[column];

                if (position > indentations.last()) {
                    // The last child of the current parent is now the new parent
                    // unless the current parent has no children.

                    if (parents.last()->childCount() > 0) {
                        parents << parents.last()->child(parents.last()->childCount()-1);
                        indentations << position;
                    }
                } else {
                    while (position < indentations.last() && parents.count() > 0) {
                        parents.pop_back();
                        indentations.pop_back();
                    }
                }

                // Append a new item to the current parent's list of children.
                parents.last()->appendChild(new PeerItem(columnData, parents.last()));
            }

            number++;
        }*/
}

