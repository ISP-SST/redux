#include "redux/gui/jobmodel.hpp"

#include "redux/gui/jobitem.hpp"

#include <QtGui>

using namespace redux::gui;
using namespace redux;

JobModel::JobModel( const Job::JobSet& j, QObject *parent ) : QAbstractItemModel( parent ), jobs(j) {
    QList<QVariant> rootData;
    rootData << "Job" << "ID" << "Label" << "Type" << "Priority" << "User" << "Status";
    rootItem = new JobItem( rootData );
    //setupModelData(data.split(QString("\n")), rootItem);
}

JobModel::~JobModel() {
    delete rootItem;
}

void JobModel::setData( const QString &data ) {
    //QList<QVariant> rootData;
    //rootData << "Name" << "Type" << "Status";
    //rootItem = new JobItem(rootData);
    setupModelData( data.split( QString( "\n" ) ), rootItem );
}

int JobModel::columnCount( const QModelIndex &parent ) const {
    if( parent.isValid() )
        return static_cast<JobItem*>( parent.internalPointer() )->columnCount();
    else
        return rootItem->columnCount();
}

QVariant JobModel::data( const QModelIndex &index, int role ) const {
    if( !index.isValid() )
        return QVariant();

    if( role != Qt::DisplayRole )
        return QVariant();

    JobItem *item = static_cast<JobItem*>( index.internalPointer() );

    return item->data( index.column() );
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

QModelIndex JobModel::index( int row, int column, const QModelIndex &parent )
const {
    if( !hasIndex( row, column, parent ) )
        return QModelIndex();

    JobItem *parentItem;

    if( !parent.isValid() )
        parentItem = rootItem;
    else
        parentItem = static_cast<JobItem*>( parent.internalPointer() );

    JobItem *childItem = parentItem->child( row );
    if( childItem )
        return createIndex( row, column, childItem );
    else
        return QModelIndex();
}

QModelIndex JobModel::parent( const QModelIndex &index ) const {
    if( !index.isValid() )
        return QModelIndex();

    JobItem *childItem = static_cast<JobItem*>( index.internalPointer() );
    JobItem *parentItem = childItem->parent();

    if( parentItem == rootItem )
        return QModelIndex();

    return createIndex( parentItem->row(), 0, parentItem );
}

int JobModel::rowCount( const QModelIndex &parent ) const {
    JobItem *parentItem;
    if( parent.column() > 0 )
        return 0;

    if( !parent.isValid() )
        parentItem = rootItem;
    else
        parentItem = static_cast<JobItem*>( parent.internalPointer() );

    return parentItem->childCount();
}

void JobModel::setupModelData( const QStringList &lines, JobItem *parent ) {
    QList<JobItem*> parents;
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
            parents.last()->appendChild( new JobItem( columnData, parents.last() ) );
        }

        number++;
    }
}

void JobModel::refreshTree( void ) {

    QList<JobItem*> parents;
    QList<int> indentations;
    parents << rootItem;
    indentations << 0;

    rootItem->clear();

    if( jobs.empty() ) return;  //|| (*rootObject).size() <1 ) return;

    QList<QVariant> rowData;

    //rootData << "#" << "ID" << "Label" << "Type" << "User" << "Status";
    int i = 1;  // index
    for( auto &it: jobs ) {
        rowData.clear();
        rowData << std::to_string( i++ ).c_str();           // rownumber+1 = jobnumber
        rowData << std::to_string( it->info.id ).c_str(); // job ID
        rowData << it->info.name.c_str();             // job-label
        rowData << it->info.typeString.c_str();        // job-type
        rowData << std::to_string( it->info.priority ).c_str(); // priority
        rowData << it->info.user.c_str();              // username
        rowData << it->info.user.c_str();     // current state
        rootItem->appendChild( new JobItem( rowData, rootItem, i ) );
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
                parents.last()->appendChild(new JobItem(columnData, parents.last()));
            }

            number++;
        }*/
}
