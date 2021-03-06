#ifndef REDUX_GUI_JOBMODEL_HPP
#define REDUX_GUI_JOBMODEL_HPP

#include "redux/gui/jobitem.hpp"

#include <vector>

#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>


namespace redux {

    namespace gui {

        class JobModel : public QAbstractItemModel {
            Q_OBJECT

        public:
            JobModel( const Job::JobSet&, QObject *parent = 0 );
            ~JobModel();


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
            const Job::JobSet& jobs;
            JobItem::Ptr rootItem;
            
        };

    }   // gui

}   // redux

#endif // REDUX_GUI_JOBMODEL_HPP
