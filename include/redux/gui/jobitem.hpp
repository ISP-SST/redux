#ifndef REDUX_GUI_JOBITEM_HPP
#define REDUX_GUI_JOBITEM_HPP

#include <QList>
#include <QVariant>

namespace redux {

    namespace gui {

        class JobItem {
        public:
            JobItem( const QList<QVariant> &data, JobItem *parent = 0, int ji = -1 );
            ~JobItem();

            void appendChild( JobItem *child );

            int getJobIndex( void ) { return jobIndex; };
            JobItem *child( int row );
            int childCount() const;
            int columnCount() const;
            QVariant data( int column ) const;
            int row() const;
            JobItem *parent();
            void clear( void ) { childItems.clear(); qDeleteAll( childItems ); };

        private:
            QList<JobItem*> childItems;
            QList<QVariant> itemData;
            JobItem *parentItem;
            int jobIndex;
        };

    }   // gui

}   // redux

#endif // REDUX_GUI_JOBITEM_HPP
