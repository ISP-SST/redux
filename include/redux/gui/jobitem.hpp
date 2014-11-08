#ifndef REDUX_GUI_JOBITEM_HPP
#define REDUX_GUI_JOBITEM_HPP

#include "redux/job.hpp"

#include <QList>
#include <QVariant>

using redux::Job;

namespace redux {

    namespace gui {

        class JobItem : public std::enable_shared_from_this<JobItem> {
                    
        public:
            typedef std::shared_ptr<JobItem> Ptr;

            static Ptr newItem( const Job::JobPtr& j, Ptr parent = 0, int ji = -1 ) {
                return Ptr( new JobItem( j, parent, ji ) );
            }
            
            JobItem( const QList<QVariant> &data );
            ~JobItem();

            void reset( Job::JobPtr j=0 );
            void appendChild( const Job::JobPtr& );
            void removeChild( const Job::JobPtr& );

            int getJobIndex( void ) { return jobIndex; };
            Ptr child( int row );
            int childCount() const;
            int columnCount() const;
            int indexOf( const Ptr& );
            QVariant data( int column ) const;

            Ptr parent();
            void clear( void ) { children.clear(); };
            Ptr shared(void) { return shared_from_this(); };

        private:
            
            JobItem( const Job::JobPtr& j, Ptr parent = 0, int ji = -1 );
            
            Job::JobPtr job;
            QList<Ptr> childItems;
            std::list<Ptr> children;
            QList<QVariant> itemData;
            Ptr parentItem;
            int jobIndex;
        };

    }   // gui

}   // redux

#endif // REDUX_GUI_JOBITEM_HPP
