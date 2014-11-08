#ifndef REDUX_GUI_HOSTITEM_HPP
#define REDUX_GUI_HOSTITEM_HPP

#include "redux/network/host.hpp"

#include <QList>
#include <QVariant>

using redux::network::Host;

namespace redux {

    namespace gui {

        class HostItem : public std::enable_shared_from_this<HostItem>  {
            
        public:
            
            typedef std::shared_ptr<HostItem> Ptr;

            static Ptr newItem( const Host::Ptr& h, Ptr parent = 0, int hi = -1 ) {
                return Ptr( new HostItem( h, parent, hi ) );
            }

            HostItem( const QList<QVariant> &data);
            ~HostItem();

            void reset( Host::Ptr h=0 );
            void appendChild( const Host::Ptr& );
            void removeChild( const Host::Ptr& );
            
            int getHostIndex( void ) { return hostIndex; };
            Ptr child( int row );
            int childCount() const;
            int columnCount() const;
            int indexOf( const Ptr& );
            QVariant data( int column ) const;

            Ptr parent();
            void clear( void ) { children.clear(); };
            Ptr shared(void) { return shared_from_this(); };

        private:
            
            HostItem( const Host::Ptr&, Ptr parent = 0, int hi = -1 );
            
            Host::Ptr host;
            std::list<Ptr> children;
            QList<QVariant> itemData;
            Ptr parentItem;
            int hostIndex;
            
        };

    }   // gui

}   // redux

#endif // REDUX_GUI_HOSTITEM_HPP
