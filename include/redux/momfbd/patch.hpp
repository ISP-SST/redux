#ifndef REDUX_MOMFBD_PATCH_HPP
#define REDUX_MOMFBD_PATCH_HPP

#include "redux/types.hpp"
#include "redux/work.hpp"
#include "redux/util/array.hpp"

#include <memory>

namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */

        struct Patch : public Part {
            typedef std::shared_ptr<Patch> Ptr;
            Point index;
            Point first,last;
            PointF residualTilts;
            std::shared_ptr<char> data;
            size_t dataSize;
            //Patch(int y, int x, uint32_t sz=1);
            void setIndex(uint32_t yid, uint32_t xid);
            size_t nPixels(void);
            size_t size( void ) const;
            char* pack( char* ) const;
            const char* unpack( const char*, bool );
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_PATCH_HPP
