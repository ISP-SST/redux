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
            Region roi;
            
            //RegionD coordinates;
            //uint32_t patchIndexX, patchIndexY,firstPixelX,lastPixelX,firstPixelY,lastPixelY;
            //double beginX, endX, beginY, endY;
            //size_t sortedID;
            //redux::util::Array<int64_t> result;      // use int64_t as temporary storage, cast to int16_t in post-processing
            Patch() : index(), roi() {};
            Patch(uint32_t y, uint32_t x, uint32_t sz);
            void setIndex(uint32_t yid, uint32_t xid);
            size_t size( void ) const;
            char* pack( char* ) const;
            const char* unpack( const char*, bool );
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_PATCH_HPP
