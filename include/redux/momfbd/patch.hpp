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

        /*! Datastructure to contain the data/info sent from the master to the slaves.
         *  Note: the configuration/settings are mostly sent via the MomfbdJob class, this structure
         *  just holds the patch-specific information and the input data.
         */
        struct PatchData : public Part {
            typedef std::shared_ptr<PatchData> Ptr;
            Point index;                        //! patch-index in full-image mozaic
            Point first,last;                   //! local offsets, in the received data-block, of the subimage to process (x,y)
            PointF residualOffset;              //! remaining part of the tilts after pixel alignment
            redux::util::Array<float> images;   //! Stack of images ([nTotalImages][yPixels][xPixels])
            //Patch(int y, int x, uint32_t sz=1);
            void setIndex(uint32_t yid, uint32_t xid);
            size_t nPixels(void);
            size_t nPixelsX(void);
            size_t nPixelsY(void);
            uint64_t size( void ) const;
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );
            bool operator==(const PatchData&);
        };

        /*! Datastructure to contain the result/info sent back from the slaves after processing.
         */
        struct PatchResult : public Part {
            PatchResult(const PatchData&);
            redux::util::Array<float> restoredObjects, modeCoefficients, PSFs;
            uint64_t size( void ) const;
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );
 
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_PATCH_HPP
