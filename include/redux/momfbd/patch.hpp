#ifndef REDUX_MOMFBD_PATCH_HPP
#define REDUX_MOMFBD_PATCH_HPP

#include "redux/types.hpp"
#include "redux/work.hpp"
#include "redux/util/array.hpp"
#include "redux/momfbd/modecache.hpp"

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
            Point index;                        //! patch-index (x,y)
            Point first,last;                   //! local offsets, in the received data-block, of the subimage to process (x,y)
            std::vector<PointF> residualOffsets;//! remaining part of the tilts after pixel alignment
            redux::util::Array<float> images;   //! Stack of images ([nTotalImages][yPixels][xPixels])
            //Patch(int y, int x, uint32_t sz=1);
            void setIndex(uint32_t yid, uint32_t xid);
            size_t nPixels(void);
            size_t nPixelsX(void);
            size_t nPixelsY(void);
            size_t size( void ) const;
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );
        };

        /*! Datastructure to contain the result/info sent back from the slaves after processing.
         */
        struct PatchResult : public Part {
            redux::util::Array<float> restoredObjects, modeCoefficients, PSFs;
            size_t size( void ) const;
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );
 
        };

        /*! Structure to group together the subimage, FT
         */
        struct Unit {
            redux::util::Array<double> img;                    //! Reference to the image.
            redux::util::Array<complex_t> ft;                  //! Reference to the Fourier-transform
        };

        /*! Structure tying together images sampling the same wavefront abberrations.
         */
        struct Simultaneous {
            // alpha -> OTF -> metric
            std::vector< redux::util::Array<double> > images;
            std::vector<double> alpha;
        };

        /*! Container used during processing. Basically temporary arrays and reorganized references to original data.
         */
        struct WorkSpace {
            
            struct eachImage {
                // alpha (sptr)
                // OTF
                // imgFT
            };
            
            struct ConstrainedGroup {
                // std::vector<eachImage>
                // alpha/OTF shared within group (sptr)
                // imgFT
            };
            
            
            struct channelSpecific {
                // phi_fixed
            };
            
            struct objectSpecific {
                // modes
                // approxObj
                objectSpecific();
                std::map<int, const ModeCache::ModePtr> modes;
                std::vector<channelSpecific> chan;
            };
            std::vector<objectSpecific> obj;
            
            WorkSpace(PatchData::Ptr);
            void init(void);
            
            PatchData::Ptr data;
            redux::util::Array<ModeCache::ModePtr> modes;
            
            // std::vector<ConstrainedGroup> groups;            // [nTotalImages - nConstraints] (some singleton groups)
           
            redux::util::Array<complex_t> imgFTs;               //! Stack of Fourier-transforms of input-images ([nTotalImages][pupilSize][pupilSize])
            redux::util::Array<complex_t> objFT;                //! Fourier-transforms of the approximate object ([pupilSize][pupilSize])
        };

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_PATCH_HPP
