#ifndef REDUX_MOMFBD_WAVEFRONT_HPP
#define REDUX_MOMFBD_WAVEFRONT_HPP

#include "redux/momfbd/subimage.hpp"

#include <memory>


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
        
        /*! @brief Structure representing a wavefront. I.e. a parametrization of the atmospheric seeing sampled by
         *         one or more (co-temporal) images.
         */
        struct WaveFront {
            
            explicit WaveFront( uint32_t );
            void setData( double* a, double* g, size_t n );
            void clear(void);
            inline size_t size(void) { return images.size(); };
            std::string print(void);
            void addImage( const std::shared_ptr<SubImage>& im );
            template <typename T> void adjustShifts( const T* a ) const;
            void adjustShifts( void ) { adjustShifts( alpha ); }
            template <typename T> void applyAlpha( const T* a ) const;
            void applyAlpha( void ) { applyAlpha( alpha ); }

            //void gradient( const bool* enabledModes, grad_tt grad );
            
            void dump( std::string tag, bool dumpImages=false ) const;

            std::vector<std::shared_ptr<SubImage>> images;              //!< List of images sampling this wavefront (i.e. co-temporal)
            
            uint32_t wfIndex;
            double *alpha, *grad_alpha;                                 //!< Shortcut to data for this wavefront.
            size_t nAlpha;
            
        };
        

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_WAVEFRONT_HPP
