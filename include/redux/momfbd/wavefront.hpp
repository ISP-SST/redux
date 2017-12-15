#ifndef REDUX_MOMFBD_WAVEFRONT_HPP
#define REDUX_MOMFBD_WAVEFRONT_HPP

#include "redux/momfbd/data.hpp"
#include "redux/momfbd/subimage.hpp"

#include <memory>


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
        
        class MomfbdJob;
        
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
        
        /*! @brief Class containing the collection of wavefronts. Basically just a container for alphas.
         * 
         */
        class WaveFronts {

        public:

            explicit WaveFronts( MomfbdJob& );
        
            MomfbdJob& myJob;
            logging::Logger& logger;
            
            void maybeInitializeStorage( void );          
            void getStorage( PatchData& );          
            void loadInit( boost::asio::io_service& service, redux::util::Array<PatchData::Ptr>& patches );
            void setLogChannel( std::string channel ) { logChannel = channel; };

            size_t nModes;
            size_t nWaveFronts;
            redux::util::Array<float> coefficients;
            std::string cacheFile;
            std::string logChannel;
            std::mutex mtx;

        };



        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_WAVEFRONT_HPP
