#ifndef REDUX_MOMFBD_WAVEFRONT_HPP
#define REDUX_MOMFBD_WAVEFRONT_HPP

#include "redux/momfbd/subimage.hpp"

#include <memory>
#include <mutex>
#include <set>


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */

        struct SubImage;

        struct WaveFront : public std::enable_shared_from_this<WaveFront> {
            
            struct ModeInfo {
                ModeInfo() : value(nullptr), weight(0), enabled(true) {}
                double* value;
                double weight;  //!< sum of inverse of atm_rms (= sqrt(covariance)*sqr(lambda)) for all images in this wavefront
                bool enabled;   // TODO: implement option to temporarily deactivate certain modes.
            };
            typedef std::map<uint16_t,ModeInfo> modeinfo_map;
            
            void clear(void);
            size_t nFreeParameters(void) const { return modes.size(); };
            void addImage(std::shared_ptr<SubImage> im);
            void addWeight(uint16_t,double);
            void zeroAlphaWeights(void);
            void computePhases(boost::asio::io_service&);
            void computePhasesAndOTF(boost::asio::io_service&);
            void setAlpha(const modeinfo_map& a);
            double coefficientMetric(void);
            void gradientTest(void);
            modeinfo_map modes;                                     //!< The coefficients defining this wavefront. (map: mode-number -> coefficient)
            std::set<std::shared_ptr<SubImage>> images;             //!< List of images sampling this wavefront (i.e. co-temporal)
            double *alpha, *grad;                                   // shortcuts to datalocations for this wavefront.
            size_t nImages;

        };
        

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_WAVEFRONT_HPP
