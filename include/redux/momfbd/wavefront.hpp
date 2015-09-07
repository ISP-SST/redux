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
        
        /*! @brief Structure representing a wavefront. Contains information about co-temporal images.
         *
         */
        struct WaveFront : public std::enable_shared_from_this<WaveFront> {
            
            struct ModeInfo {
                ModeInfo() : value(nullptr), weight(0), enabled(true) {}
                double* value;
                double weight;  //!< sum of inverse of atm_rms (= sqrt(covariance)*sqr(lambda)) for all images in this wavefront
                bool enabled;   // TODO: implement option to temporarily deactivate certain modes.
            };
            typedef std::map<uint16_t,ModeInfo> modeinfo_map;
            
            WaveFront(void);
            void init(size_t);
            void count(void) { nImages++; };
            void reset(void);
            size_t nFreeParameters(void) const { return modes.size(); };
            size_t setPointers(double* a, double* g );
            void addImage(std::shared_ptr<SubImage> im);
            void addWeight(uint16_t,double);
            void zeroAlphaWeights(void);
            void setAlpha(const modeinfo_map& a);

            void setAlphas(boost::asio::io_service&);
            void setAlphasAndUpdate(boost::asio::io_service&,bool);
            void calcGradient(boost::asio::io_service&, const grad_t&);
            double metric(boost::asio::io_service&);
            double metric(boost::asio::io_service&, const double* dAlpha);
            double coefficientMetric(void);

            modeinfo_map modes;                                     //!< The coefficients defining this wavefront. (map: mode-number -> coefficient)
            std::set<std::shared_ptr<SubImage>> images;             //!< List of images sampling this wavefront (i.e. co-temporal)
            double *alpha, *grad;                                   //!< Shortcuts to datalocations for this wavefront.
            size_t nImages;
            redux::util::Array<double> phi;                 // temporary arrays for gradient calculations
            redux::util::Array<complex_t> otf;

        };
        

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_WAVEFRONT_HPP
