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
            
            struct alpha_t {
                alpha_t() : value(0), weight(0), enabled(true) {}
                double metric(void) { std::cout << "alpha_t::metric("<<redux::util::hexString(this)<<")  v=" <<value << "  w="<< weight << std::endl; return value*value*weight; }
                double value;
                double weight;  //!< sum of inverse of atm_rms (= sqrt(covariance)*sqr(lambda)) for all images in this wavefront
                bool enabled;   // TODO implement deactivating certain modes.
            };
            typedef std::map<uint16_t,alpha_t> alpha_map;
            
            void clear(void);
            void addImage(std::shared_ptr<SubImage> im);
            void addWeight(uint16_t,double);
            void zeroWeights(void);
            void setAlpha(const alpha_map& a);
            void zeroValues(void);
            void computePhases(boost::asio::io_service&);
            void computePhasesAndOTF(boost::asio::io_service&);
            double coefficientMetric(void);
            void gradientTest(void);
            alpha_map alpha;
            std::set<std::shared_ptr<SubImage>> images;             // list of images sampling this wavefront

        };
        

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_WAVEFRONT_HPP
