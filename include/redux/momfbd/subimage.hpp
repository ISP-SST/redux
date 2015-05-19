#ifndef REDUX_MOMFBD_SUBIMAGE_HPP
#define REDUX_MOMFBD_SUBIMAGE_HPP

#include "redux/momfbd/cache.hpp"

#include "redux/types.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/statistics.hpp"
#include "redux/util/array.hpp"
#include "redux/work.hpp"

#include <memory>
#include <mutex>


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
        class Object;
        class Channel;
        class WaveFront;

        struct SubImage : public redux::util::Array<float>, public std::enable_shared_from_this<SubImage> {
            typedef std::shared_ptr<SubImage> Ptr;
            
            //SubImage(const redux::util::Array<float>& stack, uint32_t index, int firstY, int lastY, int firstX, int lastX);
            SubImage(const Object&, const Channel&, const redux::util::Array<float>& stack, uint32_t index, int firstY, int lastY, int firstX, int lastX);
            ~SubImage(void);
            
            void init(const redux::util::Array<double>&);
            //void addFT(redux::image::FourierTransform& ftsum);
            void addFT(redux::util::Array<double>& ftsum);
            void addPQ(redux::util::Array<complex_t>&,redux::util::Array<double>&) const;
            void computePhases(void);
            void oldGradientDiff(std::vector<double>&);
            void OTF(void);
            void OTF(const std::map<uint16_t,std::pair<double,bool>>&);
            redux::util::Array<double> PSF(void);
            void setWaveFront( std::shared_ptr<WaveFront> w ) { wf = w; };
            void clear(void);
            
            Point offset;
            // alpha (sptr)
            // OTF
            // imgFT
            uint32_t index;
            
            const Object& object;
            const Channel& channel;
            std::shared_ptr<WaveFront> wf;
            
            redux::util::Array<double> img;         //!< Working copy of the current subimage. Apodized & local
            redux::util::Array<double> phi;         //!< Array containing the phase of this OTF
            redux::image::FourierTransform SJ;      //!< OTF (Optical Transfer Function = FFT of the PSF)
            redux::image::FourierTransform ft;
            redux::image::Statistics stats;
            
        };


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_SUBIMAGE_HPP
