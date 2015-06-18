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
            
            SubImage(Object&, const Channel&, const redux::util::Array<double>& wind, const redux::util::Array<float>& stack,
                     uint32_t index, uint16_t firstY, uint16_t firstX, uint16_t patchSize, uint16_t pupilSize);
            ~SubImage(void);
            
            void addPQ(redux::util::Array<complex_t>&,redux::util::Array<double>&) const;
            void init(void);
            
            void addFT(redux::util::Array<double>& ftsum) const;
            void computePhases(void);
            void oldGradientDiff(std::vector<double>&);
            void calcOTF(void);
            void calcOTF(const std::map<uint16_t,std::pair<double,bool>>&);
            redux::util::Array<double> PSF(void);
            void clearModes(redux::util::Array<double>&p) const;
            void clearModes(void) { clearModes(phi); };
            
            void setWaveFront( std::shared_ptr<WaveFront> w ) { wf = w; };
            void clear(void);
            
            Point offset;
            // alpha (sptr)
            // OTF
            // imgFT
            uint32_t index;
            uint16_t imgSize, pupilSize, otfSize;
            
            Object& object;
            const Channel& channel;
            const redux::util::Array<double>& window;
            std::shared_ptr<WaveFront> wf;
            
            redux::util::Array<double> img;         //!< Working copy of the current subimage. (apodized)     size = patchSize
            redux::util::Array<double> phi;         //!< Array containing the phase of this OTF               size = pupilsize
            redux::image::FourierTransform PF;      //!< Pupil Function = pupilmask * exp(i*phi).             size = pupilsize
            redux::image::FourierTransform SJ;      //!< Optical Transfer Function = autocorrelation of PTF.  size = 2*pupilsize
            redux::image::FourierTransform ft;      //!< Fourier transform of img.                            size = 2*pupilsize
            redux::util::Array<double> vogel;       //!< used for the Vogel-method of gradient computaion     size = pupilsize
            redux::image::Statistics stats;
            

        };


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_SUBIMAGE_HPP
