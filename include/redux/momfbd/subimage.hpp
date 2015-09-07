#ifndef REDUX_MOMFBD_SUBIMAGE_HPP
#define REDUX_MOMFBD_SUBIMAGE_HPP

#include "redux/momfbd/cache.hpp"

#include "redux/types.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/util/arraystats.hpp"
#include "redux/util/array.hpp"
#include "redux/work.hpp"

#include <functional>
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
            
            void init(void);
            
            void addFT(redux::util::Array<double>& ftsum) const;
            void addPQ(complex_t* P, double* Q) const { addPQ(OTF.get(),P,Q); };
            void addPQ(const complex_t* otf, complex_t* P, double* Q) const;
            void setAlpha(double*);
            
            double metricChange(const complex_t* newOTF) const;
            double gradientFiniteDifference(uint16_t, double, complex_t*, double*) const;
            double gradientVogel(uint16_t mode, double dalpha, complex_t* tmpOTF, double* tmpPhi) const;
            void calcVogelWeight(void);
            
            void resetPhi(redux::util::Array<double>&p) const;
            void resetPhi(void) { resetPhi(phi); };
            
            void addMode(uint16_t m, double alpha) { addMode(phi.get(), m, alpha); }
            void addMode(double* phiPtr, uint16_t m, double alpha) const;
            void addModes(size_t nModes, uint16_t* modes, const double* alpha) { addModes(phi.get(), nModes, modes, alpha); }
            void addModes(double* phiPtr, size_t nModes, uint16_t* modes, const double* alpha) const;
            
            void addPhases(const double* alpha) { addPhases(phi.get(), alpha); };
            void addPhases(double* phiPtr, const double* alpha) const;
            
            void calcOTF(void) { calcPFOTF(PF.get(), OTF.get(), phi.get()); }
            void calcOTF(redux::util::Array<complex_t>& otf) const { calcOTF(otf.get(), phi.get()); }
            void calcOTF(complex_t* otfPtr, const double* phiPtr) const;
            void calcPFOTF(redux::util::Array<complex_t>& pf, redux::util::Array<complex_t>& otf) const { calcPFOTF(pf.get(), otf.get(), phi.get()); }
            void calcPFOTF(complex_t* pfPtr, complex_t* otfPtr, const double* phiPtr) const;
            
            redux::util::Array<double> getPSF(void);
            void setWaveFront( std::shared_ptr<WaveFront> w ) { wf = w; };
            
            void dump( std::string tag );

            uint32_t index;
            uint16_t imgSize, pupilSize, otfSize;
            
            Object& object;
            const Channel& channel;
            const redux::util::Array<double>& window;
            std::shared_ptr<WaveFront> wf;
            
            redux::util::Array<double> img;             //!< Working copy of the current subimage. (apodized)     size = patchSize
            redux::util::Array<double> phi,tmpPhi;      //!< Array containing the phase of this OTF               size = pupilsize
            redux::image::FourierTransform PF;          //!< Pupil Function = pupilmask * exp(i*phi).             size = pupilsize
            redux::image::FourierTransform OTF,tmpOTF;  //!< Optical Transfer Function = autocorrelation of PF.   size = 2*pupilsize
            redux::image::FourierTransform imgFT;       //!< Fourier transform of img.                            size = 2*pupilsize
            redux::util::Array<double> vogel;           //!< used for the Vogel-method of gradient computaion     size = pupilsize
            redux::util::ArrayStats stats;
            

        };
        
        
        typedef std::function<double(SubImage&, uint16_t, double, complex_t*, double*)> grad_t;


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_SUBIMAGE_HPP
