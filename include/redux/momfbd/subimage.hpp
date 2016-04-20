#ifndef REDUX_MOMFBD_SUBIMAGE_HPP
#define REDUX_MOMFBD_SUBIMAGE_HPP

#include "redux/momfbd/modes.hpp"

#include "redux/types.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/util/arraystats.hpp"
#include "redux/util/array.hpp"
#include "redux/work.hpp"

#include <functional>
#include <memory>
#include <map>
#include <mutex>


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
        class Object;
        class Channel;

        struct SubImage : public redux::util::Array<float> {
            typedef std::shared_ptr<SubImage> Ptr;
            
            SubImage(Object&, const Channel&, const redux::util::Array<double>& wind, const redux::util::Array<double>& nwind);
            /*SubImage(Object&, const Channel&, const redux::util::Array<double>& wind, const redux::util::Array<double>& nwind,
                     const redux::util::Array<float>& stack,
                     uint32_t index, const PointI& offset, uint16_t patchSize, uint16_t pupilSize);*/
            ~SubImage(void);
            
            void setPatchInfo(uint32_t, const PointI&, const PointI&, uint16_t, uint16_t, uint16_t);
            void init(void);
            void newCutout(void);
            
            void addFT(redux::util::Array<double>& ftsum) const;
            void addPQ(complex_t* P, double* Q) const { addPQ(OTF.get(),P,Q); };
            void addPQ(const complex_t* otf, complex_t* P, double* Q) const;
            void addToPQ(void) const;
            void restore(complex_t* avg_obj, double* norm) const;
            
            double metricChange(const complex_t* newOTF) const;
            double gradientFiniteDifference(uint16_t, double);
            double gradientVogel(uint16_t mode, double) const;
            void calcVogelWeight(void);
            
            //void addScaledMode( double* modePtr, double a );
            
            void addMode(double* phiPtr, const double* modePtr, double a) const;
            //void addModes(size_t nModes, uint16_t* modes, const double* a) { addModes(phi.get(), nModes, modes, a); }
            //void addModes(double* phiPtr, size_t nModes, uint16_t* modes, const double* a) const;
            
            void adjustOffset(void);

            void addPhases(const double* a) { addPhases(phi.get(), a); };
            void addPhases(double* phiPtr, const double* a) const;
            
            void addAlpha(uint16_t m, double a);
            void setAlpha(uint16_t m, double a);
            void addAlphas(const double* a);
            void setAlphas(const double* a);
            void setAlphas(const std::vector<uint16_t>& modes, const double* a);

            void getAlphas(float* alphas) const;    // copy with offset correction (results)
            void getAlphas(double* alphas) const {  // copy without offset correction
                std::copy(currentAlpha.begin(),currentAlpha.end(),alphas);
            }
            
            void resetPhi(void);
            void addPhi(const double* p, double scale=1.0);
            void calcPhi(const double* a);

            void calcOTF(complex_t* otf, const double* phiOffset, double scale);
            void calcOTF(complex_t* otf, const double* phi);
            void calcOTF(void);
            void calcPFOTF(void);
            
            template <typename T>
            void addPSF(redux::util::Array<T>& out) const {
                redux::util::Array<T> tmp;
                OTF.inv(tmp,redux::image::FT_REORDER);
                out += tmp;
            }

            template <typename T=float>
            redux::util::Array<T> getPSF (void) const {
                using namespace redux::image;
                redux::util::Array<T> tmp;
                OTF.inv(tmp,FT_REORDER|FT_NORMALIZE);
                return std::move(tmp);
            }

            
            template <typename T>
            redux::util::Array<T> convolveImage(const redux::util::Array<T>& im) {
                using namespace redux::image;
                FourierTransform imFT(im, FT_REORDER|FT_FULLCOMPLEX|FT_NORMALIZE);
                imFT *= OTF;
                imFT.reorder();
                imFT.directInverse(tmpOTF);
                FourierTransform::reorder(tmpOTF);
                return std::move(tmpOTF.copy<T>());
            }
            
            template <typename T>
            redux::util::Array<T> residual(const redux::util::Array<T>& im) {
                using namespace redux::image;
                FourierTransform imFT(im, FT_REORDER|FT_FULLCOMPLEX|FT_NORMALIZE);
                imFT *= OTF;
                imFT.reorder();
                imFT.directInverse(tmpOTF);
                FourierTransform::reorder(tmpOTF);
                redux::util::Array<T> tmp = tmpOTF.copy<T>();
                tmp -= img;
                return std::move(tmp);
            }
            template <typename T>
            redux::util::Array<T> convolvedResidual(const redux::util::Array<T>& cim) { return std::move(cim-img); }
            
            //void update(bool newVogel=false);
            void dump( std::string tag ) const;

            uint32_t index;
            PointI offset;                                      //<! Location of the current/original cutout, this typically starts at (maxLocalShift,maxLocalShift)
            PointI offsetShift;                                 //<! How the subimage has been shifted to compensate for large tip/tilt coefficients.
            uint16_t imgSize, pupilSize, nModes;
            uint32_t otfSize, pupilSize2, otfSize2;
            double oldRG;
            
            Object& object;
            const Channel& channel;
            const ModeSet& modes;
            const redux::util::Array<double>& window;
            const redux::util::Array<double>& noiseWindow;
            
            std::map<uint16_t, double> alpha;
            double** alphaRef;
            std::vector<double> currentAlpha;
            bool newPhi;
            bool newOTF;
            std::mutex mtx;
            
            redux::util::Array<complex_t> tmpC,tmpC2;   //!< Temporary arrays.                                     size = otfSize
            redux::util::Array<double> img,tmpImg;      //!< Working copy of the current subimage. (apodized)     size = patchSize
            redux::util::Array<double> phi,tmpPhi;      //!< Array containing the phase of this OTF               size = pupilsize
            redux::image::FourierTransform PF;          //!< Pupil Function = pupilmask * exp(i*phi).             size = pupilsize
            redux::image::FourierTransform OTF,tmpOTF;  //!< Optical Transfer Function = autocorrelation of PF.   size = 2*pupilsize
            redux::image::FourierTransform imgFT,tmpFT; //!< Fourier transform of img.                            size = 2*pupilsize
            redux::util::Array<double> vogel;           //!< used for the Vogel-method of gradient computaion     size = pupilsize
            redux::util::ArrayStats stats;
            

        };
        
        
        typedef std::function<double(SubImage&, uint16_t, double)> grad_t;


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_SUBIMAGE_HPP
