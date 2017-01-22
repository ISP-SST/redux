#ifndef REDUX_MOMFBD_SUBIMAGE_HPP
#define REDUX_MOMFBD_SUBIMAGE_HPP

#include "redux/momfbd/modes.hpp"

#include "redux/util/point.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/util/arraystats.hpp"
#include "redux/util/array.hpp"
#include "redux/work.hpp"

#include <functional>
#include <memory>
#include <map>
#include <mutex>


namespace redux {
    
    namespace logging {
        class Logger;
    }


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
            
            void setPatchInfo(uint32_t, const redux::util::PointI&, const redux::util::PointF&, uint16_t, uint16_t, uint16_t);
            void getWindowedImg( Array<double>&, redux::util::ArrayStats& s, bool rescaled ) const;
            void initialize( bool doReset=false );
            
            void addFT(redux::util::Array<double>& ftsum) const;
            void addPQ(complex_t* P, double* Q) const { addPQ(OTF.get(),P,Q); };
            void addPQ(const complex_t* otf, complex_t* P, double* Q) const;
            //void addToPQ(void) const;
            void restore(complex_t* avg_obj, double* norm) const;
            
            double metricChange(const complex_t* newOTF) const;
            double gradientFiniteDifference(uint16_t);
            void gradientFiniteDifference2( double* agrad, const bool* enabledModes );
            double gradientVogel(uint16_t mode);
            void gradientVogel2( double* agrad, const bool* enabledModes );
            void calcVogelWeight( complex_t*, double*, double* );
            void calcVogelWeight( void );
            
            void addToPhi(double* phiPtr, const double* modePtr, double a) const;
            

            template <typename T>
            bool adjustShifts( const T* alpha );
            void alignAgainst( const Ptr& refIm );

            void addPhases(const double* a) { addPhases(phi.get(), a); };
            void addPhases(double* phiPtr, const double* a) const;
            
            void addAlpha(uint16_t m, double a);
            void setAlpha(uint16_t m, double a);
            void addAlphas(const double* a);
            void setAlphas(const double* a);
            void setAlphas(const std::vector<uint16_t>& modes, const double* a);
            
            void resetShifts( void );
            
            void resetPhi(void);
            template <typename T> void addToPhi(const T* a, double* phiPtr) const;
            template <typename T> void addToPhi(const T* a) { addToPhi( a, phi.get() ); }
            template <typename T> void calcPhi( const T* a, double* phiPtr ) const;
            template <typename T> void calcPhi( const T* a ) { calcPhi( a, phi.get() ); }

            void calcOTF(complex_t* otf, const double* phiOffset, double scale);
            void calcOTF(complex_t* otf, const double* phi) const;
            void calcOTF(void) { calcOTF( OTF.get(), phi.get() ); }
            void calcPFOTF(void);
            
            void addPSF( double* psf ) const;
            void getPSF( double* psf ) const;
            redux::util::Array<double> getPSF( void ) const;
            
            template <typename T>
            redux::util::Array<T> convolveImage( const redux::util::Array<T>& im ) const {
                using namespace redux::image;
                redux::util::Array<T> tmp = im.copy();
                OTF.convolveInPlace( tmp, FT_FULLCOMPLEX );
                return std::move( tmp );
            }
            
            template <typename T>
            redux::util::Array<T> residual( const redux::util::Array<T>& im ) const {
                redux::util::Array<T> cim = convolveImage(im);
                redux::util::Array<double> img(imgSize,imgSize);
                redux::util::ArrayStats s;
                getWindowedImg(img,s,true);
                cim -= img;
                return std::move(cim);
            }
            template <typename T>
            redux::util::Array<T> convolvedResidual( const redux::util::Array<T>& cim ) {
                redux::util::Array<double> img(imgSize,imgSize);
                getWindowedImg(img,stats,true);
                return std::move(cim-img);
            }
            
            std::string idString( void ) const;
            void dump( std::string tag ) const;

            uint32_t index;
            redux::util::PointI initialOffset;          //<! Starting location of the patch in the datablock, this is typically (maxLocalShift, maxLocalShift)
            redux::util::PointF channelResidualOffset;  //<! Remainder after shifting the cutout integer pixels.
            redux::util::PointI imageShift;            //<! How the subimage has been shifted to compensate for large tip/tilt coefficients.
            uint16_t imgSize, pupilSize, nModes;
            uint32_t otfSize, pupilSize2, otfSize2;
            double oldRG;
            double grad_step;
            
            Object& object;
            const Channel& channel;
            logging::Logger& logger;
            const ModeSet& modes;
            const redux::util::Array<double>& window;
            const redux::util::Array<double>& noiseWindow;
            
            redux::util::PointD adjustedTilts;
            bool shifted;
            bool newPhi;
            bool newOTF;
            std::mutex mtx;
            
            std::vector<double> localAlpha;
            redux::util::Array<double> phi;         //!< Array containing the phase of this OTF               size = pupilsize
            redux::image::FourierTransform PF;      //!< Pupil Function = pupilmask * exp(i*phi).             size = pupilsize
            redux::image::FourierTransform OTF;     //!< Optical Transfer Function = autocorrelation of PF.   size = 2*pupilsize
            redux::image::FourierTransform imgFT;   //!< Fourier transform of img.                            size = 2*pupilsize
            redux::util::Array<double> vogel;       //!< used for the Vogel-method of gradient computaion     size = pupilsize
            redux::util::ArrayStats stats;
            

        };
        
        
        typedef std::function<double(SubImage&, uint16_t)> grad_t;


        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_SUBIMAGE_HPP
