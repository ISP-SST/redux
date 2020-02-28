#ifndef REDUX_MOMFBD_MODES_HPP
#define REDUX_MOMFBD_MODES_HPP

#include "redux/util/point.hpp"
#include "redux/image/pupil.hpp"

#include <memory>
#include <mutex>
#include <vector>

#define MT_NONE   0
#define MT_PUPIL  1
#define MT_OTF_BF 2

namespace redux {

    namespace momfbd {
        
        struct ModeInfo {
            
            uint16_t firstMode, lastMode, modeNumber, nPupilPixels;
            std::vector<uint16_t> modeNumbers;
            double pupilRadius, angle, cutoff;
            std::string filename;
            ModeInfo ( const std::string&, uint16_t nPixels=0 );
            ModeInfo ( uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double angle );
            ModeInfo ( uint16_t firstMode, uint16_t lastMode, uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double angle, double cutoff );
            ModeInfo ( uint16_t firstMode, uint16_t lastMode, const std::vector<uint16_t>& modeNumbers, uint16_t nPoints, double pupilRadius, double angle, double cutoff );
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            bool operator<( const ModeInfo& rhs ) const;
            operator std::string() const;
            
        };



        struct PupilMode : public redux::util::Array<double> {
            
            // Karhunen-Loeve expansion coefficients
            struct KL {
                std::vector< std::pair<uint16_t, double> > zernikeWeights; //< zernike mode numbers and corresponding weights
                double covariance;
            };
            typedef std::shared_ptr<PupilMode> Ptr;
            typedef std::shared_ptr<KL> KLPtr;

            PupilMode() : atm_rms(0) {};
            PupilMode ( uint16_t modeNumber, uint16_t nPoints, double r_c = 1.0, double angle = 0.0, int flags=0 ); // Zernike
            PupilMode ( uint16_t firstMode, uint16_t lastMode, uint16_t klModeNumber, uint16_t nPoints, double r_c = 1.0, double angle = 0.0, double cutoff=0.0, int flags=0 ); // KL

            double atm_rms;                         //!< = sqrt(covariance), used in metric computations.

        };


        struct ModeSet : public redux::util::Array<double> {
            
            ModeSet();
            ModeSet(ModeSet&& rhs);
            ModeSet(const ModeSet& rhs);
            
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            ModeSet& operator=( const ModeSet& rhs );
            ModeSet clone( void ) const;
            
            bool load( const std::string& filename, uint16_t pixels );
            void generate( uint16_t pixels, double radius, double angle, const std::vector<uint16_t>& modes, int flags );  // Zernike
            void generate( uint16_t pixels, double radius, double angle, uint16_t firstMode, uint16_t lastMode, const std::vector<uint16_t>& modes, double cutoff, int flags ); // Karhunen-Loeve
            
            void getNorms(const redux::image::Pupil&);         //!< Calculate the normalization factors so that sum(|mode*pupil|) = pupilArea
            void getNorms( void ) { getNorms( redux::image::Pupil::fetch( info.nPupilPixels, info.pupilRadius ) ); };
            
            void setPupilSize( uint16_t nPixels, double radiusInPixels, double rotation );
            
            void normalize( double scale=1.0 );
            
            ModeInfo info;
            
            redux::util::PointI tiltMode;                   //!< Contains the indices for the tilt modes. (-1 if no tilts are present)
            redux::util::PointF shiftToAlpha;               //!< Converts a shift of 1 pixel into corresponding mode-coefficient (NB. if image-size != pupil-size it has to be re-scaled accordingly)
            
            std::vector<uint16_t> modeNumbers;
            std::vector<double*> modePointers;
            std::vector<float> atm_rms;
            std::vector<float> norms;                       //!< (square of) L_{2,2} norms for the modes over the pupil
            std::mutex mtx;

            
        };

    }

}

#endif      // REDUX_MOMFBD_MODES_HPP
