#ifndef REDUX_MOMFBD_MODES_HPP
#define REDUX_MOMFBD_MODES_HPP

#include "redux/util/array.hpp"
#include "redux/image/pupil.hpp"

#include <mutex>

#define MT_NONE   0
#define MT_PUPIL  1
#define MT_OTF_BF 2

namespace redux {

    namespace momfbd {
        
        class MomfbdJob;
        class Object;
        
        struct ModeInfo {
            
            uint16_t firstMode, lastMode, modeNumber, nPupilPixels;
            double pupilRadius, angle, cutoff;
            std::string filename;
            ModeInfo ( std::string, uint16_t nPixels=0 );
            ModeInfo ( uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double angle );
            ModeInfo ( uint16_t firstMode, uint16_t lastMode, uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double angle, double cutoff );
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            bool operator<( const ModeInfo& rhs ) const;
            
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
            PupilMode ( uint16_t modeNumber, uint16_t nPoints, double r_c = 1.0, double angle = 0.0 ); // Zernike
            PupilMode ( uint16_t firstMode, uint16_t lastMode, uint16_t klModeNumber, uint16_t nPoints, double r_c = 1.0, double angle = 0.0, double cutoff=0.0 ); // KL

            operator const redux::util::Array<double>&() const { return reinterpret_cast<const redux::util::Array<double>&>(*this); }
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
            
            bool load( const std::string& filename, uint16_t pixels );
            void init( const MomfbdJob& job, const Object& ch );
            void generate( uint16_t pixels, double radius, double angle, const std::vector<uint16_t>& modes );  // Zernike
            void generate( uint16_t pixels, double radius, double angle, uint16_t firstMode, uint16_t lastMode, const std::vector<uint16_t>& modes, double cutoff ); // Karhunen-Loeve
            void normalize(const redux::image::Pupil&);         //!< Normalize modes so that sum(|mode*pupil|) = pupilArea
            
            void setPupilSize( uint16_t nPixels, double radiusInPixels, double rotation );
            
            ModeInfo info;
            
            int32_t xTiltIndex, yTiltIndex;
            std::vector<uint16_t> modeNumbers;
            std::vector<double*> modePointers;
            std::mutex mtx;

            
        };

    }

}

#endif      // REDUX_MOMFBD_MODES_HPP
