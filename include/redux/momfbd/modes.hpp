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
        
        enum ModeBase { MB_NONE=0, ZERNIKE=1, KARHUNEN_LOEVE, MB_ZKL, MB_FILE };
        
        struct ModeID {
            ModeID() : mode(0), type(ZERNIKE) {}
            ModeID( uint16_t m, ModeBase t=ZERNIKE ) : mode(m), type(t) {}
            bool operator==( const ModeID& rhs ) const { if( type==rhs.type ) return mode == rhs.mode; return type==rhs.type; }
            bool operator<( const ModeID& rhs ) const { if( type==rhs.type ) return mode < rhs.mode; return type<rhs.type; }
            uint16_t mode;
            ModeBase type;
        };
        typedef std::vector<struct ModeID> VModeID;
        
        struct ModeList : public VModeID {
            bool operator<( const ModeList& rhs ) const{
                if( defaultType==rhs.defaultType ) {
                    return (reinterpret_cast<const VModeID&>(*this) < reinterpret_cast<const VModeID&>(rhs));
                }
                return defaultType<rhs.defaultType;
            }
            void setDefaultModeType( ModeBase new_mb, ModeBase old_mb=MB_NONE ) {
                defaultType = new_mb;
                for( auto& a: *this ) {
                    if( a.type == old_mb ) a.type = new_mb;
                }
            }
            ModeBase defaultType;
        };
        
        struct DiversityValue {
            DiversityValue() : coefficient(0), physical(false) {}
            DiversityValue( double c, bool p=false ) : coefficient(c), physical(p) {}
            bool operator==( const DiversityValue& rhs ) const {
                if( physical==rhs.physical ) {
                    return coefficient == rhs.coefficient;
                }
                return physical==rhs.physical;
            }
            double coefficient;
            bool physical;
        };


        struct ModeInfo {
            
            uint16_t firstMode, lastMode, modeNumber, nPupilPixels;
            ModeList modeList;
            double pupilRadius, angle, cutoff;
            std::string filename;
            explicit ModeInfo( const std::string&, uint16_t nPixels=0 );
            ModeInfo( uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double angle );
            ModeInfo( uint16_t firstMode, uint16_t lastMode, uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double angle, double cutoff );
            ModeInfo( uint16_t firstMode, uint16_t lastMode, const ModeList& modeList, uint16_t nPoints, double pupilRadius, double angle, double cutoff );
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            bool operator<( const ModeInfo& rhs ) const;
            operator std::string() const;
        };
        inline std::ostream& operator<<( std::ostream& os, const ModeInfo& mi ) {
            os << (std::string)mi;
            return os;
        }



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
            void generate( uint16_t pixels, double radius, double angle, const ModeList& modes, int flags );  // Zernike
            void generate( uint16_t pixels, double radius, double angle, uint16_t firstMode, uint16_t lastMode, const ModeList& modes, double cutoff, int flags ); // Karhunen-Loeve
            
            void getNorms(const redux::image::Pupil&);         //!< Calculate the normalization factors so that sum(|mode*pupil|) = pupilArea
            void getNorms( void ) { getNorms( redux::image::Pupil::fetch( info.nPupilPixels, info.pupilRadius ) ); };
            
            void measureJacobian( const redux::image::Pupil&, double scale=1.0 );
            template <typename T> redux::util::PointType<T> alphaToShift( redux::util::PointType<T> alpha ) {
                return redux::util::PointType<T>( alpha.x * J_T_xy[1][0] + alpha.y * J_T_xy[1][1],
                                                  alpha.x * J_T_xy[0][0] + alpha.y * J_T_xy[0][1] );
            }
            template <typename T> redux::util::PointType<T> shiftToAlpha( redux::util::PointType<T> shift ) {
                return redux::util::PointType<T>( shift.x * J_xy_T[1][0] + shift.y * J_xy_T[1][1],
                                                  shift.x * J_xy_T[0][0] + shift.y * J_xy_T[0][1] );
            }
            
            void setPupilSize( uint16_t nPixels, double radiusInPixels, double rotation );
            
            void normalize( double scale=1.0 );
            
            ModeInfo info;
            
            double J_T_xy[2][2];          //!< Jacobian for the tilts w.r.t. the pixel coordinates. Used for converting tilts to pixels.
            double J_xy_T[2][2];          //!< Inverse of the above Jacobian. Used for converting pixels to tilts.
            
            redux::util::PointI tiltMode;                   //!< Contains the indices for the tilt modes. (-1 if no tilts are present)
            
            ModeList modeList;
            std::vector<double*> modePointers;
            std::vector<double> atm_rms;
            std::vector<double> norms;                      //!< (square of) L_{2,2} norms for the modes over the pupil
            std::mutex mtx;

            
        };

    }

}

#endif      // REDUX_MOMFBD_MODES_HPP
