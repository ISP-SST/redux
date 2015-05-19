#ifndef REDUX_MOMFBD_CACHE_HPP
#define REDUX_MOMFBD_CACHE_HPP

#include "redux/momfbd/modes.hpp"
#include "redux/image/utils.hpp"

#include <map>
#include <memory>
#include <mutex>
#include <vector>

namespace redux {

    namespace momfbd {

        class Cache {
            
        public:
            
            struct ModeID {
                uint16_t firstMode, lastMode, modeNumber, nPoints;
                double pupilRadius, wavelength, angle;
                ModeID ( uint16_t modeNumber=0, uint16_t nPoints=0, double pupilRadius=0.0, double wavelength=0.0, double angle=0.0 );
                ModeID ( uint16_t firstMode, uint16_t lastMode, uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double wavelength, double angle );
                uint64_t size(void) const;
                uint64_t pack(char*) const;
                uint64_t unpack(const char*, bool);
                bool operator<( const ModeID& rhs ) const;
            };

            static Cache& getCache( void );

            /*! @note This will invalidate all references, so do NOT clear the cache while it's in use.
             */
            void clear( void );

            const redux::image::Grid& grid( uint32_t sz, PointF=0 );
            const std::pair<redux::util::Array<double>, double>& pupil( uint32_t nPoints, float radius );
            double zernikeCovariance( int32_t m, int32_t n );
            const std::vector<double>& zernikeRadialPolynomial( int32_t m, int32_t n );
            const std::map<uint16_t, PupilMode::KLPtr>& karhunenLoeveExpansion( uint16_t first_mode, uint16_t last_mode );
            const PupilMode::Ptr addMode ( const ModeID&, PupilMode::Ptr& m);
            const PupilMode::Ptr addMode ( uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double wavelength, double angle, PupilMode::Ptr& m);
            const PupilMode::Ptr mode ( const ModeID& );
            const PupilMode::Ptr mode ( uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double wavelength, double angle );
            const PupilMode::Ptr mode ( uint16_t firstMode, uint16_t lastMode, uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double wavelength, double angle );
            
//             ~Cache() { std::cout << "~Cache(): nModes = "<< modes.size() << std::endl; }
        private:
            Cache() {}
            Cache( const Cache& ) {}

            std::mutex mtx;

            std::vector<double> factorials;
            std::map<std::pair<int32_t, int32_t>, const double> zernikeCovariances;
            std::map<std::pair<uint16_t, uint16_t>, const std::vector<double> > zernikeRadialPolynomials;
            std::map<std::pair<uint16_t, uint16_t>, const std::map<uint16_t, PupilMode::KLPtr> > karhunenLoeveExpansions;
            std::map<std::pair<uint32_t, PointF>, const redux::image::Grid> grids;
            std::map<std::pair<uint32_t, float>, const std::pair<redux::util::Array<double>, double>> pupils;
            std::map<ModeID, const PupilMode::Ptr> modes;

        };

    }   // image

}   // redux


#endif  // REDUX_MOMFBD_CACHE_HPP
