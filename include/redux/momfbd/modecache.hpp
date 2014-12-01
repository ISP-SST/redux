#ifndef REDUX_MOMFBD_MODECACHE_HPP
#define REDUX_MOMFBD_MODECACHE_HPP

#include "redux/momfbd/modes.hpp"
#include "redux/image/utils.hpp"

#include <map>
#include <memory>
#include <mutex>
#include <vector>

namespace redux {

    namespace momfbd {

        class ModeCache {
            
        public:
            
            typedef std::shared_ptr<Modes::PupilMode> ModePtr;
            
            struct index {
                index ( int modeNumber, int nPoints, double r_c, double wavelength, double angle );
                index ( int firstMode, int lastMode, int modeNumber, int nPoints, double r_c, double wavelength, double angle );
                int firstMode, lastMode, modeNumber, nPoints;
                double r_c, wavelength, angle;
                bool operator<( const index& rhs ) const;
            };

            static ModeCache& getCache( void );

            /*! @note This will invalidate all references, so do NOT clear the cache while it's in use.
             */
            void clear( void );

            const redux::image::Grid& grid( int sz );
            double zernikeCovariance( int m, int n );
            const std::vector<double>& zernikeRadialPolynomial( int m, int n );
            const std::map<int, Modes::KL_cfg>& karhunenLoeveExpansion( int first_mode, int last_mode );
            const ModePtr& mode ( int modeNumber, int nPoints, double r_c, double wavelength, double angle );
            const ModePtr& mode ( int firstMode, int lastMode, int modeNumber, int nPoints, double r_c, double wavelength, double angle );
            
        private:
            ModeCache() {}
            ModeCache( const ModeCache& ) {}

            std::mutex mtx;

            std::vector<double> factorials;
            std::map<std::pair<int, int>, const double> zernikeCovariances;
            std::map<std::pair<int, int>, const std::vector<double> > zernikeRadialPolynomials;
            std::map<std::pair<int, int>, const std::map<int, Modes::KL_cfg> > karhunenLoeveExpansions;
            std::map<int, const redux::image::Grid> grids;
            std::map<index, const ModePtr> modes;

        };

    }   // image

}   // redux


#endif  // REDUX_MOMFBD_MODECACHE_HPP
