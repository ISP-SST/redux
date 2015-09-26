#ifndef REDUX_IMAGE_ZERNIKE_HPP
#define REDUX_IMAGE_ZERNIKE_HPP

#include <map>
#include <memory>
#include <mutex>
#include <vector>

namespace redux {

    namespace image {
        
        class Zernike {
            
            struct PairID {
                PairID(int32_t a, int32_t b) : first(a), second(b) {}
                int32_t first, second;
                bool operator<(const PairID& rhs) const;
            };
            
            Zernike(){};
            
            static Zernike& get(void) { static Zernike ret; return ret; };
            
            std::vector<double> factorials;
            std::mutex mtx;
            
        public:
            
            // Karhunen-Loeve expansion coefficients
            struct KL {
                std::vector< std::pair<uint16_t, double> > zernikeWeights; //< zernike mode numbers and corresponding weights
                double covariance;
            };
            typedef std::shared_ptr<KL> KLPtr;
            
            static double covariance( int32_t m, int32_t n );
            static const std::vector<double>& radialPolynomial( int32_t m, int32_t n );
            
            static const std::map<uint16_t, KLPtr>& karhunenLoeveExpansion( uint16_t first_mode, uint16_t last_mode );
            
            static void clear(void);

        };

    }   // image

}   // redux


#endif  // REDUX_IMAGE_ZERNIKE_HPP
