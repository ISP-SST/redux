#ifndef REDUX_IMAGE_ZERNIKE_HPP
#define REDUX_IMAGE_ZERNIKE_HPP

#include <map>
#include <memory>
#include <mutex>
#include <vector>

#include <boost/multiprecision/cpp_dec_float.hpp>

 
namespace redux {

    namespace image {
        
        class Zernike {
        public:
            typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100> > mp_float;
            
            struct PairID {
                PairID(unsigned int a, unsigned int b) : first(a), second(b) {}
                unsigned int first, second;
                bool operator<(const PairID& rhs) const;
            };
            
            struct RadialID {
                RadialID( unsigned int nPixels, float rr, uint16_t nn, uint16_t mm ) : nP(nPixels), r(rr), m(mm), n(nn) {}
                unsigned int nP;
                float r;
                uint16_t m, n;
                bool operator<(const RadialID& rhs) const;
            };

            struct AngularID {
                AngularID( unsigned int nPixels, float aa, int16_t mm ) : nP(nPixels), angle(aa), m(mm) {}
                unsigned int nP;
                float angle;
                int16_t m;
                bool operator<(const AngularID& rhs) const;
            };
            
            struct PolyID {
                PolyID( uint16_t nn, uint16_t mm, int f ) : n(nn), m(mm), flags(f&0xFF) {}
                uint16_t n,m;
                int flags;
                bool operator<(const PolyID& rhs) const;
            };
            
            // Karhunen-Loeve mode
            struct KL {
                std::vector< std::pair<unsigned int, double> > zernikeWeights; //< zernike mode numbers and corresponding weights
                double covariance;
            };
            struct RadialPolynomial {
                RadialPolynomial( uint16_t nn, uint16_t mm, int f=0 ) : n(nn), m(mm), flags(f&0xFF), useRatios(false) { };
                bool empty( void ) { return poly.empty(); };
                size_t size( void ) { return values.size(); };
                void calc( int flags );
                void calcKL( int flags );
                long double eval( long double r );

                uint16_t n, m;
                int flags;
                bool useRatios;
                std::map<uint16_t,mp_float,std::greater<uint16_t>> poly;    // sort with larger keys (exponent) first.
                std::map<long double, long double> values;                  // cache calculated polynomial values given r
                
            };
        private:
            Zernike(){};
            
            static Zernike& get(void) { static Zernike ret; return ret; };
            
            std::vector<double> factorials;
            std::mutex mtx;
            
        public:
            
            enum Flags { KL_POLY         = 0x1,
                         OLD_METHOD      = 0x2,
                         GET_RADIAL      = 0x100,                          // return only radial part
                         GET_ANGULAR     = 0x200,                          // return only angluar part
                         GET_ZERNIKE     = GET_RADIAL|GET_ANGULAR,         // return radial*angular    (default)
                         NORMALIZE       = 0x10000,                        // normalize numerically
                         FORCE           = 0x20000,                        // force re-calculation (i.e. update the cached radial/angular parts)
                         VERBOSE         = 0x40000                         // print some stuff
            };
            
            typedef std::shared_ptr<KL> KLPtr;
            
            static void NollToNM( unsigned int NollIndex, uint16_t& ZernikeN , int16_t& ZernikeM );
            static unsigned int NMToNoll( const uint16_t ZernikeN, const int16_t ZernikeM );
            
            static double calcCovariance( uint16_t n1, int16_t m1, uint16_t n2, int16_t m2 );
            static double calcCovariance( unsigned int NollIndex1, unsigned int NollIndex2 );
            static double getCovariance( unsigned int NollIndex1, unsigned int NollIndex2 );
            
            static std::shared_ptr<double> getRadial( unsigned int nPixels, float radius, uint16_t n, uint16_t m, int flags );
            static std::shared_ptr<double> getAngular( unsigned int nPixels, float angle, int16_t m , int flags );
            
            static void getZernike( double* modePtr, unsigned int nPixels, float radius, float angle,
                                    uint16_t ZernikeN, int16_t ZernikeM, int flags );
            static void getZernike( double* modePtr, unsigned int nPixels, float radius, float angle,
                                    unsigned int NollIndex, int flags ) {
                uint16_t n(0);
                int16_t m(0);
                NollToNM( NollIndex, n, m );
                getZernike( modePtr, nPixels, radius, angle, n, m, flags );
            }
            
            static RadialPolynomial& getRadialPolynomial( uint16_t n, uint16_t m, int flags );
            
            static const std::map<uint16_t, KLPtr>& karhunenLoeveExpansion( uint16_t first_mode, uint16_t last_mode );
            
            static void clear(void);

        };

                
    }   // image

}   // redux


#endif  // REDUX_IMAGE_ZERNIKE_HPP
