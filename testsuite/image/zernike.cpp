#include "redux/image/zernike.hpp"

#include "catch.hpp"

#include <cmath>

using namespace redux::image;

namespace testsuite {

    namespace image {

        TEST_CASE( "Zernike Noll indexing", "[image]" ) {

            // Standard Noll (1976) ordering, cross-checked against the well-known table.
            struct { unsigned int j; uint16_t n; int16_t m; } table[] = {
                {  1, 0,  0 },      // piston
                {  2, 1,  1 },      // tip
                {  3, 1, -1 },      // tilt
                {  4, 2,  0 },      // defocus
                {  5, 2, -2 },      // astigmatism
                {  6, 2,  2 },      // astigmatism
                {  7, 3, -1 },      // coma
                {  8, 3,  1 },      // coma
                {  9, 3, -3 },
                { 10, 3,  3 },
                { 11, 4,  0 },      // spherical
                { 12, 4,  2 },
                { 13, 4, -2 },
                { 14, 4,  4 },
                { 15, 4, -4 },
            };

            for( auto& t: table ) {
                uint16_t n(0);
                int16_t m(0);
                Zernike::NollToNM( t.j, n, m );
                CHECK( n == t.n );
                CHECK( m == t.m );
                CHECK( Zernike::NMToNoll( t.n, t.m ) == t.j );
            }

            // round-trip for a wider range (through n~10)
            for( unsigned int j=1; j<=66; ++j ) {
                uint16_t n(0);
                int16_t m(0);
                Zernike::NollToNM( j, n, m );
                CHECK( Zernike::NMToNoll( n, m ) == j );
            }

        }


        TEST_CASE( "Zernike radial polynomial", "[image]" ) {

            // Classical closed-form radial polynomials R_n^m(r), used as ground truth.
            struct { uint16_t n, m; long double (*fn)(long double); } polys[] = {
                { 0, 0, []( long double r ) -> long double { return 1; } },
                { 1, 1, []( long double r ) -> long double { return r; } },
                { 2, 0, []( long double r ) -> long double { return 2*r*r - 1; } },
                { 2, 2, []( long double r ) -> long double { return r*r; } },
                { 3, 1, []( long double r ) -> long double { return 3*r*r*r - 2*r; } },
                { 3, 3, []( long double r ) -> long double { return r*r*r; } },
                { 4, 0, []( long double r ) -> long double { return 6*r*r*r*r - 6*r*r + 1; } },
                { 4, 2, []( long double r ) -> long double { return 4*r*r*r*r - 3*r*r; } },
                { 4, 4, []( long double r ) -> long double { return r*r*r*r; } },
            };

            const long double tol = 1E-9L;
            long double rvals[] = { 0.0L, 0.3L, 0.7L, 1.0L };

            for( auto& p: polys ) {
                for( int useOld=0; useOld<2; ++useOld ) {
                    int flags = useOld ? Zernike::OLD_METHOD : 0;
                    Zernike::RadialPolynomial& poly = Zernike::getRadialPolynomial( p.n, p.m, flags );
                    for( long double r: rvals ) {
                        long double expected = p.fn(r);
                        long double got = poly.eval(r);
                        CHECK( std::fabs( (double)(got-expected) ) < (double)tol );
                    }
                }
            }

        }


        TEST_CASE( "Zernike pixelized modes", "[image]" ) {

            // nPixels even -> grid origin falls exactly between pixels 4 and 5 (0-based),
            // i.e. shiftedOrigin = (5.0,5.0) for nPixels=10.
            // Gives exact, not just approximate, pixel distances/angles to check against.
            const unsigned int nPixels = 10;
            const float radius = 5.0;
            const size_t blockSize = nPixels*nPixels;
            std::vector<double> mode( blockSize );

            auto idx = [&]( unsigned y, unsigned x ) { return y*nPixels + x; };

            // Noll 1: piston (n=0,m=0) = constant 1 everywhere.
            Zernike::getZernike( mode.data(), nPixels, radius, 0, 1u, 0 );
            for( size_t i=0; i<blockSize; ++i ) {
                CHECK( std::fabs( mode[i] - 1.0 ) < 1E-9 );
            }

            // Noll 2: tip (n=1,m=1)
            // mode(y,x) = 2*(x - nPixels/2)/radius, independent of
            // y (R_1^1(r)=r, angular=cos(a), normalization=sqrt(2*(n+1))=2,
            // and r*cos(a) = xDist/radius by definition of polar coordinates).
            Zernike::getZernike( mode.data(), nPixels, radius, 0, 2u, 0 );
            for( unsigned y=0; y<nPixels; ++y ) {
                for( unsigned x=0; x<nPixels; ++x ) {
                    double xDist = (double)x - nPixels/2.0;
                    double yDist = (double)y - nPixels/2.0;
                    if( (xDist*xDist + yDist*yDist) > radius*radius ) continue;
                    double expected = 2.0*xDist / radius;
                    CHECK( std::fabs( mode[idx(y,x)] - expected ) < 1E-9 );
                }
            }

            // Noll 3: tilt (n=1,m=-1)
            // mode(y,x) = 2*(y - nPixels/2)/radius, independent of x. Same r<=1
            // restriction as Noll 2 above.
            Zernike::getZernike( mode.data(), nPixels, radius, 0, 3u, 0 );
            for( unsigned y=0; y<nPixels; ++y ) {
                for( unsigned x=0; x<nPixels; ++x ) {
                    double xDist = (double)x - nPixels/2.0;
                    double yDist = (double)y - nPixels/2.0;
                    if( (xDist*xDist + yDist*yDist) > radius*radius ) continue;
                    double expected = 2.0*yDist / radius;
                    CHECK( std::fabs( mode[idx(y,x)] - expected ) < 1E-9 );
                }
            }

            // Noll 4: defocus (n=2,m=0)
            // mode(y,x) = sqrt(3)*(2*r^2-1), r^2=(xDist^2+yDist^2)/radius^2.
            Zernike::getZernike( mode.data(), nPixels, radius, 0, 4u, 0 );
            for( unsigned y=0; y<nPixels; ++y ) {
                for( unsigned x=0; x<nPixels; ++x ) {
                    double xDist = (double)x - nPixels/2.0;
                    double yDist = (double)y - nPixels/2.0;
                    double r2 = (xDist*xDist + yDist*yDist)/(radius*radius);
                    if( r2 > 1.0 ) continue;
                    double expected = std::sqrt(3.0)*(2.0*r2 - 1.0);
                    CHECK( std::fabs( mode[idx(y,x)] - expected ) < 1E-9 );
                }
            }

        }

    }

}
