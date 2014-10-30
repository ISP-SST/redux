#include "redux/math/helpers.hpp"

#include <cstdlib>

namespace {

    /* See http://randomascii.wordpress.com/2012/01/11/tricks-with-the-floating-point-format/
     * for the potential portability problems with the union and bit-fields below.
    */
    union Float_t {
        int32_t i;
        float f;

        Float_t( float num = 0.0f ) : f( num ) {}
        bool diffSign( const Float_t& rhs ) const { return ( ( ( i ^ rhs.i ) & 0x80000000 ) != 0 ); }

    };

}

bool redux::math::almostEqual( float A, float B, int maxUlpsDiff ) {

    Float_t uA( A );
    Float_t uB( B );

    // Different signs
    if( uA.diffSign( uB ) ) {
        // Check for equality (i.e. +0 == -0)
        if( A == B ) {
            return true;
        }
        return false;
    }

    // Find the difference in ULPs.
    int ulpsDiff = abs( uA.i - uB.i );
    if( ulpsDiff <= maxUlpsDiff ) {
        return true;
    }

    return false;
}
