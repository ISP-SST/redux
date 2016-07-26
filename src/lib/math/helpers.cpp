#include "redux/math/helpers.hpp"

#include "redux/constants.hpp"

#include <cstdlib>
#include <iostream>

using namespace redux;
using namespace std;


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

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define TINY 1.0E-32

template<typename T, typename U>
void redux::math::bracket(std::function<U(T)> func, T& a, T& b, T& mid, U& fa, U& fb, U& fmid, int limit) {
    
    double u, fu, ulim;

    if( fb > fa ) {         // make sure we are going downhill
        std::swap(a,b);
        std::swap(fa,fb);
    }
    
    mid = b+GoldenRatio*(b-a);
    fmid = func(mid);
//printf ("bracket        (%.18e,%.18e,%.18e) = (%.18e,%.18e,%.18e)\n", a, mid, b, fa, fmid, fb);
    int it(0);
    while( fb > fmid ) {
        it++;
//printf ("bracket it=%d  (%.18e,%.18e,%.18e) = (%.18e,%.18e,%.18e)\n", it, a, mid, b, fa, fmid, fb);
        double r = (b-a)*(fb-fmid);
        double q = (b-mid)*(fb-fa);
        //double qr = 2.0 * (q-r);
        //double absqr = fabs(qr);
        //if( absqr < 1E-20 ) qr /= absqr*1E20;
        //u = b-((b-mid)*q-(b-a)*r)/qr;
        u = (b)-((b-mid)*q-(b-a)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
        ulim = b+limit*(mid-b);
        if ((b-u)*(u-mid) > 0.0) {      // u is between b and mid
            fu = func(u);
            if( fu < fmid ) {             // minimum is between b and mid, we're done
                a = b;
                b = u;
                fa = fb;
                fb = fu;
                break;
            } else if( fu > fb ) {      // minimum is between a and u, we're done
                mid = u;
                fmid = fu;
                break;
            }
            u = mid+GoldenRatio*(mid-b); // bad parabolic fit, try again with larger step
            fu = func(u);
        } else if( (mid-u)*(u-ulim) > 0.0 ) {        // u is between mid and ulim
            fu = func(u);
            if( fu < fmid ) {
                SHFT( b, mid, u, mid+GoldenRatio*(mid-b) )
                SHFT( fb, fmid, fu, func(u) )
            }
        } else if( (mid-ulim)*(ulim-u) >= 0.0 ) {    // u > ulim, truncate
            u = ulim;
            fu = func(u);
        } else {
            u = mid+GoldenRatio*(mid-b);
            fu = func(u);
        }
        SHFT( a, b, mid, u )
        SHFT( fa, fb, fmid, fu )
    }

    if( a > b ) {   // TODO: this can probably be optimized 
        std::swap(a,b);
        std::swap(fa,fb);
    }
    if( mid > b ) {
        std::swap(mid,b);
        std::swap(fmid,fb);
    }
    if( mid < a ) {
        std::swap(mid,a);
        std::swap(fmid,fa);
    }
//printf ("bracketE       (%.18e,%.18e,%.18e) = (%.18e,%.18e,%.18e)\n", a, mid, b, fa, fmid, fb);
    
}
template void redux::math::bracket(std::function<double(double)>, double&, double&, double&, double&, double&, double&, int);


template<typename T, typename U>
void redux::math::brent( std::function<U(T)> func, T ax, T bx, T& x, U& fx, double tol, int limit ) {
    
    static double CGOLD = 2.0-GoldenRatio;
    
    T a = std::min( ax, bx );
    T b = std::max( ax, bx );
    
    T u, v, w;
    u = v = w = x;
    
    U fu,fv,fw;
    fu = fv = fw = fx = func(x);
    
    double d(0.0);
    double e(0.0);
    
//printf ("brent0  (%e,%e,%e) = (%e)   tol=%e\n", a, x, b, fw, tol);
    for( int iter=1; iter <= limit; ++iter ) {
        double xm = 0.5*(a+b);
        double tol1 = tol*fabs(x) + TINY;
        double tol2 = 2.0*tol1;
//printf ("   iter=%d  [%e, %e]  tol1=%e  tol2=%e  x=%e  val=%e   test=%e\n", iter, b, a, tol1, tol2, x, fx, fabs(x-xm)+0.5*(b-a) );
        if( fabs(x-xm) <= (tol2-0.5*(b-a)) ) {
            return;
        }
        if( fabs(e) > tol1 ) {
            double r = (x-w)*(fx-fv);
            double q = (x-v)*(fx-fw);
            double p = (x-v)*q-(x-w)*r;
            q = 2.0*(q-r);
            if( q > 0.0 ) p = -p;
            q = fabs(q);
            double etemp = e;
            e = d;
            if( fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x) ) {
                e = (x >= xm ? a-x : b-x);
                d = CGOLD*e;
            } else {
                d = p/q;
                u = x+d;
                if( u-a < tol2 || b-u < tol2 ) {
                    //d = xm < x ? -tol1 : tol1;
                    d = SIGN(tol1,xm-x);
                }
            }
        } else {
            e = (x >= xm ? a-x : b-x);
            d = CGOLD*e;
        }
        if( fabs(d) < tol1 ) d *= (tol1/fabs(d));
        //u = x+d; //(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
        u = ((fabs(d)>=tol1)?x+d:x+SIGN(tol1,d));
        fu = func(u);
        if( fu <= fx ) {
            if( u >= x ) {
                a = x;
            } else {
                b = x;
            }
            SHFT( v, w, x, u )
            SHFT( fv, fw, fx, fu )
        } else {
            if( u < x ) {
                a = u; 
            } else {
                b = u;
            }
            if( fu <= fw || w == x ) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if( fu <= fv || v == x || v == w ) {
                v = u;
                fv = fu;
            }
        }
    }
    
    cout << "Error in brent: too many iterations!" << endl;

}
template void redux::math::brent(std::function<double(double)>, double, double, double&, double&, double, int);

#undef SHFT
