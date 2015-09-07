/*
 * Based on original DLM library by Louis Strous, Lockheed Martin Solar &
 * Astrophysics Lab, <ana@lmsal.com>. Imported into the redux library and
 * modified to fit c++/cmake by Tomas Hillberg <hillberg@astro.su.se>
 *
 */


#include "gridmatch.hpp"

#include <algorithm>
#include <math.h>

using namespace std;


#define IDL_TYP_B_BASIC (IDL_TYP_MASK(IDL_TYP_BYTE) | \
 IDL_TYP_MASK(IDL_TYP_INT) | IDL_TYP_MASK(IDL_TYP_LONG) | \
 IDL_TYP_MASK(IDL_TYP_FLOAT) | IDL_TYP_MASK(IDL_TYP_DOUBLE))


namespace {

    union pointerUnion {
        unsigned char *b;
        short         *w;
        IDL_LONG      *l;
        float         *f;
        double        *d;
    };
    typedef union pointerUnion      Pointer;
    
    IDL_LONG    badmatch, stretch_clip;

    void getmin( float *p, float *x0, float *y0 ) {

        float f11, f12, f13, f21, f22, f23, f31, f32, f33;
        float fx, fy, t, fxx, fyy, fxy;

        /* find the min, p points to a 3x3 array */
        f11 = *p++;   f21 = *p++; f31 = *p++;
        f12 = *p++;   f22 = *p++; f32 = *p++;
        f13 = *p++;   f23 = *p++; f33 = *p++;

        fx = 0.5 * ( f32 - f12 ); fy = 0.5 * ( f23 - f21 );
        t = 2.* ( f22 );
        fxx =  f32 + f12 - t; fyy = f23 + f21 - t;
        /* find in which quadrant the minimum lies */
        if( f33 < f11 ) {
            if( f33 < f31 ) {
                if( f33 < f13 )
                    fxy = f33 + f22 - f32 - f23;
                else
                    fxy = f23 + f12 - f22 - f13;
            } else {
                if( f31 < f13 )
                    fxy = f32 + f21 - f31 - f22;
                else
                    fxy = f23 + f12 - f22 - f13;
            }
        } else {
            if( f11 < f31 ) {
                if( f11 < f13 )
                    fxy = f22 + f11 - f21 - f12;
                else
                    fxy = f23 + f12 - f22 - f13;
            } else {
                if( f31 < f13 )
                    fxy = f32 + f21 - f31 - f22;
                else
                    fxy = f23 + f12 - f22 - f13;
            }
        }
        t = -1. / ( fxx * fyy - fxy * fxy );
        *x0 = t * ( fx * fyy - fy * fxy );
        *y0 = t * ( fy * fxx - fx * fxy );
        if( fabs( *x0 ) >= 0.75 || fabs( *y0 ) >= 0.75 ) {
            *x0 = -fx / fxx;
            *y0 = -fy / fyy;
        }

    }


    float averag( float *p, IDL_LONG nxa, IDL_LONG nxb, IDL_LONG nya,
                  IDL_LONG nyb, IDL_LONG nxs, IDL_LONG nys, IDL_LONG idx,
                  IDL_LONG idy, float *gx, float *gy )
    /* finds weighted average of array m over the block defined */
    {
        IDL_LONG  nxc, nxd, nyc, nyd, i, j, jj;
        float sum, sumg, sumx, sumgx;

        /* fix limits so sum doesn't run off edge of image */
        nxc = ( nxa + idx < 0 ) ? -idx : nxa;
        nyc = ( nya + idy < 0 ) ? -idy : nya;
        nxd = ( nxb + idx > nxs ) ? nxs - idx : nxb;
        nyd = ( nyb + idy > nys ) ? nys - idy : nyb;
        sum = sumg = sumgx = 0.0;
        for( i = nxc; i < nxd; i++ )
            sumgx += gx[i];
        /* weighted sum in window */
        for( j = nyc; j < nyd; j++ ) {
            sumx = 0.0;
            jj = idx + nxs * ( j + idy );
            for( i = nxc; i < nxd; i++ )
                sumx += gx[i] * p[i + jj];
            sum += gy[j] * sumx;
            sumg += gy[j] * sumgx;
        } /* end of j loop */
        return sum / sumg;
    }


    void unbias( float *m1, float *m2, IDL_LONG nxa, IDL_LONG nxb,
                 IDL_LONG nya, IDL_LONG nyb, IDL_LONG nxs, IDL_LONG nys,
                 float *gx, float *gy, float *av1, float *av2, float *cx,
                 float *cy, float *cxx, float *cxy, float *cyy,
                 IDL_LONG idelx, IDL_LONG idely ) {

        float t0, t1, t2, t3, t4, t5;


        /*  find weighted means of m1 & m2 over the window
            sets up quadratic fit to average of m2 as a fcn. of offsets */
        *av1 = averag( m1, nxa, nxb, nya, nyb, nxs, nys, 0, 0, gx, gy );
        t0 = averag( m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely, gx, gy );
        t1 = averag( m2, nxa, nxb, nya, nyb, nxs, nys, idelx + 1, idely, gx, gy );
        t2 = averag( m2, nxa, nxb, nya, nyb, nxs, nys, idelx - 1, idely, gx, gy );
        t3 = averag( m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely + 1, gx, gy );
        t4 = averag( m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely - 1, gx, gy );
        t5 = averag( m2, nxa, nxb, nya, nyb, nxs, nys, idelx + 1, idely + 1, gx, gy );
        *av2 = t0;
        *cx = 0.5 * ( t1 - t2 );
        *cy = 0.5 * ( t3 - t4 );
        *cxx = 0.5 * ( t1 - 2 * t0 + t2 );
        *cyy = 0.5 * ( t3 - 2 * t0 + t4 );
        *cxy = t5 + t0 - t1 - t3;
    }


    void gwind0( float *gwx, float *gwy, float gwid, IDL_LONG nxa, IDL_LONG nxb, IDL_LONG nya, IDL_LONG nyb ) {

        float wid, xcen, ycen, xq;
        IDL_LONG  i;
        wid = gwid * 0.6005612;
        if( wid > 0 ) {
            xcen = ( nxa + nxb ) / 2;
            ycen = ( nya + nyb ) / 2;
            for( i = nxa; i <= nxb; i++ ) {
                xq = ( i - xcen ) / wid;
                gwx[i] = exp( -( xq * xq ) );
            }
            for( i = nya; i <= nyb; i++ ) {
                xq = ( i - ycen ) / wid;
                gwy[i] = exp( -( xq * xq ) );
            }
        } else {
            for( i = nxa; i <= nxb; i++ )
                gwx[i] = 1.0;
            for( i = nya; i <= nyb; i++ )
                gwy[i] = 1.0;
        }
    }



    float resid( float *m1, float *m2, IDL_LONG idx, IDL_LONG idy, IDL_LONG nxa,
                 IDL_LONG nxb, IDL_LONG nya, IDL_LONG nyb, IDL_LONG nxs,
                 IDL_LONG nys, IDL_LONG ndmx, float *gx, float *gy, float bs ) {

        IDL_LONG  nxc, nxd, nyc, nyd, nx, ny;
        float *p1, *p2, *ps;
        float sum, sumx, t, ndmx2;
        IDL_LONG  i, j;
        float   sumg;
        static IDL_LONG   mxc, mxd, myc, myd;
        static float  gsum;

        /*set up limits */
        nxc = nxa;
        if( nxc + idx < 0 )
            nxc = -idx;
        nyc = nya;
        if( nyc + idy < 0 )
            nyc = - idy;
        nxd = nxb;
        if( nxd + idx >= nxs )
            nxd = nxs - idx - 1;
        nyd = nyb;
        if( nyd + idy >= nys )
            nyd = nys - idy - 1;
        sum = sumg = 0.0;

        nx = nxd - nxc + 1;
        p2 = gy + nyc;
        ps = gx + nxc;

        if( nxc != mxc || nxd != mxd || nyc != myc || nyd != myd ) {
            /* sum gaussians over rectangle to get normalization */
            /* (only if limits change)*/
            j = nyd - nyc + 1;
            if( j <= 0 || nxd - nxc + 1 <= 0 )
                return -1;        /* added 19feb95 LS */
            while( j ) {
                i = nx;
                p1 = ps;
                while( i ) {
                    sumg += ( *p1++ ) * ( *p2 );
                    i--;
                }
                p2++;
                j--;
            }
            gsum = sumg;
            mxc = nxc;
            mxd = nxd;
            myc = nyc;
            myd = nyd;
        } else
            sumg = gsum;

        m1 += nyc * nxs + nxc;
        m2 += ( nyc + idy ) * nxs + nxc + idx;
        ny = nxs - nx;        /*residual increment after inner loop */
        p2 = gy + nyc;

        /* now the loop to compute the residual */
        j = nyd - nyc + 1;
        ndmx2 = ndmx * ndmx;
        while( j ) {
            i = nx;
            p1 = ps;
            sumx = 0.0;
            while( i ) {
                t = *m1++ - *m2++;
                t = t + bs;
                t = t * t;
                t = min( t, ndmx2 );
                sumx += ( *p1++ ) * t;
                i--;
            }
            sum += ( *p2++ ) * sumx;
            m1 += ny;
            m2 += ny;
            j--;
        }
        /*return normalized residual */
        return sum / sumg;
    }


    void match_1( float *p1, float *p2, IDL_LONG nxa, IDL_LONG nxb,
                  IDL_LONG nya, IDL_LONG nyb, IDL_LONG nx, IDL_LONG ny,
                  float *gwx, float *gwy, float *xoffset, float *yoffset ) {

        IDL_LONG  idelx, idely, i, j, k, ndmx = 1000, done[9];
        IDL_LONG  di, dj, in, jn, iter, dd, badflag = 0;
        float av1, av2, cx, cy, cxx, cxy, cyy, avdif, t, res[9], buf[9], t1, t2;

        static IDL_LONG   itmax = 20;


        for( i = 0; i < 9; i++ )
            done[i] = 0;
        idelx = rint( *xoffset );
        idely = rint( *yoffset );

        /* look at a 3x3 matrix of residuals centered at 0 offset, find the location
           of the minimum, if not at center, then look at new 3x3 centered
        on the edge minimum; repeat until min at center */
        iter = itmax;
        badflag = 0;
        while( iter-- ) {
            for( k = 0; k < 9; k++ ) {
                if( done[k] == 0 ) {
                    i = idelx + ( k % 3 ) - 1;
                    j = idely + ( k / 3 ) - 1;
                    avdif = av2 +  i * cx + j * cy + i * i * cxx + i * j * cxy + j * j * cyy - av1;
                    res[k] = resid( p1, p2, i, j, nxa, nxb, nya, nyb, nx, ny, ndmx,
                                    gwx, gwy, avdif );
                }
            }
            t = res[0];
            i = 0;
            for( k = 1; k < 9; k++ )
                if( res[k] < t ) {
                    t = res[k];
                    i = k;
                }
            if( t < 0 ) {       /* added LS 19feb95 */
                printf( "match - ran out of data at edge\n" );
                badflag = 1;
                break;
            }
            idelx += ( i % 3 ) - 1;
            idely += ( i / 3 ) - 1;
            /* check if we have gone too far */
            if( abs( idelx ) > stretch_clip || abs( idely ) > stretch_clip ) {
                badflag++;
                break;
            }
            if( i == 4 )
                break;            /* done if it is the center one */
            /* not in center, shuffle what we have to put the edge min in center */
            di = ( i % 3 ) - 1;
            dj = ( i / 3 ) - 1;
            dd = dj * 3 + di;
            for( k = 0; k < 9; k++ ) {
                in = k % 3 + di;
                jn = k / 3 + dj;
                if( in >= 0 && jn >= 0 && in < 3 && jn < 3 ) { /* in range */
                    done[k] = 1;
                    buf[k] = res[k + dd];
                } else
                    done[k] = 0;        /* not in range, mark undone */
            }
            for( k = 0; k < 9; k++ )
                res[k] = buf[k];      /* put back in res array */
        } /* end of iter while */
        /* done or reached itmax, which ? */
        if( iter <= 0 ) {
            badflag++;
        }
        if( badflag ) {
            badmatch++;
            *xoffset = *yoffset = 0;
            return;
        }
        /* must have been OK so far */
        getmin( res, &t1, &t2 );
        *xoffset = idelx + t1;
        *yoffset = idely + t2;
    }


}


IDL_VPTR redux::gridmatch( int argc, IDL_VPTR argv[] )
/* the call is offsets = gridmatch(m1,m2,gx,gy,dx,dy,gwid,stretch_clip)
    where   m1 = reference input image
    m2 = image to compare with m1, m1 and m2 must be same size
    gx = array of x gridpoints
    gy = array of y gridpoints, gx and gy must have same size
    dx and dy are the window size, and gwid is the gaussian mask width
        stretch_clip = maximum allowed displacement before a bad match
          is declared

  Authors:  Richard A. Shine (original)
            Louis H. Strous (port from ANA to IDL)
            Lockheed-Martin Advanced Technology Center, Palo Alto, CA, USA
*/
{
    static IDL_EZ_ARG arg_struct[] = {
        {
            IDL_EZ_DIM_MASK( 2 ), IDL_TYP_B_SIMPLE, IDL_EZ_ACCESS_R, IDL_TYP_FLOAT,
            0, 0
        },           /* m1: 2D FLOAT read-only array */
        {
            IDL_EZ_DIM_MASK( 2 ), IDL_TYP_B_SIMPLE, IDL_EZ_ACCESS_R, IDL_TYP_FLOAT,
            0, 0
        },           /* m2: 2D FLOAT read-only array */
        {
            IDL_EZ_DIM_MASK( 2 ), IDL_TYP_B_SIMPLE, IDL_EZ_ACCESS_R, IDL_TYP_LONG,
            0, 0
        },           /* gx: 2D LONG read-only array */
        {
            IDL_EZ_DIM_MASK( 2 ), IDL_TYP_B_SIMPLE, IDL_EZ_ACCESS_R, IDL_TYP_LONG,
            0, 0
        },           /* gy: 2D LONG read-only array */
        {
            IDL_EZ_DIM_MASK( 0 ), IDL_TYP_B_SIMPLE, IDL_EZ_ACCESS_R, IDL_TYP_LONG,
            0, 0
        },           /* dx: FLOAT read-only scalar */
        {
            IDL_EZ_DIM_MASK( 0 ), IDL_TYP_B_SIMPLE, IDL_EZ_ACCESS_R, IDL_TYP_LONG,
            0, 0
        },           /* dy: FLOAT read-only scalar */
        {
            IDL_EZ_DIM_MASK( 0 ), IDL_TYP_B_SIMPLE, IDL_EZ_ACCESS_R, IDL_TYP_FLOAT,
            0, 0
        },           /* gwid: FLOAT read-only scalar */
        {
            IDL_EZ_DIM_MASK( 0 ), IDL_TYP_B_SIMPLE, IDL_EZ_ACCESS_R, IDL_TYP_LONG,
            0, 0
        }            /* stretch_clip: LONG read-only scalar */
    };
    IDL_VPTR  result;
    IDL_VARIABLE  zero;
    float         *out, *p1, *p2, *gwx, *gwy, xoffset, yoffset, gwid;
    IDL_LONG      *gx, *gy, nx, ny, nxg, nyg, dx2, dy2, nc, dx, dy, i1, i2, j1, j2;
    IDL_MEMINT    dims[IDL_MAX_ARRAY_DIM];

    if( argc < 7 || argc > 8 )
        IDL_Message( IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                     "GRIDMATCH needs 7 or 8 arguments:\nGRIDMATCH(image0, image1, grid_x, grid_y, window_x, window_y,\n gausswidth [, stretch_clip])" );
    IDL_EzCall( argc, argv, arg_struct );

    /* the dimensions of m1 and m2 must be equal */
    nx = arg_struct[0].value.arr->dim[0];
    ny = arg_struct[0].value.arr->dim[1];
    if( arg_struct[1].value.arr->dim[0] != nx
            || arg_struct[1].value.arr->dim[1] != ny ) {
        IDL_Message( IDL_M_NAMED_GENERIC, IDL_MSG_RET,
                     "First two arguments (images) have unequal dimensions" );
        IDL_EzCallCleanup( argc, argv, arg_struct );
        zero.type = IDL_TYP_INT;
        zero.value.i = 0;
        result = &zero;
        return result;
    }

    p1 = ( float * ) arg_struct[0].value.arr->data;
    p2 = ( float * ) arg_struct[1].value.arr->data;

    /* the dimensions of gx and gy must be equal */
    nxg = arg_struct[3].value.arr->dim[0];
    nyg = arg_struct[3].value.arr->dim[1];
    if( arg_struct[3].value.arr->dim[0] != nxg
            || arg_struct[3].value.arr->dim[1] != nyg ) {
        IDL_Message( IDL_M_NAMED_GENERIC, IDL_MSG_RET,
                     "3rd and 4th arguments (grid points) have unequal dimensions" );
        IDL_EzCallCleanup( argc, argv, arg_struct );
        zero.type = IDL_TYP_INT;
        zero.value.i = 0;
        result = &zero;
        return result;
    }

    gx = ( IDL_LONG * ) arg_struct[2].value.arr->data;
    gy = ( IDL_LONG * ) arg_struct[3].value.arr->data;

    dx = arg_struct[4].value.l;
    dy = arg_struct[5].value.l;
    gwid = arg_struct[6].value.f;
    stretch_clip = arg_struct[7].value.l;
    if( stretch_clip < 2 )
        stretch_clip = 2;
    stretch_clip--;

    /* prepare output symbol: It has the same dimensions as <gx> with an */
    /* extra dimension of 2 prepended. */
    dims[0] = 2;
    memcpy( dims + 1, arg_struct[2].value.arr->dim, 2 * sizeof( IDL_MEMINT ) );
    out = ( float * ) IDL_MakeTempArray( IDL_TYP_FLOAT, 3, dims,
                                         IDL_ARR_INI_ZERO, &result );

    /* prepare the gaussian kernel */
    gwx = ( float * ) malloc( ( nx + ny ) * sizeof( float ) );
    if( !gwx ) {
        IDL_Message( IDL_M_NAMED_GENERIC, IDL_MSG_RET,
                     "Could not allocate memory for the gaussian kernel" );
        IDL_EzCallCleanup( argc, argv, arg_struct );
        zero.type = IDL_TYP_INT;
        zero.value.i = 0;
        result = &zero;
        return result;
    }
    gwy = gwx + nx;
    nc = nxg * nyg;       /* the number of subimages */
    dx2 = dx / 2;
    dy2 = dy / 2;
    badmatch = 0;

    // i1 = 0;
    //i2 = 0;
    //#pragma omp parallel default(shared) private(i1, i2, j1, j2, xoffset, yoffset) firstprivate(gwx, gwy) num_threads(6)
    //{
    while( nc-- ) {     /* loop over all subimages */
        i1 = *gx - dx2;       /* lower x coordinate */
        if( i1 < 0 )          /* outside array? */
            i1 = 0;
        i1++;
        i2 = *gx++ + dx2;     /* upper x coordinate */
        if( i2 > nx )
            i2 = nx;
        j1 = *gy - dy2;       /* lower y coordinate */
        if( j1 < 0 )
            j1 = 0;
        j1++;
        j2 = *gy++ + dy2;     /* upper y coordinate */
        if( j2 > ny )
            j2 = ny;
        xoffset = yoffset = 0;
        i1--; i2--; j1--; j2--;
        gwind0( gwx, gwy, gwid, i1, i2, j1, j2 ); /* get gaussian kernels */
        match_1( p1, p2, i1, i2, j1, j2, nx, ny, gwx, gwy,
                 &xoffset, &yoffset ); /* get offsets */
        *out++ = xoffset;
        *out++ = yoffset;
    }
    //  }
    //if (badmatch)
    //   printf("GRIDMATCH - %d bad matches\n", badmatch);

    free( gwx );
    IDL_EzCallCleanup( argc, argv, arg_struct );
    return result;
}


IDL_VPTR redux::stretch( int argc, IDL_VPTR argv[] )
/* the call is MS = STRETCH( M2, DELTA)
   where M2 is the original array to be destretched, MS is the result, and
   DELTA is a displacement grid as generated by GRIDMATCH */
{
    static IDL_EZ_ARG arg_struct[] = {
        {
            IDL_EZ_DIM_MASK( 2 ), IDL_TYP_B_BASIC, IDL_EZ_ACCESS_R, IDL_TYP_UNDEF,
            0, 0
        },           /* <image>: 2D read-only array */
        {
            IDL_EZ_DIM_MASK( 3 ), IDL_TYP_B_SIMPLE, IDL_EZ_ACCESS_R, IDL_TYP_FLOAT,
            0, 0
        } /* <grid> 3D array */
    };
    IDL_VPTR  result_sym;
    Pointer   out, base, bb;

    int   iq, type, n, m, nxg, nyg, jy, j1, j2, nm1, nm2, mm1, mm2;
    int   nxgm, nygm, ix, iy, i1, i2, i3, i4, j3, j4, jx;
    float xd, yd, xinc, yinc, y, x, xs, dy, dy1;
    float dx, dx1, dx0, dx2, dx3, dx4, fn, fm, xq;
    float w1, w2, w3, w4, xl, yl, c1, c2, c3, c4, b1, b2, b3, b4;
    float *jpbase, *jbase, *xgbase;

    if( argc != 2 )
        IDL_Message( IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                     "STRETCH needs 2 arguments:\nSTRETCH(image, grid)" );
    IDL_EzCall( argc, argv, arg_struct );

    type = argv[0]->type;
    out.f = ( float * ) IDL_MakeTempArray( type, 2, argv[0]->value.arr->dim,
                                           IDL_BARR_INI_NOP, &result_sym );
    base.f = ( float * ) argv[0]->value.arr->data;

    n = argv[0]->value.arr->dim[0];
    m = argv[0]->value.arr->dim[1];
    nm1 = n - 1;
    mm1 = m - 1;
    nm2 = n - 2;
    mm2 = m - 2;
    fn = n;
    fm = m;

    xgbase = ( float * ) arg_struct[1].value.arr->data;
    nxg = arg_struct[1].value.arr->dim[1];
    nyg = arg_struct[1].value.arr->dim[2];
    nxgm = nxg - 1;
    nygm = nyg - 1;
    /* linearly interpolate the displacement grid values over array */
    /* similar to regrid3 in inner part */
    xd = ( float ) n / nxg;
    xinc = 1.0 / xd;
    xs = xinc + ( xd - 1.0 ) / ( 2.0 * xd );
    yd = ( float ) m / nyg;
    yinc = 1.0 / yd;
    y = yinc + ( yd - 1.0 ) / ( 2.0 * yd );
    for( iy = 0; iy < m; iy++ ) {
        x = xs;
        jy = y;
        dy = y - jy;
        dy1 = 1.0 - dy;
        if( jy < 1 )
            j1 = j2 = 0;
        else if( jy >= nyg )
            j1 = j2 = nygm;
        else {
            j1 = jy - 1;
            j2 = j1 + 1;
        }
        jbase  = xgbase + j1 * 2 * nxg;
        jpbase = xgbase + j2 * 2 * nxg;
        for( ix = 0; ix < n; ix++ ) {
            jx = x;
            dx = x - jx;
            dx1 = 1.0 - dx;
            if( jx < 1 )
                i1 = i2 = 0;
            else if( jx >= nxg )
                i1 = i2 = nxgm;
            else {
                i1 = jx - 1;
                i2 = i1 + 1;
            }
            w1 = dy1 * dx1;
            w2 = dy1 * dx;
            w3 = dy * dx1;
            w4 = dy * dx;
            i1 = 2 * i1;
            i2 = 2 * i2;
            xl = w1 * *( jbase + i1 ) + w2 * *( jbase + i2 ) + w3 * *( jpbase + i1 )
                 + w4 * *( jpbase + i2 ) + ( float ) ix;
            i1 += 1;  i2 += 1;
            yl = w1 * *( jbase + i1 ) + w2 * *( jbase + i2 ) + w3 * *( jpbase + i1 )
                 + w4 * *( jpbase + i2 ) + ( float ) iy;

            /* xl, yl is the place, now do a cubic interpolation for value */

            i2 = xl;
            j2 = yl;
            if( i2 >= 1 && i2 < nm2 ) { /* normal interior */
                i1 = i2 - 1;
                i3 = i2 + 1;
                i4 = i2 + 2;
                dx0 = xl - i2;
            } else {          /* edge cases */
                i2 = min( i2, nm1 );
                i2 = max( i2, 0 );
                xq = min<float>( xl, fn - 1.0 );
                xq = max<float>( xq, 0 );
                dx0 = xq - i2;
                i1 = min( i2 - 1, nm1 );
                i1 = max( i1, 0 );
                i3 = min( i2 + 1, nm1 );
                i3 = max( i3, 0 );
                i4 = min( i2 + 2, nm1 );
                i1 = max( i4, 0 );
            }
            dx1 = 1.0 - dx0;
            dx2 = -dx0 * 0.5;
            dx3 = dx0 * dx2;
            dx4 = 3.0 * dx0 * dx3;
            c1 = dx2 * dx1 * dx1;
            c2 = 1.0 - dx4 + 5.0 * dx3;
            c3 = dx4 - ( dx2 + 4.0 * dx3 );
            c4 = dx3 * dx1;
            if( j2 >= 1 && j2 < mm2 ) { /* normal interior */
                j1 = j2 - 1;
                j3 = j2 + 1;
                j4 = j2 + 2;
                dx0 = yl - j2;
            } else {          /* edge cases */
                j2 = min( j2, mm1 );
                j2 = max( j2, 0 );
                xq = min<float>( yl, fm - 1.0 );
                xq = max<float>( xq, 0 );
                dx0 = xq - j2;
                j1 = min( j2 - 1, mm1 );
                j1 = max( j1, 0 );
                j3 = min( j2 + 1, mm1 );
                j3 = max( j3, 0 );
                j4 = min( j2 + 2, mm1 );
                j1 = max( j4, 0 );
            }
            dx1 = 1.0 - dx0;
            dx2 = -dx0 * 0.5;
            dx3 = dx0 * dx2;
            dx4 = 3.0 * dx0 * dx3;
            b1 = dx2 * dx1 * dx1;
            b2 = 1.0 - dx4 + 5.0 * dx3;
            b3 = dx4 - ( dx2 + 4.0 * dx3 );
            b4 = dx3 * dx1;
            /* low corner offset */
            iq = i1 + j1 * n;
            i2 = i2 - i1;
            i3 = i3 - i1;
            i4 = i4 - i1;
            j4 = ( j4 - j3 ) * n;
            j3 = ( j3 - j2 ) * n;
            j2 = ( j2 - j1 ) * n;
            switch( type ) {
                case IDL_TYP_BYTE:
                    bb.b = base.b + iq;
                    xq = b1 * ( c1 * ( float ) * ( bb.b ) + c2 * ( float ) * ( bb.b + i2 )
                                + c3 * ( float ) * ( bb.b + i3 ) + c4 * ( float ) * ( bb.b + i4 ) );
                    bb.b += j2;
                    xq += b2 * ( c1 * *( bb.b ) + c2 * *( bb.b + i2 )
                                 + c3 * *( bb.b + i3 ) + c4 * *( bb.b + i4 ) );
                    bb.b += j3;
                    xq += b3 * ( c1 * *( bb.b ) + c2 * *( bb.b + i2 )
                                 + c3 * *( bb.b + i3 ) + c4 * *( bb.b + i4 ) );
                    bb.b += j4;
                    xq += b4 * ( c1 * *( bb.b ) + c2 * *( bb.b + i2 )
                                 + c3 * *( bb.b + i3 ) + c4 * *( bb.b + i4 ) );
                    *out.b++ = xq;
                    break;
                case IDL_TYP_INT:
                    bb.w = base.w + iq;
                    xq = b1 * ( c1 * *( bb.w ) + c2 * *( bb.w + i2 )
                                + c3 * *( bb.w + i3 ) + c4 * *( bb.w + i4 ) );
                    bb.w += j2;
                    xq += b2 * ( c1 * *( bb.w ) + c2 * *( bb.w + i2 )
                                 + c3 * *( bb.w + i3 ) + c4 * *( bb.w + i4 ) );
                    bb.w += j3;
                    xq += b3 * ( c1 * *( bb.w ) + c2 * *( bb.w + i2 )
                                 + c3 * *( bb.w + i3 ) + c4 * *( bb.w + i4 ) );
                    bb.w += j4;
                    xq += b4 * ( c1 * *( bb.w ) + c2 * *( bb.w + i2 )
                                 + c3 * *( bb.w + i3 ) + c4 * *( bb.w + i4 ) );
                    *out.w++ = xq;
                    break;
                case IDL_TYP_LONG:
                    bb.l = base.l + iq;
                    xq = b1 * ( c1 * *( bb.l ) + c2 * *( bb.l + i2 )
                                + c3 * *( bb.l + i3 ) + c4 * *( bb.l + i4 ) );
                    bb.l += j2;
                    xq += b2 * ( c1 * *( bb.l ) + c2 * *( bb.l + i2 )
                                 + c3 * *( bb.l + i3 ) + c4 * *( bb.l + i4 ) );
                    bb.l += j3;
                    xq += b3 * ( c1 * *( bb.l ) + c2 * *( bb.l + i2 )
                                 + c3 * *( bb.l + i3 ) + c4 * *( bb.l + i4 ) );
                    bb.l += j4;
                    xq += b4 * ( c1 * *( bb.l ) + c2 * *( bb.l + i2 )
                                 + c3 * *( bb.l + i3 ) + c4 * *( bb.l + i4 ) );
                    *out.l++ = xq;
                    break;
                case IDL_TYP_FLOAT:
                    bb.f = base.f + iq;
                    xq = b1 * ( c1 * *( bb.f ) + c2 * *( bb.f + i2 )
                                + c3 * *( bb.f + i3 ) + c4 * *( bb.f + i4 ) );
                    bb.f += j2;
                    xq += b2 * ( c1 * *( bb.f ) + c2 * *( bb.f + i2 )
                                 + c3 * *( bb.f + i3 ) + c4 * *( bb.f + i4 ) );
                    bb.f += j3;
                    xq += b3 * ( c1 * *( bb.f ) + c2 * *( bb.f + i2 )
                                 + c3 * *( bb.f + i3 ) + c4 * *( bb.f + i4 ) );
                    bb.f += j4;
                    xq += b4 * ( c1 * *( bb.f ) + c2 * *( bb.f + i2 )
                                 + c3 * *( bb.f + i3 ) + c4 * *( bb.f + i4 ) );
                    *out.f++ = xq;
                    break;
                case IDL_TYP_DOUBLE:
                    bb.d = base.d + iq;
                    xq = b1 * ( c1 * *( bb.d ) + c2 * *( bb.d + i2 )
                                + c3 * *( bb.d + i3 ) + c4 * *( bb.d + i4 ) );
                    bb.d += j2;
                    xq += b2 * ( c1 * *( bb.d ) + c2 * *( bb.d + i2 )
                                 + c3 * *( bb.d + i3 ) + c4 * *( bb.d + i4 ) );
                    bb.d += j3;
                    xq += b3 * ( c1 * *( bb.d ) + c2 * *( bb.d + i2 )
                                 + c3 * *( bb.d + i3 ) + c4 * *( bb.d + i4 ) );
                    bb.d += j4;
                    xq += b4 * ( c1 * *( bb.d ) + c2 * *( bb.d + i2 )
                                 + c3 * *( bb.d + i3 ) + c4 * *( bb.d + i4 ) );
                    *out.d++ = xq;
                    break;
            }
            x += xinc;
        }
        y += yinc;
    }
    IDL_EzCallCleanup( argc, argv, arg_struct );
    return result_sym;
}


extern "C" {

    int IDL_Load( void ) {

        static IDL_SYSFUN_DEF2 function_addr[] = {
            { {(IDL_VPTR(*)())redux::gridmatch}, "RED_GRIDMATCH", 7, 8, 0 },
            { {(IDL_VPTR(*)())redux::stretch}, "RED_STRETCH", 2, 2, 0 }
        };

        /* Register our routine. The routines must be specified exactly the same as in testmodule.dlm. */
        return IDL_SysRtnAdd( function_addr, TRUE, IDL_CARRAY_ELTS( function_addr ) );

    }

}
