#include <cmath>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <thread>
#include <vector>

#include "mathlib.hpp"
#include "idlutil.hpp"

#include "redux/util/arrayutil.hpp"
#include <redux/util/datautil.hpp>

using namespace redux::util;
using namespace redux;
using namespace std;

// ******************************************************************************************* //

template<typename T>
inline void gwind0(T* __restrict__ gwx, T* __restrict__ gwy, T const gwid, int const nxa,
		   int const nxb, int const nya, int const nyb)
{
  double const wid = gwid*0.6005612;
  
  if (wid > 0) {
    double const xcen = (nxa + nxb)/2;

    double const ycen = (nya + nyb)/2; 
    for (int i = nxa; i <= nxb; ++i) {
      double const xq = (i - xcen)/wid;
      gwx[i] = exp(-(xq*xq));
    } 
    for (int i = nya; i <= nyb; ++i) {
      double const xq = (i - ycen)/wid;
      gwy[i] = exp(-(xq * xq));
    } 
  } else {
    for (int i = nxa; i <= nxb; ++i) 
      gwx[i] = 1.0;
    for (int i = nya; i <= nyb; ++i)
      gwy[i] = 1.0; 
  }
}

// ******************************************************************************************* //
template<typename T, typename U>
inline double averag(const T* const __restrict__ p, int const nxa, int const nxb, int const nya,
	     int const nyb, int const nxs, int const nys, int const idx,
	     int idy, const U* const __restrict__ gx, const U* const __restrict__ gy)
/* finds weighted average of array m over the block defined */
{
  double sum=0.0, sumg=0.0, sumgx=0.0;

  /* fix limits so sum doesn't run off edge of image */
  int const nxc = (nxa + idx < 0)? -idx: nxa;
  int const nyc = (nya + idy < 0)? -idy: nya;
  int const nxd = (nxb + idx > nxs)? nxs - idx: nxb;
  int const nyd = (nyb + idy > nys)? nys - idy: nyb;

  for (int i = nxc; i < nxd; i++)
    sumgx += gx[i];
  /* weighted sum in window */
  for (int j = nyc; j < nyd; j++) {
    double sumx = 0.0;
    int const jj = idx + nxs*(j + idy);
    for (int i = nxc; i < nxd; i++)
      sumx += gx[i]*p[i + jj];
    sum += gy[j]*sumx;
    sumg += gy[j]*sumgx;
  } /* end of j loop */
  sum /= sumg;
  return sum;
}

// ******************************************************************************************* //
template<typename T, typename U>
inline void unbias(const T* const m1, const U* const m2, int nxa, int nxb,
		   int nya, int nyb, int nxs, int nys,
		   U *gx, U *gy, U &av1, U &av2, U &cx,
		   U &cy, U &cxx, U &cxy, U &cyy,
		   int idelx, int idely)
{ 
  /*  find weighted means of m1 & m2 over the window 
      sets up quadratic fit to average of m2 as a fcn. of offsets */
  av1 = averag(m1, nxa, nxb, nya, nyb, nxs, nys, 0, 0, gx, gy); 
  double const t0 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely, gx, gy); 
  double const t1 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx + 1, idely, gx, gy); 
  double const t2 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx - 1, idely, gx, gy); 
  double const t3 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely + 1, gx, gy); 
  double const t4 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely - 1, gx, gy); 
  double const t5 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx + 1, idely + 1, gx, gy); 
  av2 = t0; 
  cx = 0.5*(t1 - t2); 
  cy = 0.5*(t3 - t4); 
  cxx = 0.5*(t1 - 2*t0 + t2); 
  cyy = 0.5*(t3 - 2*t0 + t4); 
  cxy = t5 + t0 - t1 - t3;
}

// ******************************************************************************************* //

template<typename T>
inline int getmin(T *p, T &x0, T &y0)
{
  double	f11, f12, f13, f21, f22, f23, f31, f32, f33;
  double	fx, fy, t, fxx, fyy, fxy;

  /* find the min, p points to a 3x3 array */
  f11 = *p++;	f21 = *p++;	f31 = *p++;
  f12 = *p++;	f22 = *p++;	f32 = *p++;
  f13 = *p++;	f23 = *p++;	f33 = *p++;

  fx = 0.5 * ( f32 - f12 );	fy = 0.5 * ( f23 - f21 );
  t = 2.* ( f22 );
  fxx =  f32 + f12 - t;	fyy = f23 + f21 - t;
  /* find in which quadrant the minimum lies */
  if (f33 < f11) {
    if (f33 < f31) {
      if (f33 < f13)
	fxy = f33+f22-f32-f23;
      else
	fxy = f23+f12-f22-f13;
    } else {
      if (f31 < f13)
	fxy = f32+f21-f31-f22;
      else
	fxy = f23+f12-f22-f13;
    }
  } else { 
    if (f11 < f31) {
      if (f11 < f13) 
	fxy = f22+f11-f21-f12;
      else
	fxy = f23+f12-f22-f13;
    } else {
      if (f31 < f13)
	fxy = f32+f21-f31-f22;
      else
	fxy = f23+f12-f22-f13;
    }
  }
  t = -1./(fxx *fyy - fxy *fxy);
  x0 = t * (fx * fyy - fy * fxy);
  y0 = t * (fy * fxx - fx * fxy);
  if (std::abs(x0) >= 0.75 || std::abs(y0) >= 0.75) {
    x0 = -fx/fxx;
    y0 = -fy/fyy;
  }
  return 1;
}

// ******************************************************************************************* //

template <typename T, typename U>
inline double resid(T *m1, T *m2, int idx, int idy, int nxa,
	    int nxb, int nya, int nyb, int nxs,
	    int nys, int ndmx, U *gx, U *gy, U bs)
{
  int	nxc, nxd, nyc, nyd, nx, ny;
  U 	*p1, *p2, *ps;
  double	sum, sumx, t, ndmx2;
  int	i, j;
  double   sumg;
  static int	mxc, mxd, myc, myd;
  static double	gsum;

  /*set up limits */
  nxc = nxa;
  if (nxc + idx < 0)
    nxc = -idx;
  nyc = nya;
  if (nyc + idy < 0)
    nyc = - idy;
  nxd = nxb;
  if (nxd + idx >= nxs)
    nxd = nxs - idx - 1;
  nyd = nyb;
  if (nyd + idy >= nys)
    nyd = nys - idy - 1;
  sum = sumg = 0.0;

  nx = nxd - nxc +1;
  p2 = gy + nyc;
  ps = gx + nxc;

  if (nxc != mxc || nxd != mxd || nyc != myc || nyd != myd) {
    /* sum gaussians over rectangle to get normalization */
    /* (only if limits change)*/
    j = nyd -nyc + 1;
    if (j <= 0 || nxd - nxc + 1 <= 0)
      return -1;		/* added 19feb95 LS */
    while (j) {
      i = nx;
      p1 = ps;
      while (i) {
	sumg += (*p1++) * (*p2);
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

  m1 += nyc*nxs + nxc;
  m2 += (nyc + idy)*nxs + nxc + idx;
  ny = nxs - nx;		/*residual increment after inner loop */
  p2 = gy + nyc;
  
  /* now the loop to compute the residual */
  j = nyd - nyc +1;
  ndmx2 = ndmx*ndmx;
  while (j) {
    i = nx;
    p1 = ps;
    sumx = 0.0;
    while (i) {
      t = *m1++ - *m2++;
      t = t + bs;
      t = t*t;
      t = std::min<double>(t, ndmx2);
      sumx += (*p1++) * t;
      i--;
    }
    sum += (*p2++) * sumx;
    m1 += ny;
    m2 += ny;
    j--;
  }
  /*return normalized residual */
  sum /= sumg;
  return sum;
}

// ******************************************************************************************* //

template<typename T, typename U>
inline void match_1(T *p1, T *p2, int nxa, int nxb,
		    int nya, int nyb, int nx, int ny,
		    U *gwx, U *gwy, U &xoffset, U &yoffset, int const& stretch_clip, int &badmatch)
{
  int idelx, idely, i, j, k, ndmx = 1000, done[9]={};
  int	di, dj, in, jn, iter, dd, badflag = 0;
  double	av1, av2, cx, cy, cxx, cxy, cyy, avdif, t, res[9], buf[9]={}, t1, t2;
  static int	itmax = 40;
  
  idelx = rint(xoffset);
  idely = rint(yoffset); 
  unbias(p1, p2, nxa, nxb, nya, nyb, nx, ny, gwx, gwy, av1, av2, cx, cy,
	 cxx, cxy, cyy, idelx, idely);
  /* look at a 3x3 matrix of residuals centered at 0 offset, find the location
     of the minimum, if not at center, then look at new 3x3 centered
     on the edge minimum; repeat until min at center */
  iter = itmax;
  badflag = 0;
  while (iter--) {
    for (k = 0; k < 9; k++) {
      if (done[k] == 0) {
	i = idelx + (k % 3) - 1;
	j = idely + (k / 3) - 1;
	avdif = av2 +  i*cx + j*cy + i*i*cxx + i*j*cxy + j*j*cyy - av1;
	res[k] = resid(p1, p2, i, j, nxa, nxb, nya, nyb, nx, ny, ndmx,
		       gwx, gwy, avdif);
      }
    }
    t = res[0];
    i = 0;
    for (k = 1; k < 9; k++) 
      if (res[k] < t) {
	t = res[k];
	i = k;
      }
    if (t < 0) {		/* added LS 19feb95 */
      printf("match - ran out of data at edge\n");
      badflag = 1;
      break;
    }
    idelx += (i % 3) - 1;
    idely += (i / 3) - 1;
    /* check if we have gone too far */
    if (std::abs(idelx) > stretch_clip || std::abs(idely) > stretch_clip) {
      badflag++;
      break;
    }
    if (i == 4)
      break;			/* done if it is the center one */
    /* not in center, shuffle what we have to put the edge min in center */
    di = (i % 3) - 1;
    dj = (i / 3) - 1;
    dd = dj * 3 + di;
    for (k = 0; k < 9; k++) {
      in = k%3 + di;
      jn = k/3 + dj;
      if (in >= 0 && jn >= 0 && in < 3 && jn < 3) { /* in range */
	done[k] = 1;
	buf[k] = res[k + dd];
      } else 
	done[k] = 0;		/* not in range, mark undone */
    }
    for (k = 0; k < 9; k++)
      res[k] = buf[k];		/* put back in res array */
  } /* end of iter while */
  /* done or reached itmax, which ? */
  if (iter <= 0) {
    badflag++;
  }
  if (badflag) {
    badmatch++;
    xoffset = yoffset = 0;
    return;
  }
				/* must have been OK so far */
  getmin(res, t1, t2);
  xoffset = idelx + t1;
  yoffset = idely + t2;
}

// ******************************************************************************************* //

template<typename T, typename U>
inline void do_one(int const ii, T* __restrict__ gwx, T* __restrict__ gwy, T gwid, 
		   const U* const __restrict__ p1, const U* const __restrict__ p2, int const stretch_clip, int const dx2, int const dy2,
		   int const nyg, int const nxg, int const nx, int const ny, U* __restrict__ out,
		   const int* const __restrict__ gx, const int* const __restrict__ gy)
{
  int badmatch = 0;

  int  const i1 = std::max<int>(0,gx[ii] - dx2) ;
  int  const i2 = std::min<int>(nx,gx[ii] + dx2) - 1;

  int  const j1 = std::max<int>(gy[ii] - dy2,0);
  int  const j2 = std::min<int>(gy[ii] + dy2, ny) - 1;

  U xoffset = 0, yoffset = 0;
  
  gwind0(gwx, gwy, gwid, i1, i2, j1, j2); /* get gaussian kernels */
  
  match_1(p1, p2, i1, i2, j1, j2, nx, ny, gwx, gwy,
	  xoffset, yoffset, stretch_clip, badmatch); /* get offsets */

  out[2*ii]   = xoffset;
  out[2*ii+1] = yoffset; 
}

// ******************************************************************************************* //

template<typename T, typename U>
void gridmatch(int const ny, int const nx, int const nyg, int const nxg, const T* const __restrict__ p1,
	       const T* const __restrict__ p2, const int* const gy, const int* const gx, int const dy,
	       int const dx, U const gwid, int stretch_clip, U* __restrict__ out, int const nthreads = 6)
{
  
  /* the call is offsets = gridmatch(m1,m2,gx,gy,dx,dy,gwid,stretch_clip)
	where	m1 = reference input image
	m2 = image to compare with m1, m1 and m2 must be same size
	gx = array of x gridpoints
	gy = array of y gridpoints, gx and gy must have same size
	dx and dy are the window size, and gwid is the gaussian mask width
        stretch_clip = maximum allowed displacement before a bad match
          is declared

  Authors:  Richard A. Shine (original)
            Louis H. Strous (port from ANA to IDL)
            Lockheed-Martin Advanced Technology Center, Palo Alto, CA, USA

	    Parallel version adapted by J. de la Cruz Rodriguez (ISP-SU, 2021)

*/

  //double t0 = getTime();
  // 
  double	        /*xoffset, yoffset,*/ *gwx, *gwy;
  //int      i1, i2, j1, j2;
  //int  const dims[3] = {2,nxg,nxg};

  if (stretch_clip < 2)
    stretch_clip = 2;
  stretch_clip--;

  /* prepare the gaussian kernel */
  gwx = (double *) malloc((nx + ny)*sizeof(double));
  gwy = gwx + nx;
  
  int const nc = nxg*nyg;			/* the number of subimages */
  int const dx2 = dx/2;
  int const dy2 = dy/2;
  //int badmatch = 0;
  int ii;
#pragma omp parallel default(shared) firstprivate(ii) num_threads(nthreads)
  {
#pragma omp for schedule(static) 
    for( ii=0; ii<nc;++ii) {	
      do_one<double,T>(ii, gwx, gwy, gwid,  p1, p2, stretch_clip, dx2, dy2, nyg, nxg, nx, ny, out, gx, gy);
    }
  }
 
  
  free(gwx);

}


// ******************************************************************************************* //
 
template<typename T>
T* stretch_matrix(int const ny, int const nx,  int const npy, int const npx,
		  const T* const __restrict__ gr, int const nthreads = 2)
{  
  int   jy, j1, j2;
  int   ix, iy, i1, i2,  jx;
  double        y, x,  dy, dy1;
  double	dx, dx1;
  double	w1, w2, w3, w4;//, xl, yl;
  const T 	*jpbase, *jbase;

  int const n = nx;
  int const m = ny;

  const T* const __restrict__ xgbase = gr;
  int const nxg = npx;
  int const nyg = npy;
  int const nxgm = nxg - 1;
  int const nygm = nyg - 1;

  /* linearly interpolate the displacement grid values over array */
  /* similar to regrid3 in inner part */
  double const xd = (double) n/nxg; // number of pixels in the image per grid cell
  double const xinc = 1.0/xd; // maps the size of a pixel of the image in the matrix
  double const xs   = xinc + (xd - 1.0)/(2.0*xd);
  double const yd = (double) m/nyg;
  double const yinc = 1.0/yd;
  double const y0 = yinc + (yd - 1.0)/(2.0*yd);

  T* __restrict__ out_x = new T [2*nx*ny]();
  T* __restrict__ out_y = out_x+nx*ny;

  
#pragma omp parallel default(shared) firstprivate(ix, iy, i1, i2, jx, jy, j1, j2, w1, w2, w3, w4, jpbase, jbase, x, y, dy, dy1, dx, dx1, xgbase) num_threads(nthreads)
  {
#pragma omp for schedule(static) 
    for (iy = 0; iy < m; ++iy) {
      y = y0 + iy * yinc;
      x = xs;
      jy = y;
      dy = y - jy;
      dy1 = 1.0 - dy;
      
      if (jy < 1)
	j1 = j2 = 0;
      else if (jy >= nyg)
	j1 = j2 = nygm;
      else {
	j1 = jy - 1;
	j2 = j1 + 1;
      }
      
      jbase  = xgbase + j1 * 2 * nxg;
      jpbase = xgbase + j2 * 2 * nxg;
      
      for (ix = 0; ix < n; ++ix) {
	jx = x;
	dx = x - jx;
	dx1 = 1.0 - dx;
	
	if (jx < 1)
	  i1 = i2 = 0;
	else if (jx >= nxg)
	  i1 = i2 = nxgm;
	else {
	  i1 = jx - 1;
	  i2 = i1 + 1;
	}
	
	w1 = dy1*dx1;
	w2 = dy1*dx;
	w3 = dy*dx1;
	w4 = dy*dx;
	
	i1 = 2*i1;
	i2 = 2*i2;
	
	out_x[iy*nx+ix] = w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1)
	  + w4 * *(jpbase+i2) + (double)ix;
	
	i1 += 1;
	i2 += 1;
	
	out_y[iy*nx+ix] = w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1)
	  + w4 * *(jpbase+i2) + (double)iy;
	
	x += xinc;
      }
    } // iy    
  } // pragma

  return out_x;
}

// ******************************************************************************************* //

template<typename T>
T* stretch(int const ny, int const nx, const T* const __restrict__ im, int const npy, int const npx,
	   const T* const __restrict__ gr, int const nthreads = 4)
{  
  using fp = double;
  
  T* __restrict__ res = new T [nx*ny]();

  int   jy(0), j1(0), j2(0);
  int   ix(0), iy(0), i1(0), i2(0), jx(0); //, intx, inty;
  fp    y(0), dy(0), dy1(0), dyp(0), dyp1(0);
  fp	x(0), dx(0), dx1(0), dxp(0), dxp1(0);
  fp	w1(0), w2(0), w3(0), w4(0), xl(0), yl(0);
  const T 	*jpbase(nullptr), *jbase(nullptr);

  int const n = nx;
  int const m = ny;

  const T* const __restrict__ xgbase = gr;
  int const nxg = npx;
  int const nyg = npy;
  int const nxgm = nxg - 1;
  int const nygm = nyg - 1;
  int const nx1 = nx-1;
  int const ny1 = ny-1;
  
  /* linearly interpolate the displacement grid values over array */
  /* similar to regrid3 in inner part */
  fp const xd = (double) n/nxg; // number of pixels in the image per grid cell
  fp const xinc = 1.0/xd; // maps the size of a pixel of the image in the matrix
  fp const xs   = xinc + (xd - 1.0)/(2.0*xd);
  fp const yd = (double) m/nyg;
  fp const yinc = 1.0/yd;
  fp const y0 = yinc + (yd - 1.0)/(2.0*yd);
  y = y0;
  
#pragma omp parallel default(shared) firstprivate(ix, iy, i1, i2, jx, jy, j1, j2, w1, w2, w3, w4, jpbase, jbase, x, y, dy, dy1, dx, dx1, xgbase, xl, yl, dyp, dyp1, dxp, dxp1) num_threads(nthreads)
    {
#pragma omp for schedule(static) 
    for (iy = 0; iy < m; ++iy) {
      y = y0 + iy * yinc;
      x = xs;
      jy = y;
      dy = y - jy;
      dy1 = 1.0 - dy;
      
      if (jy < 1)
        j1 = j2 = 0;
      else if (jy >= nyg)
        j1 = j2 = nygm;
      else {
        j1 = jy - 1;
        j2 = j1 + 1;
      }
      
      jbase  = xgbase + j1 * 2 * nxg;
      jpbase = xgbase + j2 * 2 * nxg;
      
      for (ix = 0; ix < n; ++ix) {
        jx = x;
        dx = x - jx;
        dx1 = 1.0 - dx;
        
        if (jx < 1)
            i1 = i2 = 0;
        else if (jx >= nxg)
            i1 = i2 = nxgm;
        else {
            i1 = jx - 1;
            i2 = i1 + 1;
        }

        // --- bilinear interpolation of the coordinates for a pixel in the image --- //
        
        w1 = dy1*dx1;
        w2 = dy1*dx;
        w3 = dy*dx1;
        w4 = dy*dx;
        
        i1 = 2*i1;
        i2 = 2*i2;
        
        xl = std::max<double>(std::min<double>(nx1, w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1) + w4 * *(jpbase+i2) + (double)ix), 0);
        
        i1 += 1;
        i2 += 1;
        
        yl =  std::max<double>(std::min<double>(ny1, w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1) + w4 * *(jpbase+i2) + (double)iy), 0);
        

        // --- Now bilinear interpolation of the actual image --- //
        
        // -- bracket which pixels of the image must be taken --- //
        
        i2 = std::min<int>(std::max<int>(int(xl), 0), nx1-1);
        j2 = std::min<int>(std::max<int>(int(yl), 0), ny1-1);
        
        dxp  = xl - i2;
        dxp1 = 1.0 - dxp; 
        dyp  = yl - j2;
        dyp1 = 1.0  - dyp;

        w1 = *(im+j2*nx+i2);
        w2 = *(im+j2*nx+i2+1);
        w3 = *(im+(j2+1)*nx+i2);
        w4 = *(im+(j2+1)*nx+i2+1);
        
        res[iy*nx+ix] = w1*(dxp1*dyp1) + w2*(dxp*dyp1) + w3*(dxp1*dyp) + w4*(dxp*dyp);
        
        x += xinc;
      }
      // y += yinc;
    } // iy    
} // pragma
  
  return res;
}

// ******************************************************************************************* //

typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    IDL_LONG nthreads;

} KW_STRETCH;


// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR kw_stretch_pars[] = {
    IDL_KW_FAST_SCAN,
    { "HELP",             IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*)IDL_KW_OFFSETOF2(KW_STRETCH,help) },
    { "NTHREADS",         IDL_TYP_LONG,  1, 0,                      0, (char*)IDL_KW_OFFSETOF2(KW_STRETCH,nthreads) },
    { NULL }
};


string stretch_info( int lvl ) {
    string ret = "RDX_CSTRETCH";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"      ");          // newline if lvl>1
        ret += "   Syntax:   out = rdx_cstrectch( image, distorsions, /KEYWORDS )\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      NTHREADS            Number of threads (default=8).\n";
        }
    } else ret += "\n";
    return ret;
}


IDL_VPTR stretch_wrap( int nArg, IDL_VPTR argv[], char* argk ) {
    
    KW_STRETCH kw;    
    kw.nthreads = 8;
    int nPlainArgs = IDL_KWProcessByOffset( nArg, argv, argk, kw_stretch_pars, nullptr, 1, &kw );

    if( kw.help || nPlainArgs < 2 ) {
        printMessage( stretch_info(2) );
        return IDL_GettmpInt(0);
    }

    IDL_VPTR arr1 = argv[0];
    IDL_VPTR arr2   = argv[1];

    IDL_ENSURE_SIMPLE( arr1 );
    IDL_ENSURE_ARRAY( arr1 );
    IDL_ENSURE_SIMPLE( arr2 );
    IDL_ENSURE_ARRAY( arr2 );

    if( arr1->value.arr->n_dim != 2 || arr2->value.arr->n_dim != 3 ) {
        printf("rdx_cstrectch expects one 2D array, and one 3D array as input.");
    }
    
    using fp_t = float; // if you change this, make sure you also change IDL_TYP_FLOAT to double

    int const nx    = arr1->value.arr->dim[0];
    int const ny    = arr1->value.arr->dim[1];

    int const nV    = arr2->value.arr->dim[0];
    int const npx   = arr2->value.arr->dim[1];
    int const npy   = arr2->value.arr->dim[2];

    if( nV != 2 ) {
        printf("rdx_cstrectch expects the 2:nd input (3D array) to have size=2 (x/y) in the the fast dimension.");
    }
    
    int const nthreads = std::max(std::min<int>( kw.nthreads, std::thread::hardware_concurrency() ),1);

    IDL_LONG64 dims[2] = {nx, ny};

    shared_ptr<fp_t> data1 = castOrCopy<fp_t>(arr1);
    shared_ptr<fp_t> data2 = castOrCopy<fp_t>(arr2);
    
    fp_t* res = stretch<fp_t>( ny, nx, data1.get(), npy, npx, data2.get(), nthreads );

    return IDL_ImportArray( 2, dims, IDL_TYP_FLOAT, (UCHAR*)res, redux::util::castAndDelete<fp_t>, 0);

}

// ******************************************************************************************* //

typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    IDL_LONG nthreads;

} KW_DSGRIDNEST;


// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR kw_dsgridnest_pars[] = {
    IDL_KW_FAST_SCAN,
    { "HELP",             IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_DSGRIDNEST,help) },
    { "NTHREADS",         IDL_TYP_LONG,  1, 0,                      0, (char*) IDL_KW_OFFSETOF2(KW_DSGRIDNEST,nthreads) },
    { NULL }
};


string dsgridnest_info( int lvl ) {
    string ret = "RDX_CDSGRIDNEST";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"      ");          // newline if lvl>1
        ret += "   Syntax:   out = rdx_cdsgridnest( img1, img2, tiles, clips, /KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      NTHREADS            Number of threads (default=8).\n";
        }
    } else ret += "\n";
    return ret;
}


IDL_VPTR dsgridnest( int nArg, IDL_VPTR argv[], char* argk ) {

    KW_DSGRIDNEST kw;    
    kw.nthreads = 8;
    int nPlainArgs = IDL_KWProcessByOffset( nArg, argv, argk, kw_dsgridnest_pars, ( IDL_VPTR* )0, 1, &kw );

    if( kw.help || nPlainArgs < 4 ) {
        printMessage( dsgridnest_info(2) );
        return IDL_GettmpInt(0);
    }

    IDL_VPTR imgVar1 = argv[0];
    IDL_VPTR imgVar2 = argv[1];
    IDL_VPTR tileVar = argv[2];
    IDL_VPTR clipVar = argv[3];

    IDL_ENSURE_SIMPLE( imgVar1 );
    IDL_ENSURE_ARRAY( imgVar1 );

    IDL_ENSURE_SIMPLE( imgVar2 );
    IDL_ENSURE_ARRAY( imgVar2 );

    IDL_ENSURE_SIMPLE( tileVar );
    IDL_ENSURE_ARRAY( tileVar );
    IDL_ENSURE_SIMPLE( clipVar );
    IDL_ENSURE_ARRAY( clipVar );


    if( imgVar1->value.arr->n_dim != 2 || imgVar2->value.arr->n_dim != 2 || tileVar->value.arr->n_dim != 1 || clipVar->value.arr->n_dim != 1 ) {
        printf("rdx_cstrectch expects two 2D images, and two vectors (tile and clip) as input.");
    }
    
    using fp_t = double; // if changed, need to change also IDL types
    int const nthreads = std::max( std::min<int>( kw.nthreads, std::thread::hardware_concurrency() ), 1 );
  
    // --- Define variables and image sizes --- //

    int const nx    = imgVar1->value.arr->dim[0];
    int const ny    = imgVar1->value.arr->dim[1];
    int const nx1   = imgVar2->value.arr->dim[0];
    int const ny1   = imgVar2->value.arr->dim[1];
    int const nTile = tileVar->value.arr->dim[0];
    int const nClip = clipVar->value.arr->dim[0];
  
    IDL_VPTR  result;
  
  
    // --- Check dimensions --- //
  
    if( (nx != nx1) || (ny != ny1) ) {
        fprintf(stdout, "dlm::dsgridnest: ERROR, the images have different dimensions [%d,%d] != [%d,%d]\n", nx,ny,nx1,ny1);
        IDL_LONG64 dims[3] = {2,1,1};
        IDL_MakeTempArray(IDL_TYP_DOUBLE, 3, dims, IDL_ARR_INI_ZERO, &result);
        return result;
    }
  
    if( nTile != nClip ){
        fprintf(stdout, "dlm::dsgridnest: ERROR, the tile and clip arrays have different sizes: %d != %d\n", nTile, nClip);
        IDL_LONG64 dims[3] = {2,1,1};
        IDL_MakeTempArray(IDL_TYP_DOUBLE, 3, dims, IDL_ARR_INI_ZERO, &result);
        return result;
    }

    shared_ptr<fp_t> imgData1 = castOrCopy<fp_t>(imgVar1);
    shared_ptr<fp_t> imgData2 = castOrCopy<fp_t>(imgVar2);
    shared_ptr<int> tileData = castOrCopy<int>(tileVar);
    shared_ptr<int> clipData = castOrCopy<int>(clipVar);
    const fp_t* const __restrict__ im1 = imgData1.get();
    const fp_t* const __restrict__ im2 = imgData2.get();
    const int* const __restrict__ tiles = tileData.get();
    const int* const __restrict__ clips = clipData.get();

    int nprev = 0;
    std::vector<fp_t> displprev, prev, displ;
    std::vector<int> igx(nTile,0), igy(nTile,0);
  
    for(int k=0; k<nTile; ++k){
        
        int const n            = tiles[k];
        int const stretch_clip = clips[k];
        int const ngw = int((2*nx)/float(n));
        int const nw  = int(1.25f * ngw);
        
        int const ngx = ((nx > ny) ? n : int(float(n*nx)/ny +0.5));
        int const ngy = ((nx > ny) ? int(float(n*ny)/nx +0.5) : n);
        
        igx[k] = ngx;
        igy[k] = ngy;
        
        fp_t const wx = fp_t(nx) / ngx;
        fp_t const wy = fp_t(ny) / ngy;
        
        int* const __restrict__ pgx = new int [ngx*ngy];
        int* const __restrict__ pgy = new int [ngx*ngy];
        
        // --- make subfields grid --- //
        
        for(int jj = 0; jj<ngy; ++jj){
            for(int ii = 0; ii<ngx; ++ii){
                pgy[jj*ngx+ii] = int(jj*wy+wy/2-1);
                pgx[jj*ngx+ii] = int(ii*wx+wx/2-1);
            }
        }

        int const dx = nw;
        int const dy = nw;
        fp_t const gwid = ngw;

        displ.resize(2*ngx*ngy,fp_t(0));
        
        if( k == 0 ){
            gridmatch(ny, nx, ngy, ngx, im1, im2, pgy, pgx, dy, dx, gwid, stretch_clip, &displ[0], nthreads);
        } else {
            if( n != nprev ){
                int const ngtotprev = igx[k-1]*igy[k-1];
                fp_t* const __restrict__ disx = new fp_t [ngtotprev];
                fp_t* const __restrict__ disy = new fp_t [ngtotprev];

                for(int ii=0; ii<ngtotprev; ++ii){
                    disx[ii] = displprev[2*ii];
                    disy[ii] = displprev[2*ii+1];
                }
                
                const fp_t* const __restrict__ fx = redux::congrid<fp_t,fp_t>(igy[k-1],igx[k-1], disx, ngy, ngx);
                const fp_t* const __restrict__ fy = redux::congrid<fp_t,fp_t>(igy[k-1],igx[k-1], disy, ngy, ngx);

                delete [] disx;
                delete [] disy;
                
                int const ngtot = ngx*ngy;
                prev.resize(2*ngx*ngy,0);
                
                for(int ii =0; ii<ngtot; ++ii){
                    prev[2*ii] = fx[ii];
                    prev[2*ii+1] = fy[ii];
                }

                delete [] fx;
                delete [] fy;
                
            } else prev = displprev;

            {
                int const ngtot = 2*ngx*ngy;
                fp_t*  __restrict__ displnew = new fp_t [ngtot];
                
                fp_t* im3 = stretch<fp_t>(ny, nx, im2, ngy, ngx, &prev[0], nthreads);
                gridmatch<fp_t>(ny, nx, ngy, ngx, im1, im3, pgy, pgx, dy, dx, gwid, stretch_clip, displnew, nthreads);
                delete [] im3;
                
                for(int ii=0; ii<ngtot; ++ii) displ[ii] = prev[ii] + displnew[ii];
                delete [] displnew;
            }
                    
            delete [] pgx;
            delete [] pgy;
        }
        if( k < (nTile-1) ){
            nprev = n;
            displprev = displ;
        }


    
    }
  
    // --- Now store result in a real IDL array --- //
  
    {
        IDL_LONG64 dims[3] = { 2, igx[nTile-1], igy[nTile-1] };
        float* __restrict__ out = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT, 3, dims, IDL_ARR_INI_ZERO, &result);
    
        int const ngtot = igy[nTile-1]*igx[nTile-1]*2;
        for(int ii=0; ii<ngtot; ++ii){
            out[ii] = displ[ii];
        }
    }

    //return IDL_ImportArray( 2, dims, IDL_TYP_FLOAT, (UCHAR*)res, redux::util::castAndDelete<fp_t>, 0);

    return result;
    
}

// ******************************************************************************************* //

template<typename T, typename U> inline
U interpolate_one( int const ny, int const nx, const T* const __restrict__ im,
		  double const y, double const x, T const missing)
{
  if((y >= 0) && (y <= (ny-1)) && (x >= 0) && (x <= (nx-1))){

    int const iy = std::min<int>(int(y), ny-2);
    int const ix = std::min<int>(int(x), nx-2);
  
    double const dx  = x - ix;
    double const dy  = y - iy;
    double const dx1 = 1.0 - dx;
    double const dy1 = 1.0 - dy;
    
    return U(*(im+iy*nx+ix)*dx1*dy1 +  *(im+iy*nx+ix+1)*dx*dy1 +
	     *(im+(iy+1)*nx+ix)*dx1*dy + *(im+(iy+1)*nx+ix+1)*dx*dy);
  }else{
    return missing;
  }
}

// ******************************************************************************************* //

template<typename T, typename U> inline
void binterpolate_image_2D(int const ny, int const nx, const T* const __restrict__ im,
			   int const ny1, int const nx1, const T* const __restrict yin,
			   const T* const __restrict__ xin, U* const __restrict__ res,
			   int const nthreads, T const missing)
{

  int ii = 0;
  int const ntot = nx1*ny1;

#pragma omp parallel default(shared) firstprivate(ii) num_threads(nthreads)
  {
#pragma omp for schedule(static)
    for(ii = 0; ii < ntot; ++ii){
      res[ii] = interpolate_one<T,U>( ny, nx, im, yin[ii], xin[ii], missing);
    } // ii
  } //  pragma
  
}
// ******************************************************************************************* //

template<typename T, typename U> inline
void ninterpolate_image_2D(int const ny, int const nx, const T* const __restrict__ im,
			   int const ny1, int const nx1, const T* const __restrict yin,
			   const T* const __restrict__ xin, U* const __restrict__ res,
			   int const nthreads, T const missing)
{

  int ii = 0, ix=0, iy=0;
  int const ntot = nx1*ny1;
  int const nyl = ny-1;
  int const nxl = nx-1;
  
#pragma omp parallel default(shared) firstprivate(ii, ix, iy) num_threads(nthreads)
  {
#pragma omp for schedule(static)
    for(ii = 0; ii < ntot; ++ii){

      if((yin[ii] >= 0) && (yin[ii] <= nyl) && (xin[ii] >= 0) && (xin[ii] <= nxl)){
      
	iy = std::max<int>(std::min<int>(rint(yin[ii]), ny-1), 0);
	ix = std::max<int>(std::min<int>(rint(xin[ii]), nx-1), 0);
	
	res[ii] = im[iy*nx+ix];
      
      }else{
	res[ii] = missing;
      }
    } // ii
  } //  pragma
  
}

// ******************************************************************************************* //

namespace interp {
  
    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT help;
        float missing;
        int missing_set;
        IDL_INT nearest;
        IDL_INT nthreads;
    } KW_INTERP;

    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR kw_interp_pars[] = {
        IDL_KW_FAST_SCAN,
        { "HELP",             IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_INTERP,help) },
        { "MISSING",          IDL_TYP_FLOAT, 1, 0, (int*)IDL_KW_OFFSETOF2(KW_INTERP,missing_set), (char*) IDL_KW_OFFSETOF2(KW_INTERP,missing) },
        { "NEAREST",          IDL_TYP_INT,   1, 0,                      0, (char*) IDL_KW_OFFSETOF2(KW_INTERP,nearest) },
        { "NTHREADS",         IDL_TYP_INT,   1, 0,                      0, (char*) IDL_KW_OFFSETOF2(KW_INTERP,nthreads) },
        { NULL }
    };

    string interp_info( int lvl ) {
        string ret = "RDX_INTERPOL";
        if( lvl > 0 ) {
            ret += ((lvl > 1)?"\n":"      ");          // newline if lvl>1
            ret += "   Syntax:   out = rdx_interpol(file_list, /KEYWORDS)\n";
            if( lvl > 1 ) {
                ret +=  "   Accepted Keywords:\n"
                        "      HELP                Display this info.\n"
                        "      MISSING             Value to fill in.\n"
                        "      NEAREST             Use nearest method.\n"
                        "      NTHREADS            Number of threads.\n";
            }
        } else ret += "\n";
        return ret;
    }

    IDL_VPTR interp_image2D_wrap( int nArg, IDL_VPTR* argv, char* argk ) {

        // --- process optional keywords --- //

        KW_INTERP kw;    
        kw.nthreads = 1;
        kw.nearest = 0;
        int nPlainArgs = IDL_KWProcessByOffset( nArg, argv, argk, kw_interp_pars, ( IDL_VPTR* )0, 1, &kw );

        if( kw.help || nPlainArgs < 3 ) {
            printMessage( interp_info(2) );
            return IDL_GettmpInt(0);
        }

        IDL_VPTR imgVar = argv[0];
        IDL_VPTR xVar   = argv[1];
        IDL_VPTR yVar   = argv[2];

        IDL_ENSURE_SIMPLE( imgVar );
        IDL_ENSURE_ARRAY( imgVar );
        IDL_ENSURE_SIMPLE( xVar );
        IDL_ENSURE_ARRAY( xVar );
        IDL_ENSURE_SIMPLE( yVar );
        IDL_ENSURE_ARRAY( yVar );

        using fp_t = float; // if you change this, make sure you also change IDL_TYP_FLOAT to double

        int const nx    = imgVar->value.arr->dim[0];
        int const ny    = imgVar->value.arr->dim[1];

        int const nx1   = xVar->value.arr->dim[0];
        int const ny1   = xVar->value.arr->dim[1];
            
        int const nx2   = yVar->value.arr->dim[0];
        int const ny2   = yVar->value.arr->dim[1];

        IDL_VPTR result;

        if((nx1 != nx2) || (ny1 != ny2)){
            fprintf(stdout,"error, the x-y interpolation grids have different sizes x[%d,%d] != y[%d,%d]\n", nx1, ny1, nx2, ny2);
            return IDL_GettmpInt(0);
        }

        IDL_LONG64 dims[2] = {nx1, ny1};
        double* const __restrict__ out = (double *) IDL_MakeTempArray( IDL_TYP_DOUBLE, 2, dims, IDL_ARR_INI_ZERO, &result );

        shared_ptr<fp_t> imgData = castOrCopy<fp_t>(imgVar);
        shared_ptr<fp_t> xData = castOrCopy<fp_t>(xVar);
        shared_ptr<fp_t> yData = castOrCopy<fp_t>(yVar);

        const fp_t* const __restrict__ image = imgData.get();
        const fp_t* const __restrict__ x = xData.get();
        const fp_t* const __restrict__ y = yData.get();

        if( !kw.missing_set ) {     // missing was not passed, calculate and use mean of input
            kw.missing = redux::mean(nx*ny, image);
        }

        int const nthreads = std::max<int>(kw.nthreads, 4);
        int const nearest  = std::max<int>(std::abs(kw.nearest), 0);
        double const missing = kw.missing;

        if( nearest ){
            ninterpolate_image_2D<fp_t,double>( ny, nx, image, ny1, nx1, y, x, out, nthreads, missing );
        } else {
            binterpolate_image_2D<fp_t,double>( ny, nx, image, ny1, nx1, y, x, out, nthreads, missing );
        }

        IDL_KW_FREE;

        return result;

    }

}

// ******************************************************************************************* //



namespace {
    static int dummy RDX_UNUSED =
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)dsgridnest}, (char*)"RDX_CDSGRIDNEST", 4, 4, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, dsgridnest_info ) +
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)stretch_wrap}, (char*)"RDX_CSTRETCH", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, stretch_info ) +
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)interp::interp_image2D_wrap}, (char*)"RDX_INTERPOL", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, interp::interp_info );
}
