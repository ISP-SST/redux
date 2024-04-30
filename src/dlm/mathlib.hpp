#ifndef REDUX_DLM_MATHLIB_HPP
#define REDUX_DLM_MATHLIB_HPP

#include "redux/util/arrayutil.hpp"

#include <cmath>
#include <algorithm>
#include <typeinfo>
#include <omp.h>
#include <cstring>

namespace redux {

    template<typename T>
    void bilint2D ( int const ny, int const nx, const T* const __restrict__ y, const T* const __restrict__ x,
                    const T* const __restrict__ d_in, int const ny1, int const nx1, const T* const __restrict__  yy_in,
                    const T* const __restrict__  xx_in, T* __restrict__ res_in, int const nthreads )
    {

        //double* __restrict__ c_in = new double [4* ( ny-1 ) * ( nx-1 )]();

        // --- Need to pre-define these types so Clang++ swallows the "reinterpret_cast" --- //

        //using ft = T[ny][nx];
        //using ft1= T[ny1][nx1];
        //using ft2= double[ny-1][nx-1][4];

        //ft  const  &d  = *reinterpret_cast<const ft*  const> ( d_in );
        //ft1 const &xx  = *reinterpret_cast<const ft1* const> ( xx_in );
        //ft1 const &yy  = *reinterpret_cast<const ft1* const> ( yy_in );
        
        T** d = util::makePointers( d_in, ny, nx );
        T** xx = util::makePointers( xx_in, ny1, nx1 );
        T** yy = util::makePointers( yy_in, ny1, nx1 );
        
        //ft1 &res = *reinterpret_cast<ft1*> ( res_in );
        //ft2 &c   = *reinterpret_cast<ft2*> ( c_in );

        T** res = util::makePointers( res_in, ny1, nx1 );
        double*** c = util::newArray<double>( ny-1, nx-1, 4 );

        // --- compute interpolation coefficients --- //

        int const nxx = nx-1;
        int const nyy = ny-1;

        const double* __restrict__ cc = NULL;
        double dydx=0, x1=0,x2=0,dy=0,y1=0,y2=0,dx=0,ixx=0,iyy=0;
        int ii=0, jj=0,idx=0,idy=0,tt=0;

        #pragma omp parallel default(shared) firstprivate(jj,ii,dy,y1,y2,dx,x1,x2,dydx) num_threads(nthreads)
        {
            #pragma omp for schedule(static) collapse(2)
            for ( jj=0; jj<nyy; ++jj ) {
                for ( ii=0; ii<nxx; ++ii ) {

                    dy = y[jj] - y[jj+1];
                    dx = x[ii] - x[ii+1];

                    dydx = dx*dy;

                    y1 = y[jj]/dydx;
                    y2 = y[jj+1]/dydx;
                    x1 = x[ii]/dydx;
                    x2 = x[ii+1]/dydx;

                    c[jj][ii][0] = ( d[jj][ii]*x2*y2 - d[jj+1][ii]*x2*y1 - d[jj][ii+1]*x1*y2 + d[jj+1][ii+1]*x1*y1 );
                    c[jj][ii][1] = ( -d[jj][ii]*y2    + d[jj+1][ii]*y1    + d[jj][ii+1]*y2    - d[jj+1][ii+1]*y1 );
                    c[jj][ii][2] = ( -d[jj][ii]*x2    + d[jj+1][ii]*x2    + d[jj][ii+1]*x1    - d[jj+1][ii+1]*x1 );
                    c[jj][ii][3] = ( d[jj][ii]       - d[jj+1][ii]       - d[jj][ii+1]       + d[jj+1][ii+1] );
                }//ii
            }//jj
        }// pragma

        // --- Do interpolation --- //

        T const yy0 = y[0];
        T const yy1 = y[ny-1];
        T const xx0 = x[0];
        T const xx1 = x[nx-1];


        #pragma omp parallel default(shared) firstprivate(jj,ii,ixx,iyy,idx,idy,tt,cc) num_threads(nthreads)
        {
            #pragma omp for schedule(static) collapse(2)
            for ( jj=0; jj<ny1; ++jj ) {
                for ( ii=0; ii<nx1; ++ii ) {

                    // --- If the pixel is out of bounds, force coefficients at the edge --- //

                    ixx = std::min ( std::max<T> ( xx[jj][ii],xx0 ), xx1 );
                    iyy = std::min ( std::max<T> ( yy[jj][ii],yy0 ), yy1 );

                    idx = 0, idy = 0;


                    // --- bracket coordinates --- //

                    for ( tt=0; tt<nyy; ++tt ) {
                        if ( ( iyy >= y[tt] ) && ( iyy <= y[tt+1] ) ) {
                            idy = tt;
                            break;
                        } // if
                    }

                    for ( tt=0; tt<nxx; ++tt ) {
                        if ( ( ixx >= x[tt] ) && ( ixx <= x[tt+1] ) ) {
                            idx = tt;
                            break;
                        } // if
                    }


                    // --- Apply interpolation coefficients --- //

                    cc = static_cast<const double* const __restrict__> ( &c[idy][idx][0] );
                    res[jj][ii] = cc[0] + cc[1]*ixx + cc[2]*iyy + cc[3]*iyy*ixx;

                    //	} // else
                } // ii
            } // jj
        }// pragma

        util::delPointers( d );
        util::delPointers( xx );
        util::delPointers( yy );
        util::delPointers( c );
        util::delArray( c );
        
        //delete [] c_in;

    }

    // ********************************************************************************** //

    template <typename T> inline
    constexpr int sign ( T val )
    {
        return ( T ( 0 ) < val ) - ( val < T ( 0 ) );
    }

    // ********************************************************************************** //


    template<typename T, typename U>
    struct Round {
        inline constexpr static U run ( T const& val )
        {
            return ( std::abs ( val )+0.5 ) *sign<T> ( val );
        }
    };


    // ********************************************************************************** //

    template<typename T, typename U = T>
    U* congrid ( int const ny, int const nx, const T* const __restrict__ im, int const ny1, int const nx1 )
    {
        U* __restrict__ res = new U [nx1*ny1]();
        T const sclx = T ( nx-1.0 ) / T ( nx1-1.0 );
        T const scly = T ( ny-1.0 ) / T ( ny1-1.0 );

        for ( int jj=0; jj<ny1; ++jj ) {
	  T   const yl = std::min<T>(jj * scly, ny-1);
	  int const iy = std::min<int> ( int ( yl ), ny-2 );

            T const dy  = yl-iy;
            T const dy1 = 1 - dy;

            for ( int ii=0; ii<nx1; ++ii ) {
	      T const xl = std::min<T>(ii * sclx, nx-1);
	      int    const ix = std::min<int> ( int ( xl ), nx-2 );

                T const dx  = xl-ix;
                T const dx1 = 1 - dx;

                res[jj*nx1+ii] = * ( im+iy*nx+ix ) *dx1*dy1 + * ( im+iy*nx+ix+1 ) *dx*dy1 + * ( im+ ( iy+1 ) *nx+ix ) *dx1*dy + * ( im+ ( iy+1 ) *nx+ix+1 ) *dx*dy;
            }
        }

        return res;
    }

    // ********************************************************************************** //

    template<typename T> inline
    T* arange ( int const nEl )
    {
        T* __restrict__ res = new T [nEl]();
        for ( int ii=0; ii<nEl; ++ii ) {
            res[ii] = ii;
        }
        return res;
    }

    // ********************************************************************************** //

    template<typename T> inline
    T mean ( int const n, const T* const __restrict__ arr )
    {

        double sum = 0;
        #pragma omp simd reduction(+: sum) // Needed for vectorization in reduction operations
        for ( int ii=0; ii<n; ++ii ) {
            sum += arr[ii];
        }

        return T ( sum/n );
    }

    // ********************************************************************************** //


}


#endif  // REDUX_DLM_MATHLIB_HPP
