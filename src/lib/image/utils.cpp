#include "redux/image/utils.hpp"

#include "redux/image/fouriertransform.hpp"
#include "redux/image/grid.hpp"
#include "redux/image/pupil.hpp"
#include "redux/image/zernike.hpp"

#include "redux/file/fileana.hpp"
#include "redux/math/functions.hpp"
#include "redux/constants.hpp"
#include "redux/types.hpp"

#include <functional>
#include <map>
#include <set>

#include <gsl/gsl_multifit.h>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#ifdef REDUX_WITH_OPENCV
#       include <opencv2/photo/photo.hpp>
#endif
using namespace redux::image;
using namespace redux::momfbd;
using namespace redux::util;
using namespace redux;
using namespace std;

namespace {

    const int maxDistance = 64;
    const int maxDistance2 = maxDistance*maxDistance;
    const double deltasqr = 4;
    const double beta = 2;
    
    double distMap[2*maxDistance2+1];
    
    bool calculateCache( void ) {

        return true;
    }

    const double* getDistanceMap (void) {
        memset(distMap,0,(2*maxDistance2+1)*sizeof(double));
        for(int i=0; i<=2*maxDistance2; ++i) distMap[i] = pow (i + deltasqr, -beta);
        return distMap;
    }


    /*! For Zernike polynomials: convert Noll's sequential index to m,n form
    * @note: does not return the sign for negative m values (odd j).
    */
    void noll_to_mn (int j, int& m, int& n) {
        n = 0;
        int len = 1;
        for (int i = 1; len < j; ++i) {
            len += (n = i) + 1;
        }
        int dl = n + 1 - len + j;
        m = 2 * ( (dl + (n % 2)) / 2) + ! (n % 2) - 1;
    }

    void zernike_mn (int j, int &m, int &n) { // j=1...
        n = 0;
        int len = 1;
        for (int i = 1; len < j; ++i) len += (n = i) + 1;
        int dl = n + 1 - len + j;
        m = 2 * ( (dl + (n % 2)) / 2) + ! (n % 2) - 1;
    }

    double *RC (int m, int n) {
        int nmm = (n - m) / 2, npm = (n + m) / 2;
        int nmax = std::max (npm, n);
        double *f = new double [nmax + 1];
        f[0] = 1.0;
        for (int i = 1; i <= nmax; ++i) f[i] = (double) i * f[i - 1];
        double *res = new double [nmm + 1];
        for (int s = 0, pm = -1; s <= nmm; ++s)
            res[s] = (double) ( (pm *= -1) * f[n - s]) / (f[s] * f[npm - s] * f[nmm - s]);
        delete[] f;
        return res;
    }


    __inline double mvn_max (double a, double b) {
        return (a > b) ? a : b;
    }

    __inline double mvn_min (double a, double b) {
        return (a > b) ? b : a;
    }

    template <typename T>
    void maskConnected(T** data, uint8_t** mask, size_t sizeY, size_t sizeX, unsigned int yy, unsigned int xx, unsigned int y_start, unsigned int x_start, T threshold) {
        if (yy >= sizeY || xx >= sizeX) return;
        if ( !mask[yy][xx] && ((yy == y_start && xx == x_start) || data[yy][xx] > threshold)) {       // new point inserted, check neighbours
            mask[yy][xx] = 1;
            maskConnected(data, mask, sizeY, sizeX, yy + 1, xx, y_start, x_start, threshold);
            maskConnected(data, mask, sizeY, sizeX, yy - 1, xx, y_start, x_start, threshold);
            maskConnected(data, mask, sizeY, sizeX, yy, xx + 1, y_start, x_start, threshold);
            maskConnected(data, mask, sizeY, sizeX, yy, xx - 1, y_start, x_start, threshold);
        }
    }

}

template <typename T>
static __inline T sqr (T x) {
    return x * x;
}


double redux::image::makePupil_thi (double** pupil, uint32_t nPoints, double radius) {

    double area = 0.0, origin = 0.5;
    memset (*pupil, 0, nPoints * nPoints * sizeof (double));
    uint32_t mid = nPoints / 2; // N.B: Don't use odd number of points for pupils !!
    // Pupil should be centered ON pixel (mid,mid), to match the location of the origin for Fourier-transforms.

    if (nPoints % 2) {
        // TODO: warn or throw if nPoints is odd.
    }

    const Grid& grid = Grid::get(mid + 1, origin, origin);
    float** distPtr = grid.distance.get();      // distance(i,j) is the distance from the centre of the pupil, to the inner boundary of pixel (i,j)
    // i.e. dist(0,0) = dist(0,1) = dist(1,0) = dist(1,1) = sqrt(2)/2  (it is centered on that pixel)
    double val;
    for (unsigned int x = 0; x < mid; ++x) {
        for (unsigned int y = 0; y <= x; ++y) {     // We only generate the first octant, then copy.
            val = 0;
            if (distPtr[y + 1][x + 1] < radius) {
                val = 1;
            } else if (distPtr[y][x] < radius) {    // partial pixel
                if (x == 0 && y == 0) {     // central pixel = 1 for all practical cases
                    if (radius < 0.5) val = redux::PI * radius * radius;    // a pupil of size < sqrt(2) pixel is a bit absurd...
                    else val = redux::PI * radius * radius + (radius - 0.5) / (sqrt (0.5) - 0.5) * (1 - redux::PI * radius * radius);
                } else {
                    // TBD: better approximation of pixel fill-factor ??
                    val = (radius - distPtr[y][x]) / (distPtr[y + 1][x + 1] - distPtr[y][x]); // linear fill-factor from radial ratio
                }
            }
            if (val > 0) {
                pupil[mid + y][mid + x] = val;
                area += val;
                if (x != y) {
                    pupil[mid + x][mid + y] = val; // Use symmetry to fill the second octant
                    area += val;
                }
            }
        }
    }
    for (unsigned int x = 0; x < mid; ++x) {
        for (unsigned int y = 0; y < mid; ++y) {     // copy 1st quadrant to 2,3,4
            val = pupil[mid + y][mid + x];
            if (val > 0) {
                if (x) {
                    pupil[mid + y][mid - x] = val;
                    area += val;
                }
                if (y) {
                    pupil[mid - y][mid + x] = val;
                    area += val;
                }
                if (x && y) {
                    pupil[mid - y][mid - x] = val;
                    area += val;
                }
            }
        }
    }

    return area;

}


double redux::image::makePupil_mvn (double** pupil, int nph, double r_c) {

    double area = 0.0;
    memset (*pupil, 0, nph * nph * sizeof (double));
    int xo = nph / 2, yo = nph / 2;
    double dx = 0.5 / r_c, dy = 0.5 / r_c;
    for (int x = 0; x < nph; ++x) {
        double xl = fabs ( (double) (x - xo)) / r_c - dx, xh = fabs ( (double) (x - xo)) / r_c + dx;
        double xhs = sqr (xh);
        for (int y = 0; y < nph; ++y) {
            double yl = fabs ( (double) (y - yo)) / r_c - dy, yh = fabs ( (double) (y - yo)) / r_c + dy;
            double yhs = sqr (yh);
            double rsl = sqr (xl) + sqr (yl), rsh = xhs + yhs;
            if (rsl <= 1.0) {  // inside pixel
                if (rsh < 1.0) {  // full pixel
                    pupil[y][x] = 1.0;
                    //pupil[y][x] = sqrt(rsh*r_c*r_c);
                } else {           // partial pixel
                    double x2 = sqrt (mvn_max (1.0 - yhs, (double) 0.0));
                    double y3 = sqrt (mvn_max (1.0 - xhs, (double) 0.0));
                    double f = (xh > yh) ? (yh - yl) * (mvn_min (xh, mvn_max (xl, x2)) - xl) / (4 * dx * dy) :
                               (xh - xl) * (mvn_min (yh, mvn_max (yl, y3)) - yl) / (4 * dx * dy);
                    pupil[y][x] = f;
                }
                area += pupil[y][x];
            }
        }
    }

    return area;

}


void redux::image::makeZernike_thi (double** modePtr, int modeNumber, uint32_t nPoints, double r_c, double angle) {

    redux::image::Pupil pupil(nPoints, r_c);

    int m, n;
    noll_to_mn (modeNumber, m, n);

    const vector<double>& coeff = Zernike::radialPolynomial (m, n);
    const Grid& grid = Grid::get(nPoints, nPoints/2.0, nPoints/2.0);
    float** distPtr = grid.distance.get();  // distance from pixels to centre (pupil & modes are centered on pixel (mid,mid))
    float** aPtr = grid.angle.get();

    memset (*modePtr, 0, nPoints * nPoints * sizeof (double));

    double** pupPtr = makePointers (pupil.ptr(), nPoints, nPoints);

    shared_ptr<double*> r = sharedArray<double> (nPoints, nPoints);     // normalized distance ^{some order}
    shared_ptr<double*> r2 = sharedArray<double> (nPoints, nPoints);    // normalized distance squared
    double** rPtr = r.get();
    double** r2Ptr = r2.get();
    double normalization (0);

    for (unsigned int y = 0; y < nPoints; ++y) {
        for (unsigned int x = 0; x < nPoints; ++x) {
            //if(pupPtr[y][x]>0) {
            double tmp = distPtr[y][x] / r_c;                       // normalize radial distance
            //if(tmp>1) continue;
            r2Ptr[y][x] = tmp * tmp;
            if (m == 0) rPtr[y][x] = 1;
            else if (m == 1) rPtr[y][x] = tmp;
            else rPtr[y][x] = pow (tmp, m);                         // lowest order term ~ r^{m}
            modePtr[y][x] = rPtr[y][x] * coeff[0];                  // add lowest order to mode
            //}
        }
    }

    // generate polynomial part
    for (auto it = coeff.begin() + 1; it != coeff.end(); it++) {
        for (unsigned int y = 0; y < nPoints; ++y) {
            for (unsigned int x = 0; x < nPoints; ++x) {
                //if(pupPtr[y][x]>0) {
                rPtr[y][x] *= r2Ptr[y][x];                          //  next term ~ r^{m+2}
                modePtr[y][x] += rPtr[y][x] * *it;
                //}
            }
        }
    }

    // Angular component
    if (m == 0) {
        double sf = sqrt (n + 1);
        for (unsigned int y = 0; y < nPoints; ++y) {
            for (unsigned int x = 0; x < nPoints; ++x) {
                modePtr[y][x] *= sf; // * pupPtr[y][x];
                if (pupPtr[y][x] > 0) {
                    normalization +=  modePtr[y][x] * modePtr[y][x] * pupPtr[y][x];
                }
            }
        }
    } else if (modeNumber % 2) {
        double sf = sqrt ( (double) 2 * (n + 1));
        for (unsigned int y = 0; y < nPoints; ++y) {
            for (unsigned int x = 0; x < nPoints; ++x) {
                modePtr[y][x] *= sf * sin (m * (aPtr[y][x] - angle)); //* pupPtr[y][x]
                if (pupPtr[y][x] > 0) {
                    normalization +=  modePtr[y][x] * modePtr[y][x] * pupPtr[y][x];
                }
            }
        }
    } else {
        double sf = sqrt ( (double) 2 * (n + 1));
        for (unsigned int y = 0; y < nPoints; ++y) {
            for (unsigned int x = 0; x < nPoints; ++x) {
                modePtr[y][x] *= sf * cos (m * (aPtr[y][x] - angle)); //* pupPtr[y][x]
                if (pupPtr[y][x] > 0) {
                    normalization +=  modePtr[y][x] * modePtr[y][x] * pupPtr[y][x];
                }
            }
        }
    }

    // normalize
    normalization = sqrt(pupil.area/normalization);
    for (unsigned int y = 0; y < nPoints; ++y) {
        for (unsigned int x = 0; x < nPoints; ++x) {
            //if(pupPtr[y][x]>0) {
            modePtr[y][x] *= normalization;
            //}
        }
    }

    delPointers (pupPtr);

}


void redux::image::makeZernike_mvn (double** mode, int j, uint32_t nph, double r_c, double angle) {

    int xo = nph / 2, yo = nph / 2;
    if (j == 1) {
        for (unsigned int x=0; x < nph; ++x) {
            for (unsigned int y=0; y < nph; ++y) {
                mode[x][y] = 1.0;
            }
        }
    } else {                                      // j>1
        int m, n;
        zernike_mn (j, m, n);
        double *rc = RC (m, n);
        double **r = newArray<double> (nph, nph);
        double **rs = newArray<double> (nph, nph);
        for (unsigned int x=0; x < nph; ++x) {            // s=0
            for (unsigned int y=0; y < nph; ++y) {
                double rr = sqrt ( (double) sqr (x - xo) + (double) sqr (y - yo)) / r_c;
                rs[x][y] = rr * rr;
                r[x][y] = pow (rr, n);
                mode[x][y] = r[x][y] * rc[0];
            }
        }
        rs[xo][yo] = 1.0;                         // avoid division by 0
        for (int s = 1; s <= (n - m) / 2; ++s) {
            for (unsigned int x = 0; x < nph; ++x) {
                for (unsigned int y = 0; y < nph; ++y) {
                    mode[x][y] += (r[x][y] /= rs[x][y]) * rc[s];
                }
            }
            if (! (n - 2 * s)) mode[xo][yo] += rc[s]; // dividing 0 by 1 will never give 1...
        }
        delArray (rs);
        delArray (r);

        if (m) {                                  // m!=0
            double sf = sqrt ( (double) (2 * (n + 1)));
            if (j % 2)                              // odd
                for (int x = 0; x < nph; ++x) {
                    for (int y = 0; y < nph; ++y) {
                        mode[x][y] *= sf * sin ( ( (double) m) * (atan2 ( (double) (y - yo), (double) (x - xo)) + angle));
                    }
                }
            else {                                   // even
                for (int x = 0; x < nph; ++x) {
                    for (int y = 0; y < nph; ++y) {
                        mode[x][y] *= sf * cos ( ( (double) m) * (atan2 ( (double) (y - yo), (double) (x - xo)) + angle));
                    }
                }
            }
        } else {                                   // m==0
            double sf = sqrt ( (double) (n + 1));
            for (unsigned int x = 0; x < nph; ++x) {
                for (unsigned int y = 0; y < nph; ++y) {
                    mode[x][y] *= sf;
                }
            }
        }
        delete[] rc;
    }
    double sum = 0.0, N = 0.0, dx = 0.5 / r_c, dy = 0.5 / r_c;
    for (unsigned int x = 0; x < nph; ++x) {
        double xl = fabs ( ( (double) x - xo)) / r_c - dx, xh = fabs ( ( (double) x - xo)) / r_c + dx;
        double xhs = sqr (xh);
        for (unsigned int y = 0; y < nph; ++y) {
            double yl = fabs ( ( (double) y - yo)) / r_c - dy, yh = fabs ( ( (double) y - yo)) / r_c + dy;
            double yhs = sqr (yh);
            double rsl = sqr (xl) + sqr (yl), rsh = xhs + yhs;
            if (rsl <= 1.0) {   // good pixel
                if (rsh < 1.0) { // full pixel
                    sum += sqr (mode[x][y]);
                    N += 1.0;
                } else {          // partial pixel
                    double x2 = sqrt (max (1.0 - yhs, (double) 0.0));
                    double y3 = sqrt (max (1.0 - xhs, (double) 0.0));
                    double f = (xh > yh) ? (yh - yl) * (min (xh, max (xl, x2)) - xl) / (4 * dx * dy) :
                               (xh - xl) * (min (yh, max (yl, y3)) - yl) / (4 * dx * dy);
                    sum += f * sqr (mode[x][y]);
                    N += f;
                }
            }
        }
    }
    sum /= N;
    for (unsigned int x = 0; x < nph; ++x) {
        for (unsigned int y = 0; y < nph; ++y) {
            mode[x][y] /= sqrt (sum);
        }
    }
}


template <typename T>
redux::util::Array<T> redux::image::fitPlane (const redux::util::Array<T>& in, bool subtract_mean) {

    int ySize = in.dimSize(0);
    int xSize = in.dimSize(1);
    redux::util::Array<T> ret(ySize, xSize);

    int n = ySize * xSize;
    int nParams = 3;                                                        //   fit a plane as:   z = a*x + b*y + c

    gsl_vector *data = gsl_vector_alloc (n);
    gsl_vector *coeff = gsl_vector_alloc (nParams);
    gsl_matrix *X = gsl_matrix_alloc (n, nParams);
    gsl_matrix *covar = gsl_matrix_alloc (nParams, nParams);

    int yHalf = ySize/2;
    int xHalf = xSize/2;
    const T* inPtr = in.ptr();
    for (int i = 0; i < ySize; ++i) {
        double y = i-yHalf;
        for (int j = 0; j < xSize; ++j) {
            double x = j-xHalf;
            int offset = i * xSize + j;
            gsl_matrix_set (X, offset, 0, x);
            gsl_matrix_set (X, offset, 1, y);
            gsl_matrix_set (X, offset, 2, 1);
            gsl_vector_set (data, offset, inPtr[offset]);
        }
    }

    double chisq;
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, nParams);
    gsl_multifit_linear (X, data, coeff, covar, &chisq, work);
    gsl_multifit_linear_free (work);

    double a = gsl_vector_get (coeff, 0);
    double b = gsl_vector_get (coeff, 1);
    double c = gsl_vector_get (coeff, 2);

    gsl_vector_free (data);
    gsl_vector_free (coeff);
    gsl_matrix_free (X);
    gsl_matrix_free (covar);

    T* retPtr = ret.ptr();
    if (subtract_mean) {
        c = 0;                                          // ignore c-coefficient (mean), to just fit the tilts, not the offset.
    }
    for (int i = 0; i < ySize; ++i) {
        double y = static_cast<double> (i) / (ySize - 1) - 0.5;
        for (int j = 0; j < xSize; ++j) {
            double x = static_cast<double> (j) / (xSize - 1) - 0.5;
            int offset = i * xSize + j;
            retPtr[offset] = a * x + b * y + c;
        }
    }

    return ret;
}
template redux::util::Array<float> redux::image::fitPlane (const redux::util::Array<float>&, bool);
template redux::util::Array<double> redux::image::fitPlane (const redux::util::Array<double>&, bool);


template <typename T>
void redux::image::connectedRegion(T** data, uint8_t** mask, size_t sizeY, size_t sizeX, unsigned int y_start, unsigned int x_start, T threshold) {
    memset(*mask, 0, sizeY * sizeX);
    maskConnected(data, mask, sizeY, sizeX, y_start, x_start, y_start, x_start, threshold);
}
template void redux::image::connectedRegion(float**, uint8_t**, size_t, size_t, uint, uint, float);
template void redux::image::connectedRegion(double**, uint8_t**, size_t, size_t, uint, uint, double);


template <typename T>
void redux::image::smooth(T** data, size_t sizeY, size_t sizeX, size_t nY, size_t nX) {
    if (nY == 0 || nX == 0) return;
    T** tmp = newArray<T> (sizeY, sizeX);
    memset (*tmp, 0, sizeY * sizeX * sizeof (T));
    for (unsigned int y = 0; y < sizeY; ++y) {
        int yl = std::max<int>(y-nY, 0);
        int yh = std::min(y+nY, sizeY);
        for (unsigned int x = 0; x < sizeX; ++x) {
            int xl = std::max<int>(x-nX, 0);
            int xh = std::min(x+nX, sizeX);
            int cnt = 0;
            for (int yy=yl; yy < yh; ++yy) {
                for (int xx=xl; xx < xh; ++xx) {
                    tmp[y][x] += data[yy][xx];
                    cnt++;
                }
            }
            tmp[y][x] /= cnt;
        }
    }
    memcpy (*data, *tmp, sizeY * sizeX * sizeof (T));
    delArray (tmp);
}
template void redux::image::smooth (float**, size_t, size_t, size_t, size_t);
template void redux::image::smooth (double**, size_t, size_t, size_t, size_t);
template void redux::image::smooth (complex_t**, size_t, size_t, size_t, size_t);


template <typename T>
void redux::image::ScharmerFilter (T** data, double** q2_inv, size_t sizeY, size_t sizeX, double noise_power, double frequency_cutoff) {

    static const double hi = 1.0;
    static const double lo = 0.2;
    double noiseFactor = noise_power*sizeY*sizeX;
    
    double** pt = newArray<double> (sizeY, sizeX);
    double** rr = newArray<double> (sizeY, sizeX);
    
    std::transform (*data, *data + sizeY * sizeX, *q2_inv, *pt, [](const T& a, const double& b ){ return norm(a)*b; } );
    
    smooth (pt, sizeY, sizeX, 1, 1);
    
    std::transform (*pt, *pt + sizeY * sizeX, *pt, [noiseFactor](const double& a){ return noiseFactor/a;} );
    
    for (unsigned int y=0; y<sizeY; ++y) {
        rr[y][0] = pt[y][0];
        for (unsigned int x=1; x<sizeX; ++x) {
            rr[y][x] = pt[y][sizeX-x];
        }
    }
    
    double* tmp = newArray<double> (sizeY);
    for(unsigned int x=0; x<sizeX; ++x){                        // rr==filter
        rr[0][x] = max((1.0-0.5*(rr[0][x]+pt[0][x])),0.0);
        for(unsigned int y=1; y<sizeY;++y) tmp[y]=rr[y][x];       // temporary storage
        for(unsigned int y=1; y<sizeY;++y) rr[y][x]=max((1.0-0.5*(tmp[sizeY-y]+pt[y][x])), 0.0);
    }
    delArray (tmp);
    
    for(unsigned int y=0; y<sizeY; ++y) {
        for(unsigned int x=0; x<sizeX; ++x){
            if(rr[y][x]<lo) rr[y][x] = 0.0;
            if(rr[y][x]>hi) rr[y][x] = 1.0;
        }
    }
    
    unsigned int yHalf = sizeY/2;
    unsigned int xHalf = sizeX/2;
    rr[yHalf][xHalf] = 1.0;                   // ; DC gain = 1
    
    connectedRegion(rr, sizeY, sizeX, yHalf, xHalf, 0.0);
    smooth(rr, sizeY, sizeX, 4, 4);
    
    frequency_cutoff *= frequency_cutoff;
    for(unsigned int y=0; y<sizeY; ++y) {
        for(unsigned int x=0; x<sizeX; ++x) {
            if( (sqr(y-yHalf)+sqr(x-xHalf)) > frequency_cutoff){
                data[y][x] = 0.0;
            } else{
                data[y][x] *= rr[y][x];
            }
        }
    }
    delArray (pt);
    delArray (rr);
    
}
template void redux::image::ScharmerFilter (float**, double**, size_t, size_t, double, double);
template void redux::image::ScharmerFilter (double**, double**, size_t, size_t, double, double);
template void redux::image::ScharmerFilter (complex_t**, double**, size_t, size_t, double, double);


//double inv_dist_wght(double **a,int xl,int xh,int yl,int yh,int xb,int yb)
double redux::image::inv_dist_wght (float **a, size_t sizeY, size_t sizeX, size_t posY, size_t posX) {

    int xl = std::max (0L, static_cast<int64_t> (posX - maxDistance));
    int xh = std::min (sizeX, posX + maxDistance + 1);
    int yl = std::max (0L, static_cast<int64_t> (posY - maxDistance));
    int yh = std::min (sizeY, posY + maxDistance + 1);


    double weight = 0.0, res = 0.0;
    for (int y = yl; y < yh; ++y)
        for (int x = xl; x < xh; ++x)
            if (a[y][x]) {
                double c = pow ( (double) sqr (x - posX) + (double) sqr (y - posY) + deltasqr, -beta);
                res += c * a[y][x];
                weight += c;
            }
    return res / weight;
}


template <typename T>
double redux::image::inverseDistanceWeight (T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX) {

    static const double* const inverseDistanceSquared = getDistanceMap();     // will only be initialized once.

    // TODO: verify this function, results look weird
    int64_t beginX = std::max (0L, static_cast<int64_t> (posX - maxDistance));
    int64_t endX = std::min (sizeX, posX + maxDistance+1);
    int64_t beginY = std::max (0L, static_cast<int64_t> (posY - maxDistance));
    int64_t endY = std::min (sizeY, posY + maxDistance+1);

    double normalization = 0.0, weightedSum = 0.0;
    for (int y=beginY; y < endY; ++y) {
        int y2 = (y-posY)*(y-posY);
        for (int x = beginX; x < endX; ++x) {
            int x2 = (x-posX)*(x-posX);
            //if( x2+y2 > maxDistance2 ) break;
            if( array[y][x] ) {
                double tmp = inverseDistanceSquared[y2+x2];
                weightedSum += tmp * array[y][x];
                normalization += tmp;
            }
        }
    }
    if( normalization ) {
        return weightedSum / normalization;
    }
    return 0.0;
}
template double redux::image::inverseDistanceWeight (float**, size_t, size_t, size_t, size_t);
template double redux::image::inverseDistanceWeight (double**, size_t, size_t, size_t, size_t);

template <typename T>
double redux::image::horizontalInterpolation (T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX) {

    T* ptr = array[posY];

    //map the 5 pixel surrounding as bits in a byte
    int val = 0;
    if ( (posX > 1)) val |= ( (ptr[posX - 2] > 0) << 4);
    if ( (posX > 0))   val |= ( (ptr[posX - 1] > 0) << 3);
    if ( (posX+1 < sizeX))   val |= ( (ptr[posX + 1] > 0) << 1);
    if ( (posX+2 < sizeX)) val |= ( (ptr[posX + 2] > 0));
    //now select based on the number
    switch (val) {
        case (10) :     // = 0 1 x 1 0
        case (11) :     // = 0 1 x 1 1
        case (26) :     // = 1 1 x 1 0
        case (27) :     // = 1 1 x 1 1
            return (ptr[posX-1] + ptr[posX+1]) / 2;
        case (18) :     // = 1 0 x 1 0
        case (19) :     // = 1 0 x 1 1
            return (ptr[posX-2] + 2 * ptr[posX+1]) / 3;
        case (9) :      // = 0 1 x 0 1
        case (25) :     // = 1 1 x 0 1
            return (2 * ptr[posX-1] + ptr[posX+2]) / 3;
        default:
            return inverseDistanceWeight<T> (array, sizeY, sizeX, posY, posX);
    }

}
template double redux::image::horizontalInterpolation (float**, size_t, size_t, size_t, size_t);
template double redux::image::horizontalInterpolation (double**, size_t, size_t, size_t, size_t);


template <typename T>
void redux::image::apodizeInPlace (Array<T>& array, size_t blendRegion) {

    size_t sizeY = array.dimSize (0);
    size_t sizeX = array.dimSize (1);
    blendRegion = std::min (std::min (blendRegion, sizeY >> 1), sizeX >> 1);

    double* tmp = new double[blendRegion];
    tmp[0] = 0;
    redux::math::apodize (tmp, blendRegion, 1.0);
    for (size_t y = 0; y < blendRegion; ++y) {
        for (size_t x = 0; x < sizeX; ++x) {
            array (y, x) *= tmp[y];
            array (sizeY - y - 1, x)  *= tmp[y];
        }
    }
    for (size_t x = 0; x < blendRegion; ++x) {
        for (size_t y = 0; y < sizeY; ++y) {
            array (y, x)  *= tmp[x];
            array (y, sizeX - x - 1) *= tmp[x];
        }
    }
    delete[] tmp;

}
template void redux::image::apodizeInPlace (Array<int16_t>&, size_t);
template void redux::image::apodizeInPlace (Array<int32_t>&, size_t);
template void redux::image::apodizeInPlace (Array<double>&, size_t);
template void redux::image::apodizeInPlace (Array<float>&, size_t);
template void redux::image::apodizeInPlace (Array<complex_t>&, size_t);


template <typename T>
Array<T> redux::image::apodize (const Array<T>& in, size_t blendRegion) {
    if (!blendRegion) return in;       // nothing to do
    Array<T> array;
    in.copy (array);
    apodizeInPlace (array, blendRegion);
    return array;
}
template Array<int16_t> redux::image::apodize (const Array<int16_t>&, size_t);
template Array<int32_t> redux::image::apodize (const Array<int32_t>&, size_t);
template Array<double> redux::image::apodize (const Array<double>&, size_t);
template Array<float> redux::image::apodize (const Array<float>&, size_t);
template Array<complex_t> redux::image::apodize (const Array<complex_t>&, size_t);


template <typename T>
void redux::image::normalizeIfMultiFrames (redux::image::Image<T>& img) {
    if (img.meta) {
        std::string hdr = img.meta->getText();
        boost::regex re ("(\\d+)[ .]+SUM[= ]+");
        boost::smatch match;
        if (boost::regex_search (hdr, match, re)) {
            int nFrames = boost::lexical_cast<int> (match[1]);
            if (nFrames > 1) {
                img *= (1.0/nFrames);
            }
        }
    }
}
template void redux::image::normalizeIfMultiFrames (redux::image::Image<int16_t>&);
template void redux::image::normalizeIfMultiFrames (redux::image::Image<int32_t>&);
template void redux::image::normalizeIfMultiFrames (redux::image::Image<double>&);
template void redux::image::normalizeIfMultiFrames (redux::image::Image<float>&);


template <typename T, typename U>
void redux::image::inpaint( T* img, U* mask, T* out, size_t ySize, size_t xSize, double radius, int flags ) {
#ifdef REDUX_WITH_OPENCV
    cv::Mat imgMat( ySize, xSize, cv::cvType<T>(), img );
    cv::Mat maskMat( ySize, xSize, cv::cvType<U>(), mask );
    cv::Mat outMat( ySize, xSize, cv::cvType<T>(), out );
    cv::inpaint( imgMat, maskMat, outMat, radius, flags );
#else
    std::cerr << "make_mask is not yet implemented for non-OpenCV builds." << std::endl;
#endif            
}
template void redux::image::inpaint( float*, uint8_t*, float*, size_t, size_t, double, int );
template void redux::image::inpaint( double*, uint8_t*, double*, size_t, size_t, double, int );


