#include "redux/momfbd/legacy.hpp"

#include "redux/constants.hpp"
#include "redux/momfbd/modes.hpp"
#include "redux/math/linalg.hpp"
#include "redux/util/arrayutil.hpp"

#include <cmath>
#include <cstring>

using namespace redux::util;
using namespace redux::momfbd;
using namespace std;

inline double sqr(double a) {
    return a * a;
}

inline double sign(double a) {
    return (a >= 0.0) ? 1 : -1;
}

inline double sign(double a, double b) {
    return (b >= 0.0) ? fabs(a) : -fabs(a);
}

double legacy::gammln(const double& xx, double &sign) {

    static double cof[] = { 7.618009172947146E+01, -8.650532032941677E+01, 2.401409824083091E+01,
                            -1.231739572450155E+00, 1.208650973866179E-03, -5.395239384953000E-06
                          };
    /*
      double x=xx;
      double y=x;
      double tmp=x+5.5;
      tmp -= (x+0.5)*log(tmp);
      double ser=1.000000000190015;
      for (int j=0;j<=5;++j) ser+=cof[j]/++y;
      return -tmp+log(2.5066282746310005*ser/x);
    */

    double yy = xx;
    double res = 1.000000000190015;
    while(yy < 1.0) {
        res *= yy;
        yy += 1.0;
    }
    sign = (res >= 0) ? 1.0 : -1.0;
    yy -= 1.0;
    double tmp = yy + 5.5;
    tmp -= (yy + 0.5) * log(tmp);
    double ser = 1.000000000190015;
    for(int j = 0; j <= 5; ++j) ser += cof[j] / ++yy;

    return -tmp + log(2.5066282746310005 * ser) - log(fabs(res));

}


map<int, int> legacy::mappingVector(double func(int, int), int nl, int nh) {
    map<int, int> tmp;
    for(int i = nl; i <= nh; ++i) tmp[i] = i;
    for(int i = nl; i <= nh - 1; ++i) {
        double fp = func(tmp[i], tmp[nh]);
        for(int j = nh - 1; j >= i; --j) {
            double f = func(tmp[i], tmp[j]);
            if((!f) && (fp)) {
                swap(tmp[j], tmp[j + 1]);
            }
            else {
                fp = f;
            }
        }
    }
    return tmp;
}

map<int, int> legacy::reverseMappingVector(map<int, int>& mapping) {
    map<int, int> tmp;
    for(auto it : mapping) {
        tmp[it.second] = it.first;
    }
    return tmp;
}


double*** legacy::reorderedMatrix(double func(int, int), const map<int, int>& mapping, int first, int last, int &nBlocks, int *&blockFirst, int *&blockLast) {
    nBlocks = 1;
    map<int, double> offDiag;
    auto previous = mapping.begin();
    for(auto it = ++mapping.begin(); it != mapping.end(); ++it) {
        if((offDiag[it->first] = func(previous->second, it->second)) == 0) ++nBlocks;
        previous = it;
    }

    if(!nBlocks) return nullptr;

    blockFirst = new int [nBlocks];
    blockLast = new int [nBlocks];
    double ***blockMatrix = new double** [nBlocks];
    for(int b = 0, n = first; b < nBlocks; ++b) {
        blockFirst[b] = n++;
        while((n <= last) && (offDiag.at(n))) {
            ++n;
        }
        blockLast[b] = n - 1;
        size_t blockSize = blockLast[b] - blockFirst[b] + 1;
        blockMatrix[b] = newArray<double>(blockSize, blockSize);

        for(int i = blockFirst[b]; i <= blockLast[b]; ++i) {
            //          cout << "loop1:  " << b << " -> " << i << endl;
            blockMatrix[b][i - blockFirst[b]][i - blockFirst[b]] = func(mapping.at(i), mapping.at(i));
        }
        for(int i = blockFirst[b] + 1; i <= blockLast[b]; ++i) {
            //         cout << "loop2:  " << b << " -> " << i << endl;
            blockMatrix[b][i - 1 - blockFirst[b]][i - blockFirst[b]] = blockMatrix[b][i - blockFirst[b]][i - 1 - blockFirst[b]] = offDiag[i];
        }
        for(int i = blockFirst[b]; i <= blockLast[b] - 2; ++i)
            for(int j = i + 2; j <= blockLast[b]; ++j) {
                // cout << "loop3:  " << b << " -> " << i << " -> " << j << endl;
                blockMatrix[b][i - blockFirst[b]][j - blockFirst[b]] = blockMatrix[b][j - blockFirst[b]][i - blockFirst[b]] = func(mapping.at(i), mapping.at(j));
            }
    }
    return blockMatrix;
}

void legacy::blockwiseSVD(double ***blockMatrix, double *singular_values, int nBlocks, int *blockFirst, int *blockLast) {
    size_t offset = 0;
    for(int b = 0; b < nBlocks; ++b) {
        int blockSize = blockLast[b] - blockFirst[b] + 1;
        if(blockSize > 1) {
            double **v = newArray<double>(blockSize, blockSize);
            redux::math::svd(*blockMatrix[b], blockSize, blockSize, singular_values + offset, *v);
            delArray(v);
        }
        else {
            singular_values[offset] = blockMatrix[b][0][0];
            blockMatrix[b][0][0] = 1.0;
        }
        offset += blockSize;
    }
}

PupilMode::KL_cfg* legacy::klCoefficients(double ***blockMatrix, double *values, int *first_in_block, int *last_in_block, const std::map<int, int>& mapping, const std::map<int, int>& rmapping, int firstMode, int lastMode) {
    PupilMode::KL_cfg *cfs = new PupilMode::KL_cfg [lastMode - firstMode + 1];
    for(int i = firstMode; i <= lastMode; ++i) {
        int im = rmapping.at(i), s;
        for(s = 0; (first_in_block[s] > im) || (last_in_block[s] < im); ++s);
        int n = last_in_block[s] - first_in_block[s] + 1;
        cfs[i - firstMode].covariance = values[im - firstMode];
        for(int m = 0; m < n; ++m) {
            int j = m + first_in_block[s];
            cfs[i - firstMode].zernikeWeights.push_back(make_pair(mapping.at(j), blockMatrix[s][m][im - first_in_block[s]]));
        }
    }
    return cfs;
}

namespace {

    using legacy::gammln;

    /*! For Zernike polynomials: convert Noll's sequential index to m,n form
    * @note: does not return the sign for negative m values (odd j).
    */
    // void noll_to_mn( int j, int& m, int& n );

// For Zernike polynomials: convert Noll's sequential index to m,n form
// NOTE: does not return the sign for negative m values (odd j).
    void noll_to_mn(int j, int& m, int& n) {
        n = 0;
        int len = 1;
        for(int i = 1; len < j; ++i) {
            len += (n = i) + 1;
        }
        int dl = n + 1 - len + j;
        m = 2 * ((dl + (n % 2)) / 2) + !(n % 2) - 1;
    }



// TODO: cleanup
    double zernikeCovariance(int j) {
        if(j < 2) return 0.0;
        int m, n;
        noll_to_mn(j, m, n);
        int n1 = n;
//  ; Now deal with the numerical terms: Dai
        double tmp;
        double k = pow(4.8 * exp(gammln(6.0 / 5.0, tmp)), 5.0 / 6.0) * exp(gammln(14.0 / 3.0, tmp) + 2.0 * gammln(11.0 / 6.0, tmp)) / (pow(2.0, (8.0 / 3.0)) * redux::PI);
        k *= pow(-1.0, (double)((n + n1 - 2 * m) / 2)) * sqrt((double)(n + 1) * (n1 + 1));
        double g1_sgn, g2_sgn, g3_sgn, g4_sgn;
        double g1 = gammln(((double)(n + n1) - 5.0 / 3.0) / 2.0, g1_sgn);
        double g2 = gammln(((double)(n - n1) + 17.0 / 3.0) / 2.0, g2_sgn);
        double g3 = gammln(((double)(n1 - n) + 17.0 / 3.0) / 2.0, g3_sgn);
        double g4 = gammln(((double)(n + n1) + 23.0 / 3.0) / 2.0, g4_sgn);
        return k * exp(g1 - g2 - g3 - g4) * g1_sgn * g2_sgn * g3_sgn * g4_sgn;
    }

    double zernikeCovariance(int i, int j) {
        //cout << "zernikeCovariance("<<i<<","<<j<<")" << endl;
        if((i < 2) || (j < 2)) return 0.0;
        int m, n, o, p;
        noll_to_mn(i, m, n);
        noll_to_mn(j, o, p);
        if(m != o) return 0.0;
        if(m) if((i + j) % 2) return 0.0;

        //cout << "zernikeCovariance2("<<i<<","<<j<<")" << endl;
//  ; Now deal with the numerical terms: Dai
        double tmp;
        double k = pow(4.8 * exp(gammln(6.0 / 5.0, tmp)), 5.0 / 6.0) * exp(gammln(14.0 / 3.0, tmp) + 2.0 * gammln(11.0 / 6.0, tmp)) / (pow(2.0, (8.0 / 3.0)) * redux::PI);
        k *= pow(-1.0, (double)((n + p - 2 * m) / 2)) * sqrt((double)((n + 1) * (p + 1)));
        double g1_sgn, g2_sgn, g3_sgn, g4_sgn;
        double g1 = gammln(((double)(n + p) - 5.0 / 3.0) / 2.0, g1_sgn);
        double g2 = gammln(((double)(n - p) + 17.0 / 3.0) / 2.0, g2_sgn);
        double g3 = gammln(((double)(p - n) + 17.0 / 3.0) / 2.0, g3_sgn);
        double g4 = gammln(((double)(n + p) + 23.0 / 3.0) / 2.0, g4_sgn);

        // cout << "zernikeCovariance3("<<i<<","<<j<<")   v = " << k * exp( g1 - g2 - g3 - g4 ) * g1_sgn * g2_sgn * g3_sgn * g4_sgn << endl;
        return k * exp(g1 - g2 - g3 - g4) * g1_sgn * g2_sgn * g3_sgn * g4_sgn;

    }

}

PupilMode::KL_cfg* legacy::klConfig(int kl_min_mode, int kl_max_mode) {

    map<int, int> mapping = mappingVector(zernikeCovariance, kl_min_mode, kl_max_mode);
    map<int, int> rmapping = reverseMappingVector(mapping);
    int nBlocks, *first, *last;
    double ***blockMatrix = reorderedMatrix(zernikeCovariance, mapping, kl_min_mode, kl_max_mode, nBlocks, first, last);
    double *singular_values = new double [kl_max_mode - kl_min_mode + 1];
    blockwiseSVD(blockMatrix, singular_values, nBlocks, first, last);
    PupilMode::KL_cfg* ret = klCoefficients(blockMatrix, singular_values, first, last, mapping, rmapping, kl_min_mode, kl_max_mode);

    delete[] singular_values;
    for(int i = 0; i < nBlocks; ++i) {
        delArray(blockMatrix[i]);
    }
    delete[] blockMatrix;
    return ret;
}


namespace {
    double pythag(double a, double b) {
        float absa, absb;
        absa = fabs(a);
        absb = fabs(b);
        if(absa > absb)
            return absa * sqrt(1.0 + sqr(absb / absa));
        else
            return(absb == 0.0) ? 0.0 : absb * sqrt(1.0 + sqr(absa / absb));
    }
}
void legacy::svdcmp(double * const * const aaaa, int m, int n, double *wwww, double * const * const vvvv) {

    double g, scale, anorm;
    double *rv1 = new double [n];
    g = scale = anorm = 0.0;
    for(int i = 0; i < n; ++i) {
        int l = i + 1;
        rv1[i] = scale * g;
        g = scale = 0.0;
        if(i < m) {
            for(int k = i; k < m; ++k) scale += fabs(aaaa[k][i]);
            if(scale) {
                double s = 0.0;
                for(int k = i; k < m; ++k) {
                    aaaa[k][i] /= scale;
                    s += aaaa[k][i] * aaaa[k][i];
                }
                double f = aaaa[i][i];
                g = -sign(sqrt(s), f);
                double h = f * g - s;
                aaaa[i][i] = f - g;
                for(int j = l; j < n; ++j) {
                    double sum = 0.0;
                    for(int k = i; k < m; ++k) sum += aaaa[k][i] * aaaa[k][j];
                    double fct = sum / h;
                    for(int k = i; k < m; ++k) aaaa[k][j] += fct * aaaa[k][i];
                }
                for(int k = i; k < m; ++k) aaaa[k][i] *= scale;
            }
        }
        wwww[i] = scale * g;
        g = scale = 0.0;
        if((i < m) && (i != (n - 1))) {
            for(int k = l; k < n; ++k) scale += fabs(aaaa[i][k]);
            if(scale) {
                double s = 0.0;
                for(int k = l; k < n; ++k) {
                    aaaa[i][k] /= scale;
                    s += aaaa[i][k] * aaaa[i][k];
                }
                double f = aaaa[i][l];
                g = -sign(sqrt(s), f);
                double h = f * g - s;
                aaaa[i][l] = f - g;
                for(int k = l; k < n; ++k) rv1[k] = aaaa[i][k] / h;
                for(int j = l; j < m; ++j) {
                    double sum = 0.0;
                    for(int k = l; k < n; ++k) sum += aaaa[j][k] * aaaa[i][k];
                    for(int k = l; k < n; ++k) aaaa[j][k] += sum * rv1[k];
                }
                for(int k = l; k < n; ++k) aaaa[i][k] *= scale;
            }
        }
        anorm = max(anorm, (fabs(wwww[i]) + fabs(rv1[i])));
    }
    {
        double f;
        for(int i = n - 1, l; i >= 0; --i) {
            if(i < (n - 1)) {
                if(f) {
                    for(int j = l; j < n; ++j) vvvv[j][i] = (aaaa[i][j] / aaaa[i][l]) / f;
                    for(int j = l; j < n; ++j) {
                        double sum = 0.0;
                        for(int k = l; k < n; ++k) sum += aaaa[i][k] * vvvv[k][j];
                        for(int k = l; k < n; ++k) vvvv[k][j] += sum * vvvv[k][i];
                    }
                }
                for(int j = l; j < n; ++j) vvvv[i][j] = vvvv[j][i] = 0.0;
            }
            vvvv[i][i] = 1.0;
            f = rv1[i];
            l = i;
        }
    }
    for(int i = min(m, n) - 1; i >= 0; --i) {
        int l = i + 1;
        g = wwww[i];
        for(int j = l; j < n; ++j) aaaa[i][j] = 0.0;
        if(g) {
            g = 1.0 / g;
            for(int j = l; j < n; ++j) {
                double sum = 0.0;
                for(int k = l; k < m; ++k) sum += aaaa[k][i] * aaaa[k][j];
                double f = (sum / aaaa[i][i]) * g;
                for(int k = i; k < m; ++k) aaaa[k][j] += f * aaaa[k][i];
            }
            for(int j = i; j < m; ++j) aaaa[j][i] *= g;
        }
        else for(int j = i; j < m; ++j) aaaa[j][i] = 0.0;
        ++aaaa[i][i];
    }
    for(int k = n - 1; k >= 0; --k) {
        for(int its = 1; its <= 30; ++its) {
            int flag = 1, nm, l;
            for(l = k; l >= 0; --l) {
                nm = l - 1;
                if((double)(fabs(rv1[l]) + anorm) == anorm) {
                    flag = 0;
                    break;
                }
                if((double)(fabs(wwww[nm]) + anorm) == anorm) break;
            }
            if(flag) {
                double c = 0.0, s = 1.0;
                for(int i = l; i <= k; ++i) {
                    double f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if((double)(fabs(f) + anorm) == anorm) break;
                    g = wwww[i];
                    double h = pythag(f, g);
                    wwww[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for(int j = 0; j < m; ++j) {
                        double y = aaaa[j][nm];
                        double z = aaaa[j][i];
                        aaaa[j][nm] = y * c + z * s;
                        aaaa[j][i] = z * c - y * s;
                    }
                }
            }
            double z = wwww[k];
            if(l == k) {
                if(z < 0.0) {
                    wwww[k] = -z;
                    for(int j = 0; j < n; ++j) vvvv[j][k] = -vvvv[j][k];
                }
                break;
            }
            if(its == 30) exit(fprintf(stderr, "no convergence in 30 svdcmp iterations"));
            double x = wwww[l];
            double y = wwww[nm = k - 1];
            g = rv1[nm];
            double h = rv1[k];
            double f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x;
            double c = 1.0, s = 1.0;
            for(int j = l; j <= nm; ++j) {
                int i = j + 1;
                g = rv1[i];
                y = wwww[i];
                h = s * g;
                g *= c;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for(int jj = 0; jj < n; ++jj) {
                    x = vvvv[jj][j];
                    z = vvvv[jj][i];
                    vvvv[jj][j] = x * c + z * s;
                    vvvv[jj][i] = z * c - x * s;
                }
                z = pythag(f, h);
                wwww[j] = z;
                if(z) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for(int jj = 0; jj < m; ++jj) {
                    y = aaaa[jj][j];
                    z = aaaa[jj][i];
                    aaaa[jj][j] = y * c + z * s;
                    aaaa[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            wwww[k] = x;
        }
    }
    delete[] rv1;
}

namespace {

    template <typename T>
    void sort(T* x, T* y, int nx) { // inefficient algorithm, but quick and easy and nx is small anyway
        for(int i = 0; i < nx - 1; ++i)
            for(int j = 0; j < nx - i; ++j)
                if(x[j] > x[j + 1]) {
                    swap(x[j], x[j + 1]);
                    swap(y[j], y[j + 1]);
                }
    }

    int *reorder(double **mat, int nRows, int nCols) {
// associate variables
        int **asc = newArray<int>(nRows, nRows);
        memset(*asc, 0, nRows * nRows * sizeof(int));
        for(int col = 0; col < nCols; ++col) {
            for(int row = 0; row < nRows; ++row) {
                if(mat[row][col]) {
                    for(int row2 = row + 1; row2 < nRows; ++row2) {
                        if(mat[row2][col]) {
                            asc[row][row2] += 1;
                            asc[row2][row] += 1;
                        }
                    }
                }
            }
        }
//
        int *idx = new int [nRows];
        for(int i = 0; i < nRows; ++i) idx[i] = i;
//
        for(int i = nRows - 1; i >= 0; --i) {
            int nn = 0, ij = 0;
            for(int j = 0; j < nRows; ++j) {
                if(asc[i][j]) {
                    nn += asc[i][j];
                    ij = j;
                }
            }
            if(nn == 1) { // only one association: put together
                int ii = 0, ia = 0;
                for(int j = 0; j < nRows; ++j) {
                    if(ij == idx[j]) ia = j;
                    if(i == idx[j]) ii = j;
                }
                if((ia + 2) < nRows) {
                    int tmp = idx[ii];
                    memmove(idx + ia + 2, idx + ia + 1, (ii - ia - 1)*sizeof(int));
                    idx[ia + 1] = tmp;
                }
            }
        }
        delArray(asc);
// column swap to minimize distance between singly associated variables
        double **cn = newArray<double>(nRows, nCols);
        for(int i = 0; i < nRows; ++i)
            for(int k = 0; k < nCols; ++k) cn[i][k] = mat[idx[i]][k];
// sort columns in descending order
        int *len = new int [nCols];
        int *ldx = new int [nCols];
        for(int k = 0; k < nCols; ++k) {
            ldx[k] = k;
            for(int i = 0; i < nRows; ++i)
                if(cn[i][k]) { // first nonzero element
                    len[k] = i;
                    break;
                }
        }
// sort according to index of first nonzero element...
        sort(len, ldx, nCols);
// row swap to minimize Q fillup...
        for(int i = 0; i < nRows; ++i)
            for(int k = 0; k < nCols; ++k) mat[i][k] = cn[i][ldx[nCols - k + 1]];
//
        delArray(cn);
        delete[](len);
        return idx;
    }

    void reorder(double **c, int m, int n, int *idx, int dir) {

        double **cn = newArray<double>(m, n);
//
        if(dir > 0)
            for(int k = 0; k < n; ++k) for(int i = 0; i < m; ++i) cn[i][k] = c[idx[i]][k];
        else
            for(int k = 0; k < n; ++k) for(int i = 0; i < m; ++i) cn[idx[i]][k] = c[i][k];
//
        memcpy(*c, *cn, m * n * sizeof(double));
        delArray(cn);
    }
}

double** legacy::qr(double **c, int rows, int cols, uint8_t fast_QR) {
    
    int *idx = 0;

//  int ds[]={0,n,m};
//  ana_fzwrite((uint8_t*)(c[1]+1),(char*)"constr.raw.f0",ds,2,(char*)"bla",FLOAT64,io);

    if(fast_QR) {
        cout << "qr: using fast QR decomposition of " << rows << "x" << cols << endl;
        idx = reorder(c, rows, cols);
    }

//  ana_fzwrite((uint8_t*)(c[1]+1),(char*)"constr.srt.f0",ds,2,(char*)"bla",FLOAT64,io);
//
    double **q = newArray<double>(rows, rows);
    memset(*q, 0, rows * rows * sizeof(double));
    for(int row = 0; row < rows; ++row) q[row][row] = 1.0;
    double **r = c, *v = new double [rows], *u = new double [rows];
    for(int col = 0; col < cols; ++col) {
        //cout << "qr: col " << col << endl;
        double sum = 0.0;
        for(int row = col; row < rows; ++row) sum += sqr(v[row] = r[row][col]); // upper triangle + diagonal ?
        //cout << " 1"  << flush;
        v[col] += sign(r[col][col]) * sqrt(sum);            // alpha=sqrt(sum),U=X-alpha*e1 (where sign(alpha)=-sign(X[1]) for precision reasons)
        //cout << " 2"  << flush;
        sum = 0.0;
        for(int l = col; l < rows; ++l) sum += sqr(v[l]); // ||v||
        //cout << " 3"  << flush;
        sum /= 2.0;
        for(int i = 0; i < rows; ++i) {
            for(int j = col; j < rows; ++j) {               // Q.v.vT
                u[j] = 0.0;
                if(v[j]) for(int l = col; l < rows; ++l) if(v[l]) if(q[i][l]) u[j] += q[i][l] * v[l] * v[j];
            }
            for(int j = col; j < rows; ++j) if(u[j]) q[i][j] -= u[j] / sum;
        }
        //cout << " 4" << flush;
        for(int j = col + 1; j < cols; ++j) {             // only interested in r from index k+1..n
            for(int i = col + 1; i < rows; ++i) {           // vvT.R
                u[i] = 0.0;
        //cout << " 5"  << flush;
                if(v[i]) for(int l = col; l < rows; ++l) if(v[l]) if(r[l][j]) u[i] += v[i] * v[l] * r[l][j];
            }
            for(int i = col + 1; i < rows; ++i) if(u[i]) r[i][j] -= u[i] / sum;
        }
        //cout << " 6"  << flush;
        //cout << endl;
    }
    delete[](u);
    delete[](v);
//
    if(fast_QR) {
        reorder(q, rows, rows, idx, -1);
        delete[](idx);
    }

//  int ds[]={0,n,m};
//  ana_fzwrite((uint8_t*)(q[1]+1),(char*)"constr.qrd.f0",ds,2,(char*)"bla",FLOAT64,io);
//  exit(1);
    return q;
}
