#include "redux/momfbd/modes.hpp"

#include "redux/constants.hpp"
#include "redux/momfbd/modecache.hpp"
#include "redux/momfbd/legacy.hpp"
#include "redux/util/arrayutil.hpp"

#include <cmath>
#include <iostream>

using namespace redux::momfbd;
using namespace redux::momfbd::legacy;
using namespace redux::util;
using namespace std;

namespace {

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
        double k = pow(4.8 * exp(lgamma_r(6.0 / 5.0, nullptr)), 5.0 / 6.0) * exp(lgamma_r(14.0 / 3.0, nullptr) + 2.0 * lgamma_r(11.0 / 6.0, nullptr)) / (pow(2.0, (8.0 / 3.0)) * redux::PI);
        k *= pow(-1.0, (double)((n + n1 - 2 * m) / 2)) * sqrt((double)(n + 1) * (n1 + 1));
        int g1_sgn, g2_sgn, g3_sgn, g4_sgn;
        double g1 = lgamma_r(((double)(n + n1) -  5.0 / 3.0) / 2.0, &g1_sgn);
        double g2 = lgamma_r(((double)(n - n1) + 17.0 / 3.0) / 2.0, &g2_sgn);
        double g3 = lgamma_r(((double)(n1 - n) + 17.0 / 3.0) / 2.0, &g3_sgn);
        double g4 = lgamma_r(((double)(n + n1) + 23.0 / 3.0) / 2.0, &g4_sgn);
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
        double k = pow(4.8 * exp(lgamma_r(6.0 / 5.0, nullptr)), 5.0 / 6.0) * exp(lgamma_r(14.0 / 3.0, nullptr) + 2.0 * lgamma_r(11.0 / 6.0, nullptr)) / (pow(2.0, (8.0 / 3.0)) * redux::PI);
        k *= pow(-1.0, (double)((n + p - 2 * m) / 2)) * sqrt((double)((n + 1) * (p + 1)));
        int g1_sgn, g2_sgn, g3_sgn, g4_sgn;
        double g1 = lgamma_r(((double)(n + p) -  5.0 / 3.0) / 2.0, &g1_sgn);
        double g2 = lgamma_r(((double)(n - p) + 17.0 / 3.0) / 2.0, &g2_sgn);
        double g3 = lgamma_r(((double)(p - n) + 17.0 / 3.0) / 2.0, &g3_sgn);
        double g4 = lgamma_r(((double)(n + p) + 23.0 / 3.0) / 2.0, &g4_sgn);

        // cout << "zernikeCovariance3("<<i<<","<<j<<")   v = " << k * exp( g1 - g2 - g3 - g4 ) * g1_sgn * g2_sgn * g3_sgn * g4_sgn << endl;
        return k * exp(g1 - g2 - g3 - g4) * g1_sgn * g2_sgn * g3_sgn * g4_sgn;

    }

}


Modes::PupilMode::PupilMode(int modeNumber, int nPoints, double r_c, double lambda, double angle) : Array<double> (nPoints, nPoints) {      // Zernike


    static ModeCache& cache = ModeCache::getCache();

//    cout << "|Generating Z mode: (" << modeNumber << ") ..." << flush;
    double** modePtr = makePointers(ptr(), nPoints, nPoints);

    if(modeNumber == 1) {

        Array<double>::operator=(1.0);

    }
    else {

        zero();
        int m, n;
        noll_to_mn(modeNumber, m, n);

        const vector<double>& coeff = cache.zernikeRadialPolynomial(m, n);
        const ModeCache::Grid& grid = cache.grid(nPoints);
        float** distPtr = grid.distance.get();
        float** aPtr = grid.angle.get();

        shared_ptr<double*> r = sharedArray<double> (nPoints, nPoints);
        shared_ptr<double*> r2 = sharedArray<double> (nPoints, nPoints);
        size_t blockSize = nPoints * nPoints * sizeof(double);
        double** rPtr = r.get();
        double** r2Ptr = r2.get();

        for(int i = 0; i < nPoints; ++i) {
            for(int j = 0; j < nPoints; ++j) {
                double tmp = distPtr[i][j] / r_c;
                r2Ptr[i][j] = tmp * tmp;
                if(m == 0) rPtr[i][j] = 1;
                else if(m == 1) rPtr[i][j] = tmp;
                else rPtr[i][j] = pow(tmp, m);
                modePtr[i][j] = rPtr[i][j] * coeff[0];
            }
        }

        // generate polynomial part
        for(auto it = coeff.begin() + 1; it < coeff.end(); it++) {
            for(int i = 0; i < nPoints; ++i) {
                for(int j = 0; j < nPoints; ++j) {
                    rPtr[i][j] *= r2Ptr[i][j];
                    modePtr[i][j] += rPtr[i][j] * *it;
                }
            }
        }

        // Angular component
        if(m == 0) {
            double sf = sqrt(n + 1);
            for(int i = 0; i < nPoints; ++i) {
                for(int j = 0; j < nPoints; ++j) {
                    modePtr[i][j] *= sf;
                }
            }
        }
        else if(modeNumber % 2) {
            double sf = sqrt((double) 2 * (n + 1));
            for(int i = 0; i < nPoints; ++i) {
                for(int j = 0; j < nPoints; ++j) {
                    modePtr[i][j] *= sf * sin(m * (aPtr[i][j] + angle));
                }
            }
        }
        else {
            double sf = sqrt((double) 2 * (n + 1));
            for(int i = 0; i < nPoints; ++i) {
                for(int j = 0; j < nPoints; ++j) {
                    modePtr[i][j] *= sf * cos(m * (aPtr[i][j] + angle));
                }
            }
        }
    }

    // normalize
    double norm = 0.0, N = 0.0, dx = 0.5 / r_c, dy = 0.5 / r_c;
    int half = nPoints / 2;
    for(int i = 0; i < nPoints; ++i) {
        double xl = fabs(i - half) / r_c - dx;
        double xh = xl + 2 * dx;
        double xls = xl * xl;
        double xhs = xh * xh;
        for(int j = 0; j < nPoints; ++j) {
            double yl = fabs(j - half) / r_c - dy;
            double yh = yl + 2 * dy;
            double yhs = yh * yh;
            double rsl = xls + yl * yl;
            double rsh = xhs + yhs;
            if(rsl <= 1.0) {    // good pixel
                if(rsh < 1.0) {    // full pixel
                    norm += modePtr[i][j] * modePtr[i][j];
                    N += 1.0;
                }
                else {           // partial pixel
                    double x2 = sqrt(max(1.0 - yhs, (double) 0.0));
                    double y3 = sqrt(max(1.0 - xhs, (double) 0.0));
                    double f = (xh > yh) ? (yh - yl) * (min(xh, max(xl, x2)) - xl) / (4 * dx * dy) :
                               (xh - xl) * (min(yh, max(yl, y3)) - yl) / (4 * dx * dy);
                    norm += f * modePtr[i][j] * modePtr[i][j];
                    N += f;
                }
            }
        }
    }
    norm = 1 / (sqrt(norm / N) * lambda);
    for(int i = 0; i < nPoints; ++i) {
        for(int j = 0; j < nPoints; ++j) {
            modePtr[i][j] *= norm;
        }
    }

    delPointers(modePtr);

//    cout << " zdone|" << flush;

}


Modes::PupilMode::PupilMode(int firstMode, int lastMode, int klModeNumber, int nPoints, double r_c, double lambda, double angle, double cutoff) : Array<double> (nPoints, nPoints) {
//Modes::PupilMode::PupilMode( KL_cfg* kl_cfg, int klModeIndex, int nPoints, PupilMode** modeCache, double r_c, double lambda, double angle, double cutoff ) : Array<double> ( nPoints, nPoints )  {

    if(firstMode > lastMode) swap(firstMode, lastMode);

    if(klModeNumber < firstMode || klModeNumber > lastMode) {
        throw invalid_argument("klModeNumber (" + to_string(klModeNumber) +
                               ") is not in the range [ firstMode (" + to_string(firstMode) +
                               "), lastMode (" + to_string(lastMode) + ")]");
    }

    zero();
        
    static ModeCache& cache = ModeCache::getCache();
    const std::map<int, Modes::KL_cfg>& kle = cache.karhunenLoeveExpansion(firstMode, lastMode);
  //  cout << "this = "<< hexString(this) << "   ptr = " << hexString(ptr()) << endl;
  //  cout << "Generating KL mode: (" << firstMode << "," << lastMode << "," << klModeNumber << ") ..." << flush;
    double c;
  //  vector<double> bla;
    for(auto & it : kle.at(klModeNumber).zernikeWeights) {
        if(fabs(c = it.second) >= cutoff) {
     //       bla.push_back(c);
            int zernikeModeIndex = it.first;
            //if ( !modeCache[zernikeModeIndex] ) modeCache[zernikeModeIndex] = new PupilMode ( zernikeModeIndex, nPoints, r_c, lambda, angle );
            //this->add ( *modeCache[zernikeModeIndex], c );
      //      cout << "|adding " << c << " * " << zernikeModeIndex << "|" << flush;
            auto ptr = cache.mode(zernikeModeIndex, nPoints, r_c, lambda, angle);
    //        cout << endl << "BLAHA: "<< hexString(ptr.get()) << endl;
            //cout << endl << "BLAHA2: "<< hexString(ptr.get()) << printArray(this->dimensions(),"   dims") << endl;
            //cout << "||" << printArray(ptr->dimensions(),"dims") << "||" << flush;
            //cout << "VV" << ptr->at(0,0) << "VV" << flush;
     //       cout << this->operator()(nPoints/2,nPoints/2) << "|" << flush;
            this->add(*ptr, c);
    //        cout << this->operator()(nPoints/2,nPoints/2) << "|" << flush;
            //cout << "|" << flush;
        }
    }
//    cout << printArray(bla,"   coeffs") << flush;
  //  cout << " done " << printArray(bla, "   coeffs")  << endl;
}


Modes::Modes() : mode(nullptr),  pupil(nullptr),  covar(nullptr) {}

Modes::Modes(KL_cfg* kl_cfg, double lambda, double r_c, int nph_in, int basis, int nm, int *mode_num,
             int nch, int *ndo, int **dorder, int **dtype, int kl_min_mode, int kl_max_mode, double svd_reg, double angle, double **pupil_in) : Modes() {
    nph = nph_in;
    zmin = 0x7FFFFFFF;
    zmax = 0x00000000;
    kmin = 0x7FFFFFFF;
    kmax = 0x00000000;
    first_mode = 0x7FFFFFFF;
    last_mode = 0x00000000;
    double lambda2 = lambda * lambda;
    for(int c = 1; c <= nch; ++c) {
        for(int m = 1; m <= nm; ++m) {    // basis modes
            if(basis == MB_ZERNIKE) {
                zmin = min(zmin, mode_num[m]);
                zmax = max(zmax, mode_num[m]);
            }
            if(basis == MB_KARHUNEN_LOEVE) {
                kmin = min(kmin, mode_num[m]);
                kmax = max(kmax, mode_num[m]);
            }
            first_mode = min(first_mode, mode_num[m]);
            last_mode = max(last_mode, mode_num[m]);
        }
        for(int m = 1; m <= ndo[c]; ++m) {    // diversity modes
            if(dtype[c][m] == MB_ZERNIKE) {
                zmin = min(zmin, dorder[c][m]);
                zmax = max(zmax, dorder[c][m]);
            }
            if(dtype[c][m] == MB_KARHUNEN_LOEVE) {
                kmin = min(kmin, dorder[c][m]);
                kmax = max(kmax, dorder[c][m]);
            }
        }
    }
    mode = new PupilMode** [MB_END];
    mode[0] = new PupilMode* [last_mode - first_mode + 1];
    memset(mode[0], 0, (last_mode - first_mode + 1) *sizeof(PupilMode*));
    if(zmin <= zmax) {
        mode[MB_ZERNIKE] = new PupilMode* [zmax - zmin + 1];
        memset(mode[MB_ZERNIKE], 0, (zmax - zmin + 1) *sizeof(PupilMode*));
    }
    else mode[MB_ZERNIKE] = 0;
    if(kmin <= kmax) {
        mode[MB_KARHUNEN_LOEVE] = new PupilMode* [kmax - kmin + 1];
        memset(mode[MB_KARHUNEN_LOEVE], 0, (kmax - kmin + 1) *sizeof(PupilMode*));
    }
    else mode[MB_KARHUNEN_LOEVE] = 0;
//
    covar = new double* [MB_END];
    covar[0] = new double [last_mode - first_mode + 1];
    if(zmin <= zmax) covar[MB_ZERNIKE] = new double [zmax - zmin + 1];
    else covar[MB_ZERNIKE] = 0;
    if(kmin <= kmax) covar[MB_KARHUNEN_LOEVE] = new double [kmax - kmin + 1];
    else covar[MB_KARHUNEN_LOEVE] = 0;
// use temporary storage z to avoid recomputing too many Zernikes
    PupilMode **z = new PupilMode* [kl_max_mode - kl_min_mode + 1];
    memset(z, 0, (kl_max_mode - kl_min_mode + 1) *sizeof(PupilMode*));
    for(int c = 0; c < nch; ++c) {
        for(int m = 0; m < nm; ++m) {
            int mn = mode_num[m];
            if(!mode[0][mn - first_mode]) {
                if(basis == MB_ZERNIKE) {                  // Zernikes
                    mode[0][mn - first_mode] = (mode[MB_ZERNIKE][mn - zmin] = new PupilMode(mn, nph, r_c, lambda, angle));
                    covar[0][mn - first_mode] = (covar[MB_ZERNIKE][mn - zmin] = sqrt(zernikeCovariance(mn))) * lambda2;
                }
                else {                                          // Karhunen-Loeves for the rest
                    if((mn == 2) || (mn == 3)) {                  // use Zernikes for the tilts
                        mode[0][mn - first_mode] = (mode[MB_KARHUNEN_LOEVE][mn - kmin] = new PupilMode(mn, nph, r_c, lambda, angle));
                        covar[0][mn - first_mode] = (covar[MB_KARHUNEN_LOEVE][mn - kmin] = sqrt(zernikeCovariance(mn))) * lambda2;
                    }
                    else {
                        // mode[0][mn] = ( mode[MB_KARHUNEN_LOEVE][mn] = new PupilMode( lambda, r_c, nph, mn, kl_cfs, svd_reg, z, angle, io ) );
                        mode[0][mn - first_mode] = (mode[MB_KARHUNEN_LOEVE][mn - kmin] = new PupilMode(first_mode, last_mode, mn, nph, r_c, lambda, angle, svd_reg));
                        covar[0][mn - first_mode] = (covar[MB_KARHUNEN_LOEVE][mn - kmin] = sqrt(kl_cfg[mn].covariance)) * lambda2;
                    }
                }
            }
        }
    }
    for(int c = 0; c < nch; ++c) {          // compute diversity orders
        for(int m = 0; m < ndo[c]; ++m) {
            int mn = dorder[c][m];
            int type = dtype[c][m];
            if(type == MB_ZERNIKE && !mode[type][mn - zmin]) {          // use Zernikes for the tilts
                mode[MB_ZERNIKE][mn] = new PupilMode(mn, nph, r_c, lambda, angle);
                covar[MB_ZERNIKE][mn] = sqrt(zernikeCovariance(mn)) * lambda2;
            }
            else if(type == MB_KARHUNEN_LOEVE && !mode[type][mn - kmin]) {                                  // Karhunen-Loeves for the rest
                //mode[MB_KARHUNEN_LOEVE][mn] = new PupilMode( lambda, r_c, nph, mn, svd_reg, z, angle, io );
                mode[MB_KARHUNEN_LOEVE][mn] = new PupilMode(first_mode, last_mode, mn, nph, r_c, lambda, angle, svd_reg);
                covar[MB_KARHUNEN_LOEVE][mn] = sqrt(kl_cfg[mn].covariance) * lambda2;
            }
        }
    }
    for(int m = kl_max_mode - kl_min_mode; m >= 0; --m) if(z[m]) delete z[m];
    delete[] z;
// initialise pupil
    pupil = newArray<double> (nph, nph);
    if(pupil_in) {
        memcpy(*pupil, *pupil_in, nph * nph * sizeof(double));
    }
    else {
        int xo = nph / 2, yo = nph / 2;
        double dx = 0.5 / r_c, dy = 0.5 / r_c;
        for(int x = 0; x < nph; ++x) {
            double xl = fabs((double)(x - xo)) / r_c - dx, xh = fabs((double)(x - xo)) / r_c + dx;
            double xhs = xh * xh;
            for(int y = 0; y < nph; ++y) {
                double yl = fabs((double)(y - yo)) / r_c - dy, yh = fabs((double)(y - yo)) / r_c + dy;
                double yhs = yh * yh;
                double rsl = xl * xl + yl * yl, rsh = xhs + yhs;
                if(rsl <= 1.0) {    // inside pixel
                    if(rsh < 1.0)    // full pixel
                        pupil[x][y] = 1.0;
                    else {           // partial pixel
                        double x2 = sqrt(max(1.0 - yhs, (double) 0.0));
                        double y3 = sqrt(max(1.0 - xhs, (double) 0.0));
                        double f = (xh > yh) ? (yh - yl) * (min(xh, max(xl, x2)) - xl) / (4 * dx * dy) :
                                   (xh - xl) * (min(yh, max(yl, y3)) - yl) / (4 * dx * dy);
                        pupil[x][y] = f;
                    }
                }
                else pupil[x][y] = 0.0;   // outside pixel
            }
        }
    }
//
    area = 0.0;
    for(int x = 0; x < nph; ++x) {
        for(int y = 0; y < nph; ++y) {
            area += pupil[x][y];
        }
    }
}

Modes::~Modes(void) {
    if(mode) {
        if(mode[0]) delete[] mode[0];    // don't delete mode[0][m], they are co-pointers!
        for(int m = zmax - zmin; m >= 0; --m) if(mode[MB_ZERNIKE][m]) delete mode[MB_ZERNIKE][m];
        if(mode[MB_ZERNIKE]) delete[] mode[MB_ZERNIKE];
        for(int m = kmax - kmin; m >= 0; --m) if(mode[MB_KARHUNEN_LOEVE][m]) delete mode[MB_KARHUNEN_LOEVE][m];
        if(mode[MB_KARHUNEN_LOEVE]) delete[] mode[MB_KARHUNEN_LOEVE];
        delete[] mode;
    }
//
    if(covar) {
        if(covar[0]) delete[] covar[0];
        if(covar[MB_ZERNIKE]) delete[] covar[MB_ZERNIKE];
        if(covar[MB_KARHUNEN_LOEVE]) delete[] covar[MB_KARHUNEN_LOEVE];
        delete[] covar;
    }
    if(pupil) delArray(pupil);
}


void Modes::init(KL_cfg* kl_cfg, double lambda, double r_c, int nph_in, int basis, int nm, int *mode_num,
                 int nch, int *ndo, int **dorder, int **dtype, int kl_min_mode, int kl_max_mode, double svd_reg, double angle, double **pupil_in) {
    nph = nph_in;
    zmin = 0x7FFFFFFF;
    zmax = 0x00000000;
    kmin = 0x7FFFFFFF;
    kmax = 0x00000000;
    first_mode = 0x7FFFFFFF;
    last_mode = 0x00000000;
    double lambda2 = lambda * lambda;
    for(int c = 1; c <= nch; ++c) {
        for(int m = 1; m <= nm; ++m) {    // basis modes
            if(basis == MB_ZERNIKE) {
                zmin = min(zmin, mode_num[m]);
                zmax = max(zmax, mode_num[m]);
            }
            if(basis == MB_KARHUNEN_LOEVE) {
                kmin = min(kmin, mode_num[m]);
                kmax = max(kmax, mode_num[m]);
            }
            first_mode = min(first_mode, mode_num[m]);
            last_mode = max(last_mode, mode_num[m]);
        }
        for(int m = 1; m <= ndo[c]; ++m) {    // diversity modes
            if(dtype[c][m] == MB_ZERNIKE) {
                zmin = min(zmin, dorder[c][m]);
                zmax = max(zmax, dorder[c][m]);
            }
            if(dtype[c][m] == MB_KARHUNEN_LOEVE) {
                kmin = min(kmin, dorder[c][m]);
                kmax = max(kmax, dorder[c][m]);
            }
        }
    }
    mode = new PupilMode** [MB_END];
    mode[0] = new PupilMode* [last_mode - first_mode + 1];
    memset(mode[0], 0, (last_mode - first_mode + 1) *sizeof(PupilMode*));
    if(zmin <= zmax) {
        mode[MB_ZERNIKE] = new PupilMode* [zmax - zmin + 1];
        memset(mode[MB_ZERNIKE], 0, (zmax - zmin + 1) *sizeof(PupilMode*));
    }
    else mode[MB_ZERNIKE] = 0;
    if(kmin <= kmax) {
        mode[MB_KARHUNEN_LOEVE] = new PupilMode* [kmax - kmin + 1];
        memset(mode[MB_KARHUNEN_LOEVE], 0, (kmax - kmin + 1) *sizeof(PupilMode*));
    }
    else mode[MB_KARHUNEN_LOEVE] = 0;
//
    covar = new double* [MB_END];
    covar[0] = new double [last_mode - first_mode + 1];
    if(zmin <= zmax) covar[MB_ZERNIKE] = new double [zmax - zmin + 1];
    else covar[MB_ZERNIKE] = 0;
    if(kmin <= kmax) covar[MB_KARHUNEN_LOEVE] = new double [kmax - kmin + 1];
    else covar[MB_KARHUNEN_LOEVE] = 0;
// use temporary storage z to avoid recomputing too many Zernikes
    PupilMode **z = new PupilMode* [kl_max_mode - kl_min_mode + 1];
    memset(z, 0, (kl_max_mode - kl_min_mode + 1) *sizeof(PupilMode*));
    for(int c = 0; c < nch; ++c) {
        for(int m = 0; m < nm; ++m) {
            int mn = mode_num[m];
            if(!mode[0][mn - first_mode]) {
                if(basis == MB_ZERNIKE) {                  // Zernikes
                    mode[0][mn - first_mode] = (mode[MB_ZERNIKE][mn - zmin] = new PupilMode(mn, nph, r_c, lambda, angle));
                    covar[0][mn - first_mode] = (covar[MB_ZERNIKE][mn - zmin] = sqrt(zernikeCovariance(mn))) * lambda2;
                }
                else {                                          // Karhunen-Loeves for the rest
                    if((mn == 2) || (mn == 3)) {                  // use Zernikes for the tilts
                        mode[0][mn - first_mode] = (mode[MB_KARHUNEN_LOEVE][mn - kmin] = new PupilMode(mn, nph, r_c, lambda, angle));
                        covar[0][mn - first_mode] = (covar[MB_KARHUNEN_LOEVE][mn - kmin] = sqrt(zernikeCovariance(mn))) * lambda2;
                    }
                    else {
                        // mode[0][mn] = ( mode[MB_KARHUNEN_LOEVE][mn] = new PupilMode( lambda, r_c, nph, mn, kl_cfs, svd_reg, z, angle, io ) );
                        mode[0][mn - first_mode] = (mode[MB_KARHUNEN_LOEVE][mn - kmin] = new PupilMode(first_mode, last_mode, mn, nph, r_c, lambda, angle, svd_reg));
                        covar[0][mn - first_mode] = (covar[MB_KARHUNEN_LOEVE][mn - kmin] = sqrt(kl_cfg[mn].covariance)) * lambda2;
                    }
                }
            }
        }
    }
    for(int c = 0; c < nch; ++c) {          // compute diversity orders
        for(int m = 0; m < ndo[c]; ++m) {
            int mn = dorder[c][m];
            int type = dtype[c][m];
            if(type == MB_ZERNIKE && !mode[type][mn - zmin]) {          // use Zernikes for the tilts
                mode[MB_ZERNIKE][mn] = new PupilMode(mn, nph, r_c, lambda, angle);
                covar[MB_ZERNIKE][mn] = sqrt(zernikeCovariance(mn)) * lambda2;
            }
            else if(type == MB_KARHUNEN_LOEVE && !mode[type][mn - kmin]) {                                  // Karhunen-Loeves for the rest
                //mode[MB_KARHUNEN_LOEVE][mn] = new PupilMode( lambda, r_c, nph, mn, svd_reg, z, angle, io );
                mode[MB_KARHUNEN_LOEVE][mn] = new PupilMode(first_mode, last_mode, mn, nph, r_c, lambda, angle, svd_reg);
                covar[MB_KARHUNEN_LOEVE][mn] = sqrt(kl_cfg[mn].covariance) * lambda2;
            }
        }
    }
    for(int m = kl_max_mode - kl_min_mode; m >= 0; --m) if(z[m]) delete z[m];
    delete[] z;
// initialise pupil
    pupil = newArray<double> (nph, nph);
    if(pupil_in) {
        memcpy(*pupil, *pupil_in, nph * nph * sizeof(double));
    }
    else {
        int xo = nph / 2, yo = nph / 2;
        double dx = 0.5 / r_c, dy = 0.5 / r_c;
        for(int x = 0; x < nph; ++x) {
            double xl = fabs((double)(x - xo)) / r_c - dx, xh = fabs((double)(x - xo)) / r_c + dx;
            double xhs = xh * xh;
            for(int y = 0; y < nph; ++y) {
                double yl = fabs((double)(y - yo)) / r_c - dy, yh = fabs((double)(y - yo)) / r_c + dy;
                double yhs = yh * yh;
                double rsl = xl * xl + yl * yl, rsh = xhs + yhs;
                if(rsl <= 1.0) {    // inside pixel
                    if(rsh < 1.0)    // full pixel
                        pupil[x][y] = 1.0;
                    else {           // partial pixel
                        double x2 = sqrt(max(1.0 - yhs, (double) 0.0));
                        double y3 = sqrt(max(1.0 - xhs, (double) 0.0));
                        double f = (xh > yh) ? (yh - yl) * (min(xh, max(xl, x2)) - xl) / (4 * dx * dy) :
                                   (xh - xl) * (min(yh, max(yl, y3)) - yl) / (4 * dx * dy);
                        pupil[x][y] = f;
                    }
                }
                else pupil[x][y] = 0.0;   // outside pixel
            }
        }
    }
//
    area = 0.0;
    for(int x = 0; x < nph; ++x) {
        for(int y = 0; y < nph; ++y) {
            area += pupil[x][y];
        }
    }
}

