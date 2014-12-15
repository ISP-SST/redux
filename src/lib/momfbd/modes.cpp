#include "redux/momfbd/modes.hpp"

#include "redux/constants.hpp"
#include "redux/momfbd/modecache.hpp"
//#include "redux/momfbd/legacy.hpp"
#include "redux/util/arrayutil.hpp"

#include <cmath>
#include <iostream>

using namespace redux::momfbd;
//using namespace redux::momfbd::legacy;
using namespace redux::util;
using namespace std;

namespace {

    /*! For Zernike polynomials: convert Noll's sequential index to m,n form
    * @note: does not return the sign for negative m values (odd j).
    */
    void noll_to_mn(int j, int& m, int& n) {
        n = 0;
        int len = 1;
        for(int i = 1; len < j; ++i) {
            len += (n = i) + 1;
        }
        int dl = n + 1 - len + j;
        m = 2 * ((dl + (n % 2)) / 2) + !(n % 2) - 1;
    }

}


void calculate(void) {
    
}

PupilMode::PupilMode(int modeNumber, int nPoints, double r_c, double lambda, double angle) : Array<double> (nPoints, nPoints) {      // Zernike


    static ModeCache& cache = ModeCache::getCache();
    const std::pair<Array<float>,float>& pupil = cache.pupil(nPoints,r_c);

    zero();

    double** modePtr = makePointers(ptr(), nPoints, nPoints);
    const float** pupPtr = makePointers(pupil.first.ptr(), nPoints, nPoints);

    if(modeNumber == 1) {
        add(pupil.first,1.0/(pupil.second*lambda));
    } else {

        int m, n;
        noll_to_mn(modeNumber, m, n);

        const vector<double>& coeff = cache.zernikeRadialPolynomial(m, n);
        const redux::image::Grid& grid = cache.grid(nPoints,redux::PointF(nPoints/2+0.5,nPoints/2+0.5));
        float** distPtr = grid.distance.get();
        float** aPtr = grid.angle.get();

        shared_ptr<double*> r = sharedArray<double> (nPoints, nPoints);
        shared_ptr<double*> r2 = sharedArray<double> (nPoints, nPoints);
        double** rPtr = r.get();
        double** r2Ptr = r2.get();
        double squared_sum(0);

        for(int i = 0; i < nPoints; ++i) {
            for(int j = 0; j < nPoints; ++j) {
                if(pupPtr[i][j]>0) {
                    double tmp = distPtr[i][j] / r_c;
                    //if(tmp>1) continue;
                    r2Ptr[i][j] = tmp * tmp;
                    if(m == 0) rPtr[i][j] = 1;
                    else if(m == 1) rPtr[i][j] = tmp;
                    else rPtr[i][j] = pow(tmp, m);
                    modePtr[i][j] = rPtr[i][j] * coeff[0];
                }
            }
        }

        // generate polynomial part
        for(auto it = coeff.begin() + 1; it < coeff.end(); it++) {
            for(int i = 0; i < nPoints; ++i) {
                for(int j = 0; j < nPoints; ++j) {
                    if(pupPtr[i][j]>0) {
                        rPtr[i][j] *= r2Ptr[i][j];
                        modePtr[i][j] += rPtr[i][j] * *it;
                    }
                }
            }
        }

        // Angular component
        if(m == 0) {
            double sf = sqrt(n + 1);
            for(int i = 0; i < nPoints; ++i) {
                for(int j = 0; j < nPoints; ++j) {
                    if(pupPtr[i][j]>0) {
                        modePtr[i][j] *= sf * pupPtr[i][j];
                        squared_sum +=  modePtr[i][j]*modePtr[i][j];
                    }
                }
            }
        }
        else if(modeNumber % 2) {
            double sf = sqrt((double) 2 * (n + 1));
            for(int i = 0; i < nPoints; ++i) {
                for(int j = 0; j < nPoints; ++j) {
                    if(pupPtr[i][j]>0) {
                        modePtr[i][j] *= sf * pupPtr[i][j] * sin(m * (aPtr[i][j] + angle));
                        squared_sum +=  modePtr[i][j]*modePtr[i][j];
                    }
                }
            }
        }
        else {
            double sf = sqrt((double) 2 * (n + 1));
            for(int i = 0; i < nPoints; ++i) {
                for(int j = 0; j < nPoints; ++j) {
                    if(pupPtr[i][j]>0) {
                        modePtr[i][j] *= sf * pupPtr[i][j] * cos(m * (aPtr[i][j] + angle));
                        squared_sum +=  modePtr[i][j]*modePtr[i][j];
                    }
                }
            }
        }
            
        // normalize
        /*squared_sum = 1.0 / (sqrt(squared_sum / pupil.second) * lambda);
        for(int i = 0; i < nPoints; ++i) {
            for(int j = 0; j < nPoints; ++j) {
                modePtr[i][j] *= squared_sum;
            }
        }*/
        
        
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
                    double f = (xh > yh) ? 2*dy * (min(xh, max(xl, x2)) - xl) / (4 * dx * dy) :
                               2*dx* (min(yh, max(yl, y3)) - yl) / (4 * dx * dy);
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
    delPointers(pupPtr);


}


PupilMode::PupilMode(int firstMode, int lastMode, int klModeNumber, int nPoints, double r_c, double lambda, double angle, double cutoff) : Array<double> (nPoints, nPoints) {

    if(firstMode > lastMode) swap(firstMode, lastMode);

    if(klModeNumber < firstMode || klModeNumber > lastMode) {
        throw invalid_argument("klModeNumber (" + to_string(klModeNumber) +
                               ") is not in the range [ firstMode (" + to_string(firstMode) +
                               "), lastMode (" + to_string(lastMode) + ")]");
    }

    zero();
        
    static ModeCache& cache = ModeCache::getCache();
    const std::map<int, PupilMode::KL_cfg>& kle = cache.karhunenLoeveExpansion(firstMode, lastMode);

    double c;
    for(auto & it : kle.at(klModeNumber).zernikeWeights) {
        if(fabs(c = it.second) >= cutoff) {
            int zernikeModeIndex = it.first;
            auto ptr = cache.mode(zernikeModeIndex, nPoints, r_c, lambda, angle);
            this->add(*ptr, c);
        }
    }

}

