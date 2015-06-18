#include "redux/momfbd/modes.hpp"

#include "redux/constants.hpp"
#include "redux/momfbd/cache.hpp"
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

PupilMode::PupilMode(uint16_t modeNumber, uint16_t nPoints, double r_c, double angle) :
    Array<double> (nPoints, nPoints), atm_rms(0)  {      // Zernike


    static Cache& cache = Cache::getCache();
    const std::pair<Array<double>,double>& pupil = cache.pupil(nPoints,r_c);

    zero();

    if(modeNumber == 1) {
        //add(pupil.first,1.0/lambda);
        memcpy(get(),pupil.first.get(),nPoints*nPoints*sizeof(double));
    } else {

        int m, n;
        noll_to_mn(modeNumber, m, n);

        const vector<double>& coeff = cache.zernikeRadialPolynomial(m, n);
        const redux::image::Grid& grid = cache.grid(nPoints,redux::PointF(nPoints/2,nPoints/2));
        float** distPtr = grid.distance.get();  // distance from pixels to centre (pupil & modes are centered on pixel (mid,mid))
        float** aPtr = grid.angle.get();

        double** modePtr = makePointers(ptr(), nPoints, nPoints);
        const double** pupPtr = makePointers(pupil.first.ptr(), nPoints, nPoints);

        shared_ptr<double*> r = sharedArray<double> (nPoints, nPoints);     // normalized distance ^{some order}
        shared_ptr<double*> r2 = sharedArray<double> (nPoints, nPoints);    // normalized distance squared
        double** rPtr = r.get();
        double** r2Ptr = r2.get();
        double normalization(0);

        for(int y = 0; y < nPoints; ++y) {
            for(int x = 0; x < nPoints; ++x) {
                //if(pupPtr[y][x]>0) {
                    double tmp = distPtr[y][x] / r_c;                       // normalize radial distance
                    //if(tmp>1) continue;
                    r2Ptr[y][x] = tmp * tmp;
                    if(m == 0) rPtr[y][x] = 1;
                    else if(m == 1) rPtr[y][x] = tmp;
                    else rPtr[y][x] = pow(tmp, m);                          // lowest order term ~ r^{m}
                    modePtr[y][x] = rPtr[y][x] * coeff[0];                  // add lowest order to mode
                //}
            }
        }

        // generate polynomial part
        for(auto it = coeff.begin() + 1; it < coeff.end(); it++) {
            for(int y = 0; y < nPoints; ++y) {
                for(int x = 0; x < nPoints; ++x) {
                    //if(pupPtr[y][x]>0) {
                        rPtr[y][x] *= r2Ptr[y][x];                          //  next term ~ r^{m+2}
                        modePtr[y][x] += rPtr[y][x] * *it;
                    //}
                }
            }
        }

        // Angular component
        if(m == 0) {
            double sf = sqrt(n + 1);
            for(int y = 0; y < nPoints; ++y) {
                for(int x = 0; x < nPoints; ++x) {
                    modePtr[y][x] *= sf; // * pupPtr[y][x];
                    if(pupPtr[y][x]>0) {
                        normalization +=  modePtr[y][x]*modePtr[y][x] * pupPtr[y][x];
                    }
                }
            }
        }
        else if(modeNumber % 2) {
            double sf = sqrt((double) 2 * (n + 1));
            for(int y = 0; y < nPoints; ++y) {
                for(int x = 0; x < nPoints; ++x) {
                    modePtr[y][x] *= sf * sin(m * (aPtr[y][x] + angle)); //* pupPtr[y][x] 
                    if(pupPtr[y][x]>0) {
                        normalization +=  modePtr[y][x]*modePtr[y][x] * pupPtr[y][x];
                    }
                }
            }
        }
        else {
            double sf = sqrt((double) 2 * (n + 1));
            for(int y = 0; y < nPoints; ++y) {
                for(int x = 0; x < nPoints; ++x) {
                    modePtr[y][x] *= sf * cos(m * (aPtr[y][x] + angle)); //* pupPtr[y][x] 
                    if(pupPtr[y][x]>0) {
                        normalization +=  modePtr[y][x]*modePtr[y][x] * pupPtr[y][x];
                    }
                }
            }
        }
            
        // normalize
        normalization = 1.0 / (sqrt(normalization / pupil.second));
        for(int y = 0; y < nPoints; ++y) {
            for(int x = 0; x < nPoints; ++x) {
                //if(pupPtr[y][x]>0) {
                    modePtr[y][x] *= normalization;
                //}
            }
        }
        
        delPointers(modePtr);
        delPointers(pupPtr);
        
    }

    atm_rms = sqrt(cache.zernikeCovariance(modeNumber,modeNumber));

}


PupilMode::PupilMode(uint16_t firstMode, uint16_t lastMode, uint16_t klModeNumber, uint16_t nPoints, double r_c, double angle, double cutoff) :
     Array<double> (nPoints, nPoints), atm_rms(0) {

    if(firstMode > lastMode) swap(firstMode, lastMode);

    if(klModeNumber < firstMode || klModeNumber > lastMode) {
        throw invalid_argument("klModeNumber (" + to_string(klModeNumber) +
                               ") is not in the range [ firstMode (" + to_string(firstMode) +
                               "), lastMode (" + to_string(lastMode) + ")]");
    }

    zero();
        
    static Cache& cache = Cache::getCache();
    const PupilMode::KLPtr& kle = cache.karhunenLoeveExpansion(firstMode, lastMode).at(klModeNumber);
    double c;
    for(auto & it : kle->zernikeWeights) {
        if(fabs(c = it.second) >= cutoff) {
            uint16_t zernikeModeIndex = it.first;
            const PupilMode::Ptr& mode = cache.mode(zernikeModeIndex, nPoints, r_c, angle);
            this->add(*mode, c);
        }
    }
    
    atm_rms = sqrt(kle->covariance);

}

