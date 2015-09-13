#include "redux/momfbd/modes.hpp"

#include "redux/constants.hpp"
#include "redux/image/utils.hpp"
#include "redux/momfbd/cache.hpp"
//#include "redux/momfbd/legacy.hpp"
#include "redux/util/arrayutil.hpp"

#include <cmath>
#include <iostream>

using namespace redux::image;
using namespace redux::momfbd;
//using namespace redux::momfbd::legacy;
using namespace redux::util;
using namespace std;

void calculate(void) {
    
}

PupilMode::PupilMode(uint16_t modeNumber, uint16_t nPoints, double r_c, double angle) :
    Array<double> (nPoints, nPoints), atm_rms(0)  {      // Zernike

    if(modeNumber == 1) {
        Array<double>::operator=(1.0);
    } else {
        double** modePtr = makePointers(get(), nPoints, nPoints);
        //makeZernike_thi(modePtr,modeNumber,nPoints,r_c,angle);   //FIXME: using MvN's Zernike-generator for comparisons
        makeZernike_mvn(modePtr,modeNumber,nPoints,r_c,angle);
        delPointers(modePtr);
    }

    atm_rms = sqrt(Cache::getCache().zernikeCovariance(modeNumber,modeNumber));

}

namespace {
    const int signs[] = {1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1};
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
 //     cout << zernikeModeIndex << ":" << c << " " << flush;
            const PupilMode::Ptr& mode = cache.mode(zernikeModeIndex, nPoints, r_c, angle);
            this->add(*mode, c);
        }
    }

    if(klModeNumber<22) *this *= signs[klModeNumber];   // FIXME: making signs the same as MvN for testing.
    
    atm_rms = sqrt(kle->covariance);

}

