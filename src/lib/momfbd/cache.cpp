#include "redux/momfbd/cache.hpp"

#include "redux/constants.hpp"
#include "redux/math/linalg.hpp"
//#include "redux/momfbd/legacy.hpp"
#include "redux/util/arrayutil.hpp"

#include <cmath>

using namespace redux::momfbd;
//using namespace redux::momfbd::legacy;
using namespace redux::util;
using namespace redux;

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

    double calcZernikeCovariance(int i, int j) {

        if((i < 2) || (j < 2)) return 0.0;
        int m, n, o, p;
        noll_to_mn(i, m, n);
        noll_to_mn(j, o, p);
        if(m != o) return 0.0;
        if(m) if((i + j) % 2) return 0.0;

        //  ; Now deal with the numerical terms: Dai
        int g1_sgn, g2_sgn, g3_sgn, g4_sgn;
        double k = pow(4.8 * exp(lgamma_r(6.0 / 5.0, &g1_sgn)), 5.0 / 6.0) * exp(lgamma_r(14.0 / 3.0, &g1_sgn) + 2.0 * lgamma_r(11.0 / 6.0, &g1_sgn)) / (pow(2.0, (8.0 / 3.0)) * redux::PI);
        k *= pow(-1.0, (double)((n + p - 2 * m) / 2)) * sqrt((double)((n + 1) * (p + 1)));

        double g1 = lgamma_r(((double)(n + p) -  5.0 / 3.0) / 2.0, &g1_sgn);
        double g2 = lgamma_r(((double)(n - p) + 17.0 / 3.0) / 2.0, &g2_sgn);
        double g3 = lgamma_r(((double)(p - n) + 17.0 / 3.0) / 2.0, &g3_sgn);
        double g4 = lgamma_r(((double)(n + p) + 23.0 / 3.0) / 2.0, &g4_sgn);

        return k * exp(g1 - g2 - g3 - g4) * g1_sgn * g2_sgn * g3_sgn * g4_sgn;

    }

}


Cache::ModeID::ModeID(uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double wavelength, double angle)
    : firstMode(0), lastMode(0),
      modeNumber(modeNumber), nPoints(nPoints),
      pupilRadius(pupilRadius), wavelength(wavelength), angle(angle) { }


Cache::ModeID::ModeID(uint16_t firstMode, uint16_t lastMode, uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double wavelength, double angle)
    : firstMode(firstMode), lastMode(lastMode),
      modeNumber(modeNumber), nPoints(nPoints),
      pupilRadius(pupilRadius), wavelength(wavelength), angle(angle) { }

uint64_t Cache::ModeID::size( void ) const {
    static uint64_t sz = 4*sizeof(uint16_t) + 3*sizeof(double);
    return sz;
}


uint64_t Cache::ModeID::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack(ptr,firstMode);
    count += pack(ptr+count,lastMode);
    count += pack(ptr+count,modeNumber);
    count += pack(ptr+count,nPoints);
    count += pack(ptr+count,pupilRadius);
    count += pack(ptr+count,wavelength);
    count += pack(ptr+count,angle);
    return count;
}


uint64_t Cache::ModeID::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = unpack(ptr,firstMode,swap_endian);
    count += unpack(ptr+count,lastMode,swap_endian);
    count += unpack(ptr+count,modeNumber,swap_endian);
    count += unpack(ptr+count,nPoints,swap_endian);
    count += unpack(ptr+count,pupilRadius,swap_endian);
    count += unpack(ptr+count,wavelength,swap_endian);
    count += unpack(ptr+count,angle,swap_endian);
    return count;
}


bool Cache::ModeID::operator<(const ModeID& rhs) const {
    if(firstMode == rhs.firstMode) {
        if(lastMode == rhs.lastMode) {
            if(modeNumber == rhs.modeNumber) {
                if(nPoints == rhs.nPoints) {
                    if(pupilRadius == rhs.pupilRadius) {
                        if(wavelength == rhs.wavelength) {
                            return (angle < rhs.angle);
                        }
                        else return (wavelength < rhs.wavelength);
                    }
                    else return (pupilRadius < rhs.pupilRadius);
                }
                else return (nPoints < rhs.nPoints);
            }
            else return (modeNumber < rhs.modeNumber);
        }
        else return (lastMode < rhs.lastMode);
    }
    else return (firstMode < rhs.firstMode);
}


Cache& Cache::getCache(void) {
    static Cache cache;
    return cache;
}


void Cache::clear(void) {

    factorials.clear();
    grids.clear();
    zernikeRadialPolynomials.clear();

}


const redux::image::Grid& Cache::grid(uint32_t sz, PointF origin) {

    unique_lock<mutex> lock(mtx);
    auto it = grids.find(make_pair(sz,origin));
    if(it != grids.end()) return it->second;
    return grids.emplace(make_pair(sz,origin), redux::image::Grid(sz,origin)).first->second;

}


const std::pair<Array<double>, double>& Cache::pupil(uint32_t nPoints, float radius) {

    unique_lock<mutex> lock(mtx);
    auto it = pupils.find( make_pair(nPoints,radius) );
    if(it != pupils.end()) return it->second;
    Array<double> pup;
    double area = redux::image::makePupil(pup,nPoints,radius);
    return pupils.emplace(make_pair(nPoints,radius), make_pair(pup,area)).first->second;

}


double Cache::zernikeCovariance(int32_t i, int32_t j) {

    if(i > j) std::swap(i, j);      // it's symmetric, so only store 1

    unique_lock<mutex> lock(mtx);
    std::pair<int32_t, int32_t> index = make_pair(i, j);
    auto it = zernikeCovariances.find(index);
    if(it != zernikeCovariances.end()) return it->second;

    double cov = calcZernikeCovariance(i, j);
    zernikeCovariances.emplace(index, cov);

    return cov;

}


const vector<double>& Cache::zernikeRadialPolynomial(int32_t m, int32_t n) {

    unique_lock<mutex> lock(mtx);
    auto it = zernikeRadialPolynomials.find(make_pair(m, n));
    if(it != zernikeRadialPolynomials.end()) return it->second;

    int32_t nmm = (n - m) / 2;
    int32_t npm = (n + m) / 2;
    size_t nmax = max(npm, n);
    size_t fsz = factorials.size();
    if(fsz <= nmax) {
        factorials.resize(nmax + 1, 1);
        for(size_t i = max(fsz, 1UL); i <= nmax; ++i) factorials[i] = i * factorials[i - 1];
    }
    vector<double> res(nmm + 1);
    for(int32_t s = 0, pm = -1; s <= nmm; ++s) {
        res[s] = (double)((pm *= -1) * factorials[n - s]) / (factorials[s] * factorials[npm - s] * factorials[nmm - s]);
    }
    std::reverse(res.begin(), res.end());
    return zernikeRadialPolynomials.emplace(make_pair(m, n), res).first->second;

}


const std::map<uint16_t, PupilMode::KLPtr>& Cache::karhunenLoeveExpansion(uint16_t first_mode, uint16_t last_mode) {


    {
        unique_lock<mutex> lock(mtx);
        auto it = karhunenLoeveExpansions.find(make_pair(first_mode, last_mode));
        if(it != karhunenLoeveExpansions.end()) return it->second;
    }


    map<uint16_t, uint16_t> mapping, reverse_mapping;
    map<uint16_t, double> values;
    for(uint16_t i = first_mode; i <= last_mode; ++i) mapping[i] = i;
    for(uint16_t i = first_mode; i <= last_mode - 1; ++i) {
        double cov1 = zernikeCovariance(mapping[i], mapping[last_mode]);
        for(int j = last_mode - 1; j > i; --j) {
            double cov2 = zernikeCovariance(mapping[i], mapping[j]);
            if((!cov2) && (cov1)) {
                swap(mapping[j], mapping[j + 1]);
            }
            else {
                cov1 = cov2;
            }
        }
    }

    for(auto & it : mapping) {
        reverse_mapping[it.second] = it.first;
    }

    int nBlocks(1);
    auto previous = mapping.begin();
    for(auto it = ++mapping.begin(); it != mapping.end(); ++it) {
        if((values[it->first] = zernikeCovariance(previous->second, it->second)) == 0) ++nBlocks;
        previous = it;
    }

    int* first_in_block = new int [nBlocks];
    int* last_in_block = new int [nBlocks];
    double ***blockMatrix = new double** [nBlocks];
    for(int b = 0, n = first_mode; b < nBlocks; ++b) {
        first_in_block[b] = n++;
        while((n <= last_mode) && (values.at(n))) {
            ++n;
        }
        last_in_block[b] = n - 1;
        size_t blockSize = last_in_block[b] - first_in_block[b] + 1;
        blockMatrix[b] = newArray<double>(blockSize, blockSize);

        for(int i = first_in_block[b]; i <= last_in_block[b]; ++i) {        // diagonal elements (i_i)
            blockMatrix[b][i - first_in_block[b]][i - first_in_block[b]] = zernikeCovariance(mapping.at(i), mapping.at(i));
        }
        for(int i = first_in_block[b] + 1; i <= last_in_block[b]; ++i) {    // subdiagonal elements (i_i±1), already calculated and stored in values.
            blockMatrix[b][i - 1 - first_in_block[b]][i - first_in_block[b]] = blockMatrix[b][i - first_in_block[b]][i - 1 - first_in_block[b]] = values[i];
        }
        for(int i = first_in_block[b]; i <= last_in_block[b] - 2; ++i) {    // calulate the rest
            for(int j = i + 2; j <= last_in_block[b]; ++j) {
                blockMatrix[b][i - first_in_block[b]][j - first_in_block[b]] = blockMatrix[b][j - first_in_block[b]][i - first_in_block[b]] = zernikeCovariance(mapping.at(i), mapping.at(j));
            }
        }
    }

    double *singular_values = new double [last_mode - first_mode + 1];
    size_t offset = 0;
    for(int b = 0; b < nBlocks; ++b) {
        int blockSize = last_in_block[b] - first_in_block[b] + 1;
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

    std::map<uint16_t, PupilMode::KLPtr> res;
    for(uint16_t i = first_mode; i <= last_mode; ++i) {
        PupilMode::KLPtr &cfg = res[i];      // will insert a new element and return a reference if it doesn't exist
        if(!cfg) cfg.reset(new PupilMode::KL);
        int im = reverse_mapping.at(i), s;
        for(s = 0; (first_in_block[s] > im) || (last_in_block[s] < im); ++s);
        int n = last_in_block[s] - first_in_block[s] + 1;
        cfg->covariance = values[im - first_mode];
        for(int m = 0; m < n; ++m) {
            int j = m + first_in_block[s];
            // N.B. minus-sign for the coefficients below is arbitrary. Depending on the SVD solver used, the result might vary in sign.
            // This minus was added just to get the KL modes as close as possible to the Zernikes (i.e. the principal component has weight closer to 1 than -1)
            cfg->zernikeWeights.push_back(make_pair(mapping.at(j), - blockMatrix[s][m][im - first_in_block[s]]));
        }
    }

    delete[] first_in_block;
    delete[] last_in_block;
    delete[] singular_values;
    for(int i = 0; i < nBlocks; ++i) {
        delArray(blockMatrix[i]);
    }
    delete[] blockMatrix;

    unique_lock<mutex> lock(mtx);
    return karhunenLoeveExpansions.emplace(make_pair(first_mode, last_mode), res).first->second;

}


const PupilMode::Ptr& Cache::addMode(const ModeID& idx, PupilMode::Ptr& m) {

    unique_lock<mutex> lock(mtx);
    auto it = modes.find(idx);
    if(it != modes.end()) return it->second;
    return modes.emplace(idx, m).first->second;

}


const PupilMode::Ptr& Cache::addMode(uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double wavelength, double angle, PupilMode::Ptr& m) {

    ModeID idx(modeNumber, nPoints, pupilRadius, wavelength, angle);
    return addMode(idx, m);

}


const PupilMode::Ptr& Cache::mode(const ModeID& idx) {

    {
        unique_lock<mutex> lock(mtx);
        auto it = modes.find(idx);
        if(it != modes.end()) return it->second;
    }

    PupilMode::Ptr ptr;
    if( idx.firstMode == 0 || idx.lastMode == 0 ) ptr.reset(new PupilMode(idx.modeNumber, idx.nPoints, idx.pupilRadius, idx.wavelength, idx.angle));    // Zernike
    else ptr.reset(new PupilMode(idx.firstMode, idx.lastMode, idx.modeNumber, idx.nPoints, idx.pupilRadius, idx.wavelength, idx.angle));                // K-L
    
    unique_lock<mutex> lock(mtx);
    return modes.emplace(idx, ptr).first->second;

}


const PupilMode::Ptr& Cache::mode(uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double wavelength, double angle) {

    ModeID idx(modeNumber, nPoints, pupilRadius, wavelength, angle);
    return mode(idx);

}


const PupilMode::Ptr& Cache::mode(uint16_t firstMode, uint16_t lastMode, uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double wavelength, double angle) {

    ModeID idx(firstMode, lastMode, modeNumber, nPoints, pupilRadius, wavelength, angle);
    if( modeNumber == 2 || modeNumber == 3) {   // Always use Zernike modes for the tilts
        idx.firstMode = idx.lastMode == 0;
    }
    return mode(idx);

}

