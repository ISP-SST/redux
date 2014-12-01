#include "redux/momfbd/modecache.hpp"

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


ModeCache::index::index(int modeNumber, int nPoints, double r_c, double wavelength, double angle)
    : firstMode(0), lastMode(0),
      modeNumber(modeNumber), nPoints(nPoints),
      r_c(r_c), wavelength(wavelength), angle(angle) { }


ModeCache::index::index(int firstMode, int lastMode, int modeNumber, int nPoints, double r_c, double wavelength, double angle)
    : firstMode(firstMode), lastMode(lastMode),
      modeNumber(modeNumber), nPoints(nPoints),
      r_c(r_c), wavelength(wavelength), angle(angle) { }

      
bool ModeCache::index::operator<(const index& rhs) const {
    if(firstMode == rhs.firstMode) {
        if(lastMode == rhs.lastMode) {
            if(modeNumber == rhs.modeNumber) {
                if(nPoints == rhs.nPoints) {
                    if(r_c == rhs.r_c) {
                        if(wavelength == rhs.wavelength) {
                            return (angle < rhs.angle);
                        }
                        else return (wavelength < rhs.wavelength);
                    }
                    else return (r_c < rhs.r_c);
                }
                else return (nPoints < rhs.nPoints);
            }
            else return (modeNumber < rhs.modeNumber);
        }
        else return (lastMode < rhs.lastMode);
    }
    else return (firstMode < rhs.firstMode);
}


ModeCache& ModeCache::getCache(void) {
    static ModeCache cache;
    return cache;
}


void ModeCache::clear(void) {

    factorials.clear();
    grids.clear();
    zernikeRadialPolynomials.clear();

}


const redux::image::Grid& ModeCache::grid(int sz) {

    unique_lock<mutex> lock(mtx);
    auto it = grids.find(sz);
    if(it != grids.end()) return it->second;
    return grids.emplace(sz, redux::image::Grid(sz)).first->second;

}


double ModeCache::zernikeCovariance(int i, int j) {

    if(i > j) std::swap(i, j);      // it's symmetric, so only store 1

    unique_lock<mutex> lock(mtx);
    std::pair<int, int> index = make_pair(i, j);
    auto it = zernikeCovariances.find(index);
    if(it != zernikeCovariances.end()) return it->second;

    double cov = calcZernikeCovariance(i, j);
    zernikeCovariances.emplace(index, cov);

    return cov;

}


const vector<double>& ModeCache::zernikeRadialPolynomial(int m, int n) {

    unique_lock<mutex> lock(mtx);
    auto it = zernikeRadialPolynomials.find(make_pair(m, n));
    if(it != zernikeRadialPolynomials.end()) return it->second;

    int nmm = (n - m) / 2;
    int npm = (n + m) / 2;
    size_t nmax = max(npm, n);
    size_t fsz = factorials.size();
    if(fsz <= nmax) {
        factorials.resize(nmax + 1, 1);
        for(size_t i = max(fsz, 1UL); i <= nmax; ++i) factorials[i] = i * factorials[i - 1];
    }
    vector<double> res(nmm + 1);
    for(int s = 0, pm = -1; s <= nmm; ++s) {
        res[s] = (double)((pm *= -1) * factorials[n - s]) / (factorials[s] * factorials[npm - s] * factorials[nmm - s]);
    }
    std::reverse(res.begin(), res.end());
    return zernikeRadialPolynomials.emplace(make_pair(m, n), res).first->second;

}


const std::map<int, Modes::KL_cfg>& ModeCache::karhunenLoeveExpansion(int first_mode, int last_mode) {


    {
        unique_lock<mutex> lock(mtx);
        auto it = karhunenLoeveExpansions.find(make_pair(first_mode, last_mode));
        if(it != karhunenLoeveExpansions.end()) return it->second;
    }


    map<int, int> mapping, reverse_mapping;
    map<int, double> values;
    for(int i = first_mode; i <= last_mode; ++i) mapping[i] = i;
    for(int i = first_mode; i <= last_mode - 1; ++i) {
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

    std::map<int, Modes::KL_cfg> res;
    for(int i = first_mode; i <= last_mode; ++i) {
        Modes::KL_cfg &cfg = res[i];      // will insert a new element and return a reference if it doesn't exist
        int im = reverse_mapping.at(i), s;
        for(s = 0; (first_in_block[s] > im) || (last_in_block[s] < im); ++s);
        int n = last_in_block[s] - first_in_block[s] + 1;
        cfg.covariance = values[im - first_mode];
        for(int m = 0; m < n; ++m) {
            int j = m + first_in_block[s];
            cfg.zernikeWeights.push_back(make_pair(mapping.at(j), blockMatrix[s][m][im - first_in_block[s]]));
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


#include "redux/file/fileana.hpp"
const ModeCache::ModePtr& ModeCache::mode(int modeNumber, int nPoints, double r_c, double wavelength, double angle) {

    index idx(modeNumber, nPoints, r_c, wavelength, angle);
    {
        unique_lock<mutex> lock(mtx);
        auto it = modes.find(idx);
        if(it != modes.end()) return it->second;
    }

    ModePtr ptr(new Modes::PupilMode(modeNumber, nPoints, r_c, wavelength, angle));
    unique_lock<mutex> lock(mtx);
    return modes.emplace(idx, ptr).first->second;

}


const ModeCache::ModePtr& ModeCache::mode(int firstMode, int lastMode, int modeNumber, int nPoints, double r_c, double wavelength, double angle) {

    index idx(firstMode, lastMode, modeNumber, nPoints, r_c, wavelength, angle);
    {
        unique_lock<mutex> lock(mtx);
        auto it = modes.find(idx);
        if(it != modes.end()) return it->second;
    }

    ModePtr ptr(new Modes::PupilMode(firstMode, lastMode, modeNumber, nPoints, r_c, wavelength, angle));

    unique_lock<mutex> lock(mtx);
    return modes.emplace(idx, ptr).first->second;


}

