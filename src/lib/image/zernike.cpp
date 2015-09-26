#include "redux/image/zernike.hpp"

#include "redux/math/linalg.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/cache.hpp"

#include <mutex>

using namespace redux::math;
using namespace redux::image;
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

    double calcZernikeCovariance(int i, int j) {

        if((i < 2) || (j < 2)) return 0.0;
        int m, n, o, p;
        noll_to_mn(i, m, n);
        noll_to_mn(j, o, p);
        if(m != o) return 0.0;
        if(m) if((i + j) % 2) return 0.0;

        //  ; Now deal with the numerical terms: Dai
        int g1_sgn, g2_sgn, g3_sgn, g4_sgn;
        double k = pow(4.8 * exp(lgamma_r(6.0 / 5.0, &g1_sgn)), 5.0 / 6.0) * exp(lgamma_r(14.0 / 3.0, &g2_sgn) + 2.0 * lgamma_r(11.0 / 6.0, &g3_sgn)) / (pow(2.0, (8.0 / 3.0)) * M_PI);
        k *= pow(-1.0, (double)((n + p - 2 * m) / 2)) * sqrt((double)((n + 1) * (p + 1)));

        double g1 = lgamma_r(((double)(n + p) -  5.0 / 3.0) / 2.0, &g1_sgn);
        double g2 = lgamma_r(((double)(n - p) + 17.0 / 3.0) / 2.0, &g2_sgn);
        double g3 = lgamma_r(((double)(p - n) + 17.0 / 3.0) / 2.0, &g3_sgn);
        double g4 = lgamma_r(((double)(n + p) + 23.0 / 3.0) / 2.0, &g4_sgn);

        return k * exp(g1 - g2 - g3 - g4) * g1_sgn * g2_sgn * g3_sgn * g4_sgn;

    }

    mutex globalMutex;
}


bool Zernike::PairID::operator<(const PairID& rhs) const {
    if( first == rhs.first ) {
        return (second < rhs.second);
    }
    return (first < rhs.first);
}


double Zernike::covariance( int32_t i, int32_t j ) {
    
    if(i > j) std::swap(i, j);      // it's symmetric, so only store 1
    
    double& cov = Cache::get<PairID,double>( PairID(i,j), std::numeric_limits<double>::infinity() );
    unique_lock<mutex> lock(get().mtx);
    if( ! isfinite(cov) ) { // not calulated yet
        //cout << "Calulating ZernikeCov for (" << i << "," << j << ")" << endl;;
        cov = calcZernikeCovariance(i, j);
    } //else cout << "." << flush;
    return cov;

}


const vector<double>& Zernike::radialPolynomial(int32_t m, int32_t n) {

    vector<double>& poly = Cache::get<PairID,vector<double>>( PairID(m,n), vector<double>() );
    Zernike& z = get();
    unique_lock<mutex> lock(z.mtx);
    if( poly.empty() ) { // not calulated yet
        cout << "Calulating ZernikePolyfor (" << m << "," << n << ")" << endl;;
        int32_t nmm = (n - m) / 2;
        int32_t npm = (n + m) / 2;
        size_t nmax = max(npm, n);
        size_t fsz = z.factorials.size();
        if(fsz <= nmax) {
            z.factorials.resize(nmax + 1, 1);
            for(size_t i = max(fsz, 1UL); i <= nmax; ++i) z.factorials[i] = i * z.factorials[i - 1];
        }
        poly.resize(nmm + 1);
        for(int32_t s = 0, pm = -1; s <= nmm; ++s) {
            poly[s] = (double)((pm *= -1) * z.factorials[n - s]) / (z.factorials[s] * z.factorials[npm - s] * z.factorials[nmm - s]);
        }
        std::reverse(poly.begin(), poly.end());       // reverse so that the first element corresponds to the lowest-order term (r^m)
    }
    return poly;

}


const std::map<uint16_t, Zernike::KLPtr>& Zernike::karhunenLoeveExpansion(uint16_t first_mode, uint16_t last_mode) {

    std::map<uint16_t, KLPtr>& kle = Cache::get<PairID,std::map<uint16_t, KLPtr>>( PairID(first_mode,last_mode), std::map<uint16_t, KLPtr>() );
    
    unique_lock<mutex> lock(globalMutex);
    
    if( kle.empty() ) { // not calulated yet
        
        cout << "Calculating KL-expansion (" << first_mode << "," << last_mode << ")" << endl;;

        map<uint16_t, uint16_t> mapping, reverse_mapping;
        map<uint16_t, double> values;
        for(uint16_t i = first_mode; i <= last_mode; ++i) mapping[i] = i;
        for(uint16_t i = first_mode; i <= last_mode - 1; ++i) {
            double cov1 = covariance(mapping[i], mapping[last_mode]);
            for(int j = last_mode - 1; j > i; --j) {
                double cov2 = covariance(mapping[i], mapping[j]);
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
            if((values[it->first] = covariance(previous->second, it->second)) == 0) ++nBlocks;
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
                blockMatrix[b][i - first_in_block[b]][i - first_in_block[b]] = covariance(mapping.at(i), mapping.at(i));
            }
            for(int i = first_in_block[b] + 1; i <= last_in_block[b]; ++i) {    // subdiagonal elements (i_i±1), already calculated and stored in values.
                blockMatrix[b][i - 1 - first_in_block[b]][i - first_in_block[b]] = blockMatrix[b][i - first_in_block[b]][i - 1 - first_in_block[b]] = values[i];
            }
            for(int i = first_in_block[b]; i <= last_in_block[b] - 2; ++i) {    // calulate the rest
                for(int j = i + 2; j <= last_in_block[b]; ++j) {
                    blockMatrix[b][i - first_in_block[b]][j - first_in_block[b]] = blockMatrix[b][j - first_in_block[b]][i - first_in_block[b]] = covariance(mapping.at(i), mapping.at(j));
                }
            }
        }

        double *singular_values = new double [last_mode - first_mode + 1];
        size_t offset = 0;
        for(int b = 0; b < nBlocks; ++b) {
            int blockSize = last_in_block[b] - first_in_block[b] + 1;
            if(blockSize > 1) {
                double **v = newArray<double>(blockSize, blockSize);
                svd(*blockMatrix[b], blockSize, blockSize, singular_values + offset, *v);
                delArray(v);
            }
            else {
                singular_values[offset] = blockMatrix[b][0][0];
                blockMatrix[b][0][0] = 1.0;
            }
            offset += blockSize;
        }

        for(uint16_t i = first_mode; i <= last_mode; ++i) {
            KLPtr &cfg = kle[i];      // will insert a new element and return a reference if it doesn't exist
            if(!cfg) cfg.reset(new KL);
            int im = reverse_mapping.at(i), s;
            for(s = 0; (first_in_block[s] > im) || (last_in_block[s] < im); ++s);
            int n = last_in_block[s] - first_in_block[s] + 1;
            cfg->covariance = singular_values[im - first_mode];
            int sign = 1;   // Overall sign of the expansion is arbitrary. Depending on the SVD solver used, the sign varies.
                            // We choose the sign so that the KL modes as close as possible to the Zernikes (i.e. the principal component has weight closer to 1 than -1)
            double largestCoeff = 0;
            for(int m = 0; m < n; ++m) {
                double c = blockMatrix[s][m][im - first_in_block[s]];
                if ( fabs(c) > largestCoeff ) {
                    largestCoeff = fabs(c);
                    sign = (c<0)?-1:1;
                }
            }
            for(int m = 0; m < n; ++m) {
                int j = m + first_in_block[s];
                cfg->zernikeWeights.push_back(make_pair(mapping.at(j), sign * blockMatrix[s][m][im - first_in_block[s]]));
            }
        }

        delete[] first_in_block;
        delete[] last_in_block;
        delete[] singular_values;
        for(int i = 0; i < nBlocks; ++i) {
            delArray(blockMatrix[i]);
        }
        delete[] blockMatrix;
        
    }

    return kle;

}



void Zernike::clear(void) {
    Zernike& z = get();
    unique_lock<mutex> lock(z.mtx);
    Cache::clear<PairID,vector<double>>();
    Cache::clear<PairID,double>();
    z.factorials.clear();
}