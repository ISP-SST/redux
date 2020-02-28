#include "redux/image/zernike.hpp"

#include "redux/math/linalg.hpp"
#include "redux/image/grid.hpp"
#include "redux/image/pupil.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/cache.hpp"

#include <mutex>
#include <numeric>
#include <boost/multiprecision/cpp_int.hpp>

using namespace redux::math;
using namespace redux::image;
using namespace redux::util;
using namespace std;


namespace {

    mutex globalMutex;
}


bool Zernike::PairID::operator<(const PairID& rhs) const {
    if( first == rhs.first ) {
        return (second < rhs.second);
    }
    return (first < rhs.first);
}


bool Zernike::RadialID::operator<(const RadialID& rhs) const {
    if( nP == rhs.nP ) {
        if( m == rhs.m ) {
            if( n == rhs.n ) {
                return (r < rhs.r);
            }
            return (n < rhs.n);
        }
        return (m < rhs.m);
    }
    return (nP < rhs.nP);
}


bool Zernike::AngularID::operator<(const AngularID& rhs) const {
    if( nP == rhs.nP ) {
        if( m == rhs.m ) {
            return (angle < rhs.angle);
        }
        return (m < rhs.m);
    }
    return (nP < rhs.nP);
}


bool Zernike::PolyID::operator<(const PolyID& rhs) const {
    if( n == rhs.n ) {
        if( m == rhs.m ) {
            return (flags < rhs.flags);
        }
        return (m < rhs.m);
    }
    return (n < rhs.n);
}


void Zernike::NollToNM( unsigned int j, uint16_t& n, int16_t& m ) {
    
    n = 0;
    unsigned int len = 1;
    for(int i = 1; len < j; ++i) {
        len += (n = i) + 1;
    }
    int dl = n + 1 - len + j;
    m = 2 * ((dl + (n % 2)) / 2) + !(n % 2) - 1;
    if( j%2 ) m = -m;
    
}


unsigned int Zernike::NMToNoll( const uint16_t ZernikeN, const int16_t ZernikeM ) {
    
    // TODO: proper calculation instead of retarded trial and error loop!
    // note: the line below does not return corret results:
    //    return (ZernikeN*(ZernikeN+2)+ZernikeM)/2;
    
    if( (ZernikeN-abs(ZernikeM))%2 ) return 0;
    
    unsigned int j(0);
    uint16_t n(0);
    int16_t m(0);
    while( m != ZernikeM || n != ZernikeN ) {
        NollToNM( ++j, n, m );
        if( j>1E6 ) return 0;       // bailout
    }
    return j;
}


double Zernike::calcCovariance( uint16_t n1, int16_t m1, uint16_t n2, int16_t m2 ) {
    
    if( m1 != m2 ) return 0.0;
    
    //  ; Now deal with the numerical terms: Dai
    int g1_sgn, g2_sgn, g3_sgn, g4_sgn;
    double k = pow(4.8 * exp(lgamma_r(6.0 / 5.0, &g1_sgn)), 5.0 / 6.0) * exp(lgamma_r(14.0 / 3.0, &g2_sgn) + 2.0 * lgamma_r(11.0 / 6.0, &g3_sgn)) / (pow(2.0, (8.0 / 3.0)) * M_PI);
    k *= pow(-1.0, (double)(( n1 + n2 - 2 * m1 ) / 2.0)) * sqrt((double)(( n1 + 1) * ( n2 + 1)));

    double g1 = lgamma_r(((double)( n1 + n2 ) -  5.0 / 3.0) / 2.0, &g1_sgn);
    double g2 = lgamma_r(((double)( n1 - n2 ) + 17.0 / 3.0) / 2.0, &g2_sgn);
    double g3 = lgamma_r(((double)( n2 - n1 ) + 17.0 / 3.0) / 2.0, &g3_sgn);
    double g4 = lgamma_r(((double)( n1 + n2 ) + 23.0 / 3.0) / 2.0, &g4_sgn);

    return k * exp(g1 - g2 - g3 - g4) * g1_sgn * g2_sgn * g3_sgn * g4_sgn;
    
}


double Zernike::calcCovariance( unsigned int i, unsigned int j ) {
    
    if((i < 2) || (j < 2)) return 0.0;
    
    uint16_t in, jn;
    int16_t im, jm;
    
    NollToNM( i, in, im );
    if( im ) if((i + j) % 2) return 0.0;
    
    NollToNM( j, jn, jm );
    
    return calcCovariance( in, abs(im), jn, abs(jm) );
}


double Zernike::getCovariance ( unsigned int i, unsigned int j ) {
    
    if( i > j ) std::swap(i, j);      // it's symmetric, so only store 1
    
    double& cov = Cache::get<PairID,double>( PairID(i,j), std::numeric_limits<double>::infinity() );
    unique_lock<mutex> lock(get().mtx);
    if( ! isfinite(cov) ) { // not calulated yet
        cov = calcCovariance( i, j );
    }
    return cov;

}
 

shared_ptr<double> Zernike::getRadial( unsigned int nPixels, float radius, uint16_t n, uint16_t m, int flags ) {
    
    RadialID rid( nPixels, radius, n, m );
    shared_ptr<double>& rpoly = Cache::get<RadialID,shared_ptr<double>>( rid );
    
    if( !rpoly || (flags&FORCE) ) {
        
        if( flags & VERBOSE ) {
            if ( rpoly ) cout << "Re-";
            cout << "Generating Zernike radial part (" << nPixels << "x" << nPixels << " pixels, r=" << radius << " n=" << n << " |m|=" << m << endl;
            if( flags & OLD_METHOD) cout << "   (Using OLD method for evaluating radial polynomial!)" << endl;
        }

        size_t blockSize = nPixels*nPixels;
        rpoly = rdx_get_shared<double>( blockSize );
        double* polyPtr = rpoly.get();
        float midPlusHalf = nPixels/2.0 + 0.5;
        const shared_ptr<Grid> grid = Grid::get( nPixels, midPlusHalf, midPlusHalf );

        double* distPtr = grid->distance.get();
        shared_ptr<long double> r = rdx_get_shared<long double>( blockSize );
        shared_ptr<long double> r2 = rdx_get_shared<long double>( blockSize );
        shared_ptr<long double> tmpPoly = rdx_get_shared<long double>( blockSize );
        long double* rPtr = r.get();
        long double* r2Ptr = r2.get();
        long double* tmpPtr = tmpPoly.get();
        long double r_inv = 1.0/radius;

        std::transform( distPtr, distPtr+blockSize,  rPtr,
                        [r_inv](const double& d){ return std::min<double>( 1.0, r_inv*d ); });
        std::transform( rPtr, rPtr+blockSize,  r2Ptr, [](const long double& d){ return d*d; });
        memset( tmpPtr, 0, blockSize*sizeof(long double) );
        
        RadialPolynomial& poly = Zernike::getRadialPolynomial( n, m, flags );
        if( flags & OLD_METHOD ) {    // default is the old implementation
            for( auto& pc: poly.poly ) {
                double this_exp = pc.first;
                double cc = static_cast<double>(pc.second);             // N.B. Simply upgrading cc to long double, and using powl impproves edge-issues.
                std::transform( rPtr, rPtr+blockSize, tmpPtr, tmpPtr,
                                [&](const double& rr, const long double& t){ return t+pow( rr, this_exp )*cc; });
            }

        } else {
            using namespace std::placeholders;
            std::transform( rPtr, rPtr+blockSize, tmpPtr, std::bind( &RadialPolynomial::eval, poly, _1 ));
        }

        std::copy_n( tmpPtr, blockSize, polyPtr );
    }

    return rpoly;
    
}


shared_ptr<double> Zernike::getAngular( unsigned int nPixels, float angle, int16_t m , int flags ) {
    
    AngularID aid( nPixels, angle, m );
    shared_ptr<double>& ang = Cache::get<AngularID,shared_ptr<double>>( aid );
    
    if( !ang ) {    // never force re-calculation for the angular part
        if( flags & VERBOSE ) {
            cout << "Generating Zernike angular part (" << nPixels << "x" << nPixels << " pixels, m=" << m << endl;
        }
        size_t blockSize = nPixels*nPixels;
        ang = rdx_get_shared<double>( blockSize );
        double* angPtr = ang.get();
        float midPlusHalf = nPixels/2.0 + 0.5;
        const shared_ptr<Grid> grid = Grid::get( nPixels, midPlusHalf, midPlusHalf );
        double* aPtr = grid->angle.get();
        if( m < 0 ) {
            m = abs(m);
            std::transform( aPtr, aPtr+blockSize, angPtr, [&](const double& a){ return sinl(m*(a-angle)); });
        } else {
            std::transform( aPtr, aPtr+blockSize, angPtr, [&](const double& a){ return cosl(m*(a-angle)); });
        }
    }
    
    return ang;
    
}


void Zernike::getZernike( double* modePtr, unsigned int nPixels, float radius, float angle, uint16_t ZernikeN, int16_t ZernikeM, int flags ) {
    
    if( !(flags&(GET_RADIAL|GET_ANGULAR)) ) {   // neither specified, default to both
        flags |= (GET_RADIAL|GET_ANGULAR);
    }
    
    size_t blockSize = nPixels*nPixels;
    shared_ptr<long double> tmpD = rdx_get_shared<long double>( blockSize );    // temporary storage
    long double* tmpPtr = tmpD.get();
    
    if( flags&GET_RADIAL ) {
        shared_ptr<double> radial = getRadial( nPixels, radius, ZernikeN, abs(ZernikeM), flags );
        std::copy_n( radial.get(), blockSize, tmpPtr );
    } else {
        std::fill_n( tmpPtr, blockSize, 1.0 );
    }
    
    if( flags&GET_ANGULAR ) {
        shared_ptr<double> angular = getAngular( nPixels, angle, ZernikeM, flags );
        std::transform( tmpPtr, tmpPtr+blockSize, angular.get(), tmpPtr, std::multiplies<long double>());
    }
    
    long double normalization(0.0);
    if( flags&NORMALIZE ) {     // numerical normalization
        redux::image::Pupil pupil = Pupil::fetch( nPixels, radius );
        double* pupPtr = pupil.ptr();
        for( size_t i(0); i<blockSize; ++i ) {
            normalization += tmpPtr[i]*tmpPtr[i]*pupPtr[i];
        }
        normalization = sqrtl( pupil.area/normalization );
    } else {                    // analytical normalization
        if( ZernikeM ) {
            normalization = sqrtl(2*(ZernikeN + 1));
        } else {
            normalization = sqrtl(ZernikeN + 1);
        }
    }
    
    std::transform( tmpPtr, tmpPtr+blockSize, modePtr, std::bind(std::multiplies<long double>(), std::placeholders::_1, normalization));

}


Zernike::RadialPolynomial& Zernike::getRadialPolynomial( uint16_t n, uint16_t m, int flags ) {
    
    Zernike::RadialPolynomial& poly = Cache::get<PolyID,RadialPolynomial>( PolyID(n,m,flags&0xFF), RadialPolynomial(n,m,flags) );

    //unique_lock<mutex> lock(get().mtx);
    
    if( poly.empty() || flags&FORCE ) poly.calc( flags );
    
    return poly;
    
}


void Zernike::RadialPolynomial::calc( int f ) {
    
    flags = f;
    
    uint16_t bb = (n - m) / 2;
    uint16_t aa = (n + m) / 2;
    useRatios = !(flags&OLD_METHOD);
    
    if( flags&VERBOSE ) {
        if ( !poly.empty() ) cout << "Re-";
        cout << "Generating coefficient" << (useRatios?"-ratios":"s") << " for the radial polynomial [(n,m) = (" << n << "," << m << ")]: ";
    }

    poly.clear();
    values.clear();
    
    int pm = -1;
    if( useRatios ) {
        poly[m] = n_choose_k( aa, m );
        if( bb%2 ) poly[m] *= -1;
        for( uint16_t kk=0; kk < bb; ++kk ) {
            mp_float aaa = (aa-kk) * (bb-kk);
            mp_float bbb = (kk+1)*(n-kk);
            poly[ n-2*kk ] = - bbb/aaa;
        }
    } else {
        for( uint16_t s(0); s <= bb; ++s) {
            //long double v = (pm *= -1) * (long double)n_choose_k(n-s,s) * n_choose_k(n-2*s,bb-s);
            mp_float v = n_choose_k(n-s,s);
            v *= n_choose_k(n-2*s,bb-s);
            v *= (pm *= -1);
            poly[ n-2*s ] = v;
        }
    }
    
    if( flags&VERBOSE ) {
        cout << printArray( poly, "P", 8 ) << endl;
    }

}




long double Zernike::RadialPolynomial::eval( long double r ) {

    if( values.find(r) == values.end() ) {
        const mp_float rr = r;
        mp_float sum(0);
        if( useRatios ) {
            const mp_float rr2 = rr*rr;
            sum = 1;                        // this scheme is multiplicative, so start with 1.
            for( auto& p: poly ) {
                if( p.first == m ) {
                    continue;
                }
                sum *= rr2*p.second;
                sum += 1;
            }
            if( m ) sum *= poly[m]*pow( rr, m );
        } else {
            for( auto& p: poly ) {
                mp_float tmp = pow( rr, p.first );
                tmp *= p.second;
                sum += tmp;
            }
        }
        long double ldsum = static_cast<long double>(sum);
        values[ r ] = ldsum;
        return ldsum;
    }
    
    return values[ r ];
    
}


const std::map<uint16_t, Zernike::KLPtr>& Zernike::karhunenLoeveExpansion(uint16_t first_mode, uint16_t last_mode) {

    std::map<uint16_t, KLPtr>& kle = Cache::get<PairID,std::map<uint16_t, KLPtr>>( PairID(first_mode,last_mode), std::map<uint16_t, KLPtr>() );
    
    unique_lock<mutex> lock(globalMutex);
    
    if( kle.empty() ) { // not calulated yet

        map<uint16_t, uint16_t> mapping, reverse_mapping;
        map<uint16_t, double> values;
        for(uint16_t i = first_mode; i <= last_mode; ++i) mapping[i] = i;
        for(uint16_t i = first_mode; i <= last_mode - 1; ++i) {
            double cov1 = getCovariance (mapping[i], mapping[last_mode]);
            for(int j = last_mode - 1; j > i; --j) {
                double cov2 = getCovariance (mapping[i], mapping[j]);
                if((!cov2) && (cov1)) {
                    swap(mapping[j], mapping[j + 1]);
                }
                else {
                    cov1 = cov2;
                }
            }
        }

        for(auto & element : mapping) {
            reverse_mapping[element.second] = element.first;
        }

        int nBlocks(1);
        auto previous = mapping.begin();
        for(auto it = ++mapping.begin(); it != mapping.end(); ++it) {
            if((values[it->first] = getCovariance (previous->second, it->second)) == 0) ++nBlocks;
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
                blockMatrix[b][i - first_in_block[b]][i - first_in_block[b]] = getCovariance (mapping.at(i), mapping.at(i));
            }
            for(int i = first_in_block[b] + 1; i <= last_in_block[b]; ++i) {    // subdiagonal elements (i_i±1), already calculated and stored in values.
                blockMatrix[b][i - 1 - first_in_block[b]][i - first_in_block[b]] = blockMatrix[b][i - first_in_block[b]][i - 1 - first_in_block[b]] = values[i];
            }
            for(int i = first_in_block[b]; i <= last_in_block[b] - 2; ++i) {    // calulate the rest
                for(int j = i + 2; j <= last_in_block[b]; ++j) {
                    blockMatrix[b][i - first_in_block[b]][j - first_in_block[b]] = blockMatrix[b][j - first_in_block[b]][i - first_in_block[b]] = getCovariance (mapping.at(i), mapping.at(j));
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
                            // We choose the sign so that the KL modes are as close as possible to the Zernikes (i.e. the principal component has weight closer to 1 than -1)
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
    unique_lock<mutex> lock(get().mtx);
    Cache::clear<PairID,vector<double>>();          // covariances
    Cache::clear<PairID,double>();
    Cache::clear<PolyID,RadialPolynomial>();        // radial polynomials
    Cache::clear<RadialID,shared_ptr<double>>();    // radial images (i.e. polynomials evaluated for each pixel)
    Cache::clear<AngularID,shared_ptr<double>>();   // angular images (i.e. polynomials evaluated for each pixel)
    z.factorials.clear();
}
