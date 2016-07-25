#ifndef REDUX_MATH_HELPERS_HPP
#define REDUX_MATH_HELPERS_HPP

#include <functional>
#include <limits>
#include <utility>      // std::swap

#include <boost/numeric/conversion/cast.hpp>

/**
* @file
* Some general purpose functions.
*
* @author Tomas Hillberg
*/


namespace redux {

    namespace math {


        /*!
         * Almost eq. function that can be used to compare floats.
         * Ref: http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
         * maxUlps = max no. floting point possitions (tic's) away,
         *   must be in range (maxUlps > 0 && maxUlps < 4 * 1024 * 1024)
         */
        bool almostEqual( float A, float B, int maxUlps = 1 );
        
        
        template <typename T> T min(const T* data, size_t n=2) { if (n<1) return T(0); T ret = *data; while(--n) ret = std::min(*++data,ret); return ret; }
        template <typename T> T max(const T* data, size_t n=2) { if (n<1) return T(0); T ret = *data; while(--n) ret = std::max(*++data,ret); return ret; }
        template <typename T> double mean(const T* data, size_t n=2) {
            if (n<1) return NAN;
            double ret(static_cast<double>( *data ));
            for( unsigned int i=1; i<n; ++i ) ret += static_cast<double>( *++data );
            return ret/n;
        }
        template <typename T> void minMaxMean(const T* data, size_t n, T& min, T& max, double& mean) {
            if ( n > 0 ) {
                min = *data;
                max = *data;
                mean = static_cast<double>( *data );
                for( unsigned int i=1; i<n; ++i ) {
                    data++;
                    min = std::min(*data,min);
                    max = std::max(*data,max);
                    mean += static_cast<double>( *data );
                }
                mean /= n;
            }
        }
        
        template <typename T>
        void rmsStddev(const T* data, size_t n, double mean, double& rms, double& stddev) {
            if ( n > 1 ) {
                rms = stddev = 0;
                for( unsigned int i=0; i<n; ++i ) {
                    rms += static_cast<double>( *data * *data );
                    double tmp = ( static_cast<double>( *data ) - mean );
                    stddev += tmp * tmp;
                    data++;
                }
                rms = sqrt( rms / n );
                stddev = sqrt( stddev / (n-1) );
            } else if (n == 1) {
                rms = *data;
                stddev = 0;
            }
        }

        template <typename T>
        double rmsStddev(const T* data, size_t n, double& rms, double& stddev) {
            if ( n > 1 ) {
                rms = stddev = 0;
                double mn = mean(data,n);
                for( unsigned int i=0; i<n; ++i ) {
                    double tmp = static_cast<double>( *data );
                    rms += tmp * tmp;
                    tmp -= mn;
                    stddev += tmp * tmp;
                    data++;
                }
                rms = sqrt( rms / n );
                stddev = sqrt( stddev / (n-1) );
                return mn;
            } else if (n == 1) {
                rms = *data;
                stddev = 0;
                return rms;
            }
            return NAN;
        }

        

        /*!
         *  @brief      Cast with min/max boundaries.
         *  @param      src Input
         *  @param      min min
         *  @param      max max
         */
        template<typename Target, typename Source>
        Target bound_cast(Source src,
                          Target min=std::numeric_limits<Target>::lowest(),
                          Target max=std::numeric_limits<Target>::max()) {
            
            if ( min > max ) {
                std::swap( min, max );
            }

            try {
                Target tmp = boost::numeric_cast<Target>(src);
                return std::min(std::max(tmp,min),max);
            }
            catch (const boost::numeric::negative_overflow &) {
                return min;
            }
            catch (const boost::numeric::positive_overflow &) {
                return max;
            }
        }
        
        
        /*!
         *  @brief      Bracket 1D function minimum. Given a & b, find c such that f(a) > f(b) < f(c) (a & b might get modified)
         */
        template<typename T, typename U>
        void bracket(std::function<U(T)> f, T& a, T& b, T& c, U& fa, U& fb, U& fc, int limit=100);
        
        /*!
         *  @brief      Find minimum of function f on the interval [a, b], return x and f(x)
         */
        template<typename T, typename U>
        void brent(std::function<U(T)> f, T a, T b, T& x, U& fx, double tol=1E-4, int limit=100);
        
        
    }
}

#endif // REDUX_MATH_HELPERS_HPP
