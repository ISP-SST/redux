#ifndef REDUX_MATH_HELPERS_HPP
#define REDUX_MATH_HELPERS_HPP

#include <functional>
#include <limits>
#include <numeric>
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
            if( n<1 ) return NAN;
            double ret( static_cast<double>( *data ) );
            for( size_t i=1; i<n; ++i ) ret += static_cast<double>( *++data );
            //ret = std::accumulate( data+1, data+n, 0.0 );
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
        void rmsStddev( const T* data, size_t n, double mean, double& rms, double& stddev ) {
            if ( n > 1 ) {
                rms = stddev = 0;
                for( size_t i(0); i<n; ++i ) {
                    rms += static_cast<double>( *data * *data );
                    double tmp = ( static_cast<double>( *data++ ) - mean );
                    stddev += tmp * tmp;
                }
                rms = sqrt( rms / n );
                stddev = sqrt( stddev / (n-1) );
            } else if (n == 1) {
                rms = *data;
                stddev = 0;
            }
        }

        template <typename T>
        double rmsStddev( const T* data, size_t n, double& rms, double& stddev ) {
            if( n > 1 ) {
                rms = stddev = 0;
                double mn = mean( data, n );
                for( size_t i(0); i<n; ++i ) {
                    double tmp = static_cast<double>( *data++ );
                    rms += tmp * tmp;
                    tmp -= mn;
                    stddev += tmp * tmp;
                }
                rms = sqrt( rms / n );
                stddev = sqrt( stddev / (n-1) );
                return mn;
            } else if( n == 1 ) {
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
        
        template<typename T, typename U>
        void invert_4x4( const T* in, U* out ) {
            
            const double s0 = in[0] * in[5] - in[4] * in[1];
            const double s1 = in[0] * in[6] - in[4] * in[2];
            const double s2 = in[0] * in[7] - in[4] * in[3];
            const double s3 = in[1] * in[6] - in[5] * in[2];
            const double s4 = in[1] * in[7] - in[5] * in[3];
            const double s5 = in[2] * in[7] - in[6] * in[3];

            const double c5 = in[10] * in[15] - in[14] * in[11];
            const double c4 = in[ 9] * in[15] - in[13] * in[11];
            const double c3 = in[ 9] * in[14] - in[13] * in[10];
            const double c2 = in[ 8] * in[15] - in[12] * in[11];
            const double c1 = in[ 8] * in[14] - in[12] * in[10];
            const double c0 = in[ 8] * in[13] - in[12] * in[ 9];
            
            const double invdet = 1 / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0);
            
            out[0] = (in[5] * c5 - in[6] * c4 + in[7] * c3) * invdet;
            out[1] = (-in[1] * c5 + in[2] * c4 - in[3] * c3) * invdet;
            out[2] = (in[13] * s5 - in[14] * s4 + in[15] * s3) * invdet;
            out[3] = (-in[9] * s5 + in[10] * s4 - in[11] * s3) * invdet;

            out[4] = (-in[4] * c5 + in[6] * c2 - in[7] * c1) * invdet;
            out[5] = (in[0] * c5 - in[2] * c2 + in[3] * c1) * invdet;
            out[6] = (-in[12] * s5 + in[14] * s2 - in[15] * s1) * invdet;
            out[7] = (in[8] * s5 - in[10] * s2 + in[11] * s1) * invdet;

            out[8] = (in[4] * c4 - in[5] * c2 + in[7] * c0) * invdet;
            out[9] = (-in[0] * c4 + in[1] * c2 - in[3] * c0) * invdet;
            out[10] = (in[12] * s4 - in[13] * s2 + in[15] * s0) * invdet;
            out[11] = (-in[8] * s4 + in[9] * s2 - in[11] * s0) * invdet;

            out[12] = (-in[4] * c3 + in[5] * c1 - in[6] * c0) * invdet;
            out[13] = (in[0] * c3 - in[1] * c1 + in[2] * c0) * invdet;
            out[14] = (-in[12] * s3 + in[13] * s1 - in[14] * s0) * invdet;
            out[15] = (in[8] * s3 - in[9] * s1 + in[10] * s0) * invdet;

        }
        
        template<typename T, typename U>
        void invert_3x3( const T* in, U* out ) {
            
            const double i37 = in[3]*in[7];
            const double i38 = in[3]*in[8];
            const double i46 = in[4]*in[6];
            const double i48 = in[4]*in[8];
            const double i56 = in[5]*in[6];
            const double i57 = in[5]*in[7];
            
            const double invdet = 1 / (in[0]*(i48 - i57) - in[1]*(i38 - i56) + in[2]*(i37 - i46));
            
            out[0] = (i48 - i57) * invdet;
            out[1] = (in[2]*in[7] - in[1]*in[8]) * invdet;
            out[2] = (in[1]*in[5] - in[2]*in[4]) * invdet;
            out[3] = (i56 - i38) * invdet;
            out[4] = (in[0]*in[8] - in[2]*in[6]) * invdet;
            out[5] = (in[2]*in[3] - in[0]*in[5]) * invdet;
            out[6] = (i37 - i46) * invdet;
            out[7] = (in[1]*in[6] - in[0]*in[7]) * invdet;
            out[8] = (in[0]*in[4] - in[1]*in[3]) * invdet;

        }
        
        template<typename T, typename U>
        void invert_2x2( const T* in, U* out ) {
            
            const double invdet = 1 / (in[0]*in[3] - in[1]*in[2]);
            
            out[0] =  in[3] * invdet;
            out[1] = -in[2] * invdet;
            out[2] = -in[1] * invdet;
            out[3] =  in[0] * invdet;

        }
        
    }
}

#endif // REDUX_MATH_HELPERS_HPP
