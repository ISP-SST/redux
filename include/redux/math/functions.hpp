#ifndef REDUX_MATH_FUNCTIONS_HPP
#define REDUX_MATH_FUNCTIONS_HPP

#include <algorithm>
#include <complex>
#include <cstdlib>

namespace redux {

    namespace math {

        /*!  @ingroup math
         *  @{
         */

        /*!  @file      functions.hpp
         *   @brief     Collection of mathematical functions.
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         */

        /** @brief Extend the std-functors to handle different types.
         *  @{
         */
        template<typename T>
        inline T sqr( const T& a ) { return a*a; }
        
        template<typename T, typename U>
        struct assign : public std::binary_function<T, U, T> {
            T operator()( const T& a, const U& b ) const { return static_cast<T>(b); }
        };

        template<typename T>
        struct assign<T,T> {
            T operator()( const T& a, const T& b ) const { return b; }
        };

        template<typename T, typename U>
        struct assign<T,std::complex<U>> {
            T operator()( const T& a, const std::complex<U>& b ) const { return static_cast<T>(std::real(b)); }
        };
        
        template<typename T, typename U>
        struct add : public std::binary_function<T, U, T> {
            T operator()( const T& a, const U& b ) const { return static_cast<T>(a + b); }
        };

        template<typename T, typename U>
        struct add<T,std::complex<U>> {
            T operator()( const T& a, const std::complex<U>& b ) const { return static_cast<T>(a + std::real(b)); }
        };

        template<typename T, typename U>
        struct subtract : public std::binary_function<T, U, T> {
            T operator()( const T& a, const U& b ) const { return static_cast<T>(a - b); }
        };

        template<typename T, typename U>
        struct multiply : public std::binary_function<T, U, T> {
            T operator()( const T& a, const U& b ) const { return static_cast<T>(a * b); }
        };

        template<typename T, typename U>
        struct divide : public std::binary_function<T, U, T> {
            T  operator()( const T& a, const U& b ) const { return static_cast<T>(a / b); }
        };

        template<typename T, typename U>
        struct norm_multiply : public std::binary_function<T, U, U> {
            U operator()( const T& a, const U& b ) const { return static_cast<U>(std::norm(a) * b); }
        };

        template<typename T, typename U>
        struct inv_multiply : public std::binary_function<T, U, T> {
            T operator()( const T& a, const U& b ) const { return static_cast<T>(1/a * b); }
        };


        /*!  @brief         Smooth the "tail" of an array towards some specified target value.
         *   @details       Useful to make functions behave better at the edges (especially useful before applying a Fourier transform).
         *   @code
         *      #include "redux/math/functions.hpp"
         *      double data[100];
         *      data[90]=10;
         *      redux::math::apodize(data+90,9,22); // will generate a smooth transition from 10 to 22 over the last 9 values.
         *   @endcode
         *   @param data    Pointer to the last value in the array to be kept, the following \c n values will be modified.
         *   @param n       Number of points to smooth over.
         *   @param target  The desired value at point \code data[n-1] \endcode
         *   @author        Tomas Hillberg <hillberg@astro.su.se>
         */
        template <typename T> void apodize( T* data, size_t n, const T& target );


        /*!  @fn void cauchylorentz( double* data, size_t n, double fwhm, float ncav=1 )
         *   @brief  Generates a normalized Lorentzian profile of length n, and the specified fwhm.
         *   @details  The profile will be centered at the mid-point.
         *   @f[
         *    L(i) = \frac{1}{\Big(\frac{2*i}{fwhm}\Big)^{2*cav}+1}
         *   @f]
         *   @code
         *      #include "redux/math/functions.hpp"
         *      double data[100];
         *      data[90]=10;
         *      redux::math::cauchylorentz(data,100,22,2)
         *   @endcode
         *   @param data    Data array.
         *   @param n       Length of array.
         *   @param fwhm    Full width at half-maximum, in "element" units
         *   @param ncav    For multi-cavity filters the number of cavities can be specified (=1 for standard Lorentzian)
         *   @author        Tomas Hillberg <hillberg@astro.su.se>
         */
        void cauchylorentz( double* data, size_t n, double fwhm, float ncav = 1 );


        /*!  @name Gauss
         *   @brief  Generate a normalized Gauss-profile of length n, centered at the mid-point and with the specified dispersion and fwhm.
         *   @f[
         *    G(i)= \exp\Big(-\Big(\frac{4*i}{fwhm}\Big)^2\Big)
         *   @f]
         *   @code
         *      #include "redux/math/functions.hpp"
         *      double data[100];
         *      redux::math::gauss(data,100,3,55);
         *   @endcode
         *   @param data        Pointer to data array.
         *   @param n           Length of array
         *   @param dispersion  The step length in arbitrary units
         *   @param fwhm        Full width at half-maximum, in arbitrary units
         *   @author            Tomas Hillberg <hillberg@astro.su.se>
         */
        //@{
        void gauss( double* data, size_t n, double dispersion, double fwhm );
        template <typename T> void gauss( T* data, size_t sizeY, size_t sizeX, double fwhmY, double fwhmX, double centerY, double centerX );
        //@}

        /*!  @name Hann
         *   @brief         Generate a normalized Hann-function of length n in array data.
         *   @details       The templated version expects a container of type std::vector or similar.
         *   @f[
         *    H(i)= 0.5*\left(1-\cos(\frac{\pi * i}{(n-1)})\right)
         *   @f]
         *   @code
         *      #include "redux/math/functions.hpp"
         *      std::vector<double> data;
         *      redux::math::hann(data,100);    // data will be resized to 100
         *   @endcode
         *   @param data    Array where to store the output.
         *   @param n       Length of the output array (the templated version will resize the container-class)
         *   @param x       Input value in the range \f$[0, 1]\f$
         *   @returns       The scalar version returns the functional value at the point \c x
         *   @author        Tomas Hillberg <hillberg@astro.su.se>
         */
        //@{
        template <class T>
        void   hann( T& data, size_t n );
        double hann( double x );
        void   hann( double* data, size_t n );
        //@}

        /*!  @name Faraday-Voigt
         *   @brief     Generate normalized Voigt/Faraday profiles of length n, with the given dispersion/damping.
         *   @details   The profiles will centered at the mid-point.
         *   @code
         *      #include "redux/math/functions.hpp"
         *      std::vector<double> vgt(100);
         *      std::vector<double> far(100);
         *      redux::math::voigtfaraday( vgt.data(), far.data(), 100, 10.0/(99), 0.0001);   // Wavelength range will be [-5.0, 5.0] with this dispersion
         *      double vgt;
         *      double far;
         *      redux::math::voigtfaraday( &vgt, &far, 0.005, 0.0001);                        // Get the profile values at a point 0.005 "units" away from line-center
         *   @endcode
         *   @param vgt         Array where to store the Voigt profile
         *   @param far         Array where to store the Faraday profile
         *   @param n           Length of the output arrays.
         *   @param dispersion  Step-size in wavelength (arbitrary units)
         *   @param damping     Strength of wing-damping, 0 will give a Gaussian, large values a Lorentzian.
         *                      What is "large" is of course relative to the dispersion-unit.
         *   @param offset      Specifies distance from line-center
         *   @author            Tomas Hillberg <hillberg@astro.su.se>
         */
        //@{
        //! @brief     Generate normalized Voigt/Faraday profiles of length n, with the given dispersion/damping.
        void faraday_voigt( double* vgt, double* far, size_t n, double dispersion, double damping );
        //! @brief     Generate normalized Voigt/Faraday profiles of length n, with the given dispersion/damping.
        template <class T> void faraday_voigt( T& vgt, T& far, double dispersion, double damping );
        //! @brief     Generate normalized Voigt/Faraday profiles of length n, with the given dispersion/damping.
        template <class T> void faraday_voigt( T& vgt, T& far, T& lambda, double dispersion, double damping );
        //!  @brief     Calculate the Voigt/Faraday profile values at a single point "offset" from line-center, with the given damping.
        void faraday_voigt( double* vgt, double* far, double offset, double damping );
        //@}


        /*!  @fn unsigned long n_choose_k( int n, int k )
         *   @brief     Calculate the binomial coefficient, i.e. the number of ways to select \c k out of \c n objects.
         *   @code
         *      #include "redux/math/functions.hpp"
         *      unsigned long result = redux::math::n_choose_k( 52, 5 )   // 2598960
         *   @endcode
         *   @param n         Total number of objects
         *   @param k         Selection size
         *   @author          Tomas Hillberg <hillberg@astro.su.se>
         */
        unsigned long n_choose_k( int n, int k );


        /*! @} */


    }   // math


}   // redux

#endif // REDUX_MATH_FUNCTIONS_HPP

