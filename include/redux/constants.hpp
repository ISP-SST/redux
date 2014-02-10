#ifndef REDUX_CONSTANTS_HPP
#define REDUX_CONSTANTS_HPP

#include <math.h>

namespace redux {
    
    /*!  @ingroup redux
     *  @{
     */

    /*!  @file      constants.hpp
     *   @brief     Some physical/mathematical constants
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2011
     */

    const double CC = 299792458;                                              //!< Speed of light (m/s)
    const double PI = 3.141592653589793238462643383279502884197169399375105;  //!< Pi.  (128-bit precision)
    const double rtod = 180.0 / PI;                                           //!< Conversion factor: radians -> degrees
    const double dtor = PI / 180.0;                                           //!< Conversion factor: degrees -> radians
    const double PI_isqrt = 1 / sqrt( PI );                                   //!< \f$1/\sqrt(\pi)\f$
    const double planck = 4.135667516e-15;                                    //!< Plancks constant @e h (\f$eV \cdot s\f$)
    const double planckbar = 6.58211928e-16;                                  //!< Plancks constant @e \f$\hbar = h/(2 \pi)\f$ (\f$eV \cdot s\f$)
    const double BohrMagneton = 57.883817555 * planck* CC;                    //!< Bohr magneton (Gauss/Ångström)
    //  5.7883817555(79)×10−5    eV/Tesla
    //    => = 5.7883817555×10−9     eV/Gauss
    //    => = 5.7883817555×10−9 * planck*CC*1e10    /(Gauss*Ångström)

    /*! @} */


}

#endif  // REDUX_CONSTANTS_HPP
