#ifndef REDUX_TYPES_HPP
#define REDUX_TYPES_HPP

#include <complex>
#include <cstdint>

#ifdef RDX_WITH_FFTW3
#   include <fftw3.h>
#else
    typedef std::complex<double> fftw_complex;
#endif

namespace redux {


        /*!  @ingroup redux
         *  @{
         */
        
        /*!  @file      types.hpp
         *   @brief     Generic type definitions, commonly used all over the library
         *   @name      Types
         * 
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2014
         */
        
        
        typedef std::complex<double> complex_t;

        using std::int8_t;
        using std::int16_t;
        using std::int32_t;
        using std::int64_t;

        using std::uint8_t;
        using std::uint16_t;
        using std::uint32_t;
        using std::uint64_t;

#ifdef RDX_WITH_FFTW3
        inline const complex_t& operator*=( complex_t& lhs, const fftw_complex& rhs ) { return lhs *= reinterpret_cast<const complex_t&>(rhs); }
        inline const fftw_complex& operator*=( fftw_complex& lhs, const complex_t& rhs ) { reinterpret_cast<complex_t&>(lhs) *= rhs; return lhs; }
#endif
        
        /* @} */

}



#endif // REDUX_TYPES_HPP


