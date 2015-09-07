#ifndef REDUX_UTIL_GSL_HPP
#define REDUX_UTIL_GSL_HPP

#include <gsl/gsl_multimin.h>

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */

        namespace gsl {

            /*!  @file      gsl.hpp
             *   @brief     Wrappers and helpers for accessing the GSL libraries.
             *   @author    Tomas Hillberg (hillberg@astro.su.se)
             *   @date      2015
             */

            typedef double( gsl_f_t )( const gsl_vector*, void* );
            typedef void( gsl_df_t )( const gsl_vector*, void*, gsl_vector* );
            typedef void( gsl_fdf_t )( const gsl_vector*, void*, double*, gsl_vector* );

            template<typename T, typename U, typename V>
            class gsl_fdf_wrapper : public gsl_multimin_function_fdf {
            public:
                gsl_fdf_wrapper( size_t dim, const T& f_in, const U& df_in , const V& fdf_in ) :
                    f_( f_in ), df_( df_in ), fdf_( fdf_in ) {
                    n = dim;
                    f = &gsl_fdf_wrapper::call_f;
                    df = &gsl_fdf_wrapper::call_df;
                    fdf = &gsl_fdf_wrapper::call_fdf;
                    params = this;
                }
            private:
                const T& f_;
                const U& df_;
                const V& fdf_;
                static double call_f( const gsl_vector* x, void *params ) {
                    return static_cast<gsl_fdf_wrapper*>( params )->f_( x, nullptr );
                }
                static void call_df( const gsl_vector* x, void *params, gsl_vector* df ) {
                    static_cast<gsl_fdf_wrapper*>( params )->df_( x, nullptr, df );
                }
                static void call_fdf( const gsl_vector* x, void *params, double* f, gsl_vector* df ) {
                    static_cast<gsl_fdf_wrapper*>( params )->fdf_( x, nullptr, f, df );
                }
            };
        
        
        }   // gsl
        
        /*! @} */

    }   // util

}   // redux

#endif // REDUX_UTIL_GSL_HPP
