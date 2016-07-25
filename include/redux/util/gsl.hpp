#ifndef REDUX_UTIL_GSL_HPP
#define REDUX_UTIL_GSL_HPP

#include <gsl/gsl_multimin.h>

#include <functional>

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
            typedef void( gsl_precalc_t )( const gsl_vector*, const gsl_vector*, void* );
            typedef double( gsl_step_t )( const double, void* );
            
            struct rdx_fdf : public gsl_multimin_function_fdf {
            public:
                rdx_fdf( size_t dim,
                                 const std::function<gsl_f_t>& f_in,
                                 const std::function<gsl_df_t>& df_in,
                                 const std::function<gsl_fdf_t>& fdf_in,
                                 void* params_in=nullptr ) : hasPreCalc(false),
                    f_(f_in), df_(df_in), fdf_(fdf_in), realParams(params_in) {
                    n = dim;
                    f = &rdx_fdf::call_f;
                    df = &rdx_fdf::call_df;
                    fdf = &rdx_fdf::call_fdf;
                    params = this;
                }
                
                void setPreCalc( const std::function<gsl_precalc_t>& precalc_func,
                                 const std::function<gsl_step_t>& eval_func ) {
                    preCalc_ = precalc_func;
                    evaluateAt_ = eval_func;
                    if( preCalc_ && evaluateAt_ ) hasPreCalc = true;
                    else hasPreCalc = false;
                }
                void preCalc( const gsl_vector* x, const gsl_vector* df ) {
                    this->preCalc_( x, df, this->realParams );
                }
                double evaluateAt( double pos ) {
                    return this->evaluateAt_( pos, this->realParams );
                }
                
                bool hasPreCalc;
                
            private:
                std::function<gsl_f_t> f_;
                std::function<gsl_df_t> df_;
                std::function<gsl_fdf_t> fdf_;
                
                std::function<gsl_precalc_t> preCalc_;
                std::function<gsl_step_t> evaluateAt_;
                
                void* realParams;
                
                static double call_f( const gsl_vector* x, void *params ) {
                    rdx_fdf* thisPtr = static_cast<rdx_fdf*>( params );
                    return thisPtr->f_( x, thisPtr->realParams );
                }
                static void call_df( const gsl_vector* x, void *params, gsl_vector* df ) {
                    rdx_fdf* thisPtr = static_cast<rdx_fdf*>( params );
                    thisPtr->df_( x, thisPtr->realParams, df );
                }
                static void call_fdf( const gsl_vector* x, void *params, double* f, gsl_vector* df ) {
                    rdx_fdf* thisPtr = static_cast<rdx_fdf*>( params );
                    thisPtr->fdf_( x, thisPtr->realParams, f, df );
                }
            };
            
            extern const gsl_multimin_fdfminimizer_type *multimin_fdfminimizer_conjugate_rdx;

        }   // gsl
        
        /*! @} */

    }   // util

}   // redux

#endif // REDUX_UTIL_GSL_HPP
