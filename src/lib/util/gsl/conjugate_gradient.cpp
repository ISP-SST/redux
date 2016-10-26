#include "redux/util/gsl.hpp"

#include "redux/math/helpers.hpp"
#include "redux/util/stringutil.hpp"

#include <algorithm>
#include <iostream>

/* 
 * An implementation of the conjugate-gradient method for use with the GSL multimin minimizer.
 * This version is optimized for use with the specific minimization details in the MOMFBD method.
 * 
 */

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>

#include <functional>

using namespace redux::util::gsl;
using namespace redux::math;
using namespace redux::util;
using namespace std;


//#define DEBUG_


typedef struct {
    int iter;
    double step;
    double tol;
    gsl_vector *x1;     // temporary location (x)
    double pnorm;
    gsl_vector *p;      // search direction
    double g0norm;
    gsl_vector *previousGradient;     // gradient
} conjugate_rdx_state_t;


static int conjugate_rdx_alloc (void *vstate, size_t n) {
    
    conjugate_rdx_state_t *state = (conjugate_rdx_state_t *) vstate;

    state->x1 = gsl_vector_calloc (n);

    if ( state->x1 == 0 ) {
        GSL_ERROR ("failed to allocate space for x1", GSL_ENOMEM);
    }

    state->p = gsl_vector_calloc (n);

    if (state->p == 0)
    {
        gsl_vector_free (state->x1);
        GSL_ERROR ("failed to allocate space for p", GSL_ENOMEM);
    }

    state->previousGradient = gsl_vector_calloc (n);

    if (state->previousGradient == 0)
    {
        gsl_vector_free (state->p);
        gsl_vector_free (state->x1);
        GSL_ERROR ("failed to allocate space for g0", GSL_ENOMEM);
    }

    return GSL_SUCCESS;
}


static int conjugate_rdx_set (void *vstate, gsl_multimin_function_fdf * fdf,
                  const gsl_vector * x, double *f, gsl_vector * gradient,
                  double step_size, double tol) {
    
    conjugate_rdx_state_t *state = (conjugate_rdx_state_t *) vstate;

    state->iter = 0;
    state->step = step_size;
    state->tol = tol;

    *f = GSL_MULTIMIN_FN_EVAL_F( fdf, x );

    return GSL_SUCCESS;
    
}


static void conjugate_rdx_free (void *vstate) {
    
    conjugate_rdx_state_t *state = (conjugate_rdx_state_t *) vstate;

    gsl_vector_free (state->previousGradient);
    gsl_vector_free (state->p);
    gsl_vector_free (state->x1);
    
}


static int conjugate_rdx_restart (void *vstate) {
    
    conjugate_rdx_state_t *state = (conjugate_rdx_state_t *) vstate;

    state->iter = 0;
    return GSL_SUCCESS;
    
}


void take_step( const gsl_vector * x, const gsl_vector * p,
           double step, double lambda, gsl_vector * x1, gsl_vector * dx) {

    gsl_vector_set_zero (dx);
    gsl_blas_daxpy (-step * lambda, p, dx);

    gsl_vector_memcpy (x1, x);
    gsl_blas_daxpy (1.0, dx, x1);
    
}


double value_at( double step, gsl_multimin_function_fdf* fdf, gsl_vector *x, gsl_vector *p, double scale, gsl_vector *x1, gsl_vector *dx ) {
    take_step( x, p, step, scale, x1, dx );
    return GSL_MULTIMIN_FN_EVAL_F( fdf, x1 );
}


static int conjugate_rdx_iterate( void *vstate, gsl_multimin_function_fdf * fdf,
                      gsl_vector * x, double *f,
                      gsl_vector * gradient,
                      gsl_vector * dx ) {
    
    conjugate_rdx_state_t *state = (conjugate_rdx_state_t*)vstate;

    gsl_vector *x1 = state->x1;
    gsl_vector *p = state->p;
    gsl_vector *oldGradPtr = state->previousGradient;
    
    size_t n = x->size;
    double *gData = gradient->data;
    double *pData = state->p->data;

    // Get the gradient at x
    GSL_MULTIMIN_FN_EVAL_DF( fdf, x, gradient );
    
    double pnorm = state->pnorm;
    double g0norm = state->g0norm;
    double g1norm = gsl_blas_dnrm2( gradient );
    double pg;
        
    // Choose a new conjugate direction for the next step
    if( state->iter++ == 0 ) {
        gsl_vector_memcpy( p, gradient );
        pnorm = state->pnorm = g1norm;
        pg = g1norm*g1norm;
    } else {
        
        if( g0norm == 0.0 ) {
            gsl_vector_set_zero (dx);
            return GSL_ENOPROG;
        }
        
        // p' = g1 - beta * p

        double g0g1(0);
        std::transform( gData, gData+n, oldGradPtr->data, oldGradPtr->data,
                        [&g0g1]( const double& a, const double& b ) {
                            double tmp = b-a;       // g0' = g0 - g1
                            g0g1 += tmp*a;          // g1g0 = (g0-g1).g1
                            return tmp;
                        });
        double beta = g0g1 / (g0norm*g0norm);               // Polak-Ribiere:    beta = -((g1-g0).g1)/(g0.g0)
        //beta = std::min( beta, 0.0 );
        //double beta = -pow( g1norm / g0norm, 2.0 );         // Fletcher-Reeve
        std::transform( gData, gData+n, pData, pData,
                        [beta]( const double& a, const double& b ) {
                            return a - beta*b;
                        });
        
        pnorm = state->pnorm = gsl_blas_dnrm2( p );
        gsl_blas_ddot( p, gradient, &pg );
        
    }
    
    // Determine which direction is downhill, +p or -p
    double lambda = (pg >= 0.0) ? +1.0/pnorm : -1.0/pnorm;

#ifdef DEBUG_
    printf ("gnorm: %.12e   pnorm: %.12e   pg: %.12e   orth: %g\n", g1norm, pnorm, pg, fabs(pg * lambda/ g1norm));
#endif
    
    if( pnorm == 0.0 ) {
        gsl_vector_set_zero (dx);
        return GSL_ENOPROG;
    }

    double fa = *f, fb, fc;
    double stepa = 0.0, stepb, stepc = 1E-21/lambda; //1E-21/lambda; //state->step;

    rdx_fdf* rdx_fdf_ptr = nullptr;
    if( fdf == fdf->params ) {      // test if it is a rdx_fdf type.
        rdx_fdf_ptr = static_cast<rdx_fdf*>( fdf );
        if( rdx_fdf_ptr->hasPreCalc ) {
            rdx_fdf_ptr->preCalc( x, p );
        } else {
            rdx_fdf_ptr = nullptr;  // do the ordinary minimization
        }
    }

    std::function<double(double)> my_value_at;
    if( rdx_fdf_ptr ) {
        my_value_at = [rdx_fdf_ptr,lambda]( const double& val ){ return rdx_fdf_ptr->evaluateAt(-lambda*val); };
    } else {
        my_value_at = std::bind( value_at, std::placeholders::_1, fdf, x, p, lambda, x1, dx );
    }

    // Evaluate function value at stepc
    fc = my_value_at(stepc);
#ifdef DEBUG_
    printf ("got a=%.12e  f(a)=%.12e  step=%.12e f(step)=%.12e\n", stepa, fa, stepc, fc);
#endif

    // Find a region (xa,fa) (xc,fc) which contains an intermediate (xb,fb) satisifying fa > fb < fc.
    bracket( my_value_at, stepa, stepc, stepb, fa, fc, fb );

#ifdef DEBUG_
    printf ("got pos=[%.12e,%.12e,%.12e] => val=[%.12e,%.12e,%.12e]\n", stepa, stepb, stepc, fa, fb, fc);
#endif

    // Do a line minimization using Brent's method.
    brent( my_value_at, stepa, stepc, stepb, fb, 0.1 );

#ifdef DEBUG_
    printf( "got min at = %.12e f(x) = %.12e    lambda=%.12e   lambda*x=%.12e\n", stepb, fb, lambda, lambda*stepb );
#endif
  //exit(0);

    // move x to new location
    take_step( x, p, stepb, lambda, x1, dx);
    gsl_vector_memcpy( x, x1 );
    *f = fb;
    state->step = 0.5*std::max(fabs(stepb),2.0E-9);

    state->g0norm = g1norm;
    gsl_vector_memcpy( oldGradPtr, gradient );
    
    state->iter = state->iter % x->size;

#ifdef DEBUG_
    cout << "Updated conjugate directions.  iter=" << state->iter << endl;
    //cout << printArray(p->data,p->size,"p") << endl;
    //cout << printArray(gradient->data,gradient->size,"gradient") << endl;
#endif

    return GSL_SUCCESS;
}


static const gsl_multimin_fdfminimizer_type conjugate_rdx_type = {
    "conjugate_gradient_rdx",               //
    sizeof (conjugate_rdx_state_t),
    &conjugate_rdx_alloc,
    &conjugate_rdx_set,
    &conjugate_rdx_iterate,
    &conjugate_rdx_restart,
    &conjugate_rdx_free
};

const gsl_multimin_fdfminimizer_type* redux::util::gsl::multimin_fdfminimizer_conjugate_rdx = &conjugate_rdx_type;
