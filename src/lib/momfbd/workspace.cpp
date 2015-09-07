#include "redux/momfbd/workspace.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/wavefront.hpp"

#include "redux/image/utils.hpp"
#include "redux/logger.hpp"

#include "redux/file/fileana.hpp"

#include <functional>
#include <limits>

#include <gsl/gsl_blas.h>
#include <boost/timer/timer.hpp>

using namespace redux::file;
using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::util;
using namespace redux::util::gsl;
using namespace redux;
using namespace std;
namespace sp=std::placeholders;

#define lg Logger::mlg
namespace {

    const std::string thisChannel = "workspace";

    grad_t grads[] = { nullptr,
        std::bind(&SubImage::gradFiniteDifference, sp::_1, sp::_2, sp::_3),
        std::bind(&SubImage::gradVogel, sp::_1, sp::_2, sp::_3)
    };

}


WorkSpace::WorkSpace( const MomfbdJob& j ) : objects( j.getObjects() ), job(j), nFreeParameters(0) {
    //std::cout << "WorkSpace():  first = (" << d->first.x << "," << d->first.y << ")  last = (" << d->last.x << "," << d->last.y;
    //std::cout << ")   id = (" << d->index.x << "," << d->index.y << ")" << std::endl;

}


WorkSpace::~WorkSpace() {
    clear();
    //  std::cout << "~WorkSpace():  1   wfsz = " << wavefronts.size() << " osz = " << objects.size() << std::endl;

}


void WorkSpace::init( void ) {

    int nPixels = job.patchSize;//std::min( d->nPixelsY(),d->nPixelsX() );
    window.resize( nPixels, nPixels );
    window = 1.0;
    redux::image::apodizeInPlace( window, nPixels / 8);
    //redux::file::Ana::write( "window.f0", window );

    int md = std::min( 256, nPixels );                       // new size (maximum=256)
    md -= ( md % 2 );
    //int m_offs = std::max( ( nPixels-md )/2, 0 );            // centering needed?
    noiseWindow.resize(md,md);
    noiseWindow = 1.0;
    redux::image::apodizeInPlace( noiseWindow, md / 8);     // FIXME: old code spacifies md/16, but applies it after "window", so it is actually the product...
    //redux::file::Ana::write( "noisewindow.f0", noiseWindow );

    tilts.clear();
    wavefronts.clear();

    for( auto & it : job.getObjects() ) {
        it->initProcessing( shared_from_this() );
    }
 //   cout << "WorkSpace::init(1)  pP = " << job.pupilPixels << "  nFreeParameters = " << nFreeParameters << "  l: " << __LINE__ << endl;

    // Now the wavefront/tilt arrays are populated, so we can find out the total number of free parameters
    nFreeParameters = 0;
    for( auto & it : tilts ) {
        nFreeParameters += it.second->nFreeParameters();
    }
    for( auto & it : wavefronts ) {
        nFreeParameters += it.second->nFreeParameters();
    }

  //  cout << "WorkSpace::init(2)  pS = " << job.patchSize << "  nFreeParameters = " << nFreeParameters << "  l: " << __LINE__ << endl;
    alpha = sharedArray<double>( job.globalData->constraints.nParameters );
    grad_alpha = sharedArray<double>( job.globalData->constraints.nParameters );
    saved_alpha = sharedArray<double>( job.globalData->constraints.nParameters );
    
    beta = gsl_vector_alloc( job.globalData->constraints.nConstrainedParameters );
    grad_beta = gsl_vector_alloc( job.globalData->constraints.nConstrainedParameters );
    beta_init = gsl_vector_alloc( job.globalData->constraints.nConstrainedParameters );
    memset( beta_init->data, 0, job.globalData->constraints.nConstrainedParameters * sizeof( double ) );
    
    memset( saved_alpha.get(), 0, nFreeParameters * sizeof( double ) );
    alpha_init.owner = 0;
    alpha_init.stride = 1;
    alpha_init.size = nFreeParameters;

    alpha_init.data = saved_alpha.get();
    double* a = alpha.get();
    double* g = grad_alpha.get();
    size_t count(0);
    for( auto &it : tilts ) {
        count += it.second->setPointers( a+count, g+count );
    }
    for( auto &it : wavefronts ) {
        it.second->init(job.pupilPixels);
        count += it.second->setPointers( a+count, g+count );
    }


//    cout << "WorkSpace::init()  pS = " << job.patchSize << "  nFreeParameters = " << nFreeParameters << "  count = " << count << endl;

}


void WorkSpace::getMetric( boost::asio::io_service& service, uint8_t nThreads ) {

    for( auto & it : job.getObjects() ) {
        //it->addAllPQ();
        it->metric();
    }

}


void WorkSpace::resetAll( boost::asio::io_service& service ) {
    for( auto& oit: job.getObjects() ) {
        for( auto& cit: oit->channels ) {
            for( auto& it: cit->subImages ) {
                service.post( std::bind((void(SubImage::*)(void))&SubImage::resetPhi, it.get() ) );
            }
        }
        service.post( std::bind(&Object::initPQ, oit.get() ) );
        service.post( std::bind(&Object::addAllPQ, oit.get() ) );
    }
}
void WorkSpace::dumpImages( boost::asio::io_service& service, string tag ) {
    for( auto& oit: job.getObjects() ) {
        for( auto& cit: oit->channels ) {
            for( auto& it: cit->subImages ) {
                service.post( std::bind((void(SubImage::*)(string))&SubImage::dump, it.get(), tag ) );
            }
        }
    }
}

double WorkSpace::my_test( boost::asio::io_service& service, const gsl_vector* x, void* params, gsl_vector* df, string tag ) {

  //  cout << "my_test:  " << printArray(x->data,x->size,"  x") << flush;
    memcpy(alpha.get(),x->data,nFreeParameters*sizeof(double));
    //resetAll(service);
    //runThreadsAndWait( service, job.info.maxThreads );
    for( auto &it : tilts ) {
        it.second->applyTiltToImages(service);
    }
    runThreadsAndWait( service, job.info.maxThreads );
    for( auto & it : wavefronts ) {
        it.second->setAlphasAndUpdate( service, false ); //job.gradientMethod==GM_VOGEL );
    }
    runThreadsAndWait( service, job.info.maxThreads );
    //return 0;
    double m = objectMetric(service);
    //dumpImages(service, tag);
    //runThreadsAndWait( service, job.info.maxThreads );
  //  cout << "   m = " << m << endl;
    memset(grad_alpha.get(),0,nFreeParameters*sizeof(double));
    for( auto & it : wavefronts ) {
        it.second->setAlphasAndUpdate( service, true ); //job.gradientMethod==GM_VOGEL );
    }
    runThreadsAndWait( service, job.info.maxThreads );
    for( auto & it : tilts ) {
        it.second->calcGradient( service, job.info.maxThreads, gradient );
    }
    for( auto & it : wavefronts ) {
        it.second->calcGradient( service, gradient );
    }
    runThreadsAndWait( service, job.info.maxThreads );

    memcpy(df->data,grad_alpha.get(),nFreeParameters*sizeof(double));
    
 //   cout << "my_test:  " << printArray( df->data, df->size, "gradient" ) << endl;
    return m;

}

double WorkSpace::my_f( boost::asio::io_service& service, const gsl_vector* x, void* params ) {

    //cout << "my_f:  " << printArray(x->data,x->size,"  x", 6) << endl;
    //memcpy(alpha.get(),x->data,nFreeParameters*sizeof(double));
    //memset(alpha.get(),0,nFreeParameters*sizeof(double));
    //alpha.get()[0] = 8.75867e-08;
    job.globalData->constraints.reverse(x->data,alpha.get());
    //cout << "my_f:  " << printArray(x->data,x->size,"  x", 6) << endl;
    size_t nModes = job.modeNumbers.size();
    double* alphaPtr = alpha.get();
    for( const shared_ptr<Object>& o: objects ) {
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                //im->setAlphas(job.modeNumbers, alphaPtr);
                service.post ([this, &im, alphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->setAlphas(job.modeNumbers, alphaPtr);
                    im->update(false);
                });
                alphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, job.info.maxThreads );

    double m = objectMetric(service);
    //cout << "my_f:  " << printArray(x->data,x->size,"beta", 6) << endl;
    //cout << "my_f:  " << printArray(alpha.get(),job.globalData->constraints.nParameters,"alpha", 6) << endl;
    //cout << "my_f:  metric = "  << setprecision(15) << m << endl;
    return m;

}


void WorkSpace::my_df( boost::asio::io_service& service, const gsl_vector* x, void* params, gsl_vector* df ) {
    
    job.globalData->constraints.reverse(x->data,alpha.get());
    //cout << "my_df:  " << printArray(x->data,x->size,"beta", 6) << endl;
    //cout << "my_df:  " << printArray(alpha.get(),job.globalData->constraints.nParameters,"alpha", 6) << endl;
    memset(grad_alpha.get(),0,job.globalData->constraints.nParameters*sizeof(double));
    
    bool atOrigin = checkAllZero(x->data,x->size);
    if( atOrigin ) {
        //cout << "my_df:  atOrigin:" << endl;
        gradient = grads[GM_DIFF];
    } else gradient = grads[GM_VOGEL]; //grads[job.gradientMethod];

    size_t nModes = job.modeNumbers.size();
    double* alphaPtr = alpha.get();
    double* gAlphaPtr = grad_alpha.get();
    for( const shared_ptr<Object>& o: objects ) {
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                //im->setAlphas(job.modeNumbers, alphaPtr);
                service.post ([this, &im, alphaPtr, gAlphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->setAlphas(job.modeNumbers, alphaPtr);
                    im->update(true);
                    int count(0);
                    for ( auto& m: job.modeNumbers ) {
                        gAlphaPtr[count++] = gradient(*im, m, 1E-2 );
                    }
                });
                alphaPtr += nModes;
                gAlphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, job.info.maxThreads );
    
    job.globalData->constraints.apply(grad_alpha.get(),df->data);
    //cout << printArray( grad_alpha.get(), job.globalData->constraints.nParameters, "dALPHA", 8 ) << endl;
    //cout << printArray( df->data, df->size, "dBETA", 8 ) << endl;
}

void WorkSpace::my_fdf( boost::asio::io_service& service, const gsl_vector* x, void* params, double* f, gsl_vector* df ) {
    job.globalData->constraints.reverse(x->data,alpha.get());
    memset(grad_alpha.get(),0,job.globalData->constraints.nParameters*sizeof(double));
    
    bool atOrigin = checkAllZero(x->data,x->size);
    if( atOrigin ) {
      //  cout << "my_fdf:  atOrigin:" << endl;
        gradient = grads[GM_DIFF];
    } else gradient = grads[GM_VOGEL]; //grads[job.gradientMethod];

    size_t nModes = job.modeNumbers.size();
    double* alphaPtr = alpha.get();
    double* gAlphaPtr = grad_alpha.get();
    for( const shared_ptr<Object>& o: objects ) {
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                //im->setAlphas(job.modeNumbers, alphaPtr);
                service.post ([this, &im, alphaPtr, gAlphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->setAlphas(job.modeNumbers, alphaPtr);
                    im->update(true);
                    int count(0);
                    for ( auto& m: job.modeNumbers ) {
                        gAlphaPtr[count++] = gradient(*im, m, 1E-2 );
                    }
                });
                alphaPtr += nModes;
                gAlphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, job.info.maxThreads );
    
    *f = objectMetric(service);
    //cout << "   f = " << setprecision(15) << *f << endl;
    job.globalData->constraints.apply(grad_alpha.get(),df->data);
    //cout << "my_fdf:  metric = " << *f << endl;
    //cout << printArray( grad_alpha.get(), job.globalData->constraints.nParameters, "DALPHA", 8 ) << endl;
    //cout << printArray( df->data, df->size, "DBETA", 8 ) << endl;
}


void WorkSpace::run( PatchData::Ptr p, boost::asio::io_service& service, uint8_t nThreads ) {

 //   cout << "RUN:  this = " << hexString( this ) << "  serv = " << hexString( &service ) << "  dat = " << hexString( alpha_init.data ) << endl;
//    if( p->id != 1 ) return;


    //if( job.gradientMethod == GM_VOGEL) {
     //   gradient = std::bind(&SubImage::gradientVogel, sp::_1, sp::_2, sp::_3, sp::_4, sp::_5);
   // } else {
   // }
    
    
    for( auto & it : wavefronts ) {
        if( it.second ) {
            it.second->reset();     // clear image-list and set alpha-values to 0 (but keep the list/weights of the alphas)
        } else LOG << "WorkSpace::run(1)  it.second is NULL";
    }


    p->initPatch();
//dump("init");
    data = p;
//     result->id = data->id;
//     result->index = data->index;
//     result->pos = data->pos;
    LOG << "WorkSpace::run()  patch#" << data->id << "   index=" << data->index << " region=" << data->roi;

//   double gv_alpha[] = { 0, 0, 3.83e-15, 2.99e-15, 2.84e-15 };
// double gv_alpha[] = { 1.08e-07, -1.67e-07, -1.08e-07, 1.67e-07,
//     2.46e-10, 1.58e-10, 1.56e-10,
//     -1.14e-10, -5.76e-11, -6.15e-11
// };
//     double gv_alpha[] = { 7.63e-10,  4.92e-10,  4.84e-10, -2.25e-07,  1.53e-07,
//                          -3.55e-10, -1.79e-10, -1.91e-10,  2.27e-07, -1.5e-07 };

 // cout << " nFreeParams = " << nFreeParameters << endl;
//memcpy(alpha_init.data,gv_alpha,nFreeParameters*sizeof(double));
// 
//     for( auto & it : job.getObjects() ) {
//         it->addAllPQ();
//         //it->metric();
//     }
// 
// 
// return;
//     /********* gsl test ********/


   // std::function<gsl_test_t> wrapped_test = std::bind( &WorkSpace::my_test, this, std::ref( service ), sp::_1, sp::_2, sp::_3, sp::_4 );
    std::function<gsl_f_t> wrapped_f = std::bind( &WorkSpace::my_f, this, std::ref( service ), sp::_1, sp::_2 );
    std::function<gsl_df_t> wrapped_df = std::bind( &WorkSpace::my_df, this, std::ref( service ), sp::_1, sp::_2, sp::_3 );
    std::function<gsl_fdf_t> wrapped_fdf = std::bind( &WorkSpace::my_fdf, this, std::ref( service ), sp::_1, sp::_2, sp::_3 , sp::_4);

    gsl_fdf_wrapper<decltype( wrapped_f ), decltype( wrapped_df ), decltype( wrapped_fdf )>
    my_func( job.globalData->constraints.nConstrainedParameters, wrapped_f, wrapped_df, wrapped_fdf );
    
                        
    //size_t nAlpha = job.globalData->constraints.nParameters;       
    size_t nBeta = job.globalData->constraints.nConstrainedParameters;       
    //double bla;
    for(uint i=0; i<1; ++i) {
        memset(beta->data,0,nBeta*sizeof(double));
        //gsl_vector_set(beta, i, 1.0);
        wrapped_f( beta, nullptr );
        //wrapped_df( beta, nullptr, grad_beta );
       // wrapped_fdf( beta, nullptr, &bla, grad_beta );
    }
   // cout << "Metric = " << bla << endl; 
  //  return;


//    dump("test");
//    return;
//usleep(1000000);

    /***  Always run first iteration using finite-difference & steepest-descent ***/
    gradient = grads[GM_DIFF];
    //gradient = grads[job.gradientMethod];
    
    double init_step = 1e-1; //1e-18;
    double init_tol = 0.1;
    const gsl_multimin_fdfminimizer_type *minimizerType = gsl_multimin_fdfminimizer_steepest_descent;
    switch(job.getstepMethod) {
        case GSM_CNJG:
            LOG_DETAIL << "Using Conjugate-Gradient solver.";
            //minimizerType = gsl_multimin_fdfminimizer_conjugate_fr;
            minimizerType = gsl_multimin_fdfminimizer_conjugate_pr;
            //init_tol = 0; //1e-4;
            break;
        case GSM_BFGS:
            LOG_DETAIL << "Using BFGS solver.";
            minimizerType = gsl_multimin_fdfminimizer_vector_bfgs;
            break;
        case GSM_BFGS_inv:
            LOG_DETAIL << "Using BFGS-2 solver.";
            minimizerType = gsl_multimin_fdfminimizer_vector_bfgs2;
            //init_step = 1e-18; //1e-18;
            break;
        case GSM_SDSC:
            LOG_DETAIL << "Using Steepest-Descent solver.";
            minimizerType = gsl_multimin_fdfminimizer_steepest_descent;
            break;
        default:
            LOG << "No solver specified, using Steepest-Descent as default.";
            minimizerType = gsl_multimin_fdfminimizer_steepest_descent;
    }
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc( minimizerType, job.globalData->constraints.nConstrainedParameters );
    //gsl_multimin_fdfminimizer_set( s, &my_func, &beta_init, 0.01, 1e-4 );
//cout << "******************************************************************" << endl;
    gsl_multimin_fdfminimizer_set( s, static_cast<gsl_multimin_function_fdf*>(&my_func), beta_init, init_step, init_tol );

//cout << endl << endl;
//cout << "FTOL = " << job.FTOL << "   EPS = " << job.EPS << "   done_iter = " << job.targetIterations << endl;
//cout << "init_step = " << init_step << "   init_tol = " << init_tol << endl;
//cout << "******************************************************************" << endl;
wrapped_df( beta_init, nullptr, s->dx );
//cout << "******************************************************************" << endl << endl;
//return;

//return;    
   // memset(beta_init.data,0,nFreeParameters*sizeof(double));
    //beta_init.data[0] = 1.0;
    //gsl_vector_set(&beta_init, 0, 1);
  //  wrapped_f( &beta_init, nullptr );
boost::timer::auto_cpu_timer timer;
    double previousMetric(0), thisMetric;//, gradNorm;
    size_t iter = 0;
    int status, successCount(0);
    bool done(false);
    do {
        status = gsl_multimin_fdfminimizer_iterate( s );
        if( status ) {
            cout << "gsl_err1: " << gsl_strerror(status) << endl;
        }
        thisMetric = s->f;
        //gradNorm = gsl_blas_dnrm2(s->gradient)/(thisMetric*thisMetric);
        //cout << iter << ":  metric = " << thisMetric << "   gNorm = " << gradNorm << flush;
        if (iter++) {
            double relativeChange = 2.0 * fabs(thisMetric-previousMetric) / (fabs(thisMetric)+fabs(previousMetric)+job.EPS);
           // cout << "  rC = " << relativeChange << "  FTOL = " << job.FTOL << "  EPS = " << job.EPS << "  " << flush;
            if(relativeChange < job.FTOL) {      // count consecutive decreases in metric
                successCount++;
               // cout << "  successCount = " << successCount << "  " << flush;
                if( /*false &&*/ successCount >= job.targetIterations ) { // exit after targetIterations consecutive improvements.
//                    cout << "  successCount = " << successCount << "  " << flush;
                    done = true;
                }
            } else successCount = 0;
        }
        previousMetric = thisMetric;
      //  cout << endl;
        status = gsl_multimin_test_gradient( s->gradient, 1e-9 );
//         if( status ) {
//             if( status == GSL_SUCCESS ) {
//                 printf( "Minimum found at:\n" );
//             } else {
//                 cout << "gsl_err2: " << gsl_strerror(status) << endl;
//             }
//         }
//        cout << "  thisMetric = " << thisMetric << printArray(s->x->data,s->x->size,"\n  alpha")
//             << printArray(s->gradient->data,s->gradient->size,"\n  grad") << endl;

    } while( status == GSL_CONTINUE && (!done || iter < job.minIterations) && iter < job.maxIterations );

    cout << "After gsl-loop: " << gsl_strerror(status) << "  iter = " << iter << "  done=" << done << "  mI = " << job.maxIterations << endl;
    
    //memcpy(beta_init->data,s->x->data,job.globalData->constraints.nConstrainedParameters*sizeof(double));
    // do an extra call with the alphas, to make sure all data is set before getting results.
    //GSL_MULTIMIN_FN_EVAL_F( &my_func, beta_init );

    job.globalData->constraints.reverse(s->x->data,alpha.get());

//cout << "FINAL:  " << printArray(s->x->data,s->x->size,"beta", 6) << endl;    
//cout << "FINAL:  " << printArray(alpha.get(),job.globalData->constraints.nParameters,"alpha", 6) << endl;    
//    dump("final");

    gsl_multimin_fdfminimizer_free( s );
    p->collectResults();

    
    return;

    /***************************/




    {

        calcOTFs( service );
        objectMetric( service );


    }
    /*

    ModeLoop:
    if shifted/new: Cut out sub-image and apply window to cutout, get stats, compute fft and add to fftSum
    Initialize PQ, get initial metric
    IterLoop
        Calculate gradient
            OTF/psf
            Line search
            Calculate metric/diff
        Check improvement
            Keep

            Discard & Restore previous
        Keep
            Calculate new alignment from tilt-coefficients
        Discard & Restore previous


    */



    usleep( 10000 );

    // if(data->id == 1) { // DEBUG: only for 1 patch to avoid duplicate output during testing

    //service.post( std::bind( &WaveFront::computePhases, it.second.get() ) );

    for( auto & it : job.getObjects() ) {
        it->addAllPQ();
        it->metric();
    }


    coefficientMetric( service );


}


double WorkSpace::coefficientMetric( boost::asio::io_service& service ) {

    double sum( 0 );

    for( auto & it : wavefronts ) {
        if( it.second ) {
            // service.post( [it,&sum] {
            //sum.fetch_add(it.second->coefficientMetric());
            sum += it.second->coefficientMetric();
            // } );
        } else LOG << "WorkSpace::coefficientMetric()  it.second is NULL";
    }

    //runThreadsAndWait(service, job.info.maxThreads);
//    cout << "WorkSpace::coefficientMetric()   sum = " << sum << endl;
    return sum;
}


void WorkSpace::calcOTFs( boost::asio::io_service& service ) {

    for( auto & it : wavefronts ) {
        if( it.second ) {
            //it.second->setAlpha(alpha_init);
            it.second->setAlphasAndUpdate( service, true );
        } else LOG << "WorkSpace::run(2)  it.second is NULL";
    }

    runThreadsAndWait( service, job.info.maxThreads );

}


double WorkSpace::objectMetric( boost::asio::io_service& service ) {

    for( const shared_ptr<Object>& o : objects ) {
        if( o ) {
            // async: clear and accumulate P & Q for all objects (blockwise?)
            service.post( [o] {
                o->initPQ();
                o->addAllPQ();
                o->calcMetric();
                //sum.fetch_add(it.second->coefficientMetric());
            } );
            // service.post( std::bind( &Object::getPQ, o.get() ) );
        } else LOG << "WorkSpace::objectMetric()  object is NULL";
    }
    runThreadsAndWait( service, job.info.maxThreads );
//cout << "s" << flush;
    double sum( 0 );
    for( const shared_ptr<Object>& o : objects ) {
        if( o ) {
            sum += o->metric();
        }
    }
    return sum;
}



void WorkSpace::clear( void ) {
    //for( auto & it: objects ) {
    //it->clear();
    //}
    gsl_vector_free(beta);
    gsl_vector_free(grad_beta);
    gsl_vector_free(beta_init);
    
    wavefronts.clear();
}


void WorkSpace::resetAlpha( void ) {
    for( auto & it : wavefronts ) {
        if( it.second ) {
            //it.second->setAlpha( alpha_init );
        } else LOG << "WorkSpace::resetAlpha()  it.second is NULL";
    }
}

/*
PatchResult::Ptr& WorkSpace::getResult( void ) {
    return result;
    //data->images.resize();      // don't send raw data back.
}
*/


void WorkSpace::dump( string tag) {

    for( auto & it : job.getObjects() ) {
        it->dump(tag);
    }

    Ana::write (tag + "_window.f0", window);
    Ana::write (tag + "_noisewindow.f0", noiseWindow);

}



