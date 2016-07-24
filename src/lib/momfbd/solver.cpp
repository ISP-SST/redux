#include "redux/momfbd/solver.hpp"

#include "redux/momfbd/momfbdjob.hpp"

#include "redux/file/fileana.hpp"
#include "redux/image/utils.hpp"
#include "redux/logger.hpp"

#include <functional>
#include <limits>
#include <random>

#include <gsl/gsl_blas.h>

#include <boost/format.hpp>

using namespace redux::file;
using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::util;
using namespace redux::util::gsl;
using namespace redux;
using namespace std;
namespace sp=std::placeholders;

#define lg Logger::mlg
#define logChannel "solver"

//#define DEBUG_

namespace {

    grad_t gradientMethods[] = { nullptr,
        std::bind(&SubImage::gradientFiniteDifference, sp::_1, sp::_2, sp::_3),
        std::bind(&SubImage::gradientVogel, sp::_1, sp::_2, sp::_3)
    };

}


Solver::Solver( const MomfbdJob& j, boost::asio::io_service& s, uint16_t t ) : job(j), objects( j.getObjects() ),
    service(s), nThreads(t), nFreeParameters(0), nTotalImages(0), modeNumbers(nullptr), enabledModes(nullptr),alpha(nullptr), 
    grad_alpha(nullptr), init_alpha(nullptr) {

}


Solver::~Solver() {
    
    clear();
    
}


void Solver::init( void ) {
    
    clear();

    patchSize = job.patchSize;
    pupilSize = job.pupilPixels;
    nModes = job.modeNumbers.size();
    nParameters = job.globalData->constraints.nParameters;
    nFreeParameters = job.globalData->constraints.nFreeParameters;
    nTotalImages = job.nImages();
    
    window.resize( patchSize, patchSize );
    window = 1.0;
    redux::image::apodizeInPlace( window, patchSize / 8);

    int md = std::min<int>( 256, patchSize );                       // new size (maximum=256)
    md -= ( md % 2 );
    noiseWindow.resize(md,md);
    noiseWindow = 1.0;
    redux::image::apodizeInPlace( noiseWindow, md / 16);     // FIXME: old code specifies md/16, but applies it after "window", so it is actually the product...
    
    enabledModes = new uint16_t[nParameters];
    modeNumbers = new uint16_t[nParameters];
    alpha = new double[nParameters];
    grad_alpha = new double[nParameters];
    init_alpha = new double[nParameters];
    
    uint16_t* mPtr = modeNumbers;
    for( size_t n=nTotalImages; n>0; n-- ) {
        memcpy(mPtr,job.modeNumbers.data(),nModes*sizeof(uint16_t));
        mPtr += nModes;
    }

    max_mode_norm = 0;
//     Array<double> wrapper(init_alpha, nParameters/nModes, nModes);
    double* gInitPtr = init_alpha;
    memset(gInitPtr,0,nParameters*sizeof(double));
    for( auto & object : objects ) {
        object->initProcessing( *this );
        max_mode_norm = std::max( max_mode_norm, *max_element(object->modes.norms.begin(), object->modes.norms.end()));
        if(false)
        for( size_t im=object->nImages(); im; --im ) {
            for ( uint16_t m=0; m<nModes; ++m ) {
                if( m != object->modes.tiltMode.x && m != object->modes.tiltMode.y ) {
                    //gInitPtr[m] = 0.20/object->modes.norms[m];        // completely ad hoc value for starting point
                    gInitPtr[m] = 0.1*object->modes.atm_rms[m]/object->modes.norms[m];        // starting point, using atm_rms as weight
                }
            }
            gInitPtr += nModes;
        }
    }
    
}


void Solver::getMetric( boost::asio::io_service& service, uint8_t nThreads ) {

    for( auto & object : objects ) {
        object->metric();
    }

}


void Solver::reset( void ) {
    for( auto& object: objects ) {
        for( auto& ch: object->channels ) {
            for( auto& subimage: ch->subImages ) {
                service.post( std::bind((void(SubImage::*)(void))&SubImage::resetPhi, subimage.get() ) );
            }
        }
        service.post( std::bind(&Object::initPQ, object.get() ) );
        service.post( std::bind(&Object::addAllPQ, object.get() ) );
    }
}
void Solver::dumpImages( boost::asio::io_service& service, string tag ) {
    for( auto& object: objects ) {
        for( auto& ch: object->channels ) {
            for( auto& subimage: ch->subImages ) {
                service.post( std::bind((void(SubImage::*)(string))&SubImage::dump, subimage.get(), tag ) );
            }
        }
    }
}


double Solver::my_f( const gsl_vector* x, void* params ) {

    double* alphaPtr = alpha;
//    memset(alphaPtr,0,nParameters*sizeof(double));
    
    job.globalData->constraints.reverse(x->data, alphaPtr);
   
    for( const auto& o: objects ) {
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                service.post ([&im, alphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->calcPhi( alphaPtr );
                    im->calcOTF();
                });
                alphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, nThreads );
    
    for( const auto& o: objects ) {
        service.post ([&o] {    // use a lambda to ensure these calls are sequential
            o->initPQ();
            o->addAllPQ();
            o->calcMetric();
        });
    }
    runThreadsAndWait( service, nThreads );

    double sum(0);
    for( const auto& o : objects ) {
        sum += o->metric();
    }
    return sum; // /nTotalImages;
    
}


void Solver::my_df( const gsl_vector* x, void* params, gsl_vector* df ) {
    
    double* alphaPtr = alpha;
    double* gAlphaPtr = grad_alpha;
    
    job.globalData->constraints.reverse(x->data, alphaPtr);
    memset(gAlphaPtr,0,nParameters*sizeof(double));
    
    if ( checkAllSmaller( alphaPtr, nParameters, 1E-12 ) ) {
        gradientMethod = gradientMethods[GM_DIFF];
//        cout << "Forcing graddiff" << endl;
    } else {
//        cout << "Using gradvogel" << endl;
        gradientMethod = gradientMethods[GM_VOGEL];  //gradientMethods[job.gradientMethod];
    }

    for( const auto& o: objects ) {
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                service.post([this, &im, alphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->calcPhi( alphaPtr );
                    im->calcPFOTF();
                });
                alphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, nThreads );
    for( const auto& o: objects ) {
        service.post([&o] {    // use a lambda to ensure these calls are sequential
                    o->initPQ();
                    o->addAllPQ();
                });
        
    }
    runThreadsAndWait( service, nThreads );
    
    for( const auto& o: objects ) {
        const std::vector<float>& norms = o->modes.norms;
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                service.post ([this, &norms, &im, gAlphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->calcVogelWeight();
                    for ( uint16_t m=0; m<nModes; ++m ) {
                        if( enabledModes[m] ) {
                            gAlphaPtr[m] = gradientMethod(*im, m, job.graddiff_step ); // / norms[m]; ///nTotalImages;
                        }
                    }
                });
                gAlphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, nThreads );
    
    job.globalData->constraints.apply(grad_alpha,df->data);
    
    gradientMethod = gradientMethods[GM_VOGEL];     // FIXME: this is a hack to mimic the behaviour of the old code

}

void Solver::my_fdf( const gsl_vector* x, void* params, double* f, gsl_vector* df ) {
    
    double* alphaPtr = alpha;
    double* gAlphaPtr = grad_alpha;
    
    job.globalData->constraints.reverse(x->data,alphaPtr);
    memset(grad_alpha,0,nParameters*sizeof(double));
    
   if ( checkAllSmaller( alphaPtr, nParameters, 1E-12 ) ) {
 //       cout << "Forcing graddiff" << endl;
       gradientMethod = gradientMethods[GM_DIFF];
   } else {
 //       cout << "Using gradvogel" << endl;
       gradientMethod = gradientMethods[GM_VOGEL];  //gradientMethods[job.gradientMethod];
   }

    for( const auto& o: objects ) {
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                service.post([this, &im, alphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->calcPhi( alphaPtr );
                    im->calcPFOTF();
                });
                alphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, nThreads );
    for( const auto& o: objects ) {
        service.post([&o] {    // use a lambda to ensure these calls are sequential
                    o->initPQ();
                    o->addAllPQ();
                    o->calcMetric();
                });
        
    }
    runThreadsAndWait( service, nThreads );
    
    double sum=0;
    for( const auto& o: objects ) {
        const std::vector<float>& norms = o->modes.norms;
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                service.post ([this, &norms, &im, gAlphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->calcVogelWeight();
                    for ( uint16_t m=0; m<nModes; ++m ) {
                        if( enabledModes[m] ) {
                            gAlphaPtr[m] = gradientMethod(*im, m, job.graddiff_step ); // / norms[m]; ///nTotalImages;
                        }
                    }
                });
                gAlphaPtr += nModes;
            }
        }
        sum += o->metric();
    }
    runThreadsAndWait( service, nThreads );
    
    *f = sum; // /nTotalImages;

    job.globalData->constraints.apply(grad_alpha,df->data);

    gradientMethod = gradientMethods[GM_VOGEL];     // FIXME: this is a hack to mimic the behaviour of the old code
    
}


void Solver::run( PatchData::Ptr data ) {
    
    data->initPatch();

    LOG << "run()  patch#" << data->id << "   index=" << data->index << " region=" << data->roi << "   nThreads=" << nThreads;

    std::function<gsl_f_t> wrapped_f = std::bind( &Solver::my_f, this, sp::_1, sp::_2 );
    std::function<gsl_df_t> wrapped_df = std::bind( &Solver::my_df, this, sp::_1, sp::_2, sp::_3 );
    std::function<gsl_fdf_t> wrapped_fdf = std::bind( &Solver::my_fdf, this, sp::_1, sp::_2, sp::_3 , sp::_4);

    gsl_fdf_wrapper<decltype( wrapped_f ), decltype( wrapped_df ), decltype( wrapped_fdf )>
    my_func( nFreeParameters, wrapped_f, wrapped_df, wrapped_fdf );
    
    double init_step = 1E8; //10; //1.0/sqrt(max_mode_norm); //1000; //1e-18;   TODO tweak solver parameters
    double init_tol = 1E-1; //0.01; //1E-12; //0.1;
    const gsl_multimin_fdfminimizer_type *minimizerType = gsl_multimin_fdfminimizer_steepest_descent;
    switch(job.getstepMethod) {
        case GSM_BFGS:
            LOG_DETAIL << "Using BFGS solver.";
            LOG_WARN << "The BFGS method has not been properly tweaked/tested yet.";
            minimizerType = gsl_multimin_fdfminimizer_vector_bfgs;
            init_step = 1E8;            // tested for old modes only
            init_tol = 0.01;            // tested for old modes only
            break;
        case GSM_BFGS_inv:
            LOG_DETAIL << "Using BFGS-2 solver.";
            LOG_WARN << "The BFGS-2 method has not been properly tweaked/tested yet.";
            minimizerType = gsl_multimin_fdfminimizer_vector_bfgs2;
            //minimizerType = gsl_multimin_fdfminimizer_vector_bfgs;
            //init_step = 1e-18; //1e-18;
            break;
        case GSM_SDSC:
            LOG_DETAIL << "Using Steepest-Descent solver.";
            minimizerType = gsl_multimin_fdfminimizer_steepest_descent;
            break;
        case GSM_CNJG:
            LOG_DETAIL << "Using Conjugate-Gradient solver.";
            //minimizerType = gsl_multimin_fdfminimizer_conjugate_fr;
            minimizerType = gsl_multimin_fdfminimizer_conjugate_pr;
            //init_step = 1E8;            // tested for _pr, needs tweaking
            init_step = 10*max_mode_norm;            // tested for _pr, needs tweaking
            init_tol = 0.01;            // tested for _pr, needs tweaking
            //minimizerType = rdx_multimin_fdfminimizer_conjugate_pr;
            //init_tol = 0; //1e-4;
            break;
        default:
            LOG << "No solver specified, default is Conjugate-Gradient.";
            minimizerType = gsl_multimin_fdfminimizer_conjugate_pr;
            init_step = 1E8;            // tested for _pr, needs tweaking
            init_tol = 0.01;            // tested for _pr, needs tweaking
    }
    
    //init_step = 10;
    //init_tol = 0.01;

    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc( minimizerType, nFreeParameters );
    
    /***  Always run first iteration using finite-difference & steepest-descent ***/
    //gradientMethod = gradientMethods[GM_DIFF];
    gradientMethod = gradientMethods[GM_VOGEL];
    gsl_vector *beta_init = gsl_vector_alloc( nFreeParameters );
    memset(beta_init->data,0,nFreeParameters*sizeof(double));
    
    memset(s->x->data,0,nFreeParameters*sizeof(double));
    memset(s->dx->data,0,nFreeParameters*sizeof(double));
    memset(s->gradient->data,0,nFreeParameters*sizeof(double));
    s->f = 0;

    job.globalData->constraints.apply( init_alpha, beta_init->data );

    init_step /= max_mode_norm;

    double initialMetric = GSL_MULTIMIN_FN_EVAL_F( &my_func, beta_init );

    double previousMetric(0);
    double thisMetric(0);
    double gradNorm(0);
    uint16_t nModeIncrement = job.nModeIncrement;
    double* alphaPtr;

    size_t failCount(0);
    size_t totalIterations(0);
    size_t maxIterations(5);            // only 1 iteration while increasing modes, job.maxIterations for the last step.
    size_t maxFails(3);         // TODO make into a cfg parameter
    int status(0);
    uint16_t firstMode(0);

    timer.start();
    
    for( uint16_t modeCount=job.nInitialModes; modeCount; ) {
        
        modeCount = min<uint16_t>(modeCount,nModes);
        std::set<uint16_t> activeModes(job.modeNumbers.begin()+firstMode,job.modeNumbers.begin()+modeCount);
        
        transform(modeNumbers,modeNumbers+nParameters,enabledModes,
                  [&activeModes](const uint16_t& a){ return activeModes.count(a)?a:0; }
                 );
        
        if( false && (modeCount == nModes && firstMode) ) {
            firstMode = 0;                      // make a final run with all modes enabled
        } else if( modeCount == nModes ) {                                     // final loop
            modeCount = 0;                      // we use modeCount=0 to exit the external loop
            //init_tol = job.FTOL;
            maxIterations = job.maxIterations;
        } else {
            //nModeIncrement += job.nModeIncrement;                                     // first 5 modes, then 15, then 30 ... 
            //firstMode = modeCount;
            modeCount += nModeIncrement;
            failCount = 0;              // we don't worry until we have tried with all modes.
        }
        
        //cout << "Step " << totalIterations << "  tol=" << init_tol << "  step=" << init_step << "  maxIter=" << maxIterations << endl;
        gsl_multimin_fdfminimizer_set( s, static_cast<gsl_multimin_function_fdf*>(&my_func), beta_init, init_step, init_tol );
        
        size_t iter = 0;
        int successCount(0);
        bool done(false);
        
#ifdef DEBUG_
        LOG_DETAIL << "start_metric = " << s->f
                       << printArray(alpha, nParameters, "\nstart_alpha")
                       << printArray(grad_alpha, nParameters, "\nstart_grad");
            //dump("iter_"+to_string(totalIterations));
#endif

        gradientMethod = gradientMethods[GM_DIFF];
        do {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate( s );

            if( status == GSL_ENOPROG ) {
                LOG_TRACE << "iteration: " << iter << "  GSL reports no progress.";
                failCount++;
                break;
            } else if( status ) {
                LOG_WARN << "GSL error in iteration " << iter << ".  type: " << gsl_strerror(status);
            }
            thisMetric = s->f;
            if( std::isnan(thisMetric) || std::isinf(thisMetric) ) {
                failCount++;
            } 
            if( failCount > maxFails ) {
                LOG_ERR << "Giving up after " << failCount << " failures for patch#" << data->id << " (index=" << data->index << " region=" << data->roi << ")";
                status = GSL_FAILURE;
                modeCount = 0;                  // exit outer loop.
                done = true;                    // exit inner loop.
                dump("fail");
                exit(0);
            }
            gradNorm = gsl_blas_dnrm2(s->gradient)/(thisMetric);
            if (iter>1) {
                double relativeMetric = thisMetric/initialMetric; 
                double relativeChange = 2.0 * fabs(thisMetric-previousMetric) / (fabs(thisMetric)+fabs(previousMetric)+job.EPS);
//             LOG_DETAIL << "Iter: " << (totalIterations+iter) << "   thisMetric = " << thisMetric
//                         << "   relativeChange = " << relativeChange << "   relativeMetric = " << relativeMetric << "   cnt = " << successCount;
                if(relativeChange < job.FTOL) {      // count consecutive "marginal decreases" in metric
                    if( successCount++ >= job.targetIterations ) { // exit after targetIterations consecutive improvements.
                     successCount = 0;
                     done = true;
                    }
                } else {    // reset counter
                    successCount = 0;
                }
            } //else LOG_DEBUG << "Initial:   metric = " << thisMetric << "   norm(grad) = " << gradNorm;
            previousMetric = thisMetric;

            status = gsl_multimin_test_gradient( s->gradient, 1e-9 );

        } while( status == GSL_CONTINUE && (!done || iter < job.minIterations) && iter < maxIterations );

        GSL_MULTIMIN_FN_EVAL_F( &my_func, s->x );
        totalIterations += iter;
        if( status != GSL_FAILURE ) { // bad first iteration -> don't print.
            size_t nActiveModes = activeModes.size();
            LOG_DETAIL << "After " << totalIterations << " iteration" << (totalIterations>1?"s":" ") << "  metric=" << thisMetric
             << " norm(grad)=" << (gradNorm/nActiveModes) << " using " << nActiveModes << "/" << nModes << " modes.";
//            LOG_DETAIL << boost::format("After %d iteration%s  metric=%g norm(grad)=%g using %d/%d modes.") % totalIterations % (totalIterations>1?"s":" ") %
//                s->f % gradNorm % nActiveModes % nModes;
//             gradientMethod = gradientMethods[GM_VOGEL];
        }  
        //  

        alphaPtr = alpha;
        for( const auto& o: job.objects ) {
            for( const auto& c: o->channels ) {
                for( const auto& im: c->getSubImages() ) {
                    service.post( [&im,alphaPtr](){ im->adjustOffset(alphaPtr); } );
                    alphaPtr += nModes;
                }
            }
        }
        runThreadsAndWait( service, nThreads );

        job.globalData->constraints.apply( alpha, beta_init->data );
#ifdef DEBUG_
            LOG_DETAIL << printArray(alpha, nParameters, "\nalpha")
                       << printArray(grad_alpha, nParameters, "\ngrad_alpha")
                       << printArray(beta_init->data, nFreeParameters, "\nbeta")
                       << printArray(s->gradient->data, nFreeParameters, "\ngrad_beta");
            dump("iter_"+to_string(totalIterations));
#endif
    }
    
    LOG << "Elapsed: " << boost::timer::format(timer.elapsed());
    job.globalData->constraints.reverse( s->x->data, alpha );
    
    alphaPtr = alpha;
    for( auto& obj: data->objects ) {
        obj.myObject->getResults(obj,alphaPtr);
        alphaPtr += obj.myObject->nObjectImages*nModes;
        obj.setLoaded( obj.size()>6*sizeof(uint64_t) );
        for( auto& ch: obj.channels ) {
            ch.images.clear();         // don't need input data anymore.
        }
    }
    
    data->finalMetric = thisMetric;

    gsl_multimin_fdfminimizer_free( s );
    gsl_vector_free(beta_init);
    
    if( status == GSL_FAILURE ) {
        // TODO throw e.g. "part_fail" or just flag it ??
    }

}


void Solver::applyAlpha( void ) {

    double* alphaPtr = alpha;
    for( const auto& o: objects ) {
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                service.post ([&im, alphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->calcPhi( alphaPtr );
                    im->calcPFOTF();
                });
                alphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, nThreads );
    
}


void Solver::applyBeta( const gsl_vector* beta ) {

    job.globalData->constraints.reverse( beta->data, alpha );
    applyAlpha();

}


void Solver::applyBeta( const gsl_vector* beta, double scale ) {

    job.globalData->constraints.reverse( beta->data, alpha, scale );
    applyAlpha();

}




double Solver::metric(void) {

    for( const auto& o : objects ) {
        service.post( [o] {
            o->initPQ();
            o->addAllPQ();
            o->calcMetric();
        } );
    }
    runThreadsAndWait( service, nThreads );
    
    double sum(0);
    for( const auto& o : objects ) {
        sum += o->metric();
    }
    return sum;
}



void Solver::gradient(void) {

    double* gAlphaPtr = grad_alpha;
    memset(gAlphaPtr,0,nParameters*sizeof(double));

    for( const auto& o: objects ) {
        service.post([&o] {    // use a lambda to ensure these calls are sequential
                    o->initPQ();
                    o->addAllPQ();
                });
        
    }
    runThreadsAndWait( service, nThreads );
    
    for( const auto& o: objects ) {
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                service.post ([this, &im, gAlphaPtr] {    // use a lambda to ensure these calls are sequential
                    /*if ( job.gradientMethod == GM_VOGEL )*/ im->calcVogelWeight();
                    for ( uint16_t m=0; m<nModes; ++m ) {
                        if( enabledModes[m] ) {
                            gAlphaPtr[m] = gradientMethod(*im, m, job.graddiff_step );
                        }
                    }
                });
                gAlphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, nThreads );
    
}
 
 
void Solver::gradient( gsl_vector* grad ) {
    
    gradient();
    job.globalData->constraints.apply( grad_alpha, grad->data );
    
}



void Solver::clear( void ) {
    
    delete[] modeNumbers;
    delete[] enabledModes;
    delete[] alpha;
    delete[] grad_alpha;
    delete[] init_alpha;
    modeNumbers = enabledModes = nullptr;
    alpha = grad_alpha = nullptr;

}


void Solver::getAlpha(void) {

    double* alphaPtr = alpha;
    for( const auto& o: objects ) {
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                im->getAlphas(alphaPtr);
                alphaPtr += nModes;
            }
        }
    }

}


void Solver::dump( string tag ) {

    for( auto & object : objects ) {
        object->dump(tag);
    }

    Ana::write (tag + "_window.f0", window);
    Ana::write (tag + "_noisewindow.f0", noiseWindow);

}



