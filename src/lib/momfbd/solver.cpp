#include "redux/momfbd/solver.hpp"

#include "redux/momfbd/momfbdjob.hpp"

#include "redux/file/fileana.hpp"
#include "redux/image/utils.hpp"
#include "redux/logger.hpp"

#include <functional>
#include <limits>
#include <random>

#include <gsl/gsl_blas.h>

#include <boost/timer/timer.hpp>
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
namespace {

    const std::string thisChannel = "solver";

    grad_t gradientMethods[] = { nullptr,
        std::bind(&SubImage::gradientFiniteDifference, sp::_1, sp::_2, sp::_3),
        std::bind(&SubImage::gradientVogel, sp::_1, sp::_2, sp::_3)
    };

}


Solver::Solver( const MomfbdJob& j ) : objects( j.getObjects() ), job(j), nFreeParameters(0), modeNumbers(nullptr), enabledModes(nullptr), alpha(nullptr), grad_alpha(nullptr) {

}


Solver::~Solver() {
    
    clear();
    
}


void Solver::init( void ) {
    
    clear();

    int nPixels = job.patchSize;
    nModes = job.modeNumbers.size();
    nParameters = job.globalData->constraints.nParameters;
    nFreeParameters = job.globalData->constraints.nFreeParameters;
    
    window.resize( nPixels, nPixels );
    window = 1.0;
    redux::image::apodizeInPlace( window, nPixels / 8);

    int md = std::min( 256, nPixels );                       // new size (maximum=256)
    md -= ( md % 2 );
    noiseWindow.resize(md,md);
    noiseWindow = 1.0;
    redux::image::apodizeInPlace( noiseWindow, md / 8);     // FIXME: old code spacifies md/16, but applies it after "window", so it is actually the product...
    
    enabledModes = new uint16_t[nParameters];
    modeNumbers = new uint16_t[nParameters];
    alpha = new double[nParameters];
    grad_alpha = new double[nParameters];
    
    uint16_t* mPtr = modeNumbers;
    for( size_t n=job.nImages(); n>0; n-- ) {
        memcpy(mPtr,job.modeNumbers.data(),nModes*sizeof(uint16_t));
        mPtr += nModes;
    }

    for( auto & object : job.getObjects() ) {
        object->initProcessing( *this );
    }
    
}


void Solver::getMetric( boost::asio::io_service& service, uint8_t nThreads ) {

    for( auto & object : job.getObjects() ) {
        object->metric();
    }

}


void Solver::resetAll( boost::asio::io_service& service ) {
    for( auto& object: job.getObjects() ) {
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
    for( auto& object: job.getObjects() ) {
        for( auto& ch: object->channels ) {
            for( auto& subimage: ch->subImages ) {
                service.post( std::bind((void(SubImage::*)(string))&SubImage::dump, subimage.get(), tag ) );
            }
        }
    }
}


double Solver::my_f( boost::asio::io_service& service, const gsl_vector* x, void* params ) {

    double* alphaPtr = alpha;
    memset(alphaPtr,0,nParameters*sizeof(double));
    
    job.globalData->constraints.reverse(x->data, alphaPtr);
    
    size_t nModes = job.modeNumbers.size();
    for( const shared_ptr<Object>& o: objects ) {
        o->initPQ();
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                service.post ([&im, alphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->calcPhi( alphaPtr );
                    im->calcOTF();
                    im->addToPQ();
                });
                alphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, job.info.maxThreads );
    
    double sum( 0 );
    for( const shared_ptr<Object>& o : objects ) {
        if( o ) {
            service.post( [&sum,o] {
                o->calcMetric();
                sum += o->metric();
            } );
        } else LOG_ERR << "my_f()  object is NULL";
    }
    runThreadsAndWait( service, job.info.maxThreads );

    return sum;
    
}


void Solver::my_df( boost::asio::io_service& service, const gsl_vector* x, void* params, gsl_vector* df ) {
    
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

    for( const shared_ptr<Object>& o: objects ) {
        o->initPQ();
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                service.post([this, &im, alphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->calcPhi( alphaPtr );
                    im->calcPFOTF();
                    im->addToPQ();
                });
                alphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, job.info.maxThreads );
    
    for( const shared_ptr<Object>& o: objects ) {
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                service.post ([this, &im, gAlphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->calcVogelWeight();
                    for ( uint16_t m=0; m<nModes; ++m ) {
                        if( enabledModes[m] ) {
                            gAlphaPtr[m] = gradientMethod(*im, m, 1E-2 );
                        }
                    }
                });
                gAlphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, job.info.maxThreads );
    
    job.globalData->constraints.apply(grad_alpha,df->data);

}

void Solver::my_fdf( boost::asio::io_service& service, const gsl_vector* x, void* params, double* f, gsl_vector* df ) {
    
    double* alphaPtr = alpha;
    double* gAlphaPtr = grad_alpha;
    
    job.globalData->constraints.reverse(x->data,alphaPtr);
    memset(grad_alpha,0,job.globalData->constraints.nParameters*sizeof(double));
    
   if ( checkAllSmaller( alphaPtr, nParameters, 1E-12 ) ) {
 //       cout << "Forcing graddiff" << endl;
       gradientMethod = gradientMethods[GM_DIFF];
   } else {
 //       cout << "Using gradvogel" << endl;
       gradientMethod = gradientMethods[GM_VOGEL];  //gradientMethods[job.gradientMethod];
   }

    for( const shared_ptr<Object>& o: objects ) {
        o->initPQ();
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                service.post([this, &im, alphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->calcPhi( alphaPtr );
                    im->calcPFOTF();
                    im->addToPQ();
                });
                alphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, job.info.maxThreads );
    
    double sum=0;
    for( const shared_ptr<Object>& o: objects ) {
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                service.post ([this, &im, gAlphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->calcVogelWeight();
                    for ( uint16_t m=0; m<nModes; ++m ) {
                        if( enabledModes[m] ) {
                            gAlphaPtr[m] = gradientMethod(*im, m, 1E-2 );
                        }
                    }
                });
                gAlphaPtr += nModes;
            }
        }
        service.post ([this, &o, &sum] {    // use a lambda to ensure these calls are sequential
            //o->addAllPQ();
            o->calcMetric();
            sum += o->metric();
        });
        service.post(std::bind(&Object::calcMetric, o.get()));
    }
    runThreadsAndWait( service, job.info.maxThreads );
    
    *f = sum;

    job.globalData->constraints.apply(grad_alpha,df->data);

}


void Solver::run( PatchData::Ptr p, boost::asio::io_service& service, uint8_t nThreads ) {

    p->initPatch();
    data = p;

    LOG << "run()  patch#" << data->id << "   index=" << data->index << " region=" << data->roi;

    std::function<gsl_f_t> wrapped_f = std::bind( &Solver::my_f, this, std::ref( service ), sp::_1, sp::_2 );
    std::function<gsl_df_t> wrapped_df = std::bind( &Solver::my_df, this, std::ref( service ), sp::_1, sp::_2, sp::_3 );
    std::function<gsl_fdf_t> wrapped_fdf = std::bind( &Solver::my_fdf, this, std::ref( service ), sp::_1, sp::_2, sp::_3 , sp::_4);

    gsl_fdf_wrapper<decltype( wrapped_f ), decltype( wrapped_df ), decltype( wrapped_fdf )>
    my_func( nFreeParameters, wrapped_f, wrapped_df, wrapped_fdf );
    
    double init_step = 1e-1; //1e-18;   TODO tweak solver parameters
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
            LOG_WARN << "The BFGS-2 method has not been tweaked/tested yet, using BFGS instead.";
            //minimizerType = gsl_multimin_fdfminimizer_vector_bfgs2;
            minimizerType = gsl_multimin_fdfminimizer_vector_bfgs;
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

    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc( minimizerType, nFreeParameters );
    
    /***  Always run first iteration using finite-difference & steepest-descent ***/
    gradientMethod = gradientMethods[GM_DIFF];
    gsl_vector *beta_init = gsl_vector_alloc( nFreeParameters );
    memset(beta_init->data,0,nFreeParameters*sizeof(double));
    
    memset(s->x->data,0,nFreeParameters*sizeof(double));
    memset(s->dx->data,0,nFreeParameters*sizeof(double));
    memset(s->gradient->data,0,nFreeParameters*sizeof(double));
    s->f = 0;
    
    uniform_real_distribution<double> dist(-1E-8, 1E-8);
    default_random_engine re;
  /*  

    //for(uint i=0; i<nFreeParameters; ++i) s->x->data[i] = 1E-6;
    for(uint i=0; i<nParameters; ++i) {
        if((i+1)%51==3)alpha[i] = 1E-6;
        else alpha[i] = 1E-6;
    }
        std::set<uint16_t> activeModes(job.modeNumbers.begin(),job.modeNumbers.end());
        transform(modeNumbers,modeNumbers+nParameters,enabledModes,
                  [&activeModes](const uint16_t& a){ return activeModes.count(a)?a:0; }
                 );
    static int cnt(0);
    cnt++;

    GSL_MULTIMIN_FN_EVAL_F_DF( &my_func, s->x, &s->f, s->gradient);
    cout << "A  Metric: " << s->f << "   nPars" << nParameters << "  nModes=" << nModes << endl;
    Array<double> wrapper(grad_alpha, nParameters/nModes, nModes);
    Ana::write("p"+to_string(cnt)+"_graddiff.f0",wrapper);
  
    dump("p"+to_string(cnt));
    gradientMethod = gradientMethods[GM_VOGEL];
    GSL_MULTIMIN_FN_EVAL_F_DF( &my_func, s->x, &s->f, s->gradient);
    cout << "B  Metric: " << s->f << "   nPars" << nParameters << "  nModes=" << nModes << endl;
    Ana::write("p"+to_string(cnt)+"_gradvogel.f0", wrapper );
    dump("P"+to_string(cnt));
    
    cout << "C  Metric: " << GSL_MULTIMIN_FN_EVAL_F( &my_func, s->x ) << endl;
   // dump("bla");
//exit(0);
    return;*/
    
    boost::timer::auto_cpu_timer timer;

    double previousMetric(0), thisMetric(0), gradNorm(0);
    uint16_t nModeIncrement = job.nModeIncrement;

    size_t failCount(0);
    size_t totalIterations(0);
    size_t maxIterations(1);            // only 1 iteration while increasing modes, job.maxIterations for the last step.
    size_t maxFails(3);         // TODO make into a cfg parameter
    int status(0);

    for( uint16_t modeCount=job.nInitialModes; modeCount; ) {
        
        modeCount = min<uint16_t>(modeCount,job.modeNumbers.size());
        std::set<uint16_t> activeModes(job.modeNumbers.begin(),job.modeNumbers.begin()+modeCount);
        
        transform(modeNumbers,modeNumbers+nParameters,enabledModes,
                  [&activeModes](const uint16_t& a){ return activeModes.count(a)?a:0; }
                 );
        
        if( modeCount == job.modeNumbers.size() ) {                                     // final loop
            modeCount = 0;                      // we use 0 to exit the eternal loop
            init_tol = job.FTOL;
            maxIterations = job.maxIterations;
        } else {
            //nModeIncrement += job.nModeIncrement;                                     // first 5 modes, then 15, then 30 ... 
            modeCount += nModeIncrement;
            failCount = 0;              // we don't worry until we have tried with all modes.
        }
        
        gsl_multimin_fdfminimizer_set( s, static_cast<gsl_multimin_function_fdf*>(&my_func), beta_init, init_step, init_tol );
        
        size_t iter = 0;
        int successCount(0);
        bool done(false);
        
        do {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate( s );
            if( status == GSL_ENOPROG ) {
                LOG_TRACE << "iteration: " << iter << "  GSL reports no progress.";
                failCount++;
                if( iter == 1 ) {           // could not get started, perturb coefficients and try again...
                    for(uint i=0; i<nFreeParameters; ++i) {
                        s->x->data[i] += dist(re);
                    }
                    status = GSL_EFAILED;   // prevent output
                }
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
            if (iter) {
                double relativeChange = 2.0 * fabs(thisMetric-previousMetric) / (fabs(thisMetric)+fabs(previousMetric)+job.EPS);
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
            LOG_DETAIL << boost::format("After %d iteration%s  metric=%g norm(grad)=%g using %d/%d modes.") % totalIterations % (totalIterations>1?"s":" ") %
                s->f % gradNorm % activeModes.size() % job.modeNumbers.size();
        }  
        //  

/*      FIXME the shifting does not work as expected, will look at it soon.
        for( const auto& o: job.objects ) {
            for( const auto& c: o->channels ) {
                for( const auto& im: c->getSubImages() ) {
                    service.post( [&im](){ im->adjustOffset(); } );
                }
            }
        }
        runThreadsAndWait( service, job.info.maxThreads );
*/       
        getAlpha();
        job.globalData->constraints.apply( alpha, beta_init->data );

    }
  
    job.globalData->constraints.reverse(s->x->data,alpha);

    gsl_multimin_fdfminimizer_free( s );
    p->collectResults();

    gsl_vector_free(beta_init);
    
    if( status == GSL_FAILURE ) {
        // TODO throw e.g. "part_fail" or just flag it ??
    }

}



double Solver::objectMetric( boost::asio::io_service& service ) {

    for( const shared_ptr<Object>& o : objects ) {
        if( o ) {
            // async: clear and accumulate P & Q for all objects (blockwise?)
            service.post( [o] {
                o->initPQ();
                o->addAllPQ();
                o->calcMetric();
            } );
        } else LOG << "objectMetric()  object is NULL";
    }
    runThreadsAndWait( service, job.info.maxThreads );
    
    double sum( 0 );
    for( const shared_ptr<Object>& o : objects ) {
        if( o ) {
            sum += o->metric();
        }
    }
    return sum;
}



void Solver::clear( void ) {
    
    delete[] modeNumbers;
    delete[] enabledModes;
    delete[] alpha;
    delete[] grad_alpha;
    modeNumbers = enabledModes = nullptr;
    alpha = grad_alpha = nullptr;

}


void Solver::getAlpha(void) {

    double* alphaPtr = alpha;
    for( const shared_ptr<Object>& o: objects ) {
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                im->getAlphas(alphaPtr);
                alphaPtr += nModes;
            }
        }
    }

}


void Solver::dump( string tag ) {

    for( auto & object : job.getObjects() ) {
        object->dump(tag);
    }

    Ana::write (tag + "_window.f0", window);
    Ana::write (tag + "_noisewindow.f0", noiseWindow);

}



