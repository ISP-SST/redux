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

    const std::string thisChannel = "workspace";

    grad_t gradientMethods[] = { nullptr,
        std::bind(&SubImage::gradFiniteDifference, sp::_1, sp::_2, sp::_3),
        std::bind(&SubImage::gradVogel, sp::_1, sp::_2, sp::_3)
    };

}


WorkSpace::WorkSpace( const MomfbdJob& j ) : objects( j.getObjects() ), job(j), nFreeParameters(0), modeNumbers(nullptr), enabledModes(nullptr), alpha(nullptr), grad_alpha(nullptr) {
    //std::cout << "WorkSpace():  first = (" << d->first.x << "," << d->first.y << ")  last = (" << d->last.x << "," << d->last.y;
    //std::cout << ")   id = (" << d->index.x << "," << d->index.y << ")" << std::endl;

}


WorkSpace::~WorkSpace() {
    clear();
    //  std::cout << "~WorkSpace():  1   wfsz = " << wavefronts.size() << " osz = " << objects.size() << std::endl;

}


void WorkSpace::init( void ) {
    
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

    for( auto & it : job.getObjects() ) {
        it->initProcessing( shared_from_this() );
    }
    
}


void WorkSpace::getMetric( boost::asio::io_service& service, uint8_t nThreads ) {

    for( auto & it : job.getObjects() ) {
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


double WorkSpace::my_f( boost::asio::io_service& service, const gsl_vector* x, void* params ) {

    double* alphaPtr = alpha;
    job.globalData->constraints.reverse(x->data, alphaPtr);
    size_t nModes = job.modeNumbers.size();
    for( const shared_ptr<Object>& o: objects ) {
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                service.post ([this, &im, alphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->setAlphas(job.modeNumbers, alphaPtr);
                    im->update(false);
                });
                alphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, job.info.maxThreads );


    return objectMetric(service);

}


void WorkSpace::my_df( boost::asio::io_service& service, const gsl_vector* x, void* params, gsl_vector* df ) {
    
    double* alphaPtr = alpha;
    double* gAlphaPtr = grad_alpha;
    job.globalData->constraints.reverse(x->data, alphaPtr);
    
    memset(gAlphaPtr,0,nParameters*sizeof(double));
    
    gradientMethod = gradientMethods[job.gradientMethod];

    for( const shared_ptr<Object>& o: objects ) {
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                //im->setAlphas(job.modeNumbers, alphaPtr);
                service.post ([this, &im, alphaPtr, gAlphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->setAlphas(job.modeNumbers, alphaPtr);
                    im->update(true);
                    for ( uint16_t m=0; m<nModes; ++m ) {
                        if( enabledModes[m] ) {
                            gAlphaPtr[m] = gradientMethod(*im, enabledModes[m], 1E-2 );
                        }
                    }
                });
                alphaPtr += nModes;
                gAlphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, job.info.maxThreads );
    
    job.globalData->constraints.apply(grad_alpha,df->data);

}

void WorkSpace::my_fdf( boost::asio::io_service& service, const gsl_vector* x, void* params, double* f, gsl_vector* df ) {
    
    double* alphaPtr = alpha;
    double* gAlphaPtr = grad_alpha;
    
    job.globalData->constraints.reverse(x->data,alphaPtr);
    memset(gAlphaPtr,0,job.globalData->constraints.nParameters*sizeof(double));
    
    gradientMethod = gradientMethods[job.gradientMethod];


    for( const shared_ptr<Object>& o: objects ) {
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                //im->setAlphas(job.modeNumbers, alphaPtr);
                service.post ([this, &im, alphaPtr, gAlphaPtr] {    // use a lambda to ensure these calls are sequential
                    im->setAlphas(job.modeNumbers, alphaPtr);
                    im->update(true);
                    for ( uint16_t m=0; m<nModes; ++m ) {
                        if( enabledModes[m] ) {
                            gAlphaPtr[m] = gradientMethod(*im, enabledModes[m], 1E-2 );
                        }
                    }
                });
                alphaPtr += nModes;
                gAlphaPtr += nModes;
            }
        }
    }
    runThreadsAndWait( service, job.info.maxThreads );
    
    *f = objectMetric(service);
    job.globalData->constraints.apply(grad_alpha,df->data);

}


void WorkSpace::run( PatchData::Ptr p, boost::asio::io_service& service, uint8_t nThreads ) {

    p->initPatch();
    data = p;

    LOG << "WorkSpace::run()  patch#" << data->id << "   index=" << data->index << " region=" << data->roi;

    std::function<gsl_f_t> wrapped_f = std::bind( &WorkSpace::my_f, this, std::ref( service ), sp::_1, sp::_2 );
    std::function<gsl_df_t> wrapped_df = std::bind( &WorkSpace::my_df, this, std::ref( service ), sp::_1, sp::_2, sp::_3 );
    std::function<gsl_fdf_t> wrapped_fdf = std::bind( &WorkSpace::my_fdf, this, std::ref( service ), sp::_1, sp::_2, sp::_3 , sp::_4);

    gsl_fdf_wrapper<decltype( wrapped_f ), decltype( wrapped_df ), decltype( wrapped_fdf )>
    my_func( job.globalData->constraints.nFreeParameters, wrapped_f, wrapped_df, wrapped_fdf );
    
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

    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc( minimizerType, nFreeParameters );
    
    /***  Always run first iteration using finite-difference & steepest-descent ***/
    gradientMethod = gradientMethods[GM_DIFF];
    gsl_vector *beta_init = gsl_vector_alloc( nFreeParameters );
    memset(beta_init->data,0,nFreeParameters*sizeof(double));
    
    boost::timer::auto_cpu_timer timer;

    double previousMetric(0), thisMetric, gradNorm;
    uint16_t nModeIncrement = job.nModeIncrement;
    int dumpCount(0);
    size_t totalIterations(0);

    for( uint16_t modeCount=600/*job.nInitialModes*/; modeCount; ) {
        
        modeCount = min<uint16_t>(modeCount,job.modeNumbers.size());
        std::set<uint16_t> activeModes(job.modeNumbers.begin(),job.modeNumbers.begin()+modeCount);
        
        transform(modeNumbers,modeNumbers+nParameters,enabledModes,
                  [&activeModes](const uint16_t& a){ return activeModes.count(a)?a:0; }
                 );
        
        if( modeCount == job.modeNumbers.size() ) {                                     // final loop
            modeCount = 0;
            init_tol = job.FTOL;
        } else {
            //nModeIncrement += job.nModeIncrement;                                     // first 5 modes, then 15, then 30 ... 
            modeCount += nModeIncrement;
        }
        
        gsl_multimin_fdfminimizer_set( s, static_cast<gsl_multimin_function_fdf*>(&my_func), beta_init, init_step, init_tol );
        
        size_t iter = 0;
        int status, successCount(0);
        bool done(false);
        
        do {
            status = gsl_multimin_fdfminimizer_iterate( s );
            if( status == GSL_ENOPROG ) {
                //LOG_WARN << "iteration: " << iter << "  GSL reports no progress:  quitting loop.";
                break;
            } else if( status ) {
                LOG_WARN << "iteration: " << iter << "  GSL reports status: " << gsl_strerror(status);
            }
            thisMetric = s->f;
            gradNorm = gsl_blas_dnrm2(s->gradient)/(thisMetric*thisMetric);
            if (iter++) {
                double relativeChange = 2.0 * fabs(thisMetric-previousMetric) / (fabs(thisMetric)+fabs(previousMetric)+job.EPS);
                if(relativeChange < job.FTOL) {      // count consecutive "marginal decreases" in metric
                    if( successCount++ >= job.targetIterations ) { // exit after targetIterations consecutive improvements.
                        LOG_DEBUG << alignRight(to_string(iter),7) << ":   metric = " << thisMetric << "   norm(grad) = " << gradNorm
                     << "  rC = " << relativeChange << "  successCount = " << successCount << "   I would like to exit now!.";
                     successCount = 0;
                     done = true;
                    }
                } else {
                    LOG_DEBUG << alignRight(to_string(iter),7) << ":   metric = " << thisMetric << "   norm(grad) = " << gradNorm
                     << "  rC = " << relativeChange << "  successCount = " << successCount;
                    successCount = 0;
                }
            } else LOG_DEBUG << "Initial:   metric = " << thisMetric << "   norm(grad) = " << gradNorm;
            previousMetric = thisMetric;

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
        
        totalIterations += iter;
        
        LOG_DEBUG << boost::format("After %d iterations.  metric=%g  f=%g  norm(grad)=%g  %s") % totalIterations % GSL_MULTIMIN_FN_EVAL_F(&my_func,s->x)
            % s->f % gradNorm % printArray(activeModes,"modes");
        
        /*int imgIndex=0;
        dumpCount++;
        for( const ObjectData& o: data->objects ) {
            service.post(
                [&o,dumpCount](){
                    o.myObject->dump((boost::format("%d_%d")%dumpCount%o.myObject->ID).str());
                } );
            for( const ChannelData& c: o.myObject->channels ) {
                service.post(
                    [&o,&c,dumpCount](){
                        c.myChannel->dump((boost::format("%d_%d_%d")%dumpCount%o.myObject->ID%c.myChannel->ID).str());
                    } );
                for( const shared_ptr<SubImage>& im: c.myChannel->getSubImages() ) {
                    service.post(
                        [&p,&o,&c,&im,imgIndex,dumpCount](){
                            //im->adjustOffset();
                            im->dump((boost::format("%s_%d_%d_%d_%d")%p->index%dumpCount%o.myObject->ID%c.myChannel->ID%imgIndex).str());
                        } );
                    imgIndex++;
                }
            }
        }
        runThreadsAndWait( service, job.info.maxThreads );
       */
        getAlpha();
        job.globalData->constraints.apply(alpha,beta_init->data);
    }
  
    job.globalData->constraints.reverse(s->x->data,alpha);

    gsl_multimin_fdfminimizer_free( s );
    p->collectResults();

    gsl_vector_free(beta_init);
    

    
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
    
    delete[] modeNumbers;
    delete[] enabledModes;
    delete[] alpha;
    delete[] grad_alpha;
    modeNumbers = enabledModes = nullptr;
    alpha = grad_alpha = nullptr;

}


void WorkSpace::getAlpha(void) {

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



