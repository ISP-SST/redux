#include "redux/momfbd/solver.hpp"

#include "redux/momfbd/momfbdjob.hpp"

#include "redux/file/fileana.hpp"
#include "redux/image/utils.hpp"
#include "redux/logging/logger.hpp"
#include "redux/network/host.hpp"

#include <atomic>
#include <functional>
#include <limits>
#include <random>

#include <gsl/gsl_blas.h>

#include <boost/format.hpp>

using namespace redux::file;
using namespace redux::image;
using namespace redux::logging;
using namespace redux::momfbd;
using namespace redux::util;
using namespace redux::util::gsl;
using namespace redux;
using namespace std;
namespace sp=std::placeholders;



//#define DEBUG_
//RDX_DUMP_PATCHDATA

namespace {

    grad_t gradientMethods[] = { nullptr,
        std::bind(&SubImage::gradientFiniteDifference, sp::_1, sp::_2, sp::_3),
        std::bind(&SubImage::gradientVogel, sp::_1, sp::_2, sp::_3)
    };

    struct sync_counter {
        sync_counter() : counter(0) {}
        atomic<int> counter;
        mutex mtx;
        condition_variable cv;
        void reset(void) {
            unique_lock<mutex> lck(mtx);
            counter = 0;
        }
        void wait(void) {
            unique_lock<mutex> lck(mtx);
            while( counter ) cv.wait(lck);
        }
        sync_counter& operator--() {
            unique_lock<mutex> lck(mtx);
            if( --counter == 0 ) cv.notify_all();
            return *this;
        }
        sync_counter& operator++() {
            unique_lock<mutex> lck(mtx);
            counter++;
            return *this;
        }
    };
    
    
}


Solver::Solver( MomfbdJob& j, boost::asio::io_service& s, uint16_t t ) : job(j), myInfo( network::Host::myInfo() ),
    logger(j.logger), objects( j.getObjects() ), service(s), nThreads(t), nFreeParameters(0), nTotalImages(0),
    modeNumbers(nullptr), enabledModes(nullptr),alpha(nullptr), alpha_offset(nullptr), grad_alpha(nullptr),
    tmp_alpha(nullptr), beta(nullptr), grad_beta(nullptr), search_dir(nullptr), tmp_beta(nullptr),
    regAlphaWeights(nullptr) {

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
    
    tmpPhi.resize( nTotalImages, pupilSize, pupilSize );
    tmpPhiGrad.resize( nTotalImages, pupilSize, pupilSize );
    tmpOTF.resize( nTotalImages, 2*pupilSize, 2*pupilSize );

    enabledModes = new uint16_t[nParameters];
    modeNumbers = new uint16_t[nParameters];
    
    alpha = new double[nParameters];
    grad_alpha = new double[nParameters];
    alpha_offset = new double[nParameters];
    tmp_alpha = new double[nParameters];
    regAlphaWeights = new double[nParameters];
    
    beta = new double[nFreeParameters];
    grad_beta = new double[nFreeParameters];
    search_dir = new double[nFreeParameters];
    tmp_beta = new double[nFreeParameters];

    uint16_t* mPtr = modeNumbers;
    for( size_t n=nTotalImages; n; n-- ) {
        memcpy( mPtr, job.modeNumbers.data(), nModes*sizeof(uint16_t));
        mPtr += nModes;
    }

    max_wavelength = 0;
    double* raPtr = regAlphaWeights;
    for( auto & object : objects ) {
        object->initProcessing( *this );
        double scale = job.reg_alpha/object->wavelength;
        scale *= nTotalImages*patchSize*patchSize;      // FIXME the other term in the metric should be normalized instead
        for( size_t i=0; i<object->nImages(); ++i ) {
            transform(object->modes.atm_rms.begin(), object->modes.atm_rms.end(), raPtr,
               [scale](const float& a) { if(a==0) return 0.0; return scale/(a*a); } );
            raPtr += nModes;
        }
        max_wavelength = std::max( max_wavelength, object->wavelength );
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


double Solver::my_f( const gsl_vector* beta, void* params ) {

    applyBeta( beta );
    return metric();
    
}


void Solver::my_df( const gsl_vector* beta, void* params, gsl_vector* df ) {

    applyBeta( beta );
    calcPQ();
    gradient( df );

}


void Solver::my_fdf( const gsl_vector* beta, void* params, double* f, gsl_vector* df ) {
    
    applyBeta( beta );
    *f = metric();
    gradient( df );
    
}


void Solver::my_precalc( const gsl_vector* beta, const gsl_vector* beta_grad ) {
    
    sync_counter sc;
    ++sc; service.post( [this, &sc] { tmpPhiGrad.zero(); --sc; } );
    ++sc; service.post ([this, beta, &sc] { job.globalData->constraints.reverseAndAdd( beta->data, alpha_offset, alpha ); --sc; });
    ++sc; service.post ([this, beta_grad, &sc] { job.globalData->constraints.reverse( beta_grad->data, grad_alpha ); --sc; });
    sc.wait();

    size_t nElements = pupilSize*pupilSize;
    double* alphaPtr = alpha;
    double* phiPtr = tmpPhi.get();
    for( const auto& o: objects ) {
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                ++sc;
                service.post ([&im, alphaPtr, phiPtr, &sc] {    // use a lambda to ensure these calls are sequential
                    im->calcPhi( alphaPtr, phiPtr );
                    --sc;
                });
                alphaPtr += nModes;
                phiPtr += nElements;
            }
        }
    }
    
    alphaPtr = grad_alpha;
    phiPtr = tmpPhiGrad.get();
    for( const auto& o: objects ) {
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                ++sc;
                service.post ([&im, alphaPtr, phiPtr, &sc] {    // use a lambda to ensure these calls are sequential
                    im->addToPhi( alphaPtr, phiPtr );
                    --sc;
                });
                alphaPtr += nModes;
                phiPtr += nElements;
            }
        }
    }
    sc.wait();

}


double Solver::metricAt( double step ) {
    
    sync_counter sc;
    
    size_t pupilSize2 = pupilSize*pupilSize;
    size_t otfSize2 = 4*pupilSize2;
    complex_t* otfPtr = tmpOTF.get();
    double* phiPtr = tmpPhi.get();
    double* phiGradPtr = tmpPhiGrad.get();
    for( const auto& o: objects ) {
        double normalization = sqrt(1.0 / (o->pupil.area*otfSize2));
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                ++sc;
                service.post( [&im, &o, normalization, step, otfPtr, phiPtr, phiGradPtr, &sc] {
                    const double* pupilPtr = o->pupil.get();
                    for( const auto& ind: o->pupil.pupilInOTF ) {
                        otfPtr[ind.second] = polar(pupilPtr[ind.first]*normalization, phiPtr[ind.first]+step*phiGradPtr[ind.first]);
                    }
                    im->tmpOTF.ft(otfPtr);
                    im->tmpOTF.norm();
                    im->tmpOTF.getIFT(im->OTF.get());
                    FourierTransform::reorder( im->OTF.get(), im->OTF.dimSize(0), im->OTF.dimSize(1) );
                    --sc;
                });
                otfPtr += otfSize2;
                phiPtr += pupilSize2;
                phiGradPtr += pupilSize2;
            }
        }
    }
    sc.wait();

    for( const auto& o : objects ) {
        ++sc;
        service.post( [o,&sc] {
            o->initPQ();
            o->addAllPQ();
            o->calcMetric();
            --sc;
        } );
    }
    sc.wait();
    
    double sum(0);
    for( const auto& o : objects ) {
        sum += o->metric();
    }

    return sum;
    
}


void Solver::run( PatchData::Ptr data ) {

    data->initPatch();

    LOG << "Starting patch.  index=" << data->index << "  region=" << data->roi
        << "  nModes=" << nModes << "  nThreads=" << nThreads << ende;
    
    std::function<gsl_f_t> wrapped_f = std::bind( &Solver::my_f, this, sp::_1, sp::_2 );
    std::function<gsl_df_t> wrapped_df = std::bind( &Solver::my_df, this, sp::_1, sp::_2, sp::_3 );
    std::function<gsl_fdf_t> wrapped_fdf = std::bind( &Solver::my_fdf, this, sp::_1, sp::_2, sp::_3 , sp::_4);

    rdx_fdf my_func( nFreeParameters,
                     std::bind( &Solver::my_f, this, sp::_1, sp::_2 ),
                     std::bind( &Solver::my_df, this, sp::_1, sp::_2, sp::_3 ),
                     std::bind( &Solver::my_fdf, this, sp::_1, sp::_2, sp::_3 , sp::_4)
                   );
    
    // TODO: Tweak the new method calls/storage, it is only a marginal speedup at the moment
//      my_func.setPreCalc( std::bind( &Solver::my_precalc, this, sp::_1, sp::_2 ),
//                          std::bind( &Solver::metricAt, this, sp::_1 ) );
    
    double init_step = 1E-4*max_wavelength; //   TODO tweak solver parameters
    double init_tol = 1E-6;
    const gsl_multimin_fdfminimizer_type *minimizerType = gsl_multimin_fdfminimizer_steepest_descent;
    switch(job.getstepMethod) {
        case GSM_BFGS:
            LOG_DEBUG << "Using BFGS solver." << ende;
            LOG_WARN << "The BFGS method has not been properly tweaked/tested yet." << ende;
            minimizerType = gsl_multimin_fdfminimizer_vector_bfgs;
            init_step = 1E8;            // tested for old modes only
            init_tol = 0.01;            // tested for old modes only
            break;
        case GSM_BFGS_inv:
            LOG_DEBUG << "Using BFGS-2 solver." << ende;
            LOG_WARN << "The BFGS-2 method has not been properly tweaked/tested yet." << ende;
            minimizerType = gsl_multimin_fdfminimizer_vector_bfgs2;
            //minimizerType = gsl_multimin_fdfminimizer_vector_bfgs;
            //init_step = 1e-18; //1e-18;
            break;
        case GSM_SDSC:
            LOG_DEBUG << "Using Steepest-Descent solver." << ende;
            minimizerType = gsl_multimin_fdfminimizer_steepest_descent;
            break;
        case GSM_CNJG:
            LOG_DEBUG << "Using Conjugate-Gradient solver." << ende;
            //minimizerType = gsl_multimin_fdfminimizer_conjugate_fr;
            //minimizerType = gsl_multimin_fdfminimizer_conjugate_pr;
            minimizerType = multimin_fdfminimizer_conjugate_rdx;        // Use our own implementation
            break;
        default:
            LOG << "No solver specified, default is Conjugate-Gradient." << ende;
            minimizerType = multimin_fdfminimizer_conjugate_rdx;
    }

    //gradientMethod = gradientMethods[job.gradientMethod];
    gradientMethod = gradientMethods[GM_VOGEL];   // FIXME: old code just uses the cfg-setting the first iteration, then switches to Vogel.
    
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc( minimizerType, nFreeParameters );

    gsl_vector *beta_init = gsl_vector_alloc( nFreeParameters );
    memset(beta_init->data,0,nFreeParameters*sizeof(double));
    memset(alpha_offset,0,nParameters*sizeof(double));
    
    memset(s->x->data,0,nFreeParameters*sizeof(double));
    memset(s->dx->data,0,nFreeParameters*sizeof(double));
    memset(s->gradient->data,0,nFreeParameters*sizeof(double));
    s->f = 0;
        
    double initialMetric = GSL_MULTIMIN_FN_EVAL_F( &my_func, beta_init );
    LOG_TRACE << "Initial metric = " << initialMetric << ende;

    double previousMetric(0);
    double thisMetric(0);
    double gradNorm(0);
    double* alphaPtr;
    double tol = 1.0E-1;
    size_t failCount(0);
    size_t totalIterations(0);
    size_t maxIterations(10);           // fewer iterations while increasing modes, job.maxIterations for the last step.
    size_t maxFails(3);                 // TODO make into a cfg parameter
    int status(0);
    sync_counter sc;

    timer.start();
    logger.flushAll();

    memset( alpha_offset, 0, nParameters*sizeof(double) );
    
    string patchString = alignLeft(to_string(job.info.id) + ":" + to_string(data->id),8);

    for( uint16_t modeCount=job.nInitialModes; modeCount; ) {
        
        modeCount = min<uint16_t>(modeCount,nModes);
        std::set<uint16_t> activeModes(job.modeNumbers.begin(),job.modeNumbers.begin()+modeCount);
        
        string progressString = boost::str( boost::format(" (%03.1f%%)") %
        ((modeCount-min(job.nInitialModes,nModes))*100.0/nModes));
        myInfo.status.statusString = patchString + progressString;
        
        transform(modeNumbers,modeNumbers+nParameters,enabledModes,
                  [&activeModes](const uint16_t& a){ return activeModes.count(a)?a:0; }
                 );
        
        if( modeCount == nModes ) {   // final loop, use proper FTOL and maxIterations
            modeCount = 0;            // we use modeCount=0 to exit the external loop
            tol = job.FTOL;
            maxIterations = job.maxIterations;
        } else {
            modeCount += job.nModeIncrement;
            failCount = 0;           // TODO: fix/test fail-reporting before using
        }

        // The previous position is saved in in alpha_offset, beta only contains the new step, so start at 0 each iteration
        gsl_multimin_fdfminimizer_set( s, static_cast<gsl_multimin_function_fdf*>(&my_func), beta_init, init_step, init_tol );
        
        size_t iter = 0;
        int successCount(0);
        bool done(false);
        
#ifdef DEBUG_
        LOG_DETAIL << "start_metric = " << s->f
                       << printArray(alpha_offset, nParameters, "\nalpha_offset")
                       << printArray(alpha, nParameters, "\nstart_alpha")
                       << printArray(grad_alpha, nParameters, "\nstart_grad") << ende;
        Array<double> gwrapper(grad_alpha, nTotalImages, nModes);
        Ana::write("grad_alpha_"+to_string(totalIterations)+".f0",gwrapper);
        //dump("iter_"+to_string(totalIterations));
        job.globalData->constraints.reverse( s->gradient->data, grad_alpha );
        Ana::write("grad_alphac_"+to_string(totalIterations)+".f0",gwrapper);
        cout << printArray(grad_alpha, nParameters, "\nstart_gradc") << endl;
#endif

        do {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate( s );
            if( status == GSL_ENOPROG ) {
                LOG_TRACE << "iteration: " << iter << "  GSL reports no progress." << ende;
                failCount++;
                break;
            } else if( status ) {
                LOG_WARN << "GSL error in iteration " << iter << ".  type: " << gsl_strerror(status) << ende;
            }
            thisMetric = s->f;
            if( std::isnan(thisMetric) || std::isinf(thisMetric) ) {
                failCount++;
            } 
            if( failCount > maxFails ) {
                LOG_ERR << "Giving up after " << failCount << " failures for patch#" << data->id << " (index=" << data->index << " region=" << data->roi << ")" << ende;
                status = GSL_FAILURE;
                modeCount = 0;                  // exit outer loop.
                done = true;                    // exit inner loop.
                //dump("fail");
                //exit(0);
            }
            gradNorm = gsl_blas_dnrm2(s->gradient)/(thisMetric);
            if (iter>1) {
                double relativeChange = 2.0 * fabs(thisMetric-previousMetric) / (fabs(thisMetric)+fabs(previousMetric)+job.EPS);
                if( relativeChange < tol ) {
                    if( ++successCount >= job.targetIterations ) { // exit after targetIterations consecutive marginal improvements.
                     successCount = 0;
                     done = true;
                    }
                } else {
                    successCount = 0; 
                }
            }

            previousMetric = thisMetric;

            status = gsl_multimin_test_gradient( s->gradient, 1e-9 );

        } while( status == GSL_CONTINUE && (!done || iter < job.minIterations) && iter < maxIterations );

        GSL_MULTIMIN_FN_EVAL_F( &my_func, s->x );
        totalIterations += iter;
        if( status != GSL_FAILURE ) { // bad first iteration -> don't print.
            size_t nActiveModes = activeModes.size();
            if( modeCount ) {
                LOG_DETAIL << "After " << totalIterations << " iteration" << (totalIterations>1?"s":" ")
                 << "  metric=" << thisMetric << "  rmetric=" << (thisMetric/initialMetric)
                 << " norm(grad)=" << (gradNorm/nActiveModes) << " using " << nActiveModes << "/" << nModes << " modes." << ende;
            }
        }

        if( modeCount ) {
            alphaPtr = alpha;
            sc.reset();
            for( const auto& o: job.objects ) {
                for( const auto& c: o->channels ) {
                    for( const auto& im: c->getSubImages() ) {
                        ++sc;
                        service.post( [&im,alphaPtr,&sc](){
                            im->adjustOffset(alphaPtr);
                            --sc;
                        } );
                        alphaPtr += nModes;
                    }
                }
            }
            sc.wait();
            memcpy( alpha_offset, alpha, nParameters*sizeof(double) );
        }

    }       // end for-loop
  
    LOG << "After " << totalIterations << " iterations:  metric=" << thisMetric << "  (relative=" << (thisMetric/initialMetric) << ")" << ende;
    LOG << timer.print() << ende;
    myInfo.status.statusString = patchString + " completed";
    

#ifdef RDX_DUMP_PATCHDATA
    dump( "patch_"+(string)data->index );
#endif
    
    alphaPtr = alpha;
    for( auto& objData: data->objects ) {
        objData.myObject->getResults(objData,alphaPtr);
        alphaPtr += objData.myObject->nObjectImages*nModes;
        objData.setLoaded( objData.size()>6*sizeof(uint64_t) );
        for( auto& ch: objData.channels ) {
            ch.images.clear();         // don't need input data anymore.
        }
#ifdef RDX_DUMP_PATCHDATA
        Ana::write( "patch_"+(string)data->index+"_obj_"+to_string(objData.myObject->ID)+"_result.f0", objData.img );
#endif
    }
    
    data->finalMetric = thisMetric;

    gsl_multimin_fdfminimizer_free( s );
    gsl_vector_free(beta_init);
    
    if( status == GSL_FAILURE ) {
        // TODO throw e.g. "part_fail" or just flag it ??
    }

}


void Solver::applyAlpha( void ) {

    sync_counter sc;
    double* alphaPtr = alpha;
    for( const auto& o: objects ) {
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                ++sc;
                service.post ([&im, alphaPtr,&sc] {    // use a lambda to ensure these calls are sequential
                    im->calcPhi( alphaPtr );
                    im->calcPFOTF();
                    --sc;
                });
                alphaPtr += nModes;
            }
        }
    }
    sc.wait();

}


void Solver::applyBeta( const gsl_vector* beta ) {

    job.globalData->constraints.reverseAndAdd( beta->data, alpha_offset, alpha );
    applyAlpha();

}


void Solver::applyBeta( const gsl_vector* beta, double scale ) {

    job.globalData->constraints.reverse( beta->data, alpha, scale );
    applyAlpha();

}


void Solver::applyBeta( double scale ) {

    job.globalData->constraints.reverse( beta, alpha, scale );
    applyAlpha();
    
}


double Solver::metric(void) {

    sync_counter sc;
    for( const auto& o : objects ) {
        ++sc;
        service.post( [o,&sc] {
            o->initPQ();
            o->addAllPQ();
            o->calcMetric();
            --sc;
        } );
    }
    sc.wait();
    
    double sum(0);
    for( const auto& o : objects ) {
        sum += o->metric();
    }
    
    if( job.reg_alpha > 0 ) {
        for( size_t i=0; i<nParameters; ++i ) {
            sum += 0.5*alpha[i]*alpha[i]*regAlphaWeights[i];
        }
    }

    return sum;
    
}


void Solver::calcPQ(void) {

    sync_counter sc;
    for( const auto& o: objects ) {
        ++sc;
        service.post([&o,&sc] {    // use a lambda to ensure these calls are sequential
                    o->initPQ();
                    o->addAllPQ();
                    --sc;
                });
        
    }
    sc.wait();
    
}
 
 
void Solver::gradient(void) {

    if( job.reg_alpha > 0 ) {
        std::transform( alpha, alpha+nParameters, regAlphaWeights, grad_alpha,
            std::multiplies<double>() );
    } else {
        memset( grad_alpha, 0, nParameters*sizeof(double) );
    }

    double* alphaPtr = alpha;
    double* gAlphaPtr = grad_alpha;
    sync_counter sc;
    for( const auto& o: objects ) {
        double grad_step = job.graddiff_step*o->wavelength;
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                ++sc;
                service.post ([this, &im, grad_step, alphaPtr, gAlphaPtr, &sc] {    // use a lambda to ensure these calls are sequential
                    im->calcVogelWeight();
                    for ( uint16_t m=0; m<nModes; ++m ) {
                        if( enabledModes[m] ) {
                            if( fabs(alphaPtr[m]) > 1E-20 ) gAlphaPtr[m] = gradientMethod(*im, m, grad_step );
                            else gAlphaPtr[m] = gradientMethods[GM_DIFF](*im, m, grad_step );
                        }
                    }
                    --sc;
                });
                alphaPtr += nModes;
                gAlphaPtr += nModes;
            }
        }
    }
    sc.wait();

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
    delete[] alpha_offset;
    delete[] tmp_alpha;
    delete[] regAlphaWeights;
    delete[] beta;
    delete[] grad_beta;
    delete[] search_dir;
    delete[] tmp_beta;
    modeNumbers = enabledModes = nullptr;
    alpha = grad_alpha = alpha_offset = tmp_alpha = nullptr;
    regAlphaWeights = nullptr;
    beta = grad_beta = search_dir = tmp_beta = nullptr;

}


void Solver::dump( string tag ) {

    for( auto & object : objects ) {
        object->dump(tag);
    }

    Ana::write (tag + "_window.f0", window);
    Ana::write (tag + "_noisewindow.f0", noiseWindow);

    Array<double> wrapper(alpha, nTotalImages, nModes);
    Ana::write(tag + "_alpha.f0",wrapper);
    wrapper.wrap(alpha_offset, nTotalImages, nModes);
    Ana::write(tag + "_alpha_offset.f0",wrapper);
    wrapper.wrap(grad_alpha, nTotalImages, nModes);
    Ana::write(tag + "_grad_alpha.f0",wrapper);
     
}



