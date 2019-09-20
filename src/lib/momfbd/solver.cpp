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


size_t redux::momfbd::thread::TmpStorage::currentSize = 0;
uint16_t redux::momfbd::thread::TmpStorage::patchSize = 0;
uint16_t redux::momfbd::thread::TmpStorage::pupilSize = 0;

//#define DEBUG_
//#define MASSIVE_DEBUG_
//#define RDX_DUMP_PATCHDATA

namespace {

    grad_t gradientMethods[] = { nullptr,
        std::bind(&SubImage::gradientFiniteDifference, sp::_1, sp::_2),
        std::bind(&SubImage::gradientVogel, sp::_1, sp::_2)
    };
    
    thread_local shared_ptr<redux::momfbd::thread::TmpStorage> TS( new redux::momfbd::thread::TmpStorage() );
    
}

redux::momfbd::thread::TmpStorage* Solver::tmp(void) {
    return TS.get();
}

Solver::Solver( MomfbdJob& j, boost::asio::io_service& s, uint16_t t ) : job(j), myInfo( network::Host::myInfo() ),
    logger(j.logger), objects( j.getObjects() ), service(s), maxThreads(t), nFreeParameters(0), nTotalImages(0),
    beta(nullptr), grad_beta(nullptr), search_dir(nullptr), tmp_beta(nullptr),
    regAlphaWeights(nullptr), patchSize2(0), pupilSize2(0), nTotalPixels(0), otfSize(0), otfSize2(0) {

    init();

}


Solver::~Solver() {
    
    clear();
    
}


void Solver::init( void ) {
    
    LOG_TRACE << "Initializing Solver." << ende;

    clear();
    job.globalData->constraints.makeRowsCols();

    patchSize = job.patchSize;
    pupilSize = job.pupilPixels;
    patchSize2 = patchSize*patchSize;
    pupilSize2 = pupilSize*pupilSize;
    otfSize = 2*pupilSize;
    otfSize2 = otfSize*otfSize;
    nModes = job.modeNumbers.size();
    nParameters = job.globalData->constraints.nParameters;
    nFreeParameters = job.globalData->constraints.nFreeParameters;
    nTotalImages = job.nImages();
    nTotalPixels = nTotalImages*patchSize2;
    
    if( job.globalData->constraints.ns_rows.size() != nParameters ||
        job.globalData->constraints.ns_cols.size() != nFreeParameters ||
        job.globalData->constraints.ns_entries.empty() ) {
        throw std::logic_error("The constraint-set is corrupt, the nullspace has the wrong number of entries.");
    }

    window.resize( patchSize, patchSize );
    window = 1.0;
    redux::image::apodizeInPlace( window, patchSize / 8);

    noiseWindow.resize(patchSize,patchSize);
    noiseWindow = 1.0;
    redux::image::apodizeInPlace( noiseWindow, patchSize / 16);     // FIXME: old code specifies md/16, but applies it after "window", so it is actually the product...
    
    tmpPhi.resize( nTotalImages, pupilSize, pupilSize );
    tmpPhiGrad.resize( nTotalImages, pupilSize, pupilSize );
    tmpOTF.resize( nTotalImages, otfSize, otfSize );

    enabledModes = rdx_get_shared<bool>(nModes);
    
    alpha = rdx_get_shared<double>(nParameters);
    alpha_offset = rdx_get_shared<double>(nParameters);
    grad_alpha = rdx_get_shared<double>(nParameters);

    regAlphaWeights = rdx_get_shared<double>(nParameters);
    beta = rdx_get_shared<double>(nFreeParameters);
    grad_beta = rdx_get_shared<double>(nFreeParameters);
    search_dir = rdx_get_shared<double>(nFreeParameters);
    tmp_beta = rdx_get_shared<double>(nFreeParameters);

    max_wavelength = 0;
    double* raPtr = regAlphaWeights.get();
    for( auto & object : objects ) {
        object->initProcessing( *this );
        double scale = job.reg_alpha/object->wavelength;
        scale *= nTotalPixels;      // FIXME the other term in the metric should be normalized instead
        for( size_t i=0; i<object->nImages(); ++i ) {
            transform(object->modes->atm_rms.begin(), object->modes->atm_rms.end(), raPtr,
               [scale](const float& a) { if(a==0) return 0.0; return scale/(a*a); } );
            raPtr += nModes;
        }
        max_wavelength = std::max( max_wavelength, object->wavelength );
    }

    size_t offset(0);
    for( auto& it: wavefronts ) {
        it.second->setData( alpha.get()+offset, grad_alpha.get()+offset, nModes );
        offset += nModes;
    }
    
    thread::TmpStorage::setSize( patchSize, pupilSize );
    tmp()->init();     // temp-storage for the main thread.
    set<std::thread::id> initDone;
    while(true) {
        for( uint16_t i=0; i<2*maxThreads; ++i ) {
            service.post([&](){
                unique_lock<mutex> lock(mtx);
                if( initDone.count(std::this_thread::get_id()) ) return;
                tmp()->init();
                initDone.insert(std::this_thread::get_id());
            });
        }
        std::this_thread::sleep_for( std::chrono::seconds(1) );
        unique_lock<mutex> lock(mtx);
        if( initDone.size() == maxThreads ) break;
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


double Solver::my_f( const gsl_vector* b, void* params ) {

    applyBeta( b );
    
    return metric();

}


void Solver::my_df( const gsl_vector* b, void* params, gsl_vector* df ) {

    applyBeta( b );
    calcPQ();
    gradient( df );

}


void Solver::my_fdf( const gsl_vector* b, void* params, double* f, gsl_vector* df ) {
    
    applyBeta( b );
    *f = metric();
    gradient( df );
    
}


void Solver::my_precalc( const gsl_vector* b, const gsl_vector* b_grad ) {
    
    progWatch.set(3);
    service.post( [this] { tmpPhiGrad.zero(); ++progWatch; } );
    service.post ([this, b] { job.globalData->constraints.reverseAndAdd( b->data, alpha_offset.get(), alpha.get() ); ++progWatch; });
    service.post ([this, b_grad] { job.globalData->constraints.reverse( b_grad->data, grad_alpha.get() ); ++progWatch; });
    progWatch.wait();

    size_t nElements = pupilSize*pupilSize;
    double* alphaPtr = alpha.get();
    double* phiPtr = tmpPhi.get();
    progWatch.set( 2*nTotalImages );
    for( const auto& o: objects ) {
        //if( o->weight > 0 ) {
            for( const auto& c: o->getChannels() ) {
                for( const auto& im: c->getSubImages() ) {
                    service.post ([this,&im, alphaPtr, phiPtr] {    // use a lambda to ensure these calls are sequential
                        im->calcPhi( alphaPtr, phiPtr );
                        ++progWatch;
                    });
                    alphaPtr += nModes;
                    phiPtr += nElements;
                }
            }
        //}
    }
    
    alphaPtr = grad_alpha.get();
    phiPtr = tmpPhiGrad.get();
    for( const auto& o: objects ) {
        //if( o->weight > 0 ) {
            for( const auto& c: o->getChannels() ) {
                for( const auto& im: c->getSubImages() ) {
                    service.post ([this,&im, alphaPtr, phiPtr] {    // use a lambda to ensure these calls are sequential
                        im->addToPhi( alphaPtr, phiPtr );
                        ++progWatch;
                    });
                    alphaPtr += nModes;
                    phiPtr += nElements;
                }
            }
        //}
    }
    progWatch.wait();
    
}


double Solver::metricAt( double step ) {
    
    complex_t* otfPtr = tmpOTF.get();
    double* phiPtr = tmpPhi.get();
    double* phiGradPtr = tmpPhiGrad.get();
    progWatch.set( nTotalImages );
    for( const auto& o: objects ) {
        //if( o->weight > 0 ) {
            double normalization = sqrt(1.0 / (o->pupil->area*otfSize2));
            for( const auto& c: o->getChannels() ) {
                for( const auto& im RDX_UNUSED: c->getSubImages() ) {
                    service.post( [this, &o, normalization, step, otfPtr, phiPtr, phiGradPtr] {
                        const double* pupilPtr = o->pupil->get();
                        for( const auto& ind: o->pupil->pupilInOTF ) {
                            otfPtr[ind.second] = polar(pupilPtr[ind.first]*normalization, phiPtr[ind.first]+step*phiGradPtr[ind.first]);
                        }
                        //im->tmpOTF.ft(otfPtr); FIXME: broke with thread_local
                        //im->tmpOTF.norm();
                        //im->tmpOTF.getIFT(im->OTF.get());
                        //FourierTransform::reorder( im->OTF.get(), im->OTF.dimSize(0), im->OTF.dimSize(1) );
                        ++progWatch;
                    });
                    otfPtr += otfSize2;
                    phiPtr += pupilSize2;
                    phiGradPtr += pupilSize2;
                }
            }
        //}
    }
    progWatch.wait();

    return metric();
    
}


void Solver::run( PatchData::Ptr data ) {
    
    zeroAlphas();

    LOG << "Starting patch.  index=" << data->index << "  region=" << data->roi
        << "  nModes=" << nModes << "  nThreads=" << maxThreads << ende;

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
    gradientMethod = gradientMethods[GM_DIFF];   // FIXME: old code just uses the cfg-setting the first iteration, then switches to Vogel.
    
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc( minimizerType, nFreeParameters );

    gsl_vector *beta_init = gsl_vector_alloc( nFreeParameters );
    memset( beta_init->data, 0, nFreeParameters*sizeof(double) );
    memset( alpha_offset.get(), 0, nParameters*sizeof(double) );
    memset(s->x->data,0,nFreeParameters*sizeof(double));
    memset(s->dx->data,0,nFreeParameters*sizeof(double));
    memset(s->gradient->data,0,nFreeParameters*sizeof(double));
    memset( enabledModes.get(), 0, nModes );
    s->f = 0;
    

    double previousMetric(0);
    double thisMetric(0);
    double gradNorm(0);
//    double* alphaPtr = alpha.get();
    double tol = 1.0E-1;
    size_t failCount(0);
    size_t totalIterations(0);
    size_t maxIterations(10);           // fewer iterations while increasing modes, job.maxIterations for the last step.
    size_t maxFails(3);                 // TODO make into a cfg parameter
    int status(0);

    timer.start();
    logger.flushAll();
    loadInit( data, alpha_offset.get() );
    shiftAndInit( alpha_offset.get(), true );     // force initialization
    
    data->metrics.clear();
    data->metrics.reserve(10000);     // should be more than enough
    double initialMetric = GSL_MULTIMIN_FN_EVAL_F( &my_func, beta_init );
    data->metrics.push_back(initialMetric);
    LOG << "Patch" << (string)data->index << ":  Initial metric = " << initialMetric << ende;

    double gradScale = 1.0/(nTotalPixels*nTotalPixels);
    
    string patchString = alignLeft(to_string(job.info.id) + ":" + to_string(data->id),8);


    for( uint16_t modeCount=job.nInitialModes; modeCount; ) {
        
        modeCount = min<uint16_t>(modeCount,nModes);
        std::set<uint16_t> activeModes(job.modeNumbers.begin(),job.modeNumbers.begin()+modeCount);
        
        string progressString = boost::str( boost::format(" (%03.1f%%)") %
                                ((modeCount-min(job.nInitialModes,nModes))*100.0/nModes));
        myInfo.status.statusString = patchString + progressString;
        
        std::fill( enabledModes.get(), enabledModes.get()+modeCount, true );

        if( modeCount == nModes ) {   // final loop, use proper FTOL and maxIterations
            modeCount = 0;            // we use modeCount=0 to exit the external loop
            tol = job.FTOL;
            maxIterations = job.maxIterations;
        } else {
            modeCount += job.nModeIncrement;
            //failCount = 0;           // TODO: fix/test fail-reporting before using
        }

        // The previous position is saved in in alpha_offset, beta only contains the new step, so start at 0 each iteration
        gsl_multimin_fdfminimizer_set( s, static_cast<gsl_multimin_function_fdf*>(&my_func), beta_init, init_step, init_tol );

        size_t iter = 0;
        int successCount(0);
        bool done(false);
      
#ifdef MASSIVE_DEBUG_
       cout << "Iteration: " << totalIterations << "\nstart_metric = " << s->f
            << printArray(alpha_offset.get(), nParameters, "\nalpha_offset")
            << printArray(alpha.get(), nParameters, "\nstart_alpha")
            << printArray(grad_alpha.get(), nParameters, "\nstart_agrad")
            << printArray(s->gradient->data, nFreeParameters, "\nstart_bgrad") << endl;
        Array<double> gwrapper(grad_alpha.get(), nTotalImages, nModes);
        Ana::write("grad_alpha_"+to_string(totalIterations)+".f0",gwrapper);
        dump("iter_"+to_string(totalIterations));
        //job.globalData->constraints.reverse( s->gradient->data, grad_alpha );
        //Ana::write("grad_alphac_"+to_string(totalIterations)+".f0",gwrapper);
        //cout << printArray(grad_alpha.get(), nParameters, "\nreversed_agrad") << endl;
#endif

        do {
            iter++;
            boost::this_thread::interruption_point();
            status = gsl_multimin_fdfminimizer_iterate( s );
            if( failCount > maxFails ) {
                LOG_ERR << "Giving up after " << failCount << " failures for patch#" << data->id << " (index=" << data->index << " region=" << data->roi << ")" << ende;
                status = GSL_FAILURE;
                modeCount = 0;                  // exit outer loop.
                done = true;                    // exit inner loop.
                break;
            }
            
            gradientMethod = gradientMethods[GM_VOGEL];   // FIXME: old code just uses the cfg-setting the first iteration, then switches to Vogel.
            
            if( status == GSL_ENOPROG ) {
                LOG_TRACE << "iteration: " << (totalIterations+iter) << "  GSL reports no progress." << ende;
                failCount++;
                break;
            } else if( status ) {
                LOG_WARN << "GSL error in iteration " << (totalIterations+iter) << ".  type: " << gsl_strerror(status) << ende;
            }
            thisMetric = s->f;
            if( !isfinite(thisMetric) ) {
                LOG_WARN << "Error in iteration " << (totalIterations+iter) << ". Metric is not finite." << ende;
                failCount++;
                break;
            }
            data->metrics.push_back(thisMetric);
            gradNorm = gsl_blas_dnrm2(s->gradient)*gradScale;
            if( !isfinite(gradNorm) ) {
                LOG_WARN << "Error in iteration " << (totalIterations+iter) << ". Gradient is not finite." << ende;
                failCount++;
                break;
            }
            if( gradNorm < job.FTOL ) {
                LOG_DETAIL << "Norm(grad) is below threshold." << ende;
                modeCount = 0;                  // exit outer loop.
                done = true;                    // exit inner loop.
                break;
            }
            if( iter>1 ) {
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

        } while( (!done || iter < job.minIterations) && iter < maxIterations );

        GSL_MULTIMIN_FN_EVAL_F( &my_func, s->x );
        totalIterations += iter;
        if( status != GSL_FAILURE ) { // bad first iteration -> don't print.
            size_t nActiveModes = activeModes.size();
            if( modeCount ) {
                if( nActiveModes )  gradNorm /= nActiveModes;
                LOG_DETAIL << "Patch" << (string)data->index << ":  After " << totalIterations << " iteration" << (totalIterations>1?"s":" ")
                 << "  metric=" << thisMetric << "  rmetric=" << (thisMetric/initialMetric)
                 << " norm(grad)=" << gradNorm << " using " << nActiveModes << "/" << nModes << " modes." << ende;
            }
        }

        if( modeCount ) {
            shiftAndInit();
            memcpy( alpha_offset.get(), alpha.get(), nParameters*sizeof(double) );          // store current solution
        }
    }       // end for-loop
    
    LOG << "Patch" << (string)data->index << ":  After " << totalIterations << " iterations:  metric=" << thisMetric
        << "  (relative=" << (thisMetric/initialMetric) << ")  " << timer.print() << ende;
    myInfo.status.statusString = patchString + " completed";
    

#ifdef RDX_DUMP_PATCHDATA
    dump( "final_patch_"+(string)data->index );
#endif
    
//    alphaPtr = alpha.get();
    for( auto& objData: data->objects ) {
        if(!objData) continue;
        for( auto& cd: objData->channels ) {
            cd->images.clear();         // don't need input data anymore.
        }
        objData->myObject->restorePatch( *objData );
#ifdef RDX_DUMP_PATCHDATA
        Ana::write( "patch_"+(string)data->index+"_obj_"+to_string(objData->myObject->ID)+"_result.f0", objData->img );
#endif
    }
    
    data->finalMetric = thisMetric;
    data->waveFronts.ids.clear();
    data->waveFronts.alpha.resize( wavefronts.size(), nModes );
    data->storeAlpha(alpha.get());
    for( auto& it: wavefronts ) {
        data->waveFronts.ids.push_back(it.first);
    }
    
#ifdef RDX_DUMP_PATCHDATA
    data->dump( "final" );
#endif
    

    gsl_multimin_fdfminimizer_free( s );
    gsl_vector_free(beta_init);
    
    if( status == GSL_FAILURE ) {
        // TODO throw e.g. "part_fail" or just flag it ??
    }

}


template <typename T>
void Solver::shiftAndInit( const T* a, bool doReset ) {
    
    progWatch.set( nTotalImages );
    for( const auto& o: job.objects ) {
        o->progWatch.clear();
        o->imgShifted = static_cast<int>( doReset );  // force re-initialization if reset is passed
        o->progWatch.set( o->nImages() );
        o->progWatch.setHandler( std::bind( &Object::reInitialize, o.get(), std::ref(service), std::ref(progWatch), doReset ) );
        for( const auto& c: o->channels ) {
            for( const auto& im: c->getSubImages() ) {
                service.post( [&,a](){
                    o->imgShifted.fetch_or(im->adjustShifts(a));
                    ++o->progWatch;
                } );
                a += nModes;
            }
        }

    }
    progWatch.wait();

}
template void Solver::shiftAndInit( const double*, bool );
template void Solver::shiftAndInit( const float*, bool );


void Solver::alignWavefronts( void ) {

    progWatch.set( nTotalImages );
    for( const auto& o: objects ) {
        for( const auto& c: o->getChannels() ) {
            shared_ptr<SubImage> refIm;
            for( const auto& im: c->getSubImages() ) {
                service.post ([&im,refIm,this] {    // use a lambda to ensure these calls are sequential
                    im->alignAgainst( refIm );
                    ++progWatch;
                });
                refIm = im;
            }
        }
    }
    progWatch.wait();

}


void Solver::zeroAlphas( void ) {
    
    memset( alpha.get(), 0, nParameters*sizeof(double) );
    memset( grad_alpha.get(), 0, nParameters*sizeof(double) );
    memset( alpha_offset.get(), 0, nParameters*sizeof(double) );

}


template <typename T> 
void Solver::applyAlpha( T* a ) {

    progWatch.set( nTotalImages );
    for( const auto& o: objects ) {
        for( const auto& c: o->getChannels() ) {
            for( const auto& im: c->getSubImages() ) {
                service.post ([&im,a,this] {    // use a lambda to ensure these calls are sequential
                    im->calcPhi( a );
                    im->calcPFOTF();
                    ++progWatch;
                });
                a += nModes;
            }
        }
    }
    progWatch.wait();

}
template void Solver::applyAlpha( double* a );
template void Solver::applyAlpha( float* a );


void Solver::applyBeta( const gsl_vector* b ) {

    job.globalData->constraints.reverseAndAdd( b->data, alpha_offset.get(), alpha.get() );
    applyAlpha();

}


void Solver::applyBeta( const gsl_vector* b, double scale ) {

    job.globalData->constraints.reverse( b->data, alpha.get(), scale );
    applyAlpha();

}


void Solver::applyBeta( double scale ) {

    job.globalData->constraints.reverse( beta.get(), alpha.get(), scale );
    applyAlpha();
    
}


void Solver::applyConstraints( const double* a, double* b ) {

    progWatch.set( job.globalData->constraints.ns_cols.size() );
    for( auto& r: job.globalData->constraints.ns_cols ) {
        service.post( [this,&r,&a,&b]() {
            double tmp(0);
            for( auto& e: r.second ) tmp += e.second * a[e.first];
            b[r.first] += tmp;
            ++progWatch;
        });
    }
    progWatch.wait();

}


void Solver::reverseConstraints( const double* b, double* a ) {

    progWatch.set( job.globalData->constraints.ns_rows.size() );
    for( auto& r: job.globalData->constraints.ns_rows ) {
        service.post([this,&r,&a,&b]() {
            double tmpD(0);
            for( auto& e: r.second ) tmpD += e.second * b[e.first];
            a[r.first] += tmpD;
            ++progWatch;
        });
    }
    progWatch.wait();

}


void Solver::zeroAvgTilt( double* a, int m ) {
}


void Solver::zeroAvgTilts( double* a, int m1, int m2 ) {
}


void Solver::loadInit( const PatchData::Ptr pd, double* a ) const {
    
    for( auto& objData: pd->objects ) {
        if( objData ) objData->myObject->getInit( *objData, a );
        a += objData->myObject->nObjectImages*nModes;
    }
    
}


void Solver::initImages( double* a ) {

    progWatch.set( nTotalImages );
    for( const auto& o: job.objects ) {
        for( const shared_ptr<Channel>& c: o->channels ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                service.post( [&,a](){
                    im->adjustShifts( a );
                    im->initialize(true);
                    ++progWatch;
                } );
                a += nModes;
            }
        }
    }
    progWatch.wait();
    
}


namespace {
    
    const size_t nMaxImages = 30;

}


double Solver::metric(void) {
    

    double sum(0);
    progWatch.set( objects.size() );
    for( const shared_ptr<Object>& o : objects ) {
        //if( o->weight > 0 ) {
            o->initPQ();
            o->progWatch.set( o->nImages() );
            o->progWatch.setHandler( [o,&sum,this](){
                o->calcMetric();
                lock_guard<mutex> lock(mtx);
                sum += o->metric();
                ++progWatch;
            });
            for( const shared_ptr<Channel>& c: o->getChannels() ) {
                size_t nImgs = c->getSubImages().size();
                size_t begIndex(0);
                while( begIndex < nImgs ) {
                    size_t endIndex = begIndex + nMaxImages;
                    if( endIndex > nImgs ) endIndex = nImgs;
                    o->progWatch.increaseTarget();
                    service.post( [o,c,begIndex,endIndex,this] {
                        const vector< shared_ptr<SubImage> >& imgs = c->getSubImages();
                        double* tmpD = tmp()->D.get();
                        complex_t* tmpC = tmp()->C.get();
                        std::fill_n( tmpD, otfSize2, 0.0 );
                        std::fill_n( tmpC, otfSize2, complex_t(0) );
                        for( size_t i=begIndex; i<endIndex; ++i ) {
                            imgs[i]->addPQ( tmpC, tmpD );
                            ++o->progWatch;
                        }
                        o->addToPQ( tmpC, tmpD );
                        ++o->progWatch;
                    } );
                    begIndex = endIndex;
                }
            }
        //}
    }

    if( job.reg_alpha > 0 ) {
        lock_guard<mutex> lock(mtx);
        double* alphaPtr = alpha.get();
        double* rawPtr = regAlphaWeights.get();
        for( size_t i=0; i<nParameters; ++i ) {
            sum += 0.5*alphaPtr[i]*alphaPtr[i]*rawPtr[i];
        }
    }
    
    progWatch.wait();

    return sum; ///nTotalPixels;
    
}


void Solver::calcPQ(void) {

    progWatch.set( objects.size() );
    for( const shared_ptr<Object>& o : objects ) {
        o->initPQ();
        o->progWatch.set( o->nImages() );
        o->progWatch.setHandler( [o,this](){
            o->calcHelpers();
            ++progWatch;
        });
        for( const shared_ptr<Channel>& c: o->getChannels() ) {
            size_t nImgs = c->getSubImages().size();
            size_t begIndex(0);
            while( begIndex < nImgs ) {
                size_t endIndex = begIndex + nMaxImages;
                if( endIndex > nImgs ) endIndex = nImgs;
                o->progWatch.increaseTarget();
                service.post( [o,c,begIndex,endIndex,this] {
                    const vector< shared_ptr<SubImage> >& imgs = c->getSubImages();
                    double* tmpD = tmp()->D.get();
                    complex_t* tmpC = tmp()->C.get();
                    std::fill_n( tmpD, otfSize2, 0.0 );
                    std::fill_n( tmpC, otfSize2, complex_t(0) );
                    for( size_t i=begIndex; i<endIndex; ++i ) {
                        imgs[i]->addPQ( tmpC, tmpD );
                        ++o->progWatch;
                    }
                    o->addToPQ( tmpC, tmpD );
                    ++o->progWatch;
                } );
                begIndex = endIndex;
            }
        }
    }

    progWatch.wait();
    
}


void Solver::gradient(void) {

    double* alphaPtr = alpha.get();
    double* gAlphaPtr = grad_alpha.get();
    if( job.reg_alpha > 0 ) {
        std::transform( alphaPtr, alphaPtr+nParameters, regAlphaWeights.get(), gAlphaPtr,
            std::multiplies<double>() );
    } else {
        memset( gAlphaPtr, 0, nParameters*sizeof(double) );
    }

    bool* eM = enabledModes.get();
    progWatch.set(nTotalImages );
    for( const shared_ptr<Object>& o: objects ) {
        //if( o->weight > 0 ) {
            for( const shared_ptr<Channel>& c: o->getChannels() ) {
                for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                    service.post ([this, eM, o, im, gAlphaPtr] {    // use a lambda to ensure these calls are sequential
                        im->calcVogelWeight( o->PQ.get(), o->PS.get(), o->QS.get() );
                        for ( uint16_t m=0; m<nModes; ++m ) {
                            if( eM[m] ) {
                                gAlphaPtr[m] = o->weight*gradientMethod(*im, m );
                            }
                        }
                        ++progWatch;
                    });
                    alphaPtr += nModes;
                    gAlphaPtr += nModes;
                }
            }
        //}
    }
    progWatch.wait();

}
 
 
void Solver::gradient( gsl_vector* grad ) {
    
    gradient();
    job.globalData->constraints.apply( grad_alpha.get(), grad->data );
    
}


void Solver::clear( void ) {

    enabledModes.reset();
    alpha.reset();
    grad_alpha.reset();
    alpha_offset.reset();

    regAlphaWeights.reset();
    beta.reset();
    grad_beta.reset();
    search_dir.reset();
    tmp_beta.reset();
    
}


void Solver::dump( string tag ) {

    progWatch.set( objects.size()+1 );
    for( auto & o : objects ) {
        service.post( [this,tag, &o](){
            o->dump(tag);
            ++progWatch;
        });
    }
    service.post( [this,tag](){
        Ana::write (tag + "_window.f0", window);
        Ana::write (tag + "_noisewindow.f0", noiseWindow);

        Array<double> wrapper(alpha.get(), nTotalImages, nModes);
        Ana::write(tag + "_alpha.f0",wrapper);
        wrapper.wrap(alpha_offset.get(), nTotalImages, nModes);
        Ana::write(tag + "_alpha_offset.f0",wrapper);
        wrapper.wrap(grad_alpha.get(), nTotalImages, nModes);
        Ana::write(tag + "_grad_alpha.f0",wrapper);
        ++progWatch;
    });
    progWatch.wait();
}



