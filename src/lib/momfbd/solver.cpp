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


thread_local thread::TmpStorage Solver::tmp;

//#define DEBUG_
//#define RDX_DUMP_PATCHDATA

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
            while( counter > 0 ) cv.wait(lck);
        }
        sync_counter& operator--() {
            unique_lock<mutex> lck(mtx);
            if( --counter <= 0 ) cv.notify_all();
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
    redux::image::windowInPlace( window, patchSize / 8);

    int md = std::min<int>( 256, patchSize );                       // new size (maximum=256)
    md -= ( md % 2 );
    noiseWindow.resize(md,md);
    noiseWindow = 1.0;
    redux::image::windowInPlace( noiseWindow, md / 16);     // FIXME: old code specifies md/16, but applies it after "window", so it is actually the product...
    
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

    // FIXME: this is a rather kludgy way to allow for resizing the TLS, think of something neater...
    set<std::thread::id> initDone;
    progWatch.set( nThreads );
    for( uint16_t i=0; i<nThreads; ++i ) { 
        service.post([&](){
            unique_lock<mutex> lock(mtx);
            std::thread::id id = std::this_thread::get_id();
            auto ret = initDone.emplace(id);
            if( ret.second ) {
                tmp.resize( patchSize, pupilSize );
                ++progWatch;
                lock.unlock();
                progWatch.wait();
            }
        });
    }
    tmp.resize( patchSize, pupilSize );
    progWatch.wait();

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
    
    return metric2();

}


void Solver::my_df( const gsl_vector* beta, void* params, gsl_vector* df ) {

    applyBeta( beta );
    calcPQ2();
    gradient2( df );

}


void Solver::my_fdf( const gsl_vector* beta, void* params, double* f, gsl_vector* df ) {
    
    applyBeta( beta );
    *f = metric2();
    gradient2( df );
    
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
        //if( o->weight > 0 ) {
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
        //}
    }
    
    alphaPtr = grad_alpha;
    phiPtr = tmpPhiGrad.get();
    for( const auto& o: objects ) {
        //if( o->weight > 0 ) {
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
        //}
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
        //if( o->weight > 0 ) {
            double normalization = sqrt(1.0 / (o->pupil.area*otfSize2));
            for( const auto& c: o->getChannels() ) {
                for( const auto& im: c->getSubImages() ) {
                    ++sc;
                    service.post( [&im, &o, normalization, step, otfPtr, phiPtr, phiGradPtr, &sc] {
                        const double* pupilPtr = o->pupil.get();
                        for( const auto& ind: o->pupil.pupilInOTF ) {
                            otfPtr[ind.second] = polar(pupilPtr[ind.first]*normalization, phiPtr[ind.first]+step*phiGradPtr[ind.first]);
                        }
                        //im->tmpOTF.ft(otfPtr); FIXME: broke with thread_local
                        //im->tmpOTF.norm();
                        //im->tmpOTF.getIFT(im->OTF.get());
                        //FourierTransform::reorder( im->OTF.get(), im->OTF.dimSize(0), im->OTF.dimSize(1) );
                        --sc;
                    });
                    otfPtr += otfSize2;
                    phiPtr += pupilSize2;
                    phiGradPtr += pupilSize2;
                }
            }
        //}
    }
    sc.wait();

    for( const auto& o : objects ) {
        //if( o->weight > 0 ) {
            ++sc;
            service.post( [o,&sc] {
                o->initPQ();
                o->addAllPQ();
                o->calcMetric();
                --sc;
            } );
        //}
    }
    sc.wait();
    
    double sum(0);
    for( const auto& o : objects ) {
        sum += o->metric();
    }

    return sum;
    
}


void Solver::run( PatchData::Ptr data ) {

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

    gradientMethod = gradientMethods[job.gradientMethod];
    gradientMethod = gradientMethods[GM_VOGEL];   // FIXME: old code just uses the cfg-setting the first iteration, then switches to Vogel.
    gradientMethod = gradientMethods[GM_DIFF];   // FIXME: old code just uses the cfg-setting the first iteration, then switches to Vogel.
    
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc( minimizerType, nFreeParameters );

    gsl_vector *beta_init = gsl_vector_alloc( nFreeParameters );
    memset(beta_init->data,0,nFreeParameters*sizeof(double));
    memset(alpha_offset,0,nParameters*sizeof(double));
// for(uint i=0; i<nFreeParameters; ++i ) {
//  if (i%nModes==2) alpha_offset[i] = -1.68e-06;
// }
//     for(uint i=0; i<nFreeParameters; ++i ) {
//         beta_init->data[i] = 1;
//         job.globalData->constraints.reverse( beta_init->data, alpha_offset );
//         cout << i << printArray(alpha_offset, nParameters,"  alpha") << endl;
//         memset(beta_init->data,0,nFreeParameters*sizeof(double));
//         memset(alpha_offset,0,nParameters*sizeof(double));
//     }
//     
//     for(uint i=0; i<nParameters; ++i ) {
//         alpha_offset[i] = 1;
//         job.globalData->constraints.apply( alpha_offset, beta_init->data );
//         cout << i << printArray(beta_init->data, nFreeParameters,"  beta") << endl;
//         memset(beta_init->data,0,nFreeParameters*sizeof(double));
//         memset(alpha_offset,0,nParameters*sizeof(double));
//     }

    memset(s->x->data,0,nFreeParameters*sizeof(double));
    memset(s->dx->data,0,nFreeParameters*sizeof(double));
    memset(s->gradient->data,0,nFreeParameters*sizeof(double));
    s->f = 0;
    

/*
cout << "nParameters: " << nParameters << "  nFreeParameters: " << nFreeParameters << endl;
shared_ptr<Object> obj0 = *job.objects.begin();
cout << "shiftToAlpha: " << obj0->shiftToAlpha << endl;

Array<float> slaskinit;
//readFile("/scratch/tomas/momfbdsol10_5.f0",slaskinit);

memset(alpha,0,nParameters*sizeof(double));
memset(enabledModes,1,nModes*sizeof(uint16_t));

applyAlpha();
calcPQ();
gradient(beta_init);
cout << printArray(grad_alpha, nModes, "grad0") << endl;

//dump("init0");
//  Array<double> wrap(grad_alpha,nParameters);
//  Ana::write("grad1",wrap);
//cout << printArray( grad_alpha, nParameters, "\ngrad1 " ) << endl;
//cout << printArray( beta_init->data, nFreeParameters, "\ngrad_beta1 " ) << endl;
job.globalData->constraints.reverse( beta_init->data, grad_alpha );
//cout << printArray(grad_alpha, nParameters, "\ngrad2 ") << endl;
//  Ana::write("grad2",wrap);
double refm = metric();
cout << "\nMetric0: " << refm << endl;

readFile("/scratch/tomas/sol60r.f0",slaskinit);
std::transform(alpha,alpha+nParameters,slaskinit.get(),alpha,
    [](const double& a, const float&b){
    return b;
});
applyAlpha();
cout << "\nMetricR: " << (metric()/refm) << endl;

vector<float> testM,testMr,posm;
for( double pos=-2.0; pos<2; pos += 0.1 ) {
    std::transform(alpha,alpha+nParameters,slaskinit.get(),alpha,
        [pos](const double& a, const float&b){
        return pos*b;
        //return 0.5*b;
        //return -b;
    });
    applyAlpha();
    testMr.push_back(metric()/refm);
    posm.push_back(pos);
}



readFile("/scratch/tomas/sol60m.f0",slaskinit);
std::transform(alpha,alpha+nParameters,slaskinit.get(),alpha,
    [](const double& a, const float&b){
    return -b;
});
applyAlpha();
cout << "\nMetricM: " << (metric()/refm) << endl;

for( double pos=-2.0; pos<2; pos += 0.1 ) {
    std::transform(alpha,alpha+nParameters,slaskinit.get(),alpha,
        [pos](const double& a, const float&b){
        return -pos*b;
        //return 0.5*b;
        //return -b;
    });
    applyAlpha();
    testM.push_back(metric()/refm);

}
cout << printArray( testMr, "testMr" ) << endl;
cout << printArray( testM, "testMm" ) << endl;
cout << printArray( posm, "posm" ) << endl;
cout << "plot,posm,testMm" << endl;
cout << "oplot,posm,testMr" << endl;



exit(0);



std::transform(alpha,alpha+nParameters,slaskinit.get(),alpha,
    [](const double& a, const float&b){
    return -0.5*b;
    //return 0.5*b;
    //return -b;
});
applyAlpha();
calcPQ();
gradient(beta_init);
cout << printArray(grad_alpha, nModes, "gradh_diff") << endl;

gradientMethod = gradientMethods[GM_VOGEL];
gradient(beta_init);
cout << printArray(grad_alpha, nModes, "gradh_vog") << endl;

//alpha[0] = -1.81099e-07;
//alpha[1] = 6*obj0->shiftToAlpha.x;
//cout << printArray( alpha, 80, "alpha0" ) << endl;
applyAlpha();
dump("init1");
cout << "\nMetric1: " << (metric()/refm) << endl;
applyAlpha();
cout << "\nMetric2: " << (metric()/refm) << endl;
calcPQ();
gradient(beta_init);
dump("init2");
//cout << printArray( grad_alpha, nParameters, "\ngrad3 " ) << endl;
//  Ana::write("grad3",wrap);
//cout << printArray( beta_init->data, nFreeParameters, "\ngrad_beta3 " ) << endl;
job.globalData->constraints.reverse( beta_init->data, grad_alpha );
//cout << printArray(grad_alpha, nParameters, "\ngrad4 ") << endl;
//  Ana::write("grad4",wrap);
memset(beta_init->data,0,nFreeParameters*sizeof(double));

double slaskMetric = GSL_MULTIMIN_FN_EVAL_F( &my_func, beta_init );
cout << "Initial metric = " << (std::scientific) << slaskMetric << endl;
//return;
exit(0);
*/

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

    timer.start();
    logger.flushAll();

//    memset( alpha_offset, 0, nParameters*sizeof(double) );
    
    
    loadInit( data, alpha_offset );
    
    //initImages( alpha_offset );
    shiftAndInit( alpha_offset, true );     // force initialization
    
    double initialMetric = GSL_MULTIMIN_FN_EVAL_F( &my_func, beta_init );
    LOG << "Initial metric = " << initialMetric << ende;

//cout << "Solver:  " << __LINE__ << endl;
//    Array<double> wrapq(alpha_offset, nParameters);
//    Ana::write( "patch_"+(string)data->index+"_ainit.f0", wrapq );
//cout << "Solver:  " << __LINE__ << endl;
/*Array<float> slaskinit;
readFile("/scratch/tomas/sol60r.f0",slaskinit);
std::transform(alpha,alpha+nParameters,slaskinit.get(),alpha_offset,
    [](const double& a, const float&b){
    return -b;
});*/
    
    string patchString = alignLeft(to_string(job.info.id) + ":" + to_string(data->id),8);
//dump("init");
gradientMethod = gradientMethods[GM_DIFF];
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
            tol = 5E-4; //job.FTOL;
            maxIterations = job.maxIterations;
        } else {
            modeCount += job.nModeIncrement;
            //failCount = 0;           // TODO: fix/test fail-reporting before using
        }

//memset(enabledModes,1,nModes*sizeof(uint16_t));
//gradientMethod = gradientMethods[GM_DIFF];
//gradientMethod = gradientMethods[GM_VOGEL];
//cout << printArray(alpha_offset,nParameters,"alpha_offset") << endl;
//cout << printArray(enabledModes,nModes,"enabledModes") << endl;
//memset(enabledModes,1,nModes*sizeof(uint16_t));
//cout  << endl << "                                                                     Solver: " << __LINE__ << endl;
//cout  << "GradMethods: " << hexString(gradientMethods[GM_DIFF].target<double(SubImage&, uint16_t, double)>()) <<
//    "  " << hexString(gradientMethods[GM_VOGEL].target<double(SubImage&, uint16_t, double)>()) << endl;
//cout  << "GradMethod: " << hexString(gradientMethod.target<double(SubImage&, uint16_t, double)>()) << endl;
        // The previous position is saved in in alpha_offset, beta only contains the new step, so start at 0 each iteration
        gsl_multimin_fdfminimizer_set( s, static_cast<gsl_multimin_function_fdf*>(&my_func), beta_init, init_step, init_tol );
//cout << printArray(activeModes,"activeModes") << endl;
//cout << "Iteration: " << totalIterations << printArray(grad_alpha, nModes, "\nstart_agrad") << endl;
//cout  << endl << "                                                                     Solver: " << __LINE__ << endl;
//exit(0);
//gradientMethod = gradientMethods[GM_VOGEL];
        size_t iter = 0;
        int successCount(0);
        bool done(false);
//GSL_MULTIMIN_FN_EVAL_DF( &my_func, beta_init, s->gradient );        
#ifdef DEBUG_
//        cout << "Iteration: " << totalIterations << "\nstart_metric = " << s->f
//                       << printArray(alpha_offset, nParameters, "\nalpha_offset")
//                       << printArray(alpha, nParameters, "\nstart_alpha")
//                       << printArray(grad_alpha, nParameters, "\nstart_agrad")
//                       << printArray(s->gradient->data, nFreeParameters, "\nstart_bgrad") << endl;
        Array<double> gwrapper(grad_alpha, nTotalImages, nModes);
        Ana::write("grad_alpha_"+to_string(totalIterations)+".f0",gwrapper);
        //dump("iter_"+to_string(totalIterations));
        job.globalData->constraints.reverse( s->gradient->data, grad_alpha );
        Ana::write("grad_alphac_"+to_string(totalIterations)+".f0",gwrapper);
//        cout << printArray(grad_alpha, nParameters, "\nreversed_agrad") << endl;
//sleep(1);
//exit(0);
#endif
//cout << "                                                                     Solver: " << __LINE__ << endl;
        do {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate( s );
            if( failCount > maxFails ) {
                LOG_ERR << "Giving up after " << failCount << " failures for patch#" << data->id << " (index=" << data->index << " region=" << data->roi << ")" << ende;
                status = GSL_FAILURE;
                modeCount = 0;                  // exit outer loop.
                done = true;                    // exit inner loop.
                //dump("fail");
                //exit(0);
            }
gradientMethod = gradientMethods[GM_VOGEL];
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
            gradNorm = gsl_blas_dnrm2(s->gradient)/thisMetric;
            if( !isfinite(gradNorm) ) {
                LOG_WARN << "Error in iteration " << (totalIterations+iter) << ". Gradient is not finite." << ende;
//                 LOG_ERR << "Gradient is infinite (or NaN) for patch#" << data->id << " (index=" << data->index << " region=" << data->roi << ")" << ende;
//                 status = GSL_FAILURE;
//                 modeCount = 0;                  // exit outer loop.
//                 done = true;                    // exit inner loop.
//                dump("gradfail_"+to_string(failCount));
                failCount++;
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

            status = gsl_multimin_test_gradient( s->gradient, 1e-9 );

        } while( status == GSL_CONTINUE && (!done || iter < job.minIterations) && iter < maxIterations );

        GSL_MULTIMIN_FN_EVAL_F( &my_func, s->x );
        totalIterations += iter;
        if( status != GSL_FAILURE ) { // bad first iteration -> don't print.
            size_t nActiveModes = activeModes.size();
            if( modeCount ) {
                if( nActiveModes )  gradNorm /= nActiveModes;
                LOG_DETAIL << "After " << totalIterations << " iteration" << (totalIterations>1?"s":" ")
                 << "  metric=" << thisMetric << "  rmetric=" << (thisMetric/initialMetric)
                 << " norm(grad)=" << gradNorm << " using " << nActiveModes << "/" << nModes << " modes." << ende;
            }
        }

        if( modeCount ) {
            shiftAndInit( alpha );
            memcpy( alpha_offset, alpha, nParameters*sizeof(double) );
        }
    }       // end for-loop
    
    //dump("final");
    //double finalMetric = GSL_MULTIMIN_FN_EVAL_F( &my_func, s->x );
    //LOG << "Final metric = " << finalMetric << ende;
  
//cout << "                                                                     Solver: " << __LINE__ << endl;
    LOG << "After " << totalIterations << " iterations:  metric=" << thisMetric << "  (relative=" << (thisMetric/initialMetric) << ")" << ende;
    LOG << timer.print() << ende;
    myInfo.status.statusString = patchString + " completed";
    

#ifdef RDX_DUMP_PATCHDATA
    dump( "patch_"+(string)data->index );
#endif
    
    alphaPtr = alpha;
    for( auto& objData: data->objects ) {
        if(!objData) continue;
        objData->myObject->getResults(*objData,alphaPtr);
        alphaPtr += objData->myObject->nObjectImages*nModes;
        objData->setLoaded( objData->size()>6*sizeof(uint64_t) );
        for( auto& cd: objData->channels ) {
            cd->images.clear();         // don't need input data anymore.
        }
#ifdef RDX_DUMP_PATCHDATA
        Ana::write( "patch_"+(string)data->index+"_obj_"+to_string(objData->myObject->ID)+"_result.f0", objData->img );
#endif
    }
    
    data->finalMetric = thisMetric;

    gsl_multimin_fdfminimizer_free( s );
    gsl_vector_free(beta_init);
    
    if( status == GSL_FAILURE ) {
        // TODO throw e.g. "part_fail" or just flag it ??
    }

//cout << "                                                                     Solver: " << __LINE__ << endl;
}


void Solver::shiftAndInit( double* a, bool doReset ) {
    
    job.progWatch.set( job.nImages() );
    for( const auto& o: job.objects ) {
        o->progWatch.clear();
        o->imgShifted = static_cast<int>( doReset );  // force re-initialization if reset is passed
        o->progWatch.set( o->nImages() );
        o->progWatch.setHandler( std::bind( &Object::reInitialize, o.get(), std::ref(service), doReset ) );
        for( const auto& c: o->channels ) {
            for( const auto& im: c->getSubImages() ) {
                service.post( [&,a](){
                    o->imgShifted.fetch_or(im->adjustOffset(a));
                    ++o->progWatch;
                } );
                a += nModes;
            }
        }

    }
    job.progWatch.wait();

}


void Solver::applyAlpha( double* a ) {

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

/*
            template <typename T> void apply(const T* in, T* out) const {
                memset(out,0,nFreeParameters*sizeof(T));
                for( auto& entry: ns_entries ) {
                    out[entry.first.x] += entry.second * in[entry.first.y];
                }
            }

*/
void Solver::applyConstraints( const double* alpha, double* beta ) {

    progWatch.set( job.globalData->constraints.ns_cols.size() );
    for( auto& r: job.globalData->constraints.ns_cols ) {
        service.post( [this,&r,&alpha,&beta]() {
            double tmp(0);
            for( auto& e: r.second ) tmp += e.second * alpha[e.first];
            beta[r.first] += tmp;
            ++progWatch;
        });
    }
    progWatch.wait();

}


void Solver::reverseConstraints( const double* beta, double* a ) {

    progWatch.set( job.globalData->constraints.ns_rows.size() );
    for( auto& r: job.globalData->constraints.ns_rows ) {
        service.post([this,&r,&a,&beta]() {
            double tmp(0);
            for( auto& e: r.second ) tmp += e.second * beta[e.first];
            a[r.first] += tmp;
            ++progWatch;
        });
    }
    progWatch.wait();

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
                    im->adjustOffset( a );
                    im->initialize(true);
                    ++progWatch;
                } );
                a += nModes;
            }
        }
    }
    progWatch.wait();
    
}


double Solver::metric(void) {


    sync_counter sc;
    for( const auto& o : objects ) {
        //if( o->weight > 0 ) {
            ++sc;
            service.post( [o,&sc] {
                o->initPQ();
                o->addAllPQ();
                o->calcMetric();
                --sc;
            } );
        //}
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

namespace {
    
    const size_t nMaxImages = 30;

}


double Solver::metric2(void) {
    

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
                        tmp.C.zero();
                        tmp.D.zero();
                        for( size_t i=begIndex; i<endIndex; ++i ) {
                            imgs[i]->addPQ( tmp.C.get(), tmp.D.get() );
                            ++o->progWatch;
                        }
                        o->addToPQ( tmp.C.get(), tmp.D.get() );
                        ++o->progWatch;
                    } );
                    begIndex = endIndex;
                }
            }
        //}
    }

    if( job.reg_alpha > 0 ) {
        lock_guard<mutex> lock(mtx);
        for( size_t i=0; i<nParameters; ++i ) {
            sum += 0.5*alpha[i]*alpha[i]*regAlphaWeights[i];
        }
    }
    
    progWatch.wait();

    return sum;
    
}


void Solver::calcPQ(void) {

    sync_counter sc;
    for( const auto& o: objects ) {
        //if( o->weight > 0 ) {
            ++sc;
            service.post([&o,&sc] {    // use a lambda to ensure these calls are sequential
                    o->initPQ();
                    o->addAllPQ();
                    o->calcHelpers();
                    --sc;
                });
        //}
    }
    sc.wait();
    
}


void Solver::calcPQ2(void) {

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
                    tmp.C.zero();
                    tmp.D.zero();
                    int cnt(0);
                    for( size_t i=begIndex; i<endIndex; ++i ) {
                        cnt++;
                        imgs[i]->addPQ( tmp.C.get(), tmp.D.get() );
                        ++o->progWatch;
                    }
                    o->addToPQ( tmp.C.get(), tmp.D.get() );
                    ++o->progWatch;
                } );
                begIndex = endIndex;
            }
        }
    }

    progWatch.wait();
    
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
    for( const shared_ptr<Object>& o: objects ) {
        //if( o->weight > 0 ) {
            double grad_step = job.graddiff_step*o->wavelength;
            for( const shared_ptr<Channel>& c: o->getChannels() ) {
                for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                    ++sc;
                    service.post ([this, o, im, grad_step, alphaPtr, gAlphaPtr, &sc] {    // use a lambda to ensure these calls are sequential
                        im->calcVogelWeight( o->PQ.get(), o->PS.get(), o->QS.get() );
                        for ( uint16_t m=0; m<nModes; ++m ) {
                            if( enabledModes[m] ) {
                                //if( fabs(alphaPtr[m]) < 1E-16 ) gAlphaPtr[m] = o->weight*gradientMethod(*im, m, grad_step );
                                //if( fabs(alphaPtr[m]) < 1E-16 ) gAlphaPtr[m] = o->weight*gradientMethods[GM_DIFF](*im, m, grad_step );
                                //else gAlphaPtr[m] = o->weight*gradientMethods[GM_VOGEL](*im, m, grad_step );
                                gAlphaPtr[m] = o->weight*gradientMethod(*im, m, grad_step );
                            }
                        }
                        --sc;
                    });
                    alphaPtr += nModes;
                    gAlphaPtr += nModes;
                }
            }
        //}
    }
    sc.wait();

}
 
 
void Solver::gradient2(void) {

    if( job.reg_alpha > 0 ) {
        std::transform( alpha, alpha+nParameters, regAlphaWeights, grad_alpha,
            std::multiplies<double>() );
    } else {
        memset( grad_alpha, 0, nParameters*sizeof(double) );
    }

    double* alphaPtr = alpha;
    double* gAlphaPtr = grad_alpha;
    progWatch.set( nTotalImages );
    for( const shared_ptr<Object>& o: objects ) {
        //if( o->weight > 0 ) {
            double grad_step = job.graddiff_step*o->wavelength;
            for( const shared_ptr<Channel>& c: o->getChannels() ) {
                for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                    service.post ([this, o, im, grad_step, alphaPtr, gAlphaPtr] {    // use a lambda to ensure these calls are sequential
                        im->calcVogelWeight( o->PQ.get(), o->PS.get(), o->QS.get() );
                        for ( uint16_t m=0; m<nModes; ++m ) {
                            if( enabledModes[m] ) {
                                //if( fabs(alphaPtr[m]) < 1E-16 ) gAlphaPtr[m] = o->weight*gradientMethod(*im, m, grad_step );
                                //if( fabs(alphaPtr[m]) < 1E-16 ) gAlphaPtr[m] = o->weight*gradientMethods[GM_DIFF](*im, m, grad_step );
                                //else gAlphaPtr[m] = o->weight*gradientMethods[GM_VOGEL](*im, m, grad_step );
                                gAlphaPtr[m] = o->weight*gradientMethod(*im, m, grad_step );
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
    job.globalData->constraints.apply( grad_alpha, grad->data );
    
}


void Solver::gradient2( gsl_vector* grad ) {
    
    gradient2();
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



