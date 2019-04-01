#include "redux/momfbd/subimage.hpp"

#include "redux/momfbd/momfbdjob.hpp"

#include "redux/file/fileana.hpp"
#include "redux/logging/logger.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <algorithm>

using namespace redux::file;
using namespace redux::logging;
using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

//#define DEBUG_SIMG_

//#define USE_LUT
#define ALPHA_CUTOFF 0.0

namespace {

#ifdef USE_LUT
#define QUADRANT_SAMPLES    1000000       // angular sampling points (per quadrant) for the sine/cosine lookup table
    const size_t total_samples(5*QUADRANT_SAMPLES);             // extend by pi/2 to fit sine/cosine in one table (cos(x) = sin(x+pi/2))
    double SineLUT[total_samples+1];                            // ...and add 1 to make room for the case of exactly 2*PI when truncating
    double* const CosineLUT = SineLUT + QUADRANT_SAMPLES;
    const size_t period_samples( QUADRANT_SAMPLES<<2 );
    const double pi_x_2 = 2 * M_PI;
    const double angular_step( pi_x_2/(period_samples) );
    const double angular_step_inv( 1.0/angular_step );
    
    int initSineLUT (void) {
        double previousVal = 0, phi = 0;
        for ( int i = 0; i < total_samples + 1; i++ ) {
            phi += angular_step;
            double newVal = sin( phi );
            SineLUT[i] = ( previousVal + newVal ) * 0.5;
            previousVal = newVal;
        }
        return 0;
    }

    inline complex_t getPolar(double magnitude, double phase) {
        double idxD = phase - pi_x_2 * floor( phase / pi_x_2 );
        size_t idx = static_cast<size_t> (idxD*angular_step_inv);
        return std::move( complex_t(magnitude*CosineLUT[idx], magnitude*SineLUT[idx]) );
    }
#undef QUADRANT_SAMPLES
#endif

}

SubImage::SubImage (Object& obj, const Channel& ch, const Array<double>& wind, const Array<double>& nwind)
    : imgSize(0), pupilSize(0), nModes(0), oldRG(0), grad_step(0), object (obj), channel(ch), logger(ch.logger), modes(obj.modes),
      window (wind), noiseWindow(nwind), shifted(false) {

#ifdef USE_LUT
    static int dummy RDX_UNUSED = initSineLUT(); 
#endif
}

SubImage::~SubImage (void) {

}


void SubImage::setPatchInfo( uint32_t i, const PointI& pos, const PointF& resOffs, uint16_t patchSize, uint16_t pupSz, uint16_t nM ) {

    index = i;
    initialOffset = pos;
    channelResidualOffset = resOffs;
    imageShift = 0;
    adjustedTilts = 0;
    
    if( patchSize != imgSize ) {
        imgSize = patchSize;
    }
    
    if( pupSz != pupilSize ) {
        pupilSize = pupSz;
        pupilSize2 = pupilSize*pupilSize;
        otfSize = 2*pupilSize;
        otfSize2 = otfSize*otfSize;
        phi.resize( pupilSize, pupilSize );
        PF.resize( pupilSize, pupilSize );
        vogel.resize( pupilSize, pupilSize );
        OTF.resize( otfSize, otfSize, FT_REORDER|FT_FULLCOMPLEX );
        imgFT.resize( otfSize, otfSize );
        imgFT.zero();
        vogel.zero();
    }
    
    nModes = nM;
    localAlpha.assign( nM, 0.0 );
    grad_step = object.myJob.graddiff_step * object.wavelength;

}


void SubImage::getWindowedImg( Array<double>& img, redux::util::Array<float>& surf, ArrayStats& imgStats, bool rescaled ) const {
    
    if( img.dimSize(0) < imgSize || img.dimSize(1) < imgSize  ) {
        string msg = "SubImage::getWindowedImg: Output container is smaller than image! " +
               printArray( dimensions(), "\n\timgdims" ) +
               printArray( img.dimensions(), "\n\toutdims" );
        throw std::logic_error(msg);
    }

    double* imgPtr;
    for( size_t y=0; y<imgSize; ++y ) {
        const float* rowPtr = ptr(0,y,0);
        std::copy( rowPtr, rowPtr+imgSize, img.ptr(y,0) );
    }
    
    size_t nEl = img.nElements();
    imgStats.getStats( img.get(), nEl, ST_VALUES );
    double avg = imgStats.mean;
    double newAvg(0);
    imgPtr = img.get();
    const double* winPtr = window.get();
    if( surf.nElements() == nEl ) {
        const float* surfPtr = surf.get();
        for( size_t i=0; i<nEl; ++i) {
            // Note: mean(plane) = 0, so subtracting it will not affect the mean value.
//            imgPtr[i] = (imgPtr[i]-avg-surfPtr[i])*winPtr[i]+avg;
            double tmp = (imgPtr[i]-avg-surfPtr[i])*winPtr[i]+avg;
            newAvg += tmp;
            imgPtr[i] = tmp;
        }
    } else {
        transform(imgPtr, imgPtr+nEl, window.get(), imgPtr,
                [avg,&newAvg](const double& a, const double& b) {
//                    return (a-avg)*b+avg;
                    double tmp = (a-avg)*b+avg;
                    newAvg += tmp;
                    return tmp;
                    
                });
    }
    newAvg /= nEl;
    if( rescaled ) {
        avg = object.objMaxMean;
    } /*else {
        avg = 1.0/newAvg;
        s.mean = newAvg;
    }*/
    imgStats.getStats( img.get(), imgSize*imgSize, ST_VALUES|ST_RMS );     // TODO test if this is necessary or if the stats for the larger area is sufficient
    
    transpose( img.get(), imgSize, imgSize );                      // to match MvN
    
}


void SubImage::getWindowedImg( Array<double>& im, redux::util::ArrayStats& s, bool rescaled ) const {
    getWindowedImg( im, object.fittedPlane, s, rescaled);
}


void SubImage::initialize( Object& o, bool doReset ) {
    
    if( doReset ) {
        Solver::tmp->FT.zero();
    } else {
        memcpy( Solver::tmp->FT.get(), imgFT.get(), otfSize*otfSize*sizeof(complex_t));                          // make a temporary copy to pass to addDifftoFT below
    }
    
    getWindowedImg( Solver::tmp->D, o.fittedPlane, stats, false );
    double* imgPtr = Solver::tmp->D.get();
    
    size_t noiseSize = noiseWindow.dimSize(0);
    size_t nEl = Solver::tmp->D.nElements();
    string msg;
    if( doReset ) {
        msg = "Initializing image " + to_string(o.ID) + ":" + to_string(channel.ID) + ":" + to_string(index)
            + "   mean=" + to_string (stats.mean) + " stddev=" + to_string (stats.stddev);
        if( imgSize == noiseSize ) {
            transform( imgPtr, imgPtr+nEl, noiseWindow.get(), Solver::tmp->D2.get(),
                    [&](const double& a, const double& b) { return (a-stats.mean)*b; }
                    );
            imgFT.reset( Solver::tmp->D2.get(), imgSize, imgSize, FT_FULLCOMPLEX);
            stats.noise = imgFT.noise(-1,-1);
        } else {
            size_t offset = (imgSize-noiseSize) / 2;
            Array<double> tmpD(Solver::tmp->D, offset, offset+noiseSize-1, offset ,offset+noiseSize-1);
            tmpD.trim();
            transform(tmpD.get(), tmpD.get()+tmpD.nElements(), noiseWindow.get(), tmpD.get(),
                    [&](const double& a, const double& b) { return (a-stats.mean)*b; }
                    );
            imgFT.reset(tmpD.get(), imgSize, imgSize, FT_FULLCOMPLEX);
            stats.noise = imgFT.noise(-1,-1);
        }

        stats.noise *= channel.noiseFudge;
        double rg = stats.noise/stats.stddev;
        msg += " noise=" + to_string (stats.noise) + " rg=" + to_string(rg);
        msg += "  initial shift=" + (string)imageShift;
        LOG_TRACE << msg << ende;
        o.addRegGamma( rg );
    }
    
    if( imgSize == otfSize ) {                                                                   // imgSize = 2*pupilSize
        imgFT.reset( imgPtr, otfSize, otfSize, FT_FULLCOMPLEX ); //|FT_NORMALIZE );        // full-complex for now, perhaps half-complex later for performance
    } else {                                                                                    // imgSize > 2*pupilSize should never happen (cf. calculatePupilSize)
        int offset = (otfSize - imgSize) / 2;
        Array<double> tmpD (otfSize, otfSize);
        tmpD = stats.mean;
        double* tmpPtr = tmpD.get();
        for (int i = 0; i < imgSize; ++i) {
            memcpy (tmpPtr + offset * (otfSize + 1) + i * otfSize, imgPtr + i * imgSize, imgSize * sizeof (double));
        }
        imgFT.reset (tmpD.get(), otfSize, otfSize, FT_FULLCOMPLEX); //|FT_NORMALIZE );               // full-complex for now, perhaps half-complex later for performance
    }
    
    FourierTransform::reorder(imgFT);                                                          // keep FT in centered form
    
    if( doReset ) {
        o.addToFT( imgFT );
    } else {
        o.addDiffToFT( imgFT, Solver::tmp->FT );
    }
    
}


void SubImage::initialize( bool doReset ) {
    initialize( object, doReset );
}


void SubImage::addFT(Array<double>& ftsum) const {
    const complex_t* ftPtr = imgFT.get();
    double* ftsPtr = ftsum.get();
    for (size_t ind = 0; ind < imgFT.nElements(); ++ind) {
        ftsPtr[ind] += norm (ftPtr[ind]);
    }
}


void SubImage::addPQ (const complex_t* otf, complex_t* P, double* Q) const {

    const complex_t* ftPtr = imgFT.get();
    for( const size_t& ind: object.pupil->otfSupport ) {
        Q[ind] += norm(otf[ind]);                    // Q += sj.re^2 + sj.im^2 = norm(sj)
        P[ind] += conj(ftPtr[ind]) * otf[ind];       // P += conj(ft)*sj            c.f. Vogel
    }

}


// void SubImage::addToPQ(void) const {
//     object.addToPQ( imgFT.get(), OTF.get() );
// }


void SubImage::restore( complex_t* obj, double* obj_norm ) const {

//     bool no_restore(false);         // TODO: implement NO_RESTORE cfg flag
//     if( !no_restore ) {
//         addPQ( OTF.get(), obj, obj_norm );      // Note: should really return obj = ft*conj(sj) for the deconvolution,
//                                                 // and addPQ returns the conjugate.
//     }

    const complex_t* ftPtr = imgFT.get();
    const complex_t* otfPtr = OTF.get();
    for( const size_t& ind: object.pupil->otfSupport ) {
        obj_norm[ind] += norm(otfPtr[ind]);
        obj[ind] += conj(ftPtr[ind]) * (otfPtr[ind]);
    }
    
    
}


double SubImage::metricChange( const complex_t* newOTF ) const {

    const complex_t* oldOTF = OTF.get();
    const complex_t* p = object.P.get();
    const complex_t* ftPtr = imgFT.get();
    const double* q = object.Q.get();
    
    complex_t dp, dsj;
    double dl = 0.0, dq, dn;
    for( size_t& ind: object.pupil->otfSupport ) {
        dsj = newOTF[ind] - oldOTF[ind];            // change in sj
        dp = conj(ftPtr[ind]) * dsj;                // change p and q
        dq = 2.0 * (oldOTF[ind].real() * dsj.real() + oldOTF[ind].imag() * dsj.imag()) + norm (dsj);
        dn = 2.0 * (dp.real() * p[ind].real() + dp.imag() * p[ind].imag()) + norm (dp);
        dl -= (q[ind] * dn - dq * (norm (p[ind]))) / (q[ind] * (q[ind]+dq));
    }
    return dl / otfSize2;
}


double SubImage::gradientFiniteDifference( uint16_t modeIndex ) {

    if( object.weight == 0 ) return 0;
        
    complex_t* otfPtr = Solver::tmp->C.get();
    double* phiPtr = Solver::tmp->D.get();
    memcpy( phiPtr, phi.get(), pupilSize2*sizeof(double) );
    addToPhi( phiPtr, modes->modePointers[modeIndex], grad_step );
    calcOTF( otfPtr, phiPtr );
    return metricChange( otfPtr )/grad_step*object.weight;
    
}


void SubImage::gradientFiniteDifference2( double* agrad, const bool* enabledModes ) {

    if( object.weight == 0 ) return;
    
    complex_t* otfPtr = Solver::tmp->C.get();
    double* phiPtr = Solver::tmp->D.get();
    for( uint16_t m=0; m<nModes; ++m ) {
        if( enabledModes[m] ) {
            memcpy( phiPtr, phi.get(), pupilSize2*sizeof(double) );
            addToPhi( phiPtr, modes->modePointers[m], grad_step );
            calcOTF( otfPtr, phiPtr );
            agrad[m] += metricChange( otfPtr )/grad_step*object.weight;
        }
    }
}


double SubImage::gradientVogel(uint16_t modeIndex ) {

    if( object.weight == 0 ) return 0;
    
    double ret = 0;
    const double* modePtr = modes->modePointers[modeIndex];
    const double* vogPtr = vogel.get();
    double scale = -2.0 * object.pupil->area / otfSize2;
    for( auto & ind : object.pupil->pupilSupport ) {
        ret += scale * vogPtr[ind] * modePtr[ind];
    }
    
    return ret*object.weight;
    
}


void SubImage::gradientVogel2( double* agrad, const bool* enabledModes ) {

    if( object.weight == 0 ) return;
    
    const double* vogPtr = vogel.get();
    double scale = -2.0 * object.pupil->area / otfSize2 * object.weight;
    for( uint16_t m=0; m<nModes; ++m ) {
        if( enabledModes[m] ) {
            const double* modePtr = modes->modePointers[m];
            double tmp(0);
            for( auto & ind : object.pupil->pupilSupport ) {
                tmp += vogPtr[ind] * modePtr[ind];
            }
            agrad[m] += tmp*scale;
        }
    }
}


void SubImage::calcVogelWeight( complex_t* pq, double* ps, double* qs ) {

#ifdef DEBUG_SIMG_
    LOG_TRACE << "SubImage::calcVogelWeight(" << hexString(this) << ")   indexSize=" << object.pupil->pupilInOTF.size() << ende;
#endif
  
    const complex_t* otfPtr = OTF.get();
    const complex_t* ftPtr = imgFT.get();
    const complex_t* pfPtr = PF.get();
    
    complex_t* tmpOtfPtr = Solver::tmp->OTF.get();
    complex_t* glPtr = Solver::tmp->C.get();
    complex_t* hjPtr = Solver::tmp->C2.get();
    
    Solver::tmp->OTF.zero();
    FourierTransform::reorderInto( pfPtr, pupilSize, pupilSize, tmpOtfPtr, otfSize, otfSize );

    Solver::tmp->OTF.getIFT(hjPtr);       // normalize by otfSize2 below
    Solver::tmp->OTF.zero();
    for( const size_t& ind: object.pupil->otfSupport ) {
        tmpOtfPtr[ind] = (pq[ind]*ftPtr[ind] - ps[ind]*otfPtr[ind]) / qs[ind];
    }

    FourierTransform::reorder( tmpOtfPtr, otfSize, otfSize );
    Solver::tmp->OTF.getIFT(glPtr);
    
    double normalization = 1.0/(otfSize2*otfSize2);
    transform( glPtr, glPtr+otfSize2, hjPtr, glPtr,
              [normalization](const complex_t& g,const complex_t& h) {
                  return normalization * h * g.real();
              });

    Solver::tmp->OTF.ft( glPtr );

    FourierTransform::reorder( tmpOtfPtr, otfSize, otfSize );

    double* vogPtr = vogel.get();
    const double* pupilPtr = object.pupil->get();
    for( auto & ind : object.pupil->pupilInOTF ) {
        vogPtr[ind.first] = imag(conj(pfPtr[ind.first])*tmpOtfPtr[ind.second])*pupilPtr[ind.first];
    }


}


void SubImage::calcVogelWeight( void ) {
    calcVogelWeight( object.PQ.get(), object.PS.get(), object.QS.get() );
}


void SubImage::resetPhi(void) {
    //memset(phi.get(), 0, pupilSize2*sizeof (double));    // FIXME: use phi_fixed
    memcpy( phi.get(), channel.phi_fixed.get(), channel.phi_fixed.nElements()*sizeof (double));
}


void SubImage::zeroPhi(void) {
    memset( phi.get(), 0, pupilSize2*sizeof(double) );
}


void SubImage::alignAgainst( const Ptr& refIm ) {

    if( !refIm ) {
        return;
    }
    const complex_t* ftPtr = imgFT.get();
    const complex_t* refPtr = refIm->imgFT.get();
    complex_t* tmpPtr = Solver::tmp->FT.get();
    std::transform( ftPtr, ftPtr+imgFT.nElements(), refPtr, tmpPtr,
        []( const complex_t& a, const complex_t& b ){
            return b*std::conj(a);
    });
    FourierTransform::reorder( Solver::tmp->FT );
    Solver::tmp->FT.getIFT( Solver::tmp->D.get() );
    FourierTransform::reorder( Solver::tmp->D.get(), imgSize, imgSize );
    
    auto arr = reshapeArray( Solver::tmp->D.get(), imgSize, imgSize );
    double** arrPtr = arr.get();
    int mid = imgSize/2;
    imageShift = 0;
    double maxVal = arrPtr[mid][mid];
    for( int x=-19; x<20; ++x ) {
        for( int y=-19; y<20; ++y ) {
            if( arrPtr[mid+y][mid+x] > maxVal ) {
                maxVal = arrPtr[mid+y][mid+x];
                imageShift.x = x;
                imageShift.y = y;
            }
        }
    }
    
}

            
template <typename T>
bool SubImage::adjustShifts( const T* alpha ) {

    PointI oldShift = imageShift;
    PointD oldVal(0,0),newVal(0,0);
    bool ret(false);
    
    int32_t mIndex = object.modes->tiltMode.y;   // FIXME: should be x, but image is transposed
    if( mIndex >= 0 ) {
        double shiftToAlpha = object.shiftToAlpha.y;
        double alphaToShift = 1.0/shiftToAlpha;
        newVal.x = oldVal.x = alpha[mIndex]+localAlpha[mIndex]; // + channelResidualOffset.x;
        int adjust = -lround( oldVal.x*alphaToShift );
        if( adjust && (adjust = shift(1/*FIXME 2*/,adjust)) ) {    // will return the "actual" shift. (the cube-edge might restrict it)
            imageShift.x += adjust;
            localAlpha[mIndex] += adjust*shiftToAlpha;
            newVal.x += adjust*shiftToAlpha;
        }
    }

    mIndex = object.modes->tiltMode.x;   // FIXME: should be y, but image is transposed
    if( mIndex >= 0 ) {
        double shiftToAlpha = object.shiftToAlpha.x;
        double alphaToShift = 1.0/shiftToAlpha;
        newVal.y = oldVal.y = alpha[mIndex]+localAlpha[mIndex]; // + channelResidualOffset.y;
        int adjust = -lround( oldVal.y*alphaToShift );
        if( adjust && (adjust = shift(2/*FIXME 1*/,adjust)) ) {        // will return the "actual" shift. (the cube-edge might restrict it)
            imageShift.y += adjust;
            localAlpha[mIndex] += adjust*shiftToAlpha;
            newVal.y += adjust*shiftToAlpha;
        }
    }

    if( oldShift != imageShift ) {
        //LOG_DEBUG << "SubImage Shifting:  pix2cf=" << object.shiftToAlpha << "   tiltMode=" << object.modes->tiltMode << printArray(alpha,2,"  tilts") << ende;
        LOG_TRACE << "SubImage " << to_string(object.ID) << ":" << to_string(channel.ID) << ":" << to_string(index)
                  << ":  cutout was shifted, from " << oldShift << " to " << imageShift
                  << std::scientific << " oldVal=" << oldVal << "  newVal=" << newVal << ende;
        //newCutout();
        ret = shifted = true;
        //LOG_TRACE << "SubImage Shifting:  " << printArray(first(),"\nfirst") << printArray(last(),"\nlast") << ende;
//LOG << printArray(localAlpha,"newLocalAlpha") << ende;
    }
    return ret;
}
template bool SubImage::adjustShifts( const double* );
template bool SubImage::adjustShifts( const float* );


void SubImage::resetShifts( void ) {
    
    shift( 2/*FIXME 1*/, -imageShift.y );
    shift( 1/*FIXME 2*/, -imageShift.x );
    imageShift = 0;
    int32_t mIndex = object.modes->tiltMode.x;
    if( mIndex >= 0 ) {
        localAlpha[mIndex] = 0;
    }
    
    mIndex = object.modes->tiltMode.y;
    if( mIndex >= 0 ) {
        localAlpha[mIndex] = 0;
    }
    
}


void SubImage::addToPhi( double* phiPtr, const double* modePtr, double w ) const {

    transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
              [w](const double& p, const double& m) {
                  return p + w*m;
              });

}

template <typename T>
void SubImage::addToPhi( const T* a, double* phiPtr ) const {

#ifdef DEBUG_SIMG_
    LOG_TRACE << "SubImage(" << object.ID << ":" << index << ")::addToPhi()" << printArray( a, nModes, "  newAlpha" ) << ende;
#endif

    for( unsigned int i=0; i<nModes; ++i ) {
        double alpha = a[i] + localAlpha[i];
        if( fabs(alpha) > ALPHA_CUTOFF ) {
            const double* modePtr = modes->modePointers[i];
            transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
                [alpha](const double& p, const double& m) {
                    return p + alpha*m;
                });
        }
    }

}
template void SubImage::addToPhi( const double* a, double* phiPtr ) const;
template void SubImage::addToPhi( const float* a, double* phiPtr ) const;


template <typename T>
void SubImage::calcPhi( const T* a, double* phiPtr ) const {

    memcpy( phiPtr, channel.phi_channel.get(), pupilSize2*sizeof(double) );
    //memset( phiPtr, 0, pupilSize2*sizeof(double) );    // FIXME: use phi_fixed
    addToPhi( a, phiPtr );

}
template void SubImage::calcPhi( const double* a, double* phiPtr ) const;
template void SubImage::calcPhi( const float* a, double* phiPtr ) const;


void SubImage::calcOTF(complex_t* otfPtr, const double* phiOffset, double scale) {

    const double* phiPtr = phi.get();
    const double* pupilPtr = object.pupil->get();
    
    for (auto & ind : object.pupil->pupilInOTF) {
#ifdef USE_LUT
        otfPtr[ind.second] = getPolar( pupilPtr[ind.first]*channel.otfNormalization, phiPtr[ind.first]+scale*phiOffset[ind.first]);
#else
        otfPtr[ind.second] = polar(pupilPtr[ind.first]*channel.otfNormalization, phiPtr[ind.first]+scale*phiOffset[ind.first]);
#endif
    }

    Solver::tmp->OTF.ft(otfPtr);
    Solver::tmp->OTF.norm();
    Solver::tmp->OTF.getIFT(otfPtr);
    FourierTransform::reorder(otfPtr, otfSize, otfSize);

}


void SubImage::calcOTF( complex_t* otfPtr, const double* phiPtr ) const {

    std::fill_n( otfPtr, otfSize2, complex_t(0) );

    const double* pupilPtr = object.pupil->get();

    for( auto& ind : object.pupil->pupilInOTF ) {
#ifdef USE_LUT
        otfPtr[ind.second] = getPolar( pupilPtr[ind.first]*channel.otfNormalization, phiPtr[ind.first]);
#else
        otfPtr[ind.second] = polar(pupilPtr[ind.first]*channel.otfNormalization, phiPtr[ind.first]);
#endif
    }
    
    Solver::tmp->OTF.ft(otfPtr);
    Solver::tmp->OTF.norm();
    Solver::tmp->OTF.getIFT(otfPtr);
    FourierTransform::reorder(otfPtr, otfSize, otfSize);

}


void SubImage::calcPFOTF(void) {
    
#ifdef DEBUG_SIMG_
    LOG_TRACE << "SubImage::calcPFOTF(" << hexString(this) << ")   indexSize=" << object.pupil->pupilInOTF.size() << ende;
#endif

    complex_t* pfPtr = PF.get();
    complex_t* otfPtr = OTF.get();
    const double* phiPtr = phi.get();

    std::fill_n( pfPtr, pupilSize2, complex_t(0) );
    std::fill_n( otfPtr, otfSize2, complex_t(0) );
    
    const double* pupilPtr = object.pupil->get();
    
    for( const auto& ind: object.pupil->pupilInOTF ) {
#ifdef USE_LUT
        pfPtr[ind.first] = getPolar( pupilPtr[ind.first], phiPtr[ind.first]);
#else
        pfPtr[ind.first] = polar(pupilPtr[ind.first], phiPtr[ind.first]);
#endif
        otfPtr[ind.second] = channel.otfNormalization*pfPtr[ind.first];
    }

    Solver::tmp->OTF.ft( otfPtr );
    Solver::tmp->OTF.norm();
    Solver::tmp->OTF.getIFT( otfPtr );
    FourierTransform::reorder( otfPtr, otfSize, otfSize );
    
}


void SubImage::addPSF( double* outPSF ) const {
    
    complex_t* cPtr = Solver::tmp->C.get();
    OTF.getIFT( cPtr );     // FIXME: Using forward transform to match MvN (who conjugates the transform)
    FourierTransform::reorder( cPtr, otfSize, otfSize );
//    double normalization(0.0);        // TBD: should the PSF be auto-scaled or filtered ??
//     std::for_each( tmp.get(), tmp.get()+otfSize2,
//         [&normalization]( const complex_t& c ) {
//             normalization += std::abs(std::real(c));
//         });
    double normalization = 1.0 / otfSize2;
    if( normalization > 0.0 ) {
//         normalization = 1.0/normalization;
        std::transform( outPSF, outPSF+otfSize2, cPtr, outPSF,
            [normalization]( const double& p, const complex_t& c ) {
                //double val = std::max( std::real(c), 0.0 );
                //double val = std::abs( std::real(c) );
                //double val = std::real(c);
                return p+normalization*std::real(c);
            });
    }

}


void SubImage::getPSF( double* psf ) const {
    memset( psf, 0, otfSize2*sizeof(double) );
    addPSF( psf );
}


redux::util::Array<double> SubImage::getPSF( void ) const {
    using namespace redux::image;
    redux::util::Array<double> tmp( otfSize, otfSize );
    getPSF(tmp.get());
    return tmp;
}


string SubImage::idString( void ) const {
    return to_string(object.ID) + ":" + to_string(channel.ID) + ":" + to_string(index);
}


void SubImage::dump( std::string tag ) const {

    tag += "_" + idString();
 //   cout << "    dumping image:  this=" << hexString(this) << " with tag=" << tag << endl;
 //   cout << "                  phiPtr=" << hexString(phi.get()) << endl;
    Array<double> img( imgSize, imgSize );
    ArrayStats s;
    getWindowedImg(img,s,true);
    Ana::write (tag + "_rimg.f0", *this);
    Ana::write (tag + "_img.f0", img);
    Ana::write (tag + "_phi.f0", phi);
    Ana::write (tag + "_pf.f0", PF);
    Ana::write (tag + "_otf.f0", OTF);
    Ana::write (tag + "_psf.f0", getPSF());
    Ana::write (tag + "_imgFT.f0", imgFT);
    Ana::write (tag + "_window.f0", window);
    Ana::write (tag + "_vogel.f0", vogel);
    Ana::write (tag + "_alpha.f0", localAlpha);
    


}

#undef ALPHA_CUTOFF

