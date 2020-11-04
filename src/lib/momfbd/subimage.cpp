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
        return complex_t(magnitude*CosineLUT[idx], magnitude*SineLUT[idx]);
    }
#undef QUADRANT_SAMPLES
#endif

}

SubImage::SubImage (Object& obj, const Channel& ch, const Array<double>& wind, const Array<double>& nwind)
    : imgSize(0), pupilSize(0), nModes(0), rowStride(0), imgSize2(0), oldRG(0), grad_step(0), object (obj), channel(ch), logger(ch.logger), modes(obj.modes),
      window (wind), noiseWindow(nwind) {
#ifdef USE_LUT
    static int dummy RDX_UNUSED = initSineLUT(); 
#endif
}

SubImage::~SubImage (void) {

}


void SubImage::setPatchInfo( uint32_t i, const PointI& pos, const PointF& resOffs, uint16_t patchSize, size_t bStride, uint16_t pupSz, uint16_t nM ) {

    index = i;
    initialOffset = pos;
    channelResidualOffset = resOffs;
    currentShift = 0;
    rowStride = bStride;
    imgSize = patchSize;
    imgSize2 = imgSize*imgSize;
    
    if( pupSz != pupilSize ) {
        pupilSize = pupSz;
        pupilSize2 = pupilSize*pupilSize;
        otfSize = 2*pupilSize;
        otfSize2 = otfSize*otfSize;
        phi.resize( pupilSize, pupilSize );
        PF = redux::util::rdx_get_shared<complex_t>(pupilSize2);
        OTF.init( otfSize, otfSize, REORDER_IMG|REORDER_FT|FULLCOMPLEX|NORMALIZE_FT );
        OTF.setCentered();  // OTF is nat called through "ft", so we need to set "centered" manually.
        imgFT.init( otfSize, otfSize, REORDER_FT|FULLCOMPLEX );
        vogel.resize( pupilSize, pupilSize );
        phi.zero();
        std::fill_n( PF.get(), pupilSize2, complex_t(0) );
        OTF.zero();
        imgFT.zero();
        vogel.zero();
    }
    
    nModes = nM;
    grad_step = object.myJob.graddiff_step * object.wavelength;

}


void SubImage::getWindowedImg( double* out, float* plane, ArrayStats& imgStats, bool rescaled ) const {

    const float* rowPtr = ptr(0,0,0);
    const double* winPtr = window.get();
    
    for( size_t y=0; y<imgSize; ++y ) {
        std::copy_n( rowPtr+y*rowStride, imgSize, out+y*imgSize );
    }
    
    imgStats.getStats( out, imgSize2, ST_VALUES );
    double avg = imgStats.mean;
    double newAvg(0);
    if( plane ) {
        for( size_t i=0; i<imgSize2; ++i) {
            // Note: mean(plane) = 0, so subtracting it will not affect the mean value.
            double tmp = (out[i]-avg-plane[i])*winPtr[i]+avg;
            newAvg += tmp;
            out[i] = tmp;
        }
    } else {
        transform(out, out+imgSize2, winPtr, out,
                [avg,&newAvg](const double& a, const double& b) {
                    double tmp = (a-avg)*b+avg;
                    newAvg += tmp;
                    return tmp;
                    
                });
    }
    newAvg /= imgSize2;
    if( rescaled ) {
        avg = object.objMaxMean;
    } /*else {
        avg = 1.0/newAvg;
        s.mean = newAvg;
    }*/
    imgStats.getStats( out, imgSize2, ST_VALUES|ST_RMS );     // TODO test if this is necessary or if the stats for the larger area is sufficient
    
#ifdef RDX_DO_TRANSPOSE
    transpose( out, imgSize, imgSize );                      // to match MvN
#endif

}


void SubImage::getWindowedImg( Array<double>& im, redux::util::ArrayStats& s, bool rescaled ) const {
    getWindowedImg( im.get(), object.fittedPlane.get(), s, rescaled);
}


void SubImage::initialize( Object& o, bool doReset ) {
    
    auto tmp = Solver::tmp();
    double* imgPtr = tmp->D.get();
    double* d2Ptr = tmp->D2.get();
    complex_t* oldFT = tmp->C.get();
    complex_t* c2Ptr = tmp->C2.get();
    
    if( doReset ) {
        std::fill_n( oldFT, otfSize2, complex_t(0) );
    } else {
        std::copy_n( imgFT.get(), otfSize2, oldFT );
    }
    
    getWindowedImg( imgPtr, o.fittedPlane.get(), stats, false );

    string msg;
    if( doReset ) {
        msg = "Initializing image " + to_string(o.ID) + ":" + to_string(channel.ID) + ":" + to_string(index)
            + "   mean=" + to_string (stats.mean) + " stddev=" + to_string (stats.stddev);

        transform( imgPtr, imgPtr+imgSize2, noiseWindow.get(), d2Ptr,
                [&](const double& a, const double& b) { return (a-stats.mean)*b; });
        FourierTransform& tmpFT = tmp->FT;
        tmpFT.init( d2Ptr, imgSize, imgSize, FULLCOMPLEX );
        stats.noise = tmpFT.noise(-1,-1);
        stats.noise *= channel.noiseFudge;
        double rg = stats.noise/stats.stddev;
        msg += " noise=" + to_string (stats.noise) + " rg=" + to_string(rg);
        msg += "  initial shift=" + (string) currentShift;
        LOG_TRACE << msg << ende;
        o.addRegGamma( rg );
    }
    
    if( imgSize != otfSize ) {                      // Generate an image padded to the size of the OTF.
        int offset = (otfSize - imgSize) / 2;
        std::fill_n( c2Ptr, otfSize2, complex_t(stats.mean) );
        complex_t* tmpPtr = c2Ptr + offset*(otfSize+1);
        for( uint16_t i(0); i < imgSize; ++i) {
            std::copy_n( imgPtr+i*imgSize, imgSize, tmpPtr+i*otfSize );
        }
    } else {
        d2Ptr = imgPtr;
        std::copy_n( imgPtr, otfSize2, c2Ptr );
    }

    imgFT.ft( c2Ptr );

    if( doReset ) {
        o.addToFT( imgFT.get() );
    } else {
        o.addDiffToFT( imgFT.get(), oldFT );
    }

}


void SubImage::initialize( bool doReset ) {
    initialize( object, doReset );
}


void SubImage::addFT(Array<double>& ftsum) const {
    const complex_t* ftPtr = imgFT.get();
    double* ftsPtr = ftsum.get();
    for (size_t ind = 0; ind < otfSize2; ++ind) {
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
    double dl(0.0);
    for( size_t& ind: object.pupil->otfSupport ) {
        dsj = newOTF[ind] - oldOTF[ind];            // change in sj
        dp = conj(ftPtr[ind]) * dsj;                // change p and q
        double dq = 2.0 * (oldOTF[ind].real() * dsj.real() + oldOTF[ind].imag() * dsj.imag()) + norm (dsj);
        double dn = 2.0 * (dp.real() * p[ind].real() + dp.imag() * p[ind].imag()) + norm (dp);
        dl -= (q[ind] * dn - dq * (norm (p[ind]))) / (q[ind] * (q[ind]+dq));
    }
    return dl / otfSize2;
}


double SubImage::gradientFiniteDifference( uint16_t modeIndex ) {

    if( object.weight == 0 ) return 0;
        
    complex_t* tmpOTF = Solver::tmp()->C.get();
    double* tmpPhi = Solver::tmp()->D.get();
    std::copy_n( phi.get(), pupilSize2, tmpPhi );
    addToPhi( tmpPhi, modes->modePointers[modeIndex], grad_step );
    calcOTF( tmpOTF, tmpPhi );
    return metricChange( tmpOTF )/grad_step*object.weight;
    
}


void SubImage::gradientFiniteDifference2( double* agrad, const bool* enabledModes ) {

    if( object.weight == 0 ) return;
    
    complex_t* tmpOTF = Solver::tmp()->C.get();
    double* tmpPhi = Solver::tmp()->D.get();
    for( uint16_t m=0; m<nModes; ++m ) {
        if( enabledModes[m] ) {
            std::copy_n( phi.get(), pupilSize2, tmpPhi );
            addToPhi( tmpPhi, modes->modePointers[m], grad_step );
            calcOTF( tmpOTF, tmpPhi );
            agrad[m] += metricChange( tmpOTF )/grad_step*object.weight;
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

    auto tmp = Solver::tmp();
    complex_t* tmpOtfPtr = tmp->OTF.get();
    complex_t* glPtr = tmp->C.get();
    complex_t* hjPtr = tmp->C2.get();

    tmp->OTF.zero();
    FourierTransform::reorderInto( pfPtr, pupilSize, pupilSize, tmpOtfPtr, otfSize, otfSize );

    tmp->OTF.ift(hjPtr);       // normalize by otfSize2 below
    tmp->OTF.zero();
    for( const size_t& ind: object.pupil->otfSupport ) {
        tmpOtfPtr[ind] = (pq[ind]*ftPtr[ind] - ps[ind]*otfPtr[ind]) / qs[ind];
    }

    FourierTransform::reorder( tmpOtfPtr, otfSize, otfSize );
    tmp->OTF.ift(glPtr);
    
    transform( glPtr, glPtr+otfSize2, hjPtr, glPtr,
              [](const complex_t& g,const complex_t& h) {
                  return h * g.real();
              });

    tmp->OTF.ft( glPtr );

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
    complex_t* tmpPtrFT = Solver::tmp()->FT.get();
    double* tmpPtrD = Solver::tmp()->D.get();
    std::transform( ftPtr, ftPtr+imgFT.nElements(), refPtr, tmpPtrFT,
        []( const complex_t& a, const complex_t& b ){
            return b*std::conj(a);
    });
    FourierTransform::reorder( Solver::tmp()->FT );
    Solver::tmp()->FT.ift( tmpPtrD );
    FourierTransform::reorder( tmpPtrD, imgSize, imgSize );
    
    auto arr = reshapeArray( tmpPtrD, imgSize, imgSize );
    double** arrPtr = arr.get();
    int mid = imgSize/2;
    currentShift = 0;
    double maxVal = arrPtr[mid][mid];
    for( int x=-19; x<20; ++x ) {
        for( int y=-19; y<20; ++y ) {
            if( arrPtr[mid+y][mid+x] > maxVal ) {
                maxVal = arrPtr[mid+y][mid+x];
                currentShift.x = x;
                currentShift.y = y;
            }
        }
    }
    
}

            
template <typename T>
bool SubImage::adjustShifts( const T* alpha ) {

    bool ret(false);
    
    PointI xyDims(1,2);                     // N.B. dim=0 is the index of the image-stack
#ifdef RDX_DO_TRANSPOSE
    std::swap( xyDims.x, xyDims.y );
#endif
    const PointI& tm = modes->tiltMode;
    const PointD thisAlpha = PointD( alpha[tm.y], alpha[tm.x] ) + currentTiltOffset;

    const PointI thisShift = -modes->alphaToShift( thisAlpha ).round();
    
    PointI actualShift;
    if( thisShift.x ) {
        actualShift.x = shift( xyDims.x, thisShift.x );
    }
    if( thisShift.y ) {
        actualShift.y = shift( xyDims.y, thisShift.y );
    }
    
    if( actualShift != 0 ) {
        //cout << index << "  thisShift = " << thisShift << "  actualShift = " << actualShift
        //     << "  aToS = " << modes->alphaToShift( PointD(1,1) ) << endl;
        currentShift += actualShift;
        currentTiltOffset += modes->shiftToAlpha( PointD(actualShift) );
        ret = true;
    }

    if( actualShift != thisShift ) {
        LOG_NOTICE << "SubImage " << to_string(object.ID) << ":" << to_string(channel.ID) << ":" << to_string(index)
                   << ": Patch is too close to edge, could only shift " << actualShift << "/" << thisShift << " pixels!" << ende;
    }
    return ret;
    
}
template bool SubImage::adjustShifts( const double* );
template bool SubImage::adjustShifts( const float* );


void SubImage::resetShifts( void ) {
    
    PointI xyDims(1,2);                     // N.B. dim=0 is the index of the image-stack
#ifdef RDX_DO_TRANSPOSE
    std::swap( xyDims.x, xyDims.y );
#endif
    
    shift( xyDims.x, -currentShift.x );
    shift( xyDims.y, -currentShift.y );
    currentShift = 0;
    currentTiltOffset = 0;
    
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

    for( int mi(0); mi <nModes; ++mi ) {
        double alpha = a[mi];
        if( mi == modes->tiltMode.x ) {
            alpha += currentTiltOffset.x;
        } else if( mi == modes->tiltMode.y ) {
            alpha += currentTiltOffset.y;
        }
        if( fabs(alpha) > ALPHA_CUTOFF ) {
            const double* modePtr = modes->modePointers[mi];
            transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
                [alpha]( const double& p, const double& m ) {
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

   Solver::tmp()->OTF.autocorrelate( otfPtr, true );

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

    Solver::tmp()->OTF.autocorrelate( otfPtr, true );

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

    Solver::tmp()->OTF.autocorrelate( otfPtr, true );

}


void SubImage::addPSF( double* outPSF ) const {
    
    complex_t* cPtr = Solver::tmp()->C.get();
    OTF.ift( cPtr );

    std::transform( outPSF, outPSF+otfSize2, cPtr, outPSF,
        []( const double& p, const complex_t& c ) {
            return p+std::real(c);
        });


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

    Array<double> img( imgSize, imgSize );
    Array<complex_t> tmp;
    ArrayStats s;
    getWindowedImg( img,s, true );
    Ana::write( tag + "_rimg.f0", *this );
    Ana::write( tag + "_img.f0", img );
    Ana::write( tag + "_phi.f0", phi );
    tmp.wrap( PF.get(), pupilSize, pupilSize );
    Ana::write( tag + "_pf.f0", tmp );
    Ana::write( tag + "_otf.f0", OTF );
    Ana::write( tag + "_psf.f0", getPSF() );
    Ana::write( tag + "_imgFT.f0", imgFT );
    Ana::write( tag + "_window.f0", window );
    Ana::write( tag + "_vogel.f0", vogel );


}

#undef ALPHA_CUTOFF

