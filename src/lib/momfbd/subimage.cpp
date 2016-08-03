#include "redux/momfbd/subimage.hpp"

#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/object.hpp"

#include "redux/file/fileana.hpp"
#include "redux/logger.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <algorithm>

using namespace redux::file;
using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
#define logChannel "subimage"

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
    : imgSize(0), pupilSize(0), nModes(0), oldRG(0), object (obj), channel(ch), modes(obj.modes),
      window (wind), noiseWindow(nwind) {

#ifdef USE_LUT
    static int dummy RDX_UNUSED = initSineLUT(); 
#endif
}

SubImage::~SubImage (void) {

}


void SubImage::setPatchInfo( uint32_t i, const PointI& offs, const PointI& shift, uint16_t patchSize, uint16_t pupSz, uint16_t nM ) {

    index = i;
    offset = offs;
    offsetShift = shift;
    adjustedTilts = 0;
    
    if( patchSize != imgSize ) {
        imgSize = patchSize;
        img.resize( imgSize, imgSize );
        tmpImg.resize( imgSize, imgSize );
        imgFT.resize( imgSize, imgSize );
        tmpFT.resize( imgSize, imgSize );
        imgFT.zero();
    }
    
    if( pupSz != pupilSize ) {
        pupilSize = pupSz;
        pupilSize2 = pupilSize*pupilSize;
        otfSize = 2*pupilSize;
        otfSize2 = otfSize*otfSize;
        phi.resize( pupilSize, pupilSize );
        tmpPhi.resize( pupilSize, pupilSize );
        PF.resize( pupilSize, pupilSize );
        vogel.resize( pupilSize, pupilSize );
        OTF.resize( otfSize, otfSize, FT_REORDER|FT_FULLCOMPLEX );
        tmpOTF.resize( otfSize, otfSize, FT_FULLCOMPLEX );
        tmpC.resize( otfSize, otfSize );
        tmpC2.resize( otfSize, otfSize );
        vogel.zero();
    }
    
    nModes = nM;
    
}


void SubImage::init (void) {

    tmpFT.zero();
    
    copy(img);                    // copy current cut-out to (double) working copy.
    
    stats.getStats( img.get(), imgSize*imgSize, ST_VALUES );

    transform(img.get(), img.get()+img.nElements(), window.get(), img.get(),
              [&](const double& a, const double& b) { return (a-stats.mean)*b+stats.mean; }
             );
    stats.getStats( img.get(), imgSize*imgSize, ST_VALUES|ST_RMS );     // TODO test if this is necessary or if the stats for the larger area is sufficient
    string str = "Initializing image " + to_string(object.ID) + ":" + to_string(channel.ID) + ":" + to_string(index)
    + "   mean=" + to_string (stats.mean) + " stddev=" + to_string (stats.stddev);
    
    transpose(img.get(),imgSize,imgSize);                                                       // to match MvN
    
    size_t noiseSize = noiseWindow.dimSize(0);
    if (imgSize == noiseSize) {
        img.copy(tmpImg);
        transform(tmpImg.get(), tmpImg.get()+tmpImg.nElements(), noiseWindow.get(), tmpImg.get(),
                [&](const double& a, const double& b) { return (a-stats.mean)*b; }
                );
        imgFT.reset(tmpImg.get(), imgSize, imgSize, FT_FULLCOMPLEX);
        stats.noise = imgFT.noise(-1,-1);
    } else {
        size_t offset = (imgSize-noiseSize) / 2;
        Array<double> tmp(img, offset, offset+noiseSize-1, offset ,offset+noiseSize-1);
        tmp.trim();
        transform(tmp.get(), tmp.get()+tmp.nElements(), noiseWindow.get(), tmp.get(),
                [&](const double& a, const double& b) { return (a-stats.mean)*b; }
                );
        imgFT.reset(tmp.get(), imgSize, imgSize, FT_FULLCOMPLEX);
        stats.noise = imgFT.noise(-1,-1);
    }

    if (imgSize == otfSize) {                                                                   // imgSize = 2*pupilSize
        imgFT.reset(img.get(), imgSize, imgSize, FT_FULLCOMPLEX); //|FT_NORMALIZE );                 // full-complex for now, perhaps half-complex later for performance
    } else {                                                                                    // imgSize > 2*pupilSize should never happen (cf. calculatePupilSize)
        int offset = (otfSize - imgSize) / 2;
        Array<double> tmp (otfSize, otfSize);
        tmp.zero();
        double* imgPtr = img.get();
        double* tmpPtr = tmp.get();
        for (int i = 0; i < imgSize; ++i) {
            memcpy (tmpPtr + offset * (otfSize + 1) + i * otfSize, imgPtr + i * imgSize, imgSize * sizeof (double));
        }
        imgFT.reset (tmp.get(), otfSize, otfSize, FT_FULLCOMPLEX); //|FT_NORMALIZE );               // full-complex for now, perhaps half-complex later for performance
    }

    FourierTransform::reorder(imgFT);                                                          // keep FT in centered form
    stats.noise *= channel.noiseFudge;
    double rg = stats.noise/stats.stddev;
    str += " noise=" + to_string (stats.noise) + " rg=" + to_string(rg);
    object.addRegGamma( rg );
    object.addToFT( imgFT );

    LOG_TRACE << str << "  initial shift=" << (string)offsetShift;
    
}


void SubImage::newCutout(void) {

    memcpy(tmpFT.get(),imgFT.get(),imgSize*imgSize*sizeof(complex_t));                          // make a temporary copy to pass to addDifftoFT below
    
    copy(img);                                                                                  // copy current cut-out to (double) working copy.

    stats.getStats( img.get(), imgSize*imgSize, ST_VALUES );     // TODO test if this is necessary or if the stats for the larger area is sufficient

    double avg = stats.mean;
    transform(img.get(), img.get()+img.nElements(), window.get(), img.get(),
              [avg](const double& a, const double& b) { return (a-avg)*b+avg; }
             );
    
    transpose(img.get(),imgSize,imgSize);                                                       // to match MvN
    
    if (imgSize == otfSize) {                                                                   // imgSize = 2*pupilSize
        imgFT.reset(img.get(), imgSize, imgSize, FT_FULLCOMPLEX );                              // full-complex for now, perhaps half-complex later for performance
    } else {                                                                                    // imgSize > 2*pupilSize should never happen (cf. calculatePupilSize)
        int offset = (otfSize - imgSize) / 2;
        Array<double> tmp (otfSize, otfSize);
        tmp.zero();
        double* imgPtr = img.get();
        double* tmpPtr = tmp.get();
        for (int i = 0; i < imgSize; ++i) {
            memcpy (tmpPtr + offset * (otfSize + 1) + i * otfSize, imgPtr + i * imgSize, imgSize * sizeof (double));
        }
        imgFT.reset (tmp.get(), otfSize, otfSize, FT_FULLCOMPLEX );                             // full-complex for now, perhaps half-complex later for performance
    }
   
  
    FourierTransform::reorder(imgFT);                                                           // keep FT in centered form
    object.addDiffToFT( imgFT, tmpFT );
    
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
    for( const size_t& ind: object.pupil.otfSupport ) {
        Q[ind] += norm(otf[ind]);                    // Q += sj.re^2 + sj.im^2 = norm(sj)
        P[ind] += conj(ftPtr[ind]) * otf[ind];       // P += conj(ft)*sj            c.f. Vogel
    }

}


void SubImage::addToPQ(void) const {
    object.addToPQ( imgFT.get(), OTF.get() );
}


void SubImage::restore( complex_t* obj, double* obj_norm ) const {

    bool no_restore(false);         // TODO: implement NO_RESTORE cfg flag
    if( !no_restore ) {
        addPQ( OTF.get(), obj, obj_norm );      // Note: should really return obj = ft*conj(sj) for the deconvolution,
                                                // and addPQ returns the conjugate.
    }

}


double SubImage::metricChange(const complex_t* newOTF) const {

    const complex_t* oldOTF = OTF.get();
    const complex_t* p = object.P.get();
    const complex_t* ftPtr = imgFT.get();
    const double* q = object.Q.get();
    
    complex_t dp, dsj;
    double dl = 0.0, dq, dn;
    for( size_t& ind: object.pupil.otfSupport ) {
        dsj = newOTF[ind] - oldOTF[ind];            // change in sj
        dp = conj(ftPtr[ind]) * dsj;                // change p and q
        dq = 2.0 * (oldOTF[ind].real() * dsj.real() + oldOTF[ind].imag() * dsj.imag()) + norm (dsj);
        dn = 2.0 * (dp.real() * p[ind].real() + dp.imag() * p[ind].imag()) + norm (dp);
        dl -= (q[ind] * dn - dq * (norm (p[ind]))) / (q[ind] * (q[ind]+dq));
    }
    return dl / otfSize2;
}


double SubImage::gradientFiniteDifference( uint16_t modeIndex, double step ) {
    
    complex_t* otfPtr = tmpC.get();
    double* phiPtr = tmpPhi.get();
    memcpy(phiPtr, phi.get(), pupilSize2*sizeof(double));
    addToPhi(phiPtr, modes.modePointers[modeIndex], step );
    calcOTF(otfPtr, phiPtr);
    return metricChange(otfPtr)/step;
    
}


double SubImage::gradientVogel(uint16_t modeIndex, double ) const {

    double ret = 0;
    const double* modePtr = modes.modePointers[modeIndex];
    const double* vogPtr = vogel.get();
    double scale = -2.0 * object.pupil.area / otfSize2;
    for (auto & ind : object.pupil.pupilSupport) {
        ret += scale * vogPtr[ind] * modePtr[ind];
    }
    
    return ret;
    
}


void SubImage::calcVogelWeight(void) {

#ifdef DEBUG_
    LOG_TRACE << "SubImage::calcVogelWeight(" << hexString(this) << ")   indexSize=" << object.pupil.pupilInOTF.size();
#endif
  
    const complex_t* pPtr = object.P.get();
    const double* qPtr = object.Q.get();
    const complex_t* otfPtr = OTF.get();
    const complex_t* ftPtr = imgFT.get();
    const complex_t* pfPtr = PF.get();
    
    complex_t* tmpOtfPtr = tmpOTF.get();
    complex_t* glPtr = tmpC.get();
    complex_t* hjPtr = tmpC2.get();
    
    tmpOTF.zero();
    FourierTransform::reorderInto( pfPtr, pupilSize, pupilSize, tmpOtfPtr, otfSize, otfSize );

    tmpOTF.getIFT(hjPtr);       // normalize by otfSize2 below
    tmpOTF.zero();
    for( const size_t& ind: object.pupil.otfSupport ) {
        complex_t pq = pPtr[ind] * qPtr[ind];
        double ps = norm (pPtr[ind]);
        double qs = qPtr[ind] * qPtr[ind];
        tmpOtfPtr[ind] = (pq*ftPtr[ind] - ps*otfPtr[ind]) / qs;
    }

    FourierTransform::reorder( tmpOtfPtr, otfSize, otfSize );
    tmpOTF.getIFT(glPtr);
    
    double normalization = 1.0/(otfSize2*otfSize2);
    transform( glPtr, glPtr+otfSize2, hjPtr, glPtr,
              [normalization](const complex_t& g,const complex_t& h) {
                  return normalization * h * g.real();
              });

    tmpOTF.ft( glPtr );

    FourierTransform::reorder( tmpOtfPtr, otfSize, otfSize );

    double* vogPtr = vogel.get();
    const double* pupilPtr = object.pupil.get();
    for( auto & ind : object.pupil.pupilInOTF ) {
        vogPtr[ind.first] = imag(conj(pfPtr[ind.first])*tmpOtfPtr[ind.second])*pupilPtr[ind.first];
    }


}


void SubImage::resetPhi(void) {
    memset(phi.get(), 0, pupilSize2*sizeof (double));    // FIXME: use phi_fixed
    //memcpy( phi.get(), channel.phi_fixed.get(), channel.phi_fixed.nElements()*sizeof (double));
}


void SubImage::adjustOffset( double* alpha ) {

    PointI oldOffset = offsetShift;
    PointD oldVal(0,0),newVal(0,0);
    
    int32_t mIndex = object.modes.tiltMode.x;
    if( mIndex >= 0 ) {
        double shiftToAlpha = object.shiftToAlpha.x;
        double alphaToShift = 1.0/shiftToAlpha;
        newVal.x = oldVal.x = alpha[mIndex];
        int adjust = -lround( oldVal.x*alphaToShift );
        if( adjust ) {  // FIXME: should be 2, this is just because subimage is transposed!!
            adjust = shift(1,adjust);           // will return the "actual" shift. (the cube-edge might restrict it)
            if( adjust ) {
                offsetShift.x += adjust;
                alpha[mIndex] += adjust*shiftToAlpha;
                adjustedTilts.x = offsetShift.x*shiftToAlpha;
                newVal.x = alpha[mIndex];
            }
        }
    }

    mIndex = object.modes.tiltMode.y;
    if( mIndex >= 0 ) {
        double shiftToAlpha = object.shiftToAlpha.y;
        double alphaToShift = 1.0/shiftToAlpha;
        newVal.y = oldVal.y = alpha[mIndex];
        int adjust = -lround( oldVal.y*alphaToShift );
        if( adjust ) {  // FIXME: should be 1, this is just because subimage is transposed!!
            adjust = shift(2,adjust);           // will return the "actual" shift. (the cube-edge might restrict it)
            if(adjust) {
                offsetShift.y += adjust;
                alpha[mIndex] += adjust*shiftToAlpha;
                adjustedTilts.y = offsetShift.y*shiftToAlpha;
                newVal.y = alpha[mIndex];
            }
        }
    }

    if( oldOffset != offsetShift ) {
        //LOG_DEBUG << "SubImage Shifting:  pix2cf=" << object.shiftToAlpha;
        LOG_DEBUG << "SubImage " << to_string(object.ID) << ":" << to_string(channel.ID) << ":" << to_string(index)
                   << ":  cutout was shifted, from " << oldOffset << " to " << offsetShift
                   << " oldVal=" << oldVal << "  newVal=" << newVal << "  adj=" << adjustedTilts;
        newCutout();
        //LOG_TRACE << "SubImage Shifting:  " << printArray(first(),"\nfirst") << printArray(last(),"\nlast");
    }
    
}


void SubImage::addAlphaOffsets(double* alphas, float* alphaOut) const {
    
    std::copy( alphas, alphas+nModes, alphaOut );
    
    int32_t mIndex = object.modes.tiltMode.x;
    if( mIndex >= 0 && offsetShift.x ) {
        alphaOut[mIndex] -= adjustedTilts.x;
    }
    
    mIndex = object.modes.tiltMode.y;
    if( mIndex >= 0 && offsetShift.y ) {
        alphaOut[mIndex] -= adjustedTilts.y;
    }
    
}


void SubImage::addToPhi( double* phiPtr, const double* modePtr, double a ) const {

    transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
              [a](const double& p, const double& m) {
                  return p + a*m;
              });

}


void SubImage::addToPhi( const double* a, double* phiPtr ) const {
    
#ifdef DEBUG_
    LOG_TRACE << "SubImage(" << index << ")::addToPhi(" << hexString(this) << ")   nModes=" << nModes << "  pupilSize2=" << pupilSize2
    << printArray( a, nModes, "  newAlpha" );
#endif

    for( unsigned int i=0; i<nModes; ++i ) {
        if( fabs(a[i]) > ALPHA_CUTOFF ) {
            double scaledAlpha = a[i]; ///object.wavelength;
            const double* modePtr = modes.modePointers[i];
            transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
                [scaledAlpha](const double& p, const double& m) {
                    return p + scaledAlpha*m;
                });
        }
    }
    
}


void SubImage::calcPhi( const double* a, double* phiPtr ) const {

    memcpy( phiPtr, channel.phi_channel.get(), pupilSize2*sizeof(double) );
    //memset( phiPtr, 0, pupilSize2*sizeof(double) );    // FIXME: use phi_fixed
    addToPhi( a, phiPtr );

}


void SubImage::calcOTF(complex_t* otfPtr, const double* phiOffset, double scale) {

    const double* phiPtr = phi.get();
    const double* pupilPtr = object.pupil.get();
    
    double normalization = sqrt(1.0 / object.pupil.area);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)

    for (auto & ind : object.pupil.pupilInOTF) {
#ifdef USE_LUT
        otfPtr[ind.second] = getPolar( pupilPtr[ind.first]*normalization, phiPtr[ind.first]+scale*phiOffset[ind.first]);
#else
        otfPtr[ind.second] = polar(pupilPtr[ind.first]*normalization, phiPtr[ind.first]+scale*phiOffset[ind.first]);
#endif
    }
    
    tmpOTF.ft(otfPtr);
    tmpOTF.norm();
    tmpOTF.getIFT(otfPtr);
    FourierTransform::reorder(otfPtr, otfSize, otfSize);

}


void SubImage::calcOTF(complex_t* otfPtr, const double* phiPtr) {

    memset (otfPtr, 0, otfSize2*sizeof (complex_t));

    const double* pupilPtr = object.pupil.get();
    double normalization = sqrt(1.0 / object.pupil.area / otfSize2);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)
    //double normalization = 1.0; // / object.pupil.area;   // normalize OTF by pupil-area

    for (auto & ind : object.pupil.pupilInOTF) {
#ifdef USE_LUT
        otfPtr[ind.second] = getPolar( pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#else
        otfPtr[ind.second] = polar(pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#endif
    }
    
    tmpOTF.ft(otfPtr);
    tmpOTF.norm();
    tmpOTF.getIFT(otfPtr);
    FourierTransform::reorder(otfPtr, otfSize, otfSize);
    //FourierTransform::autocorrelate(otfPtr, otfSize, otfSize);

    newPhi = false;
    newOTF = true;

}


void SubImage::calcOTF(void) {
    
    //if( !newPhi ) return;
   //     LOG_TRACE << "SubImage::calcOTF(" << hexString(this) << ")";
    
    complex_t* otfPtr = OTF.get();
    const double* phiPtr = phi.get();

    memset(otfPtr, 0, otfSize2*sizeof (complex_t));
    
    const double* pupilPtr = object.pupil.get();
    double normalization = sqrt(1.0 / object.pupil.area / otfSize2);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)
    //double normalization = 1.0; // / object.pupil.area;   // normalize OTF by pupil-area
    
    for (auto & ind : object.pupil.pupilInOTF) {
#ifdef USE_LUT
        otfPtr[ind.second] = getPolar( pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#else
        otfPtr[ind.second] = polar(pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#endif
    }

    tmpOTF.ft(otfPtr);
    tmpOTF.norm();
    tmpOTF.getIFT(otfPtr);
    FourierTransform::reorder(otfPtr, otfSize, otfSize);
    //FourierTransform::autocorrelate(otfPtr, otfSize, otfSize);

    newPhi = false;
    newOTF = true;

}


void SubImage::calcPFOTF(void) {
    
#ifdef DEBUG_
    LOG_TRACE << "SubImage::calcPFOTF(" << hexString(this) << ")   indexSize=" << object.pupil.pupilInOTF.size();
#endif

    complex_t* pfPtr = PF.get();
    complex_t* otfPtr = OTF.get();
    const double* phiPtr = phi.get();

    memset (pfPtr, 0, pupilSize2*sizeof (complex_t));
    memset (otfPtr, 0, otfSize2*sizeof (complex_t));
    
    const double* pupilPtr = object.pupil.get();
    double normalization = sqrt(1.0 / object.pupil.area / otfSize2);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)
    //double normalization = 1.0; // / object.pupil.area;   // normalize OTF by pupil-area

    for( const auto& ind: object.pupil.pupilInOTF ) {
#ifdef USE_LUT
        pfPtr[ind.first] = getPolar( pupilPtr[ind.first], phiPtr[ind.first]);
#else
        pfPtr[ind.first] = polar(pupilPtr[ind.first], phiPtr[ind.first]);
#endif
        otfPtr[ind.second] = normalization*pfPtr[ind.first];
    }


    tmpOTF.ft(otfPtr);
    tmpOTF.norm();
    tmpOTF.getIFT(otfPtr);
    FourierTransform::reorder(otfPtr, otfSize, otfSize);
    
}


void SubImage::addPSF( double* outPSF ) const {
    
    redux::util::Array<complex_t> tmpC( otfSize, otfSize );
    OTF.getIFT( tmpC.get() );     // FIXME: Using forward transform to match MvN (who conjugates the transform)
    FourierTransform::reorder( tmpC.get(), otfSize, otfSize );
//    double normalization(0.0);        // TBD: should the PSF be auto-scaled or filtered ??
//     std::for_each( tmp.get(), tmp.get()+otfSize2,
//         [&normalization]( const complex_t& c ) {
//             normalization += std::abs(std::real(c));
//         });
    double normalization = 1.0 / otfSize2;
    if( normalization > 0.0 ) {
//         normalization = 1.0/normalization;
        std::transform( outPSF, outPSF+otfSize2, tmpC.get(), outPSF,
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
    return std::move(tmp);
}



void SubImage::dump (std::string tag) const {

 //   cout << "    dumping image:  this=" << hexString(this) << " with tag=" << tag << endl;
 //   cout << "                  phiPtr=" << hexString(phi.get()) << endl;
    Ana::write (tag + "_rimg.f0", *this);
    Ana::write (tag + "_img.f0", img);
    Ana::write (tag + "_phi.f0", phi);
    Ana::write (tag + "_pf.f0", PF);
    Ana::write (tag + "_otf.f0", OTF);
    Ana::write (tag + "_psf.f0", getPSF());
    Ana::write (tag + "_imgFT.f0", imgFT);
    Ana::write (tag + "_window.f0", window);
    Ana::write (tag + "_vogel.f0", vogel);
    
    Array<double> tmp( imgSize, imgSize );
    tmp = 0;
    for(int x=0;x<imgSize/2;++x)
      for(int y=0;y<imgSize/2;++y) {
          tmp(x,y) = 1;
          tmp(imgSize/2+x,y) = 2;
          tmp(x,imgSize/2+y) = 3;
          tmp(imgSize/2+x,imgSize/2+y) = 4;
      }
    Ana::write (tag + "_g.f0", tmp);
    FourierTransform::reorder( tmp.get(), imgSize, imgSize );
    Ana::write (tag + "_g2.f0", tmp);

}

#undef ALPHA_CUTOFF

