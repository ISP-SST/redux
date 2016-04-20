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

//#define USE_LUT
#define ALPHA_CUTOFF 1E-12

namespace {
    
    const std::string thisChannel = "subimage";

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
    static int dummy UNUSED = initSineLUT(); 
#endif
}

SubImage::~SubImage (void) {

}


void SubImage::setPatchInfo( uint32_t i, const PointI& offs, const PointI& shift, uint16_t patchSize, uint16_t pupSz, uint16_t nM ) {

    index = i;
    offset = offs;
    offsetShift = shift;
    
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
    
    if( nM != nModes ) {
        nModes = nM;
        currentAlpha.resize(nM);
    }

}


void SubImage::init (void) {

    tmpFT.zero();
    
    copy(img);                    // copy current cut-out to (double) working copy.
    
    stats.getStats( img.get(), imgSize*imgSize, ST_VALUES|ST_RMS );     // TODO test if this is necessary or if the stats for the larger area is sufficient

    double avg = stats.mean;
    string str = "Initializing image " + to_string(object.ID) + ":" + to_string(channel.ID) + ":" + to_string(index)
    + "   mean=" + to_string (stats.mean) + " stddev=" + to_string (stats.stddev);
    
    transform(img.get(), img.get()+img.nElements(), window.get(), img.get(),
              [avg](const double& a, const double& b) { return (a-avg)*b+avg; }
             );
    
    transpose(img.get(),imgSize,imgSize);                                                       // to match MvN
    
    size_t noiseSize = noiseWindow.dimSize(0);
    if (imgSize == noiseSize) {
        img.copy(tmpImg);
        transform(tmpImg.get(), tmpImg.get()+tmpImg.nElements(), noiseWindow.get(), tmpImg.get(),
                [avg](const double& a, const double& b) { return (a-avg)*b; }
                );
        stats.getNoise(tmpImg);
    } else {
        size_t offset = (imgSize-noiseSize) / 2;
        Array<double> tmp(img, offset, offset+noiseSize-1, offset ,offset+noiseSize-1);
        tmp.trim();
        transform(tmp.get(), tmp.get()+tmp.nElements(), noiseWindow.get(), tmp.get(),
                [avg](const double& a, const double& b) { return (a-avg)*b; }
                );
        stats.getNoise(tmp);
    }
    
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


void SubImage::restore(complex_t* obj, double* obj_norm) const {

    const complex_t* ftPtr = imgFT.get();
    const complex_t* otfPtr = OTF.get();
    for ( const size_t& ind: object.pupil.otfSupport ) {
        obj_norm[ind] += norm(otfPtr[ind]);
        obj[ind] += ftPtr[ind] * conj(otfPtr[ind]);
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


double SubImage::gradientFiniteDifference( uint16_t modeIndex, double dalpha ) {
    complex_t* otfPtr = tmpC.get();
    double* phiPtr = tmpPhi.get();
    double scale = 1.0/(dalpha * object.wavelength);
    memcpy(phiPtr, phi.get(), pupilSize2*sizeof(double));
    addMode(phiPtr, modes.modePointers[modeIndex], dalpha);
    calcOTF(otfPtr, phiPtr);
    return scale*metricChange(otfPtr);
    
}


double SubImage::gradientVogel(uint16_t modeIndex, double ) const {
    double ret = 0;
    const double* modePtr = modes.modePointers[modeIndex];
    const double* vogPtr = vogel.get();
    double scale = -2.0 / (pupilSize2 * object.wavelength);
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
    tmpOTF.getIFT(hjPtr);
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
    for (auto & ind : object.pupil.pupilInOTF) {
        vogPtr[ind.first] = imag(conj(pfPtr[ind.first])*tmpOtfPtr[ind.second]);//*pupPtr[ind.first];
    }

}


void SubImage::resetPhi(void) {
    memset(phi.get(), 0, pupilSize2*sizeof (double));    // FIXME: use phi_fixed
    //memcpy( phi.get(), channel.phi_fixed.get(), channel.phi_fixed.nElements()*sizeof (double));
}

/*
void SubImage::addScaledMode( double* modePtr, double a ) {
    addMode(phi.get(), modePtr, a/object.wavelength);
}*/

void SubImage::addMode( double* phiPtr, const double* modePtr, double a ) const {

    //if(fabs(a) < ALPHA_CUTOFF) return;
    //if(index==0 && channel.ID==0 && object.ID == 0 ) {
    //    cout << "this=" << hexString(this) << "  modePtr=" << hexString(modePtr) << "  phiPtr=" << hexString(phiPtr) << "  a=" << a << endl;
    //}
        
    transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
              [a](const double& p, const double& m) {
                  return p + a*m;
              });
//     for (auto & ind : channel.pupilIndices ) {
//         phiPtr[ind] += a * modePtr[ind];
//     }

}

/*
void SubImage::addModes (double* phiPtr, size_t nModes, uint16_t* modes, const double* alphas) const {
    for (size_t m = 0; m < nModes; ++m) {
        addMode (phiPtr, modes[m], alphas[m]);
    }
}
*/

void SubImage::adjustOffset(void) {

    PointI oldOffset = offsetShift;
    PointD oldVal(0,0),newVal(0,0);
    int32_t mIndex = object.modes.tiltMode.x;  // NOTE: should be x, this is just because subimage is transposed!!
    if( mIndex >= 0 ) {
        double shiftToAlpha = object.shiftToAlpha.x;
        double alphaToShift = 1.0/shiftToAlpha;
        oldVal.x = currentAlpha[mIndex]*alphaToShift;
        int adjust = -lround(currentAlpha[mIndex]*alphaToShift);
        if( adjust ) {
            adjust = shift(2,adjust);           // will return the "actual" shift. (the cube-edge might restrict it)
            if(adjust) {
                offsetShift.x += adjust;
                currentAlpha[mIndex] += adjust*shiftToAlpha;
                newVal.x = currentAlpha[mIndex]*alphaToShift;
            }
        }
    }

    mIndex = object.modes.tiltMode.y;  // NOTE: should be y, this is just because subimage is transposed!!
    if( mIndex >= 0 ) {
        double shiftToAlpha = object.shiftToAlpha.y;
        double alphaToShift = 1.0/shiftToAlpha;
        oldVal.y = currentAlpha[mIndex]*alphaToShift;
        int adjust = -lround(currentAlpha[mIndex]*alphaToShift);
        if( adjust ) {
            adjust = shift(1,adjust);           // will return the "actual" shift. (the cube-edge might restrict it)
            if(adjust) {
                offsetShift.y += adjust;
                currentAlpha[mIndex] += adjust*shiftToAlpha;
                newVal.y = currentAlpha[mIndex]*alphaToShift;
            }
        }
    }
    
    if(oldOffset != offsetShift) {
        LOG_TRACE << "SubImage " << to_string(object.ID) << ":" << to_string(channel.ID) << ":" << to_string(index)
                  << ":  cutout was shifted, from " << oldOffset << " to " << offsetShift
                  << " oldVal=" << oldVal << "  newVal=" << newVal;
        newCutout();
    }
    
}


void SubImage::addAlpha(uint16_t m, double a) {
    
    unique_lock<mutex> lock (mtx);
    alpha[m] += a;
    
}


void SubImage::setAlpha(uint16_t m, double a) {
    
    unique_lock<mutex> lock (mtx);
    alpha[m] = a;
    
}


void SubImage::addAlphas(const double* alphas) {
    unique_lock<mutex> lock (mtx);
    for(unsigned int m=0; m<nModes; ++m) {
        alpha[m] += alphas[m];
    }
}


void SubImage::setAlphas(const double* alphas) {
    unique_lock<mutex> lock (mtx);
    for(unsigned int m=0; m<nModes; ++m) {
        alpha[m] = alphas[m];
    }
}


void SubImage::setAlphas(const std::vector<uint16_t>& modes, const double* a) {

    int cnt=0;
    for (auto & m : modes) {
        //if( m > 2 ) {
            alpha[m] = a[cnt];
        //}
        cnt++;
    }

}


void SubImage::getAlphas(float* alphas) const {
    std::copy( currentAlpha.begin(), currentAlpha.end(), alphas );
    int32_t mIndex = object.modes.tiltMode.x;  // NOTE: should be x, this is just because subimage is transposed!!
    if( mIndex >= 0 ) {
        alphas[mIndex] -= offsetShift.x*object.shiftToAlpha.x;
    }
    mIndex = object.modes.tiltMode.y;  // NOTE: should be y, this is just because subimage is transposed!!
    if( mIndex >= 0 ) {
        alphas[mIndex] -= offsetShift.y*object.shiftToAlpha.y;
    }
}


void SubImage::addPhi( const double* p, double scale ) {
    
    double* phiPtr = phi.get();
    transform( phiPtr, phiPtr+pupilSize2, p, phiPtr,
            [scale](const double& a, const double& b) {
                return a + scale*b;
            });
    
    newPhi = true;
    
}


void SubImage::calcPhi( const double* a ) {
#ifdef DEBUG_
    LOG_TRACE << "SubImage::calcPhi(" << hexString(this) << ")   nModes=" << nModes << "  pupilSize2=" << pupilSize2;
#endif

    resetPhi();

    double* phiPtr = phi.get();
    for( unsigned int i=0; i<nModes; ++i) {
        double scaledAlpha = a[i]/object.wavelength;
        const double* modePtr = modes.modePointers[i];
        transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
            [scaledAlpha](const double& p, const double& m) {
                return p + scaledAlpha*m;
            });

    }
    
    memcpy( currentAlpha.data(), a, nModes*sizeof(double) );
    
    newPhi = true;
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
    
    tmpOTF.ift(otfPtr);
    tmpOTF.autocorrelate();
    tmpOTF.getFT(otfPtr);
    FourierTransform::reorder(otfPtr, otfSize, otfSize);

}


void SubImage::calcOTF(complex_t* otfPtr, const double* phiPtr) {

    memset (otfPtr, 0, otfSize2*sizeof (complex_t));

    const double* pupilPtr = object.pupil.get();
    double normalization = sqrt(1.0 / object.pupil.area);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)

    for (auto & ind : object.pupil.pupilInOTF) {
#ifdef USE_LUT
        otfPtr[ind.second] = getPolar( pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#else
        otfPtr[ind.second] = polar(pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#endif
    }
    
    tmpOTF.ift(otfPtr);
    tmpOTF.autocorrelate();
    tmpOTF.getFT(otfPtr);
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
    double normalization = sqrt(1.0 / object.pupil.area);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)
    
    for (auto & ind : object.pupil.pupilInOTF) {
#ifdef USE_LUT
        otfPtr[ind.second] = getPolar( pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#else
        otfPtr[ind.second] = polar(pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#endif
    }

    tmpOTF.ift(otfPtr);
    tmpOTF.autocorrelate();
    tmpOTF.getFT(otfPtr);
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
    double normalization = sqrt(1.0 / object.pupil.area);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)

    for( const auto& ind: object.pupil.pupilInOTF ) {
#ifdef USE_LUT
        pfPtr[ind.first] = getPolar( pupilPtr[ind.first], phiPtr[ind.first]);
#else
        pfPtr[ind.first] = polar(pupilPtr[ind.first], phiPtr[ind.first]);
#endif
        otfPtr[ind.second] = normalization*pfPtr[ind.first];
    }

    tmpOTF.ift(otfPtr);
    tmpOTF.autocorrelate();
    tmpOTF.getFT(otfPtr);
    FourierTransform::reorder(otfPtr, otfSize, otfSize);
    //FourierTransform::autocorrelate(otfPtr, otfSize, otfSize);
    //object.addToPQ(imgFT,OTF,tmpC);
    
}

/*
void SubImage::update(bool newVogel) {
 //   cout << "SubImage::update(" << hexString(this) << ")  nV = " << newVogel << endl;
    //calcPhi();
    if( newPhi ) calcOTF();
    //if( newOTF ) object.addDiffToPQ(imgFT, OTF, oldOTF);
    if( newVogel ) calcVogelWeight(); 
//    cout << "SubImage::update(" << hexString(this) << ")   L:" << __LINE__ << endl;
}*/

void SubImage::dump (std::string tag) const {

 //   cout << "    dumping image:  this=" << hexString(this) << " with tag=" << tag << endl;
 //   cout << "                  phiPtr=" << hexString(phi.get()) << endl;
    Ana::write (tag + "_rimg.f0", *this);
    Ana::write (tag + "_img.f0", img);
    Ana::write (tag + "_phi.f0", phi);
    Ana::write (tag + "_pf.f0", PF);
    Ana::write (tag + "_otf.f0", OTF);
    Ana::write (tag + "_imgFT.f0", imgFT);
    Ana::write (tag + "_window.f0", window);
    Ana::write (tag + "_vogel.f0", vogel);

}

#undef ALPHA_CUTOFF

