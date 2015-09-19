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

namespace {
    
    const std::string thisChannel = "subimage";

#ifdef USE_LUT
#define QUADRANT_SAMPLES    1000000       // angular sampling points (per quadrant) for the sine/cosine lookup table
    const size_t total_samples(5*QUADRANT_SAMPLES);             // extend by pi/2 to fit sine/cosine in one table (cos(x) = sin(x+pi/2))
    double SineLUT[total_samples+1];                            // ...and add 1 to make room for the case of exactly 2*PI when truncating
    double* const CosineLUT = SineLUT + QUADRANT_SAMPLES;
    const double pi_x_2(2*M_PI);
    const size_t period_samples( QUADRANT_SAMPLES<<2 );
    const double angular_step( pi_x_2/(period_samples) );
    const double angular_step_inv( 1.0/angular_step );
    
    int initSineLUT (void) {
        memset(SineLUT,0,total_samples*sizeof(double));
        long double previousPhi(0), thisPhi(0);
        transform( SineLUT, SineLUT+period_samples, SineLUT,    // initialize LUT at i as average of values at i & i+1
            [&previousPhi,&thisPhi](const double&){             // so we can just truncate instead of round. 
                thisPhi += angular_step;
                long double val = 0.5L*previousPhi;
                previousPhi = sin(thisPhi);
                val += 0.5L*previousPhi;
                return val;
            }
        );
        memcpy( SineLUT+period_samples, SineLUT, (QUADRANT_SAMPLES+1)*sizeof(double));
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

SubImage::SubImage (Object& obj, const Channel& ch, const Array<double>& wind, const Array<double>& nwind, const Array<float>& stack,
                    uint32_t index, const PointI& offset, uint16_t patchSize, uint16_t pupilSize)
    : Array<float>(stack, index, index, offset.y, offset.y+patchSize-1, offset.x, offset.x+patchSize-1),
      index(index), offset(offset), offsetShift(0,0), imgSize(patchSize), pupilSize(pupilSize), otfSize(2*pupilSize), oldRG(0), object (obj),
      channel(ch), window (wind), noiseWwindow(nwind), newPhi(false), newOTF(false),
      OTF(otfSize, otfSize, FT_REORDER|FT_FULLCOMPLEX), tmpOTF(otfSize, otfSize, FT_FULLCOMPLEX) {

    tmpC.resize(otfSize, otfSize);
    tmpC2.resize(otfSize, otfSize);
    img.resize (imgSize, imgSize);
    tmpImg.resize (imgSize, imgSize);
    imgFT.resize (imgSize, imgSize);
    tmpFT.resize (imgSize, imgSize);
    phi.resize (pupilSize, pupilSize);
    tmpPhi.resize (pupilSize, pupilSize);
    PF.resize (pupilSize, pupilSize);
    vogel.resize (pupilSize, pupilSize);

    vogel.zero();
    imgFT.zero();
    resetPhi();

#ifdef USE_LUT
    static int dummy UNUSED = initSineLUT(); 
#endif
}


SubImage::~SubImage (void) {

}


void SubImage::init (void) {

    memcpy(tmpFT.get(),imgFT.get(),imgSize*imgSize*sizeof(complex_t));                          // make a temporary copy to pass to addDifftoFT below
    
    copy(img);                                                                                  // copy current cut-out to (double) working copy.
    stats.getStats(img, ST_VALUES);                                                             // get statistics before windowing
    double avg = stats.mean;
    string str = "init(" + to_string(index) + "): mean1=" + to_string (stats.mean);
    
    transform(img.get(), img.get()+img.nElements(), window.get(), img.get(),
              [avg](const double& a, const double& b) { return (a-avg)*b+avg; }
             );
    
    stats.getStats(img, ST_VALUES | ST_RMS);                                                    // get statistics after windowing
    str += " mean2=" + to_string (stats.mean) + " std=" + to_string (stats.stddev);
    avg = stats.mean;
    transpose(img.get(),imgSize,imgSize);                                                       // to match MvN
    
    size_t noiseSize = noiseWwindow.dimSize(0);
    if (imgSize == noiseSize) {
        img.copy(tmpImg);
        transform(tmpImg.get(), tmpImg.get()+tmpImg.nElements(), noiseWwindow.get(), tmpImg.get(),
                [avg](const double& a, const double& b) { return (a-avg)*b; }
                );
        stats.getNoise(tmpImg);
    } else {
        size_t offset = (imgSize-noiseSize) / 2;
        Array<double> tmp(img, offset, offset+noiseSize-1, offset ,offset+noiseSize-1);
        tmp.trim();
        transform(tmp.get(), tmp.get()+tmp.nElements(), noiseWwindow.get(), tmp.get(),
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
   
  
    FourierTransform::reorder (imgFT);                                                          // keep FT in centered form
    stats.noise *= channel.noiseFudge;
    str += " noise=" + to_string (stats.noise) + " rg=" + to_string (stats.noise/stats.stddev);
    
    object.addDiffToFT( imgFT, tmpFT, stats.noise/stats.stddev-oldRG );
    oldRG = stats.noise/stats.stddev;
    LOG_TRACE << str << (string)offsetShift;
    
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
    //for (auto & ind : channel.otfIndices) {
    size_t N = imgSize*imgSize;
    for(size_t ind=0; ind<N; ++ind) {
        Q[ind] += norm(otf[ind]);                    // Q += sj.re^2 + sj.im^2 = norm(sj)
        P[ind] += conj(ftPtr[ind]) * otf[ind];       // P += conj(ft)*sj            c.f. Vogel
        //P[ind] += ftPtr[ind] * otf[ind];       // to make P match MvN
    }

}

void SubImage::restore(complex_t* obj, double* obj_norm) const {

    const complex_t* ftPtr = imgFT.get();
    const complex_t* otfPtr = OTF.get();
    //for (auto & ind : channel.otfIndices) {
    size_t N = imgSize*imgSize;
    for(size_t ind=0; ind<N; ++ind) {
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
    //for (auto & ind : channel.otfIndices) {
    size_t N = imgSize*imgSize;
    for(size_t y=0; y<imgSize; ++y) {
        for(size_t x=0; x<imgSize; ++x) {
            size_t ind = y*imgSize+x;
            dsj = newOTF[ind] - oldOTF[ind];            // change in sj
            dp = conj (ftPtr[ind]) * dsj;               // change p and q
            //dp = ftPtr[ind] * dsj;         // to match MvN
            dq = 2.0 * (oldOTF[ind].real() * dsj.real() + oldOTF[ind].imag() * dsj.imag()) + norm (dsj);
            dn = 2.0 * (dp.real() * p[ind].real() + dp.imag() * p[ind].imag()) + norm (dp);
            dl -= (q[ind] * dn - dq * (norm (p[ind]))) / (q[ind] * (q[ind]+dq));
        }
    }
    return dl / N; //channel.otfIndices.size();
}


double SubImage::gradientFiniteDifference( uint16_t mode, double dalpha ) {
    
    complex_t* otfPtr = tmpC.get();
    double* phiPtr = tmpPhi.get();
    
    memcpy(phiPtr, phi.get(), pupilSize*pupilSize*sizeof(double));
    addMode(phiPtr, mode, dalpha);
    calcOTF(otfPtr, phiPtr);
    double ret = metricChange(otfPtr) / object.wavelength / dalpha;     // FIXME: multiply by lambda to match old MvN code.

    return ret;
}


double SubImage::gradientVogel(uint16_t mode, double dalpha) const {
    double ret = 0;
    const double* modePtr = channel.modes.at(mode)->get();
    const double* vogPtr = vogel.get();
    double scale = -2.0/object.wavelength * channel.pupil.second / (imgSize*imgSize) ;
    for (auto & ind : channel.pupilIndices) {
        ret += scale * vogPtr[ind] * modePtr[ind];
    }
    return ret;
}


void SubImage::calcVogelWeight(void) {
  
    const complex_t* pPtr = object.P.get();
    const double* qPtr = object.Q.get();
    const complex_t* otfPtr = OTF.get();
    const complex_t* ftPtr = imgFT.get();
    const complex_t* pfPtr = PF.get();
    
    complex_t* tmpOtfPtr = tmpOTF.get();
    complex_t* glPtr = tmpC.get();
    complex_t* hjPtr = tmpC2.get();

    size_t otf2 = otfSize*otfSize;
    
    tmpOTF.zero();
    FourierTransform::reorderInto( pfPtr, pupilSize, pupilSize, tmpOtfPtr, otfSize, otfSize );
    tmpOTF.ift(hjPtr);
    tmpOTF.zero();
    for (size_t ind=0; ind < otf2; ++ind) {
    //for( const auto& ind : channel.otfIndices) {
        complex_t pq = pPtr[ind] * qPtr[ind];
        double ps = norm (pPtr[ind]);
        double qs = qPtr[ind] * qPtr[ind];
        tmpOtfPtr[ind] = (pq*ftPtr[ind] - ps*otfPtr[ind]) / qs;
    }
    FourierTransform::reorder( tmpOtfPtr, otfSize, otfSize );
    tmpOTF.ift(glPtr);
    
    tmpOTF.zero();
    double normalization = 1.0/(otf2*otf2);
    transform( glPtr, glPtr+otf2, hjPtr, glPtr,
              [normalization](const complex_t& g,const complex_t& h) {
                  return normalization * h * g.real();
              });

    tmpOTF.ft( glPtr );
    FourierTransform::reorder( tmpOtfPtr, otfSize, otfSize );

    double* vogPtr = vogel.get();

    for (auto & ind : channel.pupilInOTF) {
        vogPtr[ind.first] = imag(conj(pfPtr[ind.first])*tmpOtfPtr[ind.second]);//*pupPtr[ind.first];
    }

}


void SubImage::resetPhi (Array<double>&p) const {
    memset(p.get(), 0, pupilSize*pupilSize*sizeof (double));    // FIXME: use phi_fixed
    //memcpy(p.get(), channel.phi_fixed.get(), channel.phi_fixed.nElements() *sizeof (double));
}


#define ALPHA_CUTOFF 1E-12
void SubImage::addScaledMode(double* phiPtr, uint16_t m, double a) const {
 //   if( fabs(a) > 0 ) cout << "SubImage::addScaledMode(" << hexString(this) << ")  m = " << m << "  a = " << a << "\n";
    addMode(phiPtr, m, a/object.wavelength);
}


void SubImage::addMode (double* phiPtr, uint16_t m, double a) const {

    //if(fabs(a) < ALPHA_CUTOFF) return;
    const double* modePtr = channel.modes.at(m)->get();
    for (auto & ind : channel.pupilIndices ) {
        phiPtr[ind] += a * modePtr[ind];
    }

}


void SubImage::addModes (double* phiPtr, size_t nModes, uint16_t* modes, const double* alphas) const {
    for (size_t m = 0; m < nModes; ++m) {
        addMode (phiPtr, modes[m], alphas[m]);
    }
}


void SubImage::adjustOffset(void) {
    bool adjusted(false);
    if( alpha.count(2) ) {
        int adjust = lround(alpha[2]*channel.alphaToPixels);
        if( adjust ) {
            adjust = shift(2,adjust);           // will return the "actual" shift. (the cube-edge might restrict it)
            if(adjust) {
                offsetShift.x += adjust;        // Noll-index = 2 corresponds to x-tilt.
                double oldval = alpha[2];
                alpha[2] -= adjust*channel.pixelsToAlpha;
                LOG_TRACE << "adjustOffset:  mode 2 was adjusted.  from " << oldval << " to " << alpha[2] << "  (" << adjust<< " pixels)";
                adjusted = true;
            }
        }
    }
    if( alpha.count(3) ) {
        int adjust = lround(alpha[3]*channel.alphaToPixels);
        if( adjust ) {
            adjust = shift(1,adjust);           // will return the "actual" shift. (the cube-edge might restrict it)
            if(adjust) {
                offsetShift.y += adjust;        // Noll-index = 3 corresponds to y-tilt.
                double oldval = alpha[3];
                alpha[3] -= adjust*channel.pixelsToAlpha;
                LOG_TRACE << "adjustOffset:  mode 3 was adjusted.  from " << oldval << " to " << alpha[3] << "  (" << adjust<< " pixels)";
                adjusted = true;
            }
        }
    }
    if(adjusted) init();
    
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
    if (wf) {
        size_t cnt (0);
        for (auto &it : wf->modes) {
            addAlpha(it.first, alphas[cnt++]);
        }
    }
}


void SubImage::setAlphas(const double* alphas) {
    if (wf) {
        size_t cnt (0);
        for (auto & it : wf->modes) {
            setAlpha(it.first, alphas[cnt++]);
        }
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




void SubImage::resetPhi(void) {
    resetPhi(phi);
    for( auto& it: alpha ) {
        it.second = 0;
    }
    currentAlpha = alpha;
    calcPFOTF();
    newPhi = false;
    newOTF = true;
}

void SubImage::calcPhi(void) {

    memset(phi.get(), 0, pupilSize*pupilSize*sizeof (double));
    for( auto& it: alpha ) {
        addScaledMode( it.first, it.second );  // only add difference
    }
    newPhi = true;
}


void SubImage::calcOTF(complex_t* otfPtr, const double* phiPtr) {

    memset (otfPtr, 0, otfSize * otfSize * sizeof (complex_t));

    const double* pupilPtr = channel.pupil.first.get();
    double normalization = sqrt(1.0 / channel.pupil.second);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)

    for (auto & ind : channel.pupilInOTF) {
#ifdef USE_LUT
        otfPtr[ind.second] = getPolar( pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#else
        otfPtr[ind.second] = polar(pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#endif
    }
    
    //tmpOTF.ft(otfPtr);
    //tmpOTF.autocorrelate(1.0/(otfSize*otfSize));
    //tmpOTF.ift(otfPtr);
    //FourierTransform::reorder(otfPtr, otfSize, otfSize);
    FourierTransform::autocorrelate(otfPtr, otfSize, otfSize);

    newPhi = false;
    newOTF = true;

}


void SubImage::calcOTF(void) {
    
    //if( !newPhi ) return;
   //     LOG_TRACE << "SubImage::calcOTF(" << hexString(this) << ")";
    
    complex_t* otfPtr = OTF.get();
    const double* phiPtr = phi.get();

    memset (otfPtr, 0, otfSize * otfSize * sizeof (complex_t));
    
    const double* pupilPtr = channel.pupil.first.get();
    double normalization = sqrt(1.0 / channel.pupil.second);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)

    for (auto & ind : channel.pupilInOTF) {
#ifdef USE_LUT
        otfPtr[ind.second] = getPolar( pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#else
        otfPtr[ind.second] = polar(pupilPtr[ind.first]*normalization, phiPtr[ind.first]);
#endif
    }
    
    //tmpOTF.ft(otfPtr);
    //tmpOTF.autocorrelate(1.0/(otfSize*otfSize));
    //tmpOTF.ift(otfPtr);
    //FourierTransform::reorder(otfPtr, otfSize, otfSize);
    FourierTransform::autocorrelate(otfPtr, otfSize, otfSize);

    newPhi = false;
    newOTF = true;

}


void SubImage::calcPFOTF(void) {
    
    //    LOG_TRACE << "SubImage::calcPFOTF(" << hexString(this) << ")";
    complex_t* pfPtr = PF.get();
    complex_t* otfPtr = OTF.get();
    const double* phiPtr = phi.get();
   // memcpy(tmpC.get(),otfPtr,otfSize*otfSize*sizeof(complex_t));
    
    memset (pfPtr, 0, pupilSize * pupilSize * sizeof (complex_t));
    memset (otfPtr, 0, otfSize * otfSize * sizeof (complex_t));
    
    const double* pupilPtr = channel.pupil.first.get();
    double normalization = sqrt(1.0 / channel.pupil.second);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)

    for (auto & ind : channel.pupilInOTF) {
#ifdef USE_LUT
        pfPtr[ind.first] = getPolar( pupilPtr[ind.first], phiPtr[ind.first]);
#else
        pfPtr[ind.first] = polar(pupilPtr[ind.first], phiPtr[ind.first]);
#endif
        otfPtr[ind.second] = normalization*pfPtr[ind.first];
    }

    tmpOTF.ft(otfPtr);
    tmpOTF.autocorrelate(1.0/(otfSize*otfSize));
    tmpOTF.ift(otfPtr);
    FourierTransform::reorder(otfPtr, otfSize, otfSize);
    //FourierTransform::autocorrelate(otfPtr, otfSize, otfSize);
    //object.addToPQ(imgFT,OTF,tmpC);
    
}


void SubImage::update(bool newVogel) {
 //   cout << "SubImage::update(" << hexString(this) << ")  nV = " << newVogel << endl;
    calcPhi();
    if( newPhi ) calcOTF();
    //if( newOTF ) object.addDiffToPQ(imgFT, OTF, oldOTF);
    if( newVogel ) calcVogelWeight(); 
//    cout << "SubImage::update(" << hexString(this) << ")   L:" << __LINE__ << endl;
}

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
