#include "redux/momfbd/subimage.hpp"

#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/object.hpp"

#include "redux/file/fileana.hpp"
#include "redux/logger.hpp"
#include "redux/util/datautil.hpp"

#include <algorithm>

using namespace redux::file;

using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

#define ANGULAR_SAMPLES    1000000                          // angular sampling points for sine/cosine lookup table
#define lg Logger::mlg
namespace {

    const std::string thisChannel = "subimage";

    const uint32_t total_samples(5*ANGULAR_SAMPLES/4);     // extend by pi/2 to fit sine/cosine in one table (cos(x) = sin(x+pi/2))
    double SineLUT[total_samples];
    double* const CosineLUT = SineLUT + ANGULAR_SAMPLES/4;
    const double pi_x_2(2*M_PI);
    const double angular_step(2 * M_PI/(ANGULAR_SAMPLES-1));
    const double angular_step_inv(1.0/angular_step);
    int initSineLUT (void) {
        for (uint i = 0; i < total_samples; ++i) {
            SineLUT[i] = sin(i * angular_step);
        }
    Array<double> slaskD(SineLUT, total_samples);
    Ana::write("sine_lut.f0", slaskD);
    slaskD.wrap(CosineLUT,ANGULAR_SAMPLES);
    Ana::write("cosine_lut.f0", slaskD);
        return 0;
    }

}


SubImage::SubImage (Object& obj, const Channel& ch, const Array<double>& wind, const Array<double>& nwind, const Array<float>& stack,
                    uint32_t index, const PointI& offset, uint16_t patchSize, uint16_t pupilSize)
    : Array<float>(stack, index, index, offset.y, offset.y+patchSize-1, offset.x, offset.x+patchSize-1),
      index(index), offset(offset), offsetShift(0,0), imgSize(patchSize), pupilSize(pupilSize), otfSize(2*pupilSize), oldRG(0), object (obj),
      channel(ch), window (wind), noiseWwindow(nwind), newPhi(false), newOTF(false), OTF(otfSize, otfSize, FT_REORDER|FT_FULLCOMPLEX) {

    img.resize (imgSize, imgSize);
    tmpImg.resize (imgSize, imgSize);
    imgFT.resize (imgSize, imgSize);
    tmpFT.resize (imgSize, imgSize);
    phi.resize (pupilSize, pupilSize);
    tmpPhi.resize (pupilSize, pupilSize);
    PF.resize (pupilSize, pupilSize);
    //OTF.resize (otfSize, otfSize);                  // big enough for autocorrelation of the PF
    tmpOTF.resize (otfSize, otfSize);
    vogel.resize (pupilSize, pupilSize);

    vogel.zero();
    imgFT.zero();
    resetPhi();

    static int dummy UNUSED = initSineLUT(); 
}


SubImage::~SubImage (void) {

}


void SubImage::init (void) {
    
    
   // stats.getStats (*this, ST_VALUES );               // get statistics

    string str = "SubImage::init(" + to_string (index) + "):  mean=" + to_string (stats.mean);

    auto wit = window.begin();
    auto nwit = noiseWwindow.begin();
    Array<float>::const_iterator dit = this->begin();
    size_t noiseSize = noiseWwindow.dimSize(0);

    for (auto & iit : img) {                    // windowing: subtract and re-add mean afterwards
        iit = static_cast<int16_t>(*dit++);     // FIXME: this cast is just to mimic old momfbd code for testing.
    }

     stats.getStats(img, ST_VALUES);               // get statistics before windowing
     double avg = stats.mean;
    dit = this->begin();
    for (auto & iit : img) {                    // windowing: subtract and re-add mean afterwards
        iit = (static_cast<int16_t>(*dit++) - avg) * *wit++ + avg;
    }
    stats.getStats(img, ST_VALUES | ST_RMS);               // get statistics after windowing

    transpose(img.get(),imgSize,imgSize);     // to match MvN

 //   Ana::write("img_windowed"+to_string(cnt)+".f0", img);
    if (imgSize == noiseSize) {                                                       // imgSize = 2*pupilSize
        Array<double> tmp = img.copy();
        tmp -= stats.mean;
        tmp *= noiseWwindow;
        //imgFT.reset(tmp.get(), imgSize, imgSize, FT_FULLCOMPLEX|FT_NORMALIZE );   // full-complex for now, perhaps half-complex later for performance
  //  Ana::write("img_noise"+to_string(cnt)+".f0", tmp);

        stats.getStats(tmp,ST_NOISE);
    } else {                                                                        // imgSize > 2*pupilSize should never happen (cf. calculatePupilSize)
        size_t offset = (imgSize-noiseSize) / 2;
        Array<double> tmp(img, offset, offset+noiseSize-1, offset ,offset+noiseSize-1);
        tmp.trim();
        tmp -= stats.mean;
        tmp *= noiseWwindow;
        //imgFT.reset (tmp.get(), noiseSize, noiseSize, FT_FULLCOMPLEX|FT_NORMALIZE );   // full-complex for now, perhaps half-complex later for performance
        stats.getStats(tmp,ST_NOISE);
    }

   // FourierTransform::reorder(imgFT);                             // keep FT in centered form
   // Ana::write("img_noiseft"+to_string(cnt)+".f0", imgFT);
   // stats.noise = channel.noiseFudge * imgFT.noise (-1, -1);       // mask/cutoff < 0 will revert to hardcoded values used by MvN

//     dit = this->begin();
//     for (auto & iit : img) {                    // windowing: subtract and re-add mean afterwards
//         iit = (*dit++ - avg) * *wit++ + avg;
//     }
//     stats.getStats (img, ST_VALUES | ST_RMS);               // get statistics after windowing

    str += "  avg=" + to_string (avg) + "  mean2=" + to_string (stats.mean) + "  std=" + to_string (stats.stddev) + "  sum=" + to_string (stats.sum);
    str += "  noise=" + to_string (stats.noise) + "  rg=" + to_string (stats.noise/stats.stddev);
    if (imgSize == otfSize) {                                                       // imgSize = 2*pupilSize
        imgFT.reset (img.get(), imgSize, imgSize, FT_FULLCOMPLEX );   // full-complex for now, perhaps half-complex later for performance
    } else {                                                                        // imgSize > 2*pupilSize should never happen (cf. calculatePupilSize)
        int offset = (otfSize - imgSize) / 2;
        Array<double> tmp (otfSize, otfSize);
        tmp.zero();
        double* imgPtr = img.get();
        double* tmpPtr = tmp.get();
        for (int i = 0; i < imgSize; ++i) {
            memcpy (tmpPtr + offset * (otfSize + 1) + i * otfSize, imgPtr + i * imgSize, imgSize * sizeof (double));
        }
        imgFT.reset (tmp.get(), otfSize, otfSize, FT_FULLCOMPLEX );   // full-complex for now, perhaps half-complex later for performance
    }

    FourierTransform::reorder (imgFT);                             // keep FT in centered form
    stats.noise *= channel.noiseFudge;
//    Ana::write("img_ft"+to_string(cnt)+".f0", imgFT);

   // str += "  noise=" + to_string (stats.noise) + "  rg=" + to_string (stats.noise/stats.stddev);
    object.addToFT( imgFT, stats.noise/stats.stddev );
vogel.zero();
    
resetPhi();

 //   cout << str << printArray (dimensions(), "  dims") << endl;

}


void SubImage::addFT(Array<double>& ftsum) const {
    const complex_t* ftPtr = imgFT.get();
    double* ftsPtr = ftsum.get();
    for (size_t ind = 0; ind < imgFT.nElements(); ++ind) {
        ftsPtr[ind] += norm (ftPtr[ind]);
    }
}


void SubImage::addPQ (const complex_t* otf, complex_t* P, double* Q) const {

  //  cout << "-" << flush;
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


double SubImage::gradientFiniteDifference(uint16_t mode, double dalpha, complex_t* myOtf, double* myPhi) const {
    memcpy(myPhi, phi.get(), pupilSize*pupilSize*sizeof(double));
    //cout << "d" << flush;
    //memset(myPhi, 0, pupilSize*pupilSize*sizeof(double));
    addMode(myPhi, mode, dalpha);
    calcOTF(myOtf, myPhi);
    double ret = metricChange(myOtf) / object.wavelength / dalpha;     // FIXME: multiply by lambda to match old MvN code.
//    Array<complex_t> cwrapper (tmpOTF, otfSize, otfSize);
//Ana::write("gm_"+to_string(mode)+"_tmpotf.f0", cwrapper);
//dump("gm"+to_string(mode));
// cout << "SubImage::gradientFiniteDifference  dalpha = " <<  dalpha << "  mode = " << mode << "  part_der = " << ret << endl;
    return ret;
}


double SubImage::gradientVogel(uint16_t mode, double dalpha, complex_t* tmpOTF, double* tmpPhi) const {
    double ret = 0;
    //cout << "v" << flush;
    const double* modePtr = channel.modes.at(mode)->get();
    const double* vogPtr = vogel.get();
    double scale = -2.0/object.wavelength * channel.pupil.second / (imgSize*imgSize) ;
    for (auto & ind : channel.pupilIndices) {
        ret += scale * vogPtr[ind] * modePtr[ind];
    }
    //cout << "SubImage::gradientVogel(" << hexString(this) << ")   g = " << ret << endl;
    return ret;
}


void SubImage::calcVogelWeight(void) {
    
   // cout << "SubImage::calcVogelWeight(" << hexString(this) << ")"  << endl;
    //cout << "V" << flush;
    int offset = max ( (otfSize - pupilSize) / 2, 0);

    FourierTransform pj (otfSize, otfSize, FT_FULLCOMPLEX | FT_REORDER );    // input is already centered, so set the flag
    FourierTransform glFT (otfSize, otfSize, FT_FULLCOMPLEX | FT_REORDER);
    Array<complex_t> hj (otfSize, otfSize);
    Array<complex_t> gl (otfSize, otfSize);
    pj.zero();

    const complex_t* pfPtr = PF.get();
    complex_t* pjPtr = pj.ptr (offset, offset);
    for (uint i = 0; i < pupilSize; ++i) {
        memcpy (pjPtr + i * otfSize, pfPtr + i * pupilSize, pupilSize * sizeof (complex_t));
    }
//   transform(pjPtr,pjPtr+otfSize*otfSize, pjPtr,
//             [](complex_t a){
//                 return conj(a);
//             });

//  Ana::write("vogel_pj.f0", pj);
    pj.directInverse(hj);      // note: this will first reorder pj
//  Ana::write("vogel_hj.f0", pj);
//  Ana::write("vogel_hj2.f0", hj);

    const complex_t* pPtr = object.P.get();
    const double* qPtr = object.Q.get();
    const complex_t* otfPtr = OTF.get();
    const complex_t* ftPtr = imgFT.get();
    complex_t* glFTPtr = glFT.get();
    complex_t* glPtr = gl.get();
    complex_t* hjPtr = hj.get();

    glFT.zero();
    size_t N = otfSize*otfSize;
    for (size_t ind=0; ind < N; ++ind) {
    //for (auto & ind : channel.otfIndices) {
        complex_t pq = pPtr[ind] * qPtr[ind];
        //complex_t pq = conj(pPtr[ind]) * qPtr[ind];     // to match MvN
        double ps = norm (pPtr[ind]);
        double qs = qPtr[ind] * qPtr[ind];
        //glFTPtr[ind] = (pq*conj(ftPtr[ind]) - ps*otfPtr[ind]) / qs;
        glFTPtr[ind] = (pq*ftPtr[ind] - ps*otfPtr[ind]) / qs;
        //glFTPtr[ind] = (conj(pq*conj(ftPtr[ind])) - ps*otfPtr[ind]) / qs;     // to match MvN
    }
//  Ana::write("vogel_glft.f0", glFT);

    glFT.directInverse(gl);      // note: this will first reorder glFT
 // Ana::write("vogel_gl.f0", gl);
    glFT.zero();
   for (size_t ind=0; ind < N; ++ind) {
    //for (auto & ind : channel.otfIndices) {
        //glFTPtr[ind] = real(hjPtr[ind] * glPtr[ind]);     // eq. 22 in Vogel
        glFTPtr[ind] = hjPtr[ind] * glPtr[ind].real();     //  from appendix (eq 34.5) matches MvN
    }
    memcpy (glPtr, glFTPtr, glFT.nElements() *sizeof (complex_t));
 // Ana::write("vogel_gl2.f0", gl);

    glFT.reset (gl.get(), otfSize, otfSize, FT_FULLCOMPLEX );
    glFT.reorder();
  //Ana::write("vogel_glft2.f0", glFT);
    glFTPtr = glFT.get();
    double* vogPtr = vogel.get();
    //const double* pupPtr = channel.pupil.first.get();
    //double scale = -2/object.wavelength;
    for (auto & ind : channel.pupilInOTF) {
        //vogPtr[ind.first] = -2.0*imag(pfPtr[ind.first]*conj(glFTPtr[ind.second]))*pupPtr[ind.first];
        vogPtr[ind.first] = imag(conj(pfPtr[ind.first])*glFTPtr[ind.second]);//*pupPtr[ind.first];
        //vogPtr[ind.first] = imag(pfPtr[ind.first]*glFTPtr[ind.second]);
    }
  
 // Ana::write("vogel_p.f0", object.P);
 // Ana::write("vogel_q.f0", object.Q);
 // Ana::write("vogel_raw.f0", vogel);
//     std::vector<double> v;
//     for( auto& it : currentAlpha) v.push_back(it.second);
//     cout << printArray(v,"alpha") << endl;
//   cout << "BLA " << __LINE__ << endl;


}


void SubImage::resetPhi (Array<double>&p) const {
 //   cout << "SubImage::resetPhi(" << hexString(this) << ")"  << endl;
    memset(p.get(), 0, pupilSize*pupilSize*sizeof (double));    // FIXME: use phi_fixed
    //memcpy(p.get(), channel.phi_fixed.get(), channel.phi_fixed.nElements() *sizeof (double));
}


#define ALPHA_CUTOFF 1E-12
void SubImage::addScaledMode(double* phiPtr, uint16_t m, double a) const {
 //   cout << "a" << flush;
 //   if( fabs(a) > 0 ) cout << "SubImage::addScaledMode(" << hexString(this) << ")  m = " << m << "  a = " << a << "\n";
    addMode(phiPtr, m, a/object.wavelength);
}


void SubImage::addMode (double* phiPtr, uint16_t m, double a) const {
  //  cout << "SubImage::addMode(" << hexString(this) << ")  m = " << m << "   a = " << a << endl;
    //if(fabs(a) < ALPHA_CUTOFF) return;
    //memset(phiPtr,0,pupilSize*pupilSize*sizeof(double));
    const double* modePtr = channel.modes.at(m)->get();
    for (auto & ind : channel.pupilIndices) {
        phiPtr[ind] += a * modePtr[ind];
    }
//     Array<double> awrapper(phiPtr, pupilSize, pupilSize);
//     Ana::write("phi_"+hexString(this)+"_"+to_string(m)+".f0", awrapper);
//     awrapper.wrap(const_cast<double*>(modePtr), pupilSize, pupilSize);
//     Ana::write("mode_"+hexString(this)+"_"+to_string(m)+".f0", awrapper);
}


void SubImage::addModes (double* phiPtr, size_t nModes, uint16_t* modes, const double* alphas) const {
//    cout << "SubImage::addModes(" << hexString(this) << ")  " << printArray(alphas,nModes,"alphas") << endl;
    for (size_t m = 0; m < nModes; ++m) {
        addMode (phiPtr, modes[m], alphas[m]);
    }
}


void SubImage::adjustOffset(void) {
    bool adjusted(false);
    if( alpha.count(2) ) {
        double dadjust = alpha[2]*channel.alphaToPixels;
        int adjust = lround(alpha[2]*channel.alphaToPixels);
        int oadjust = adjust;
        if( adjust ) {
            adjust = shift(2,adjust);           // will return the "actual" shift. (the cube-edge might restrict it)
            if(adjust) {
                offsetShift.x += adjust;        // Noll-index = 2 corresponds to x-tilt.
                double oldval = alpha[2];
                double vshift = adjust*channel.pixelsToAlpha;
                alpha[2] -= vshift;
                LOG_TRACE << "adjustOffset:  mode 2 was adjusted.  from " << oldval << " to " << alpha[2] << "  (" << adjust<< " pixels)";
                adjusted = true;
            }
        }
    }
    if( alpha.count(3) ) {
        double dadjust = alpha[2]*channel.alphaToPixels;
        int adjust = lround(alpha[3]*channel.alphaToPixels);
        int oadjust = adjust;
        if( adjust ) {
            adjust = shift(1,adjust);           // will return the "actual" shift. (the cube-edge might restrict it)
            if(adjust) {
                offsetShift.y += adjust;        // Noll-index = 3 corresponds to y-tilt.
                double oldval = alpha[3];
                double vshift = adjust*channel.pixelsToAlpha;
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
   // std::cout << "addA: m="<<m<<"  a=" << a << std::endl;
    alpha[m] += a;
    
}


void SubImage::setAlpha(uint16_t m, double a) {
    
    unique_lock<mutex> lock (mtx);
   // std::cout << "setA: m="<<m<<"  a=" << a << std::endl;
    alpha[m] = a;
    
}


void SubImage::addAlphas(const double* alphas) {
    if (wf) {
        size_t cnt (0);
  //  cout << "SubImage::addAlphas(" << hexString(this) << ")  " << printArray(alphas,wf->modes.size(),"alphas") << endl;
        for (auto &it : wf->modes) {
            addAlpha(it.first, alphas[cnt++]);
        }
    }
}


void SubImage::setAlphas(const double* alphas) {
    if (wf) {
        size_t cnt (0);
   // cout << "SubImage::setAlphas(" << hexString(this) << ")  " << printArray(alphas,wf->modes.size(),"alphas") << endl;
        for (auto & it : wf->modes) {
            setAlpha(it.first, alphas[cnt++]);
        }
    }
}


namespace {
    mutex gmtx;
}

void SubImage::setAlphas(const std::vector<uint16_t>& modes, const double* a) {
  //  unique_lock<mutex> lock (gmtx);
    //LOG_TRACE << "SubImage::setAlphas(" << hexString(this) << ")  " << printArray(a,modes.size(),"alphas");
    int cnt=0;
    for (auto & m : modes) {
        //if( m > 2 ) {
            alpha[m] = a[cnt];
        //}
        cnt++;
    }

}


void SubImage::getAlphas(float* alphas) const {
    size_t cnt (0);
    for (auto& it: alpha) {
        alphas[cnt++] = it.second;
    }
    //LOG_TRACE << "SubImage::getAlphas(" << hexString(this) << ")  " << printArray(alphas,alpha.size(),"alphas");
}


void SubImage::resetPhi(void) {
    resetPhi(phi);
    for( auto& it: alpha ) {
        it.second = 0;
    }
    currentAlpha = alpha;
    calcPFOTF(PF.get(), OTF.get(), phi.get());
    newPhi = false;
    newOTF = true;
}

void SubImage::calcPhi(void) {
  //  unique_lock<mutex> lock (gmtx);
    memset(phi.get(), 0, pupilSize*pupilSize*sizeof (double));
    for( auto& it: alpha ) {
     //   if( fabs(it.second) > 0 ) cout << "SubImage::calcPhi():  adding mode " <<  it.first << " - " << it.second << endl;
        addScaledMode( it.first, it.second );  // only add difference
    }
    newPhi = true;
}


void SubImage::calcOTF(void) {
    if(newPhi) {
        calcPFOTF(PF.get(), OTF.get(), phi.get());
        newPhi = false;
        newOTF = true;
    }
}


void SubImage::calcOTF (complex_t* otfPtr, const double* phiPtr) const {

 //   cout << "SubImage::calcOTF(" << hexString(this) << ")" << endl;
// Array<double> awrapper(const_cast<double*>(phiPtr), pupilSize, pupilSize);
//  Ana::write("phisum.f0", awrapper);
    memset (otfPtr, 0, otfSize * otfSize * sizeof (complex_t));
    const double* pupilPtr = channel.pupil.first.ptr();
    double norm = sqrt(1.0 / channel.pupil.second);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)
    //double norm = 1.0 / channel.pupil.second;   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)
    for (auto & ind : channel.pupilInOTF) {
        double tmp = fmod(phiPtr[ind.first],pi_x_2);
        if( tmp < 0 ) tmp += pi_x_2;
        //uint32_t idx = static_cast<uint32_t> (tmp*angular_step_inv + 0.5);
        otfPtr[ind.second] = polar(pupilPtr[ind.first]*norm, phiPtr[ind.first]);
//         otfPtr[ind.second] = complex_t( pupilPtr[ind.first]*cos(phiPtr[ind.first]),
//                                         pupilPtr[ind.first]*sin(phiPtr[ind.first]) );
//         otfPtr[ind.second] = complex_t( pupilPtr[ind.first]*CosineLUT[idx],
//                                         pupilPtr[ind.first]*SineLUT[idx] );
    }
    FourierTransform::autocorrelate (otfPtr, otfSize, otfSize);
//     for (auto & ind : channel.pupilInOTF) {
//         otfPtr[ind.second] *= norm;
//     }

//   Ana::write("otf.f0", OTF);
//   Ana::write("otfP.f0", object.P);
//   Ana::write("otfQ.f0", object.Q);

}


void SubImage::calcPFOTF(complex_t* pfPtr, complex_t* otfPtr, const double* phiPtr) const {
    
 //   cout << "o" << flush;
//    cout << "SubImage::calcPFOTF(" << hexString(this) << ")" << endl;
//  Array<double> awrapper(const_cast<double*>(phiPtr), pupilSize, pupilSize);
//  Ana::write("phisum2.f0", awrapper);
    memset (pfPtr, 0, pupilSize * pupilSize * sizeof (complex_t));
    memset (otfPtr, 0, otfSize * otfSize * sizeof (complex_t));
//    redux::file::Ana::write("blu_pup.f0", channel.pupil.first);
    const double* pupilPtr = channel.pupil.first.ptr();
    double normalization = sqrt(1.0 / channel.pupil.second);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)
    //double norm = 1.0 / channel.pupil.second;   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)
    for (auto & ind : channel.pupilInOTF) {
 //       double tmp = fmod(phiPtr[ind.first],pi_x_2);
       // string tmpStr = "phi = " + to_string(phiPtr[ind.first]) + "  tmp = " + to_string(tmp);
 //       if( tmp < 0 ) {
//            tmp += pi_x_2;
          //  tmpStr += "  ntmp = " + to_string(tmp);
 //       }
 //       uint32_t idx = static_cast<uint32_t> (tmp*angular_step_inv + 0.5);
      //  tmpStr += "  idx = " + to_string(idx);
       pfPtr[ind.first] = polar(pupilPtr[ind.first], phiPtr[ind.first]);
        //pfPtr[ind.first].real(pupilPtr[ind.first]*cos(phiPtr[ind.first]));
        //pfPtr[ind.first].imag(pupilPtr[ind.first]*sin(phiPtr[ind.first]));
//         pfPtr[ind.first] = complex_t( pupilPtr[ind.first]*cos(phiPtr[ind.first]),
//                                       pupilPtr[ind.first]*sin(phiPtr[ind.first]) );
//        pfPtr[ind.first] = complex_t( pupilPtr[ind.first]*CosineLUT[idx],
//                                        pupilPtr[ind.first]*SineLUT[idx] );
        otfPtr[ind.second] = normalization*pfPtr[ind.first];
      //  if(ind.first%1000 == 0) cout << tmpStr << endl;
    }
//    redux::file::Ana::write("blu_pf.f0", PF);
//    redux::file::Ana::write("blu_phi.f0", phi);
    FourierTransform::autocorrelate(otfPtr, otfSize, otfSize);
 //   cout << "SubImage::calcPFOTF(" << hexString(this) << ")  norm = " << normalization << endl;
//     for (auto & ind : channel.pupilInOTF) {
//         otfPtr[ind.second] *= norm;
//     }
//Array<complex_t> bwrapper(otfPtr, otfSize, otfSize);
//Ana::write("potf.f0", OTF);
//   Ana::write("potfP.f0", object.P);
//   Ana::write("potfQ.f0", object.Q);
}


void SubImage::addPSF(Array<float>& out) {
    Array<complex_t> tmp(OTF.dimensions());
    OTF.directInverse (tmp);
    out += tmp;
}


Array<double> SubImage::getPSF (void) {
    Array<complex_t> tmp(OTF.dimensions());
    OTF.directInverse(tmp);
    return tmp.copy<double>();
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
