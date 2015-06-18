#include "redux/momfbd/subimage.hpp"

#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/object.hpp"

#include "redux/file/fileana.hpp"
using namespace redux::file;

using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

    const std::string thisChannel = "subimage";

}


SubImage::SubImage (Object& obj, const Channel& ch, const redux::util::Array<double>& wind, const redux::util::Array<float>& stack,
                    uint32_t index, uint16_t firstY, uint16_t firstX, uint16_t patchSize, uint16_t pupilSize)
    : Array<float> (stack, index, index, firstY, firstY + patchSize - 1, firstX, firstX + patchSize - 1),
      index (index), imgSize (patchSize), pupilSize (pupilSize), otfSize (2 * pupilSize), object (obj), channel (ch), window (wind) {

    img.resize (imgSize, imgSize);
    phi.resize (pupilSize, pupilSize);
    PF.resize (pupilSize, pupilSize);
    OTF.resize (otfSize, otfSize);                  // big enough for autocorrelation of the PF
    vogel.resize (pupilSize, pupilSize);

}


SubImage::~SubImage (void) {

}


void SubImage::init (void) {

    stats.getStats (*this, ST_VALUES);                      // get mean of un-windowed image
    string str = "    SubImage::init(" + to_string (index) + "):  mean=" + to_string (stats.mean);

    auto wit = window.begin();
    Array<float>::const_iterator dit = this->begin();

    for (auto & iit : img) {                    // windowing: subtract and re-add mean afterwards
        iit = (*dit++ - stats.mean) * *wit++ + stats.mean;
    }

    stats.getStats (img, ST_VALUES | ST_RMS);               // get statistics
    str += "  mean2=" + to_string (stats.mean) + "  std=" + to_string (stats.stddev) + "  sum=" + to_string (stats.sum);

    if (imgSize == otfSize) {                                                       // imgSize = 2*pupilSize
        imgFT.reset (img.get(), imgSize, imgSize, FT_FULLCOMPLEX | FT_NORMALIZE);   // full-complex for now, perhaps half-complex later for performance
    } else {                                                                        // imgSize > 2*pupilSize should never happen (cf. calculatePupilSize)
        int offset = (otfSize - imgSize) / 2;
        Array<double> tmp (otfSize, otfSize);
        tmp.zero();
        double* imgPtr = img.get();
        double* tmpPtr = tmp.get();
        for (int i = 0; i < imgSize; ++i) {
            memcpy (tmpPtr + offset * (otfSize + 1) + i * otfSize, imgPtr + i * imgSize, imgSize * sizeof (double));
        }
        imgFT.reset (tmp.get(), otfSize, otfSize, FT_FULLCOMPLEX | FT_NORMALIZE);   // full-complex for now, perhaps half-complex later for performance
    }

    FourierTransform::reorder (imgFT);                             // keep FT in centered form
    stats.noise = channel.noiseFudge * imgFT.noise (-1, -1);       // mask/cutoff < 0 will revert to hardcoded values used by MvN

    str += "  noise=" + to_string (stats.noise);
    object.addToFT( imgFT, stats.noise/stats.stddev );



    cout << str << printArray (dimensions(), "  dims") << endl;

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
    for (auto & ind : channel.otfIndices) {
        Q[ind] += norm (otf[ind]);                    // Q += sj.re^2 + sj.im^2 = norm(sj)
        P[ind] += ftPtr[ind] * conj (otf[ind]);       // P += ft * conj(sj)
    }

}


double SubImage::metricChange(const complex_t* newOTF) const {

    const complex_t* oldOTF = OTF.get();
    const complex_t* p = object.P.get();
    const complex_t* ftPtr = imgFT.get();
    const double* q = object.Q.get();
    
    complex_t dp, dsj;
    double dl = 0.0, dq, dn;
    for (auto & ind : channel.otfIndices) {
        dsj = newOTF[ind] - oldOTF[ind];           // change in sj
        dp = conj (ftPtr[ind]) * dsj;         // change p and q
        dq = 2.0 * (oldOTF[ind].real() * dsj.real() + oldOTF[ind].imag() * dsj.imag()) + norm (dsj);
        dn = 2.0 * (dp.real() * p[ind].real() + dp.imag() * p[ind].imag()) + norm (dp);
        dl -= (q[ind] * dn - dq * (norm (p[ind]))) / (q[ind] * (q[ind] + dq));
    }
    return dl / channel.otfIndices.size();
}


double SubImage::gradientFiniteDifference(uint16_t mode, double dalpha, complex_t* tmpOTF, double* tmpPhi) const {
    memcpy(tmpPhi, phi.get(), pupilSize*pupilSize*sizeof(double));
    addMode(tmpPhi, mode, dalpha);
    calcOTF(tmpOTF, tmpPhi);
    return metricChange(tmpOTF) / dalpha;
}


double SubImage::gradientVogel(uint16_t mode, double dalpha, complex_t* tmpOTF, double* tmpPhi) const {
    double ret = 0;
    memset(tmpPhi, 0, pupilSize*pupilSize*sizeof(double));
    addMode(tmpPhi, mode, -2);
    const double* vogPtr = vogel.get();
    const double* pupPtr = channel.pupil.first.get();
    for (auto & ind : channel.pupilIndices) {
        ret += vogPtr[ind] * tmpPhi[ind] * pupPtr[ind];
    }
    return ret;
}


void SubImage::calcVogelWeight(void) {

    int offset = max ( (otfSize - pupilSize) / 2, 0);

    FourierTransform pj (otfSize, otfSize, FT_FULLCOMPLEX | FT_REORDER);    // input is already centered, so set the flag
    FourierTransform glFT (otfSize, otfSize, FT_FULLCOMPLEX | FT_REORDER);
    Array<complex_t> hj (otfSize, otfSize);
    Array<complex_t> gl (otfSize, otfSize);
    pj.zero();

    const complex_t* ptfPtr = PF.get();
    complex_t* pjPtr = pj.ptr (offset, offset);

    for (uint i = 0; i < pupilSize; ++i) {
        memcpy (pjPtr + i * otfSize, ptfPtr + i * pupilSize, pupilSize * sizeof (complex_t));
    }

    pj.directInverse (hj);      // note: this will first reorder pj

    const complex_t* pPtr = object.P.get();
    const double* qPtr = object.Q.get();
    const complex_t* otfPtr = OTF.get();
    const complex_t* ftPtr = imgFT.get();
    complex_t* glFTPtr = glFT.get();
    complex_t* glPtr = gl.get();
    complex_t* hjPtr = hj.get();

    glFT.zero();
    for (auto & ind : channel.otfIndices) {
        complex_t pq = pPtr[ind] * qPtr[ind];
        double ps = norm (pPtr[ind]);
        double qs = qPtr[ind] * qPtr[ind];
        glFTPtr[ind] = (pq * ftPtr[ind] - ps * otfPtr[ind]) / qs;
    }

    glFT.directInverse (gl);      // note: this will first reorder glFT
    glFT.zero();
    for (auto & ind : channel.otfIndices) {
        glFTPtr[ind] = hjPtr[ind] * glPtr[ind].real();
    }
    memcpy (glPtr, glFTPtr, glFT.nElements() *sizeof (complex_t));
    glFT.reset (gl.get(), otfSize, otfSize, FT_FULLCOMPLEX | FT_REORDER);
    glFTPtr = glFT.get();
    
    double* vogPtr = vogel.get();
    for (auto & ind : channel.pupilInOTF) {
        vogPtr[ind.first] = imag(conj(ptfPtr[ind.first])*glFTPtr[ind.second]);
    }

}


void SubImage::clearModes (redux::util::Array<double>&p) const {
    memcpy(p.get(), channel.phi_fixed.get(), channel.phi_fixed.nElements() *sizeof (double));
}


#define ALPHA_CUTOFF 1E-12
void SubImage::addMode (double* phiPtr, uint16_t m, double alpha) const {
    //alpha *= object.wavelength;
    //if(fabs(alpha) < ALPHA_CUTOFF) return;
    const double* modePtr = channel.modes.at(m)->get();
    for (auto & ind : channel.pupilIndices) {
        phiPtr[ind] += alpha * modePtr[ind];
    }
}


void SubImage::addModes (double* phiPtr, size_t nModes, uint16_t* modes, const double* alphas) const {
    for (size_t m = 0; m < nModes; ++m) {
        addMode (phiPtr, modes[m], alphas[m]);
    }
}


void SubImage::addPhases(double* phiPtr, const double* alpha) const {
    if (wf) {
        size_t cnt (0);
        for (auto & it : wf->modes) {
            addMode(phiPtr, it.first, alpha[cnt++]);
        }
    }
}


void SubImage::calcOTF (complex_t* otfPtr, const double* phiPtr) const {

    memset (otfPtr, 0, otfSize * otfSize * sizeof (complex_t));
    const double* pupilPtr = channel.pupil.first.ptr();
    double norm = sqrt(1.0 / channel.pupil.second);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)
    for (auto & ind : channel.pupilInOTF) {
        otfPtr[ind.second] = polar(pupilPtr[ind.first] * norm, phiPtr[ind.first]);
    }
    FourierTransform::autocorrelate (otfPtr, otfSize, otfSize);

}


void SubImage::calcPFOTF(complex_t* pfPtr, complex_t* otfPtr, const double* phiPtr) const {
    
    memset (otfPtr, 0, otfSize * otfSize * sizeof (complex_t));
    const double* pupilPtr = channel.pupil.first.ptr();
    double norm = sqrt(1.0 / channel.pupil.second);   // normalize OTF by pupil-area (sqrt since the autocorrelation squares the OTF)
    for (auto & ind : channel.pupilInOTF) {
        pfPtr[ind.first] = polar(pupilPtr[ind.first], phiPtr[ind.first]);
        otfPtr[ind.second] = pfPtr[ind.first]*norm;
    }
    FourierTransform::autocorrelate(otfPtr, otfSize, otfSize);

}


redux::util::Array<double> SubImage::getPSF (void) {
    redux::util::Array<double> tmp;
    OTF.directInverse (tmp);
    return std::move (tmp);
}


void SubImage::dump (std::string tag) {

    cout << "    dumping image:  this=" << hexString(this) << " with tag=" << tag << endl;
    cout << "                  phiPtr=" << hexString(phi.get()) << endl;
    Ana::write (tag + "_img.f0", img);
    Ana::write (tag + "_phi.f0", phi);
    Ana::write (tag + "_ptf.f0", PF);
    Ana::write (tag + "_otf.f0", OTF);
    Ana::write (tag + "_imgFT.f0", imgFT);
    Ana::write (tag + "_window.f0", window);
    Ana::write (tag + "_pupil.f0", channel.pupil.first);
    Ana::write (tag + "_P.f0", object.P);
    Ana::write (tag + "_Q.f0", object.Q);
    Ana::write (tag + "_vogel.f0", vogel);

    /*if( wf && !wf->alphas.empty() ) {
        vector<float> alphas;
        for( auto& it: wf->alphas ) {
            if( it.first >= alphas.size() ) alphas.resize(it.first,0.0);
            alphas[it.first] = it.second.value;
        }
        Ana::write ("alpha_" + tag + ".f0", alphas);
    }*/

}
