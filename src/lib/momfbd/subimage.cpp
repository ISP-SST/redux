#include "redux/momfbd/subimage.hpp"

#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/object.hpp"

#include "redux/file/fileana.hpp"

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
    SJ.resize (otfSize, otfSize);                  // big enough for autocorrelation of the PTF
    vogel.resize (pupilSize, pupilSize);

}


SubImage::~SubImage (void) {

    //  std::cout << "~SubImage():  1   " << hexString(this) << std::endl;
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
        ft.reset (img.get(), imgSize, imgSize, FT_FULLCOMPLEX | FT_NORMALIZE);   // full-complex for now, perhaps half-complex later for performance
    } else {                                                                        // imgSize > 2*pupilSize should never happen (cf. calculatePupilSize)
        int offset = (otfSize - imgSize) / 2;
        Array<double> tmp (otfSize, otfSize);
        tmp.zero();
        double* imgPtr = img.get();
        double* tmpPtr = tmp.get();
        for (int i = 0; i < imgSize; ++i) {
            memcpy (tmpPtr + offset * (otfSize + 1) + i * otfSize, imgPtr + i * imgSize, imgSize * sizeof (double));
        }
        ft.reset (tmp.get(), otfSize, otfSize, FT_FULLCOMPLEX | FT_NORMALIZE);   // full-complex for now, perhaps half-complex later for performance
    }

    FourierTransform::reorder (ft);                             // keep FT in centered form
    stats.noise = channel.noiseFudge * ft.noise (-1, -1);       // mask/cutoff < 0 will revert to hardcoded values used by MvN

    str += "  noise=" + to_string (stats.noise);
    object.addToFT( ft );


//      static int bla(0);
//      dump("init_"+to_string(++bla));


    cout << str << printArray (dimensions(), "  dims") << endl;

}


void SubImage::addFT(Array<double>& ftsum) const {
//    cout << "SubImage::addFT:  " << __LINE__ << "   this=" << hexString(this) << endl;
    const complex_t* ftPtr = ft.get();
    double* ftsPtr = ftsum.get();
    for (size_t ind = 0; ind < ft.nElements(); ++ind) {
        ftsPtr[ind] += norm (ftPtr[ind]);
    }
}

void SubImage::addPQ (redux::util::Array<complex_t>& P, redux::util::Array<double>& Q) const {
//    cout << "SubImage::addPQ(" << hexString(this) << ")  " << printArray(P.dimensions(), "pdims") << printArray(Q.dimensions(), "  qdims") << endl;
//    cout << "SubImage::addPQ(" << hexString(this) << ")  " << printArray(ft.dimensions(), "ftdims") << printArray(SJ.dimensions(), "  sjdims") << endl;

//     static int bla(0);
//     redux::file::Ana::write("sj_" + to_string(bla++) + ".f0", SJ);
    complex_t mx = SJ (0, 0);
    auto ftit = ft.begin();
    auto pit = P.begin();
    auto qit = Q.begin();

    for (auto & sjit : SJ) {
        mx.real (max (mx.real(), abs (sjit.real())));
        mx.imag (max (mx.imag(), abs (sjit.imag())));
        *qit++ += norm (sjit);                  // Q += sj.re^2 + sj.im^2 = norm(sj)
        *pit++ += *ftit++ * conj (sjit);        // P += ft * conj(sj)
    }

//    cout << "SubImage::addPQ(" << hexString(this) << ")  max = " << mx << endl;
}


void SubImage::computePhases (void) {
    cout << "SubImage::computePhases(" << hexString (this) << ")   wf = " << hexString (wf.get()) << endl;

    if (wf) {
        channel.getPhi (phi, *wf);
    }

//    cout << "SubImage::computePhases(" << hexString(this) << ")2" << printArray(phi.dimensions(),"dims") << endl;
}


namespace {

    double old_dmetric (int nPixels, const double **q, const complex_t **p, const complex_t **sj, const complex_t **osj, const complex_t **imgFT) {
        complex_t dp, dsj, max_diff;
        double dl = 0.0, dq, dn;
        //nPixels = 2;
        static int cnt (0);
        Array<complex_t> diff (nPixels, nPixels);

        //cout << "dmetric(): nPixels = " << nPixels << "   dl = " << dl << endl;
        for (int y = 0; y < nPixels; ++y) {
            for (int x = 0; x < nPixels; ++x) {                       // combination of make_qp and metric
                diff (y, x) = dsj = sj[y][x] - osj[y][x];           // change in sj

                if (norm (dsj) > norm (max_diff)) max_diff = dsj;

                dp = imgFT[y][x] * conj (dsj);          // change p and q
                dq = 2.0 * (osj[y][x].real() * dsj.real() + osj[y][x].imag() * dsj.imag()) + norm (dsj);
                dn = 2.0 * (dp.real() * p[y][x].real() + dp.imag() * p[y][x].imag()) + norm (dp);
                dl -= (q[y][x] * dn - dq * (norm (p[y][x]))) / (q[y][x] * (q[y][x] + dq));
//         if(x == nPixels/3 && y == nPixels/4)
//             cout << "dmetric(): x = " << x << "  y = " << y << "  dp = " << dp << "  dq = " << dq << "  sj = " << sj[y][x] << "  osj = " << osj[y][x] << "  dsj = " << dsj << endl;
            }
        }

        cout << "  dmetric(): max_diff = " << max_diff << endl;
        redux::file::Ana::write ("diff_" + to_string (cnt++) + ".f0", diff);
        dl /= (double) (nPixels * nPixels);
        //if(reg_alpha) dl+=0.5*reg_alpha*(2.0*alpha[m]*dalpha+sqr(dalpha))/atm_rms[m]; // only one term contributes
        return dl;
    }

}



// note: q and p are not modified!
void SubImage::oldGradientDiff (vector<double>& grad) {

    int pupilPixels = channel.pupilPixels;
    int imgPixels = channel.patchSize;
    int offset = max ( (2 * pupilPixels - imgPixels) / 2, 0);
    cout << "oldGradientDiff()  pupilPixels=" << pupilPixels << "  imgPixels=" << imgPixels << "  offset=" << offset << endl;
    redux::file::Ana::write ("p.f0", object.P);
    redux::file::Ana::write ("q.f0", object.Q);
    double dalpha = (double) 1.0E-2 * object.wavelength;
    int cnt (0);
    static int obj_cnt (0);
    obj_cnt++;
    SJ.zero();
    Array<complex_t> OSJ;
    SJ.copy (OSJ);
    // SJ += 1.2E-3;
    const complex_t** p = makePointers<const complex_t> (object.P.ptr(), object.P.dimSize (0), object.P.dimSize (1));
    const double** q = makePointers<const double> (object.Q.ptr(), object.Q.dimSize (0), object.Q.dimSize (1));
    const complex_t** sj = makePointers ( (const complex_t*) SJ.ptr(), SJ.dimSize (0), SJ.dimSize (1));
    const complex_t** osj = makePointers ( (const complex_t*) OSJ.ptr(), OSJ.dimSize (0), OSJ.dimSize (1));
    const complex_t** ftp = makePointers ( (const complex_t*) ft.ptr(), ft.dimSize (0), ft.dimSize (1));

    //complex_t **osj = ct2dim(1, nPixels, 1, nPixels);                    // original sj: only np x np since not used in autocorrelation
    //for(int x=0; x < nPixels; ++x) memcpy(osj[x], sj[x], nPixels*sizeof(complex_t)); // save original sj: beware: osj=npxnp but sj=2*nphx2*nph
    for (auto & it : wf->alpha) {
        cout << "oldGradientDiff(2)  cnt = " << cnt << "   first=" << it.first << "   sj=" << sj[imgPixels / 4][imgPixels / 3] << flush;
        redux::file::Ana::write ("phia_" + to_string (obj_cnt) + "_" + to_string (it.first) + ".f0", phi);
        redux::file::Ana::write ("osj_" + to_string (obj_cnt) + "_" + to_string (it.first) + ".f0", OSJ);
        channel.addMode (phi, it.first, -dalpha);              // new phi
        redux::file::Ana::write ("phib_" + to_string (obj_cnt) + "_" + to_string (it.first) + ".f0", phi);
        calcOTF();                                    // new sj
        //SJ += (cnt+1);
        cout << "   sj2=" << sj[imgPixels / 4][imgPixels / 3] << flush;
        redux::file::Ana::write ("sj_" + to_string (obj_cnt) + "_" + to_string (it.first) + ".f0", SJ);
        redux::file::Ana::write ("sjdiff_" + to_string (obj_cnt) + "_" + to_string (it.first) + ".f0", SJ - OSJ);
        double g = old_dmetric (imgPixels, q, p, sj, osj, ftp); // / dalpha; // calculate partial derivative
        cout << "  g = " << g << endl;
        grad[cnt++] += g;
        channel.addMode (phi, it.first, dalpha);              // old phi
    }

    delPointers (p);
    delPointers (q);
    delPointers (sj);
    delPointers (osj);
    delPointers (ftp);
    OSJ.copy (SJ);
    //for(int x=0; x < nPixels; ++x) memcpy(sj[x], osj[x], nPixels*sizeof(complex_t)); // restore original sj
    //del_ct2dim(osj, 1, nPixels, 1, nPixels);

}
/*
mn=3
pupil=f0('pupil.f0')
;mode2=f0('mode_2.f0')*pupil
;mode3=f0('mode_3.f0')*pupil
mode4=f0('mode_4.f0')*pupil
;mode5=f0('mode_5.f0')*pupil
phia=f0('phia_4.f0')
phib=f0('phib_4.f0')
sj=f0('sj_4.f0')
osj=f0('osj_4.f0')
diff0=f0('diff_0.f0')
tvscl,phia
tvscl,phib
tvscl,phib
print,min(phia),max(phia),mean(phia)
print,min(phib),max(phib),mean(phib)
;print,min(mode3),max(mode3),mean(mode3)
print,min(mode4),max(mode4),mean(mode4)
print,min(sj),max(sj),mean(sj)
*/


void SubImage::clearModes (redux::util::Array<double>&p) const {
    //cout << "SubImage::clearModes()  this=" << hexString (this) << "  phiPtr=" << hexString(p.get()) << endl;
    memcpy(p.get(), channel.phi_fixed.get(), channel.phi_fixed.nElements() *sizeof (double));
}

#define ALPHA_CUTOFF 1E-12
void SubImage::calcOTF (void) {
    //    cout << "SubImage::OTF()  channel.pupilPixels = " << channel.pupilPixels << endl;
    SJ.zero();
    complex_t** sj = makePointers (SJ.ptr(), SJ.dimSize (0), SJ.dimSize (1));
    const double** pupil = makePointers (channel.pupil.first.ptr(), channel.pupil.first.dimSize (0), channel.pupil.first.dimSize (1));
    double pval;

    for (int y = 0; y < channel.pupilPixels; ++y) {
        for (int x = 0; x < channel.pupilPixels; ++x) {
            pval = pupil[y][x];
            //if( (pval=pupil[y][x]) > 0 ) {
            //sj[y][x] = pval * phi(y, x);
            sj[y][x].real (pval * cos (phi (y, x)));
            sj[y][x].imag (pval * sin (phi (y, x)));
            //}
        }
    }

    delPointers (sj);
    delPointers (pupil);
//     if(channel.diversityModes.size()) {
//         redux::file::Ana::write("otfphi.f0", phi);
//     }
//     if(channel.diversityModes.size()) {
//         redux::file::Ana::write("newrawotf.f0", SJ);
//     }

    FourierTransform::autocorrelate (SJ);
    SJ /= channel.pupil.second;     // normalize by pupil-area

//     if(channel.diversityModes.size()) {
//         redux::file::Ana::write("newotf.f0", SJ);
//     }
}


void SubImage::calcOTF (const std::map<uint16_t, std::pair<double, bool>>& alpha) {

}


redux::util::Array<double> SubImage::PSF (void) {
    redux::util::Array<double> tmp;
    SJ.inv (tmp);

    if (channel.diversityModes.size()) {
        redux::file::Ana::write ("newpsf.f0", tmp);
    }

    return tmp;
}


void SubImage::clear (void) {
    wf.reset();
}
