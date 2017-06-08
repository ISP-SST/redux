#include "redux/image/descatter.hpp"


#include "redux/image/fouriertransform.hpp"
#include "redux/image/utils.hpp"

#include "redux/util/arrayutil.hpp"

using namespace redux;
using namespace redux::image;
using namespace redux::util;

using namespace std;


void redux::image::descatter( double* img, size_t imgY, size_t imgX, double* gain, fftw_complex* otf, size_t paddedY, size_t paddedX, int nthreads, int niter, float epsilon ) {
    
    size_t posX = (paddedX-imgX)/2;
    size_t posY = (paddedY-imgY)/2;
    
    size_t ftSize = paddedY*(paddedX/2+1);
    size_t nPixels = imgY*imgX;
    size_t nPaddedPixels = paddedY*paddedX;
    
    FourierTransform::Plan::Ptr plan = FourierTransform::Plan::get( paddedY, paddedX, FourierTransform::Plan::R2C, nthreads );
    
    std::shared_ptr<double> ft( (double*)fftw_malloc(ftSize*sizeof(fftw_complex)), fftw_free );
    fftw_complex* ftPtr = reinterpret_cast<fftw_complex*>(ft.get());
    std::shared_ptr<double> tmp( (double*)fftw_malloc(nPaddedPixels*sizeof(double)), fftw_free );
    double* tmpPtr = tmp.get();
    std::shared_ptr<double> tmpC( (double*)fftw_malloc(nPaddedPixels*sizeof(double)), fftw_free );
    double* tmpPtrC = tmpC.get();
    
    if( (paddedX==imgX) && (paddedY==imgY) ) {
        memcpy( tmpPtr, img, nPixels*sizeof(double) );
    } else {
        copyInto( img, imgY, imgX, tmpPtr, paddedY, paddedX, posY, posX );
    }
        
    double chiSq;
    double limit = nPixels*epsilon;                // multiply with image-size and skip dividing chiSq with it in every iteration

    
    int i(0);
    do {
        std::transform( tmpPtr, tmpPtr+nPaddedPixels, gain, tmpPtrC, []( const double& im, const double& g ){ return im*g; });

        plan->forward( tmpPtrC, ftPtr );
        for( size_t n=0; n<ftSize; ++n ) reinterpret_cast<complex_t&>(ftPtr[n]) *= reinterpret_cast<complex_t&>(otf[n]);
        plan->backward( ftPtr, tmpPtrC );
        
        std::transform( tmpPtrC, tmpPtrC+nPaddedPixels, gain, tmpPtrC, []( const double& im, const double& g ){ return im*g; });

        // Compute ChiSq
        chiSq = 0.0;
        double* tmpRowC = tmpPtrC + posY*paddedX + posX;
        double* tmpRow = tmpPtr + posY*paddedX + posX;
        double* imgRow = img;
        for( size_t y = 0; y < imgY; ++y ) {
            for( size_t x = 0; x < imgX; ++x ) {
                double newValue = imgRow[x] - tmpRowC[x];
                chiSq += norm(tmpRow[x] - newValue);   // norm() works for complex numbers too.
                tmpRow[x] = newValue;
            }
            tmpRowC += paddedX;
            tmpRow  += paddedX;
            imgRow  += imgX;
        }
        
        //cout << "cdescatter:  iteration #" << i << "  ChiSq = " << chiSq << "  ChiSq/ nPixels = " << (chiSq / nPixels) << endl;

    } while ( (chiSq > limit) && (++i < niter) );
        
    copyPart( tmpPtr, paddedY, paddedX, img, imgY, imgX, posY, posX );
 
}


void redux::image::descatter( double* img, size_t nRows, size_t nCols, double* gain, double* psf, size_t paddedY, size_t paddedX, int nthreads, int niter, float epsilon ) {

    size_t ftSize = nRows*(nCols/2+1);
    FourierTransform::Plan::Ptr plan = FourierTransform::Plan::get( nRows, nCols, FourierTransform::Plan::R2C, 1 );
    std::shared_ptr<double> otf( (double*)fftw_malloc(ftSize*sizeof(fftw_complex)), fftw_free );
    fftw_complex* otfPtr = reinterpret_cast<fftw_complex*>(otf.get());
    plan->forward( psf, otfPtr );
    descatter( img, nRows, nCols, gain, otfPtr, paddedY, paddedX, nthreads, niter, epsilon );

}


template <typename T, typename U>
void redux::image::descatter (Array<T>& data, const Array<U>& ccdgain, const Array<U>& psf_in, int maxIterations, double minImprovement, double epsilon) {

    vector<size_t> dims = data.dimensions (true);

    if (dims.size() != 2 || dims != ccdgain.dimensions() || dims != psf_in.dimensions()) {
        cout << "descatter(): dimensions of gain/psf does not match image." << endl;
        return;
    }

    for (auto & dim : dims) dim *= 2;

    Array<double> img (dims);                   // Twice the size of input
    Array<double> img_center (img, dims[0] / 4, 3 * dims[0] / 4 - 1, dims[1] / 4, 3 * dims[1] / 4 - 1); // centered subimage of half (i.e. original) size. N.B: shares data with img.
    img.zero();                                 // whole array = 0
    img_center.assign(data);                    // central subimage = img_in

    Array<double> psf (dims);                  // Twice the size of input
    Array<double> psf_center (psf, dims[0] / 4, 3 * dims[0] / 4 - 1, dims[1] / 4, 3 * dims[1] / 4 - 1); // centered subimage of half (i.e. original) size. N.B: shares data with psf.
    psf.zero();                                 // whole array = 0
    psf_center.assign(psf_in);                  // central subimage = psf_in
    double sum = total (psf_in);
    if (sum) {                                  // normalize
        psf_center *= 1.0 / sum;
    }

    Array<double> gain;
    ccdgain.copy(gain);

    redux::image::FourierTransform otf (psf, FT_REORDER | FT_NORMALIZE);

    Array<double>::const_iterator itg = gain.begin();
    for (auto & value : img_center) {
        double g = *itg++;
        value /= (1.0 + g * g);
    }

    Array<double> tmp (dims);                  // Twice the size of input
    Array<double> tmp_center (tmp, dims[0] / 4, 3 * dims[0] / 4 - 1, dims[1] / 4, 3 * dims[1] / 4 - 1); // centered subimage of half (i.e. original) size. N.B: shares data with tmp.

    double metric = std::numeric_limits< double >::max();
    double delta;
    int i = 0;
    do {
        img.copy (tmp);
        tmp_center *= gain;
        otf.convolveInPlace (tmp);
        tmp_center *= gain;
        delta = metric;
        metric = 0.0;
        typename Array<T>::const_iterator in_it = data.begin();
        Array<double>::const_iterator tmp_it = tmp_center.begin();
        for (auto & value : img_center) {
            double new_value = (*in_it++ - *tmp_it++);
            metric += (value - new_value) * (value - new_value);
            value = new_value;
        }
        metric /= data.nElements();
        delta /= metric;
    } while ( (delta > minImprovement) && (metric > epsilon) && (++i < maxIterations));

    img_center.copy (data);

}
template void redux::image::descatter (Array<float>&, const Array<float>&, const Array<float>&, int, double, double);
template void redux::image::descatter (Array<double>&, const Array<float>&, const Array<float>&, int, double, double);
template void redux::image::descatter (Array<float>&, const Array<double>&, const Array<double>&, int, double, double);
