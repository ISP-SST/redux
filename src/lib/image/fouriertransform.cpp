#include "redux/image/fouriertransform.hpp"

#include <algorithm>
#include <map>
#include <mutex>
#include <set>

using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;


FourierTransform::Plan::Index::Index (const std::vector<size_t>& dims, TYPE t, uint8_t nt) : tp (t), nThreads (nt), sizes (dims) {
    sizes.erase( remove_if( sizes.begin(), sizes.end(), [](size_t i){ return i <= 1;}), sizes.end());
    if ( sizes.empty() ) {
        throw std::logic_error ("FT::Plan constructed with no non-trivial dimensions:  " + printArray (dims, "in"));
    }
}

bool FourierTransform::Plan::Index::operator< (const Index& rhs) const {

    if (tp == rhs.tp) {
        if (nThreads == rhs.nThreads) {
            return (sizes < rhs.sizes);
        } else return (nThreads < rhs.nThreads);
    } else return tp < rhs.tp;

}



FourierTransform::Plan::Plan ( const Index& i ) : id(i) {
    init();
}


FourierTransform::Plan::~Plan() {
    if (! id.sizes.empty()) {
        fftw_destroy_plan (forward_plan);
        fftw_destroy_plan (backward_plan);
    }
}


void FourierTransform::Plan::init (void) {
    
    fftw_plan_with_nthreads(id.nThreads);
    
    if (id.tp == R2C) {
        if (id.sizes.size() == 2) {
            auto in = sharedArray<double> (id.sizes[0], id.sizes[1]);
            auto out = sharedArray<fftw_complex> (id.sizes[0], id.sizes[1] / 2 + 1);
            forward_plan = fftw_plan_dft_r2c_2d (id.sizes[0], id.sizes[1], *in.get(), *out.get(), FFTW_MEASURE);
            backward_plan = fftw_plan_dft_c2r_2d (id.sizes[0], id.sizes[1], *out.get(), *in.get(), FFTW_MEASURE);
        } else
            if (id.sizes.size() == 1) {
                auto in = sharedArray<double> (id.sizes[0]);
                auto out = sharedArray<fftw_complex> (id.sizes[0] / 2 + 1);
                forward_plan = fftw_plan_dft_r2c_1d (id.sizes[0], in.get(), out.get(), FFTW_MEASURE);
                backward_plan = fftw_plan_dft_c2r_1d (id.sizes[0] / 2 + 1, out.get(), in.get(), FFTW_MEASURE);
            } else {
                throw std::logic_error ("FT::Plan::init() is only implemented for 1/2 dimensions, add more when/if needed: " + printArray (id.sizes, "dims"));
            }
    } else {
        if (id.tp == C2C) {
            if (id.sizes.size() == 2) {
                auto in = sharedArray<fftw_complex> (id.sizes[0], id.sizes[1]);
                auto out = sharedArray<fftw_complex> (id.sizes[0], id.sizes[1]);
                forward_plan = fftw_plan_dft_2d (id.sizes[0], id.sizes[1], *in.get(), *out.get(), FFTW_FORWARD, FFTW_MEASURE);
                backward_plan = fftw_plan_dft_2d (id.sizes[0], id.sizes[1], *in.get(), *out.get(), FFTW_BACKWARD, FFTW_MEASURE);
            } else {
                if (id.sizes.size() == 1) {
                    auto in = sharedArray<fftw_complex> (id.sizes[0]);
                    auto out = sharedArray<fftw_complex> (id.sizes[0]);
                    forward_plan = fftw_plan_dft_1d (id.sizes[0], in.get(), out.get(), FFTW_FORWARD, FFTW_MEASURE);
                    backward_plan =  fftw_plan_dft_1d (id.sizes[0], in.get(), out.get(), FFTW_BACKWARD, FFTW_MEASURE);
                } else {
                    throw std::logic_error ("FT::Plan::init() is only implemented for 1/2 dimensions, add more when/if needed: " + printArray (id.sizes, "dims"));
                }
            }
        } else {
            throw std::logic_error ("FT::Plan::init() unknown tp.");
        }
    }
}


namespace {
    
    struct PlansContainer {  // using a static instance below ensures that fftw_init_threads is called only once, and cleanup is done when program exits.
        PlansContainer() { fftw_init_threads(); };
        ~PlansContainer(){
            for( auto& it : plans) {
                if( it.second ) {
                    it.second.reset();
                }
            }
            fftw_cleanup_threads();
        };
        std::map<FourierTransform::Plan::Index, FourierTransform::Plan::Ptr> plans;
    };
    
}


FourierTransform::Plan::Ptr FourierTransform::getPlan (const std::vector<size_t>& dims, Plan::TYPE tp, uint8_t nThreads) {

    static PlansContainer pc;
    
    static mutex mtx;
    unique_lock<mutex> lock (mtx);
    auto plan = pc.plans.emplace( Plan::Index(dims, tp, nThreads), nullptr );

    if (plan.second) {  // insertion successful, i.e. new plan -> initialize it
        plan.first->second.reset( new Plan( plan.first->first ) );
    }

    return plan.first->second;
}


template <typename T>
void FourierTransform::reset (const T* rhs, size_t ySize, size_t xSize, int flags, uint8_t nT) {

    inputSize = ySize*xSize;
    normalized = false;
    centered = false;       // output from FFTW is in re-ordered form, i.e. not centered
    nThreads = nT;

    if (flags & FT_FULLCOMPLEX) {
        halfComplex = false;
        Array<complex_t> tmp(ySize,xSize);
        tmp.copyFrom<T>(rhs);
        if (flags & FT_REORDER) {
            reorder(tmp.get(),ySize,xSize);
        }
        init(tmp.get(),ySize,xSize);
    } else {
        halfComplex = true;     // transform of real data, let's save half the space.
        Array<double> tmp(ySize,xSize);
        tmp.copyFrom<T>(rhs);
        if (flags & FT_REORDER) {
            reorder(tmp.get(),ySize,xSize);
        }
        init(tmp.get(),ySize,xSize);
    }

    if (flags & FT_NORMALIZE) {
        normalize();
    }

}


namespace redux {
    namespace image {
        template <>
        void FourierTransform::reset (const complex_t* rhs, size_t ySize, size_t xSize, int flags, uint8_t nT) {

            inputSize = ySize*xSize;
            normalized = false;
            centered = false;       // output from FFTW is in re-ordered form, i.e. not centered
            halfComplex = false;
            nThreads = nT;

            Array<complex_t> tmp(ySize,xSize);
            memcpy(tmp.get(),rhs,ySize*xSize*sizeof(complex_t));

            if (flags & FT_REORDER) {
                reorder(tmp.get(),ySize,xSize);
            }

            init(tmp.get(),ySize,xSize);

            if (flags & FT_NORMALIZE) {
                normalize();
            }

        }
    }
}


FourierTransform::FourierTransform() : centered(false), halfComplex(false), normalized(true), nThreads(1), inputSize(0) {

}


FourierTransform::FourierTransform (size_t ySize, size_t xSize, int flags, uint8_t nT) :
    centered (flags & FT_REORDER), halfComplex (! (flags & FT_FULLCOMPLEX)), normalized(flags & FT_NORMALIZE), nThreads(nT), inputSize(ySize*xSize) {

    if (halfComplex) {
        Array<complex_t>::resize (ySize, xSize / 2 + 1);
    } else {
        Array<complex_t>::resize (ySize, xSize);
    }

    init();

}


FourierTransform::FourierTransform (const FourierTransform& rhs) : plan (rhs.plan), centered (rhs.centered),
    halfComplex (rhs.halfComplex), normalized (rhs.normalized), nThreads (rhs.nThreads), inputSize (rhs.inputSize) {
    rhs.copy (*this);
}


template <typename T>
FourierTransform::FourierTransform (const Array<T>& rhs, int flags, uint8_t nT) :
    centered (false), halfComplex (false), normalized (false), nThreads(nT), inputSize(0) {

    vector<size_t> dims = rhs.dimensions(true);
    if (dims.size() != 2 ) {
        throw logic_error ("FourierTransform only supports 2 dimensions at the moment: " + printArray (dimensions(), "dims"));
    }

    reset(rhs.ptr(), dims[0], dims[1], flags, nT);

}


template <typename T>
FourierTransform::FourierTransform (const T* rhs, size_t ySize, size_t xSize, int flags, uint8_t nT) :
    centered (false), halfComplex (false), normalized (false) {
        
    reset(rhs, ySize, xSize, flags, nT);

}
template FourierTransform::FourierTransform (const float*, size_t, size_t, int, uint8_t);
template FourierTransform::FourierTransform (const double*, size_t, size_t, int, uint8_t);
template FourierTransform::FourierTransform (const complex_t*, size_t, size_t, int, uint8_t);
template FourierTransform::FourierTransform (const int32_t*, size_t, size_t, int, uint8_t);


void FourierTransform::ft( double* data ) {

    fftw_execute_dft_r2c( plan->forward_plan, data, reinterpret_cast<fftw_complex*>(get()) );

}


void FourierTransform::ft( complex_t* data ) {

    fftw_execute_dft( plan->forward_plan, reinterpret_cast<fftw_complex*>(data), reinterpret_cast<fftw_complex*>(get()) );

}


void FourierTransform::ift( double* out ) {

    fftw_execute_dft_c2r( plan->backward_plan, reinterpret_cast<fftw_complex*>(get()), out);

}


void FourierTransform::ift( complex_t* out ) {

    fftw_execute_dft( plan->backward_plan, reinterpret_cast<fftw_complex*>(get()), reinterpret_cast<fftw_complex*>(out) );

}


void FourierTransform::init (void) {

    if( halfComplex ) { 
        plan = getPlan(dimensions(), Plan::R2C, nThreads);
    } else {
        plan = getPlan(dimensions(), Plan::C2C, nThreads);
    }

}


void FourierTransform::init (const double* in, size_t ySize, size_t xSize) {

    // for r2c transforms, the last dimension has size = n/2+1
    if( (ySize != dimSize(0)) || (xSize != 2*(dimSize(1)-1))) {
        Array<complex_t>::resize (ySize,xSize/2+1);
    }

    plan = getPlan({ySize, xSize}, Plan::R2C, nThreads);
    
    double* tmpData = new double[ySize*xSize];
    memcpy (tmpData, in, ySize*xSize*sizeof (double));  // copy data because r2c-transforms modifies input
    fftw_execute_dft_r2c (plan->forward_plan, tmpData, reinterpret_cast<fftw_complex*> (ptr()));
    delete[] tmpData;

}


void FourierTransform::init (const complex_t* in, size_t ySize, size_t xSize) {
    if( (ySize != dimSize(0)) || (xSize != dimSize(1))) {
        Array<complex_t>::resize(ySize,xSize);
    }

    plan = getPlan({ySize, xSize}, Plan::C2C, nThreads);
    complex_t* dataPtr = const_cast<complex_t*>(in);    // fftw takes non-const, even if input is not modified.
    fftw_execute_dft (plan->forward_plan, reinterpret_cast<fftw_complex*> (dataPtr), reinterpret_cast<fftw_complex*> (ptr()));
}


template <typename T>
void FourierTransform::directInverse(T* out) {

    if (centered) {
        reorder();
    }
    
    if (!normalized) {
        normalize();
    }

    if (halfComplex) {
        fftw_execute_dft_c2r (plan->backward_plan, reinterpret_cast<fftw_complex*> (ptr()), reinterpret_cast<double*> (out));
    } else {
        fftw_execute_dft (plan->backward_plan, reinterpret_cast<fftw_complex*> (ptr()), reinterpret_cast<fftw_complex*> (out));
    }

}
template void FourierTransform::directInverse (double* out);
template void FourierTransform::directInverse (complex_t* out);


template <typename T>
void FourierTransform::inv (Array<T>& out, int flags) const {
    vector<size_t> dims = dimensions();
    if (halfComplex) {
        dims.back() -= 1;
        dims.back() <<= 1;      // for r2c transforms, the last dimension has half size (n/2+1)
        Array<double> tmp(dims);
        size_t sz = nElements();
        complex_t* dataPtr = new complex_t[sz];
        memcpy (dataPtr, ptr(), sz * sizeof (complex_t));  // copy since r2c transforms modifies input
        fftw_execute_dft_c2r (plan->backward_plan, reinterpret_cast<fftw_complex*> (dataPtr), tmp.ptr());
        delete[] dataPtr;
        tmp.copy (out);
    } else {
        if (centered) {
            FourierTransform tmp(*this);
            tmp.reorder();
            tmp.inv(out,flags);
            return;
        }
        Array<complex_t> tmp(dims);
        complex_t* dataPtr = const_cast<complex_t*> (ptr());
        fftw_execute_dft (plan->backward_plan, reinterpret_cast<fftw_complex*> (dataPtr), reinterpret_cast<fftw_complex*> (tmp.ptr()));
        tmp.copy( out );
    }

    if (flags & FT_REORDER) {
        reorder(out);
    }

    if (flags & FT_NORMALIZE) {
        out *= (1.0/out.nElements());
    }

}
namespace redux {
    namespace image {
        template <>
        void FourierTransform::inv (Array<double>& out, int flags) const {
           const vector<size_t>& dims = dimensions();
            if (halfComplex) {
               if( (dims[0] != out.dimSize(0)) || (2*(dims[1]-1) != out.dimSize(1))) {
                    out.resize(dims[0],2*(dims[1]-1));
                }
                size_t sz = nElements();
                complex_t* dataPtr = new complex_t[sz];
                memcpy (dataPtr, ptr(), sz * sizeof (complex_t));  // copy since r2c transforms modifies input
                fftw_execute_dft_c2r (plan->backward_plan, reinterpret_cast<fftw_complex*> (dataPtr), out.ptr());
                delete[] dataPtr;
            } else {
                if (centered) {
                    FourierTransform tmp(*this);
                    tmp.reorder();
                    tmp.inv(out,flags);
                    return;
                }
                Array<complex_t> tmp(dims);
                complex_t* dataPtr = const_cast<complex_t*> (ptr());
                fftw_execute_dft (plan->backward_plan, reinterpret_cast<fftw_complex*> (dataPtr), reinterpret_cast<fftw_complex*> (tmp.ptr()));
                out = tmp;
            }

            if (flags & FT_REORDER) {
                reorder (out);
            }

            if (flags & FT_NORMALIZE) {
                out *= (1.0/out.nElements());
            }
        
        }
    }
}
template void FourierTransform::inv (Array<float>& out, int flags) const;
template void FourierTransform::inv (Array<double>& out, int flags) const;
template void FourierTransform::inv (Array<complex_t>& out, int flags) const;


template <typename T>
Array<T> FourierTransform::correlate (const Array<T>& in) const {
    
    int f = FT_REORDER;
    if (centered) {
        f = 0;
    }

    FourierTransform inFT (in, f);
    
    transform( get(), get()+nElements(), inFT.get(), inFT.get(),
        [](const complex_t& a, const complex_t& b) { return b*conj(a); } );

    Array<T> out (in.dimensions());
    inFT.inv (out, FT_REORDER);
    return out;

}
template Array<float> FourierTransform::correlate (const Array<float>&) const;
template Array<double> FourierTransform::correlate (const Array<double>&) const;
template Array<complex_t> FourierTransform::correlate (const Array<complex_t>&) const;


void FourierTransform::autocorrelate(void) {

    // FourierTransform is dense by construction, so this should always be ok
    transform( get(), get()+nElements(), get(), [](const complex_t&a){ return norm(a); } );

}


void FourierTransform::autocorrelate( double scale ) {

    // FourierTransform is dense by construction, so this should always be ok
    transform( get(), get()+nElements(), get(), [scale](const complex_t&a){ return scale*norm(a); } );

}


template <typename T>
void FourierTransform::autocorrelate(T* data, size_t ySize, size_t xSize ) {
    FourierTransform ft(data, ySize, xSize);
    ft.autocorrelate();
    ft.directInverse(data);
    reorder(data, ySize, xSize);
}
template void FourierTransform::autocorrelate(double*, size_t, size_t);
template void FourierTransform::autocorrelate(complex_t*, size_t, size_t);


template <typename T>
void FourierTransform::autocorrelate (const Array<T>& in, Array<T>& out) {
    FourierTransform ft (in);
    ft.autocorrelate();
    out.resize (in.dimensions());
    ft.directInverse(out.get());
    reorder(out);
}
template void FourierTransform::autocorrelate (const Array<double>&, Array<double>&);
template void FourierTransform::autocorrelate (const Array<complex_t>&, Array<complex_t>&);


Array<double> FourierTransform::power (void) const {
    Array<double> tmp (dimensions());
    transform( get(), get()+nElements(), tmp.get(),[](const complex_t&a){ return norm(a); } );
    return std::move(tmp);
}


double FourierTransform::noise (int mask, double limit) const {

    double noise_power = 0.0;

    int N = 0;
    double xx, yy;
    int npY = dimSize(0);
    int npX = dimSize(1);
    
    int64_t nElements = npY*npX;
    if(halfComplex) {
        nElements = npY*(npX/2+1);
    }
    
    const complex_t* dataPtr = get();

    if (mask < 0) {         // specifying mask < 0 will match MvN's noise_level computation   ( mask = (np/6) )
        mask = std::max (npY, npX) / 6 + 1;
    }

    if (limit < 0) {        // specifying limit < 0 will match MvN's noise_level computation  ( sqrt((npY*npX)/4.0) )
        if (halfComplex) {
            limit = sqrt ( (npY * (npX - 1)) / 2.0);
        } else  {
            limit = sqrt ( (npY * npX) / 4.0);
        }
    }

    double limits = limit * limit;
    int endX = npX;
    int endY = npY;

    if (mask > 0) {
        endY = npY - mask + 1;
        if (!halfComplex) {
            endX = npX - mask + 1;
        }
    }

    for (int x = mask; x < endX; ++x) {
        if (halfComplex) {
            xx = x;
        } else {
            xx = std::min (x, npX - x);
        }

        xx *= xx;
        int count = 1 + (halfComplex && x > 0 && x < npX - 1); // if halfComplex and x is not first/last column, value is counted twice.

        for (int y = mask; y < endY; ++y) {
            yy = std::min (y, npY - y);
            if ( (xx + yy * yy) < limits) continue;         // inside region to ignore
            noise_power += count * norm(dataPtr[y*npX+x]);
            N += count;
        }
    }
    if( N && nElements ) noise_power = std::sqrt (noise_power / (N * nElements));
    return noise_power;

}



template <typename T>
void FourierTransform::convolveInPlace (Array<T>& inout, int flags) const {

    auto dims = inout.dimensions();

    if (halfComplex) {            // for half-complex, the last dimension has half size (n/2+1)
        dims.back() >>= 1;
        dims.back() += 1;
    }

    if (dims != dimensions()) {
        cout << "FourierTransform::convolveInPlace(): input dimensions does not match this FourierTransform.\n\t"
             << printArray (dims, "input-dims") << printArray (dimensions(), "   FT-dims") << endl;
        return;
    }

    if (! normalized) {             // normalize only 1 of the inputs
        flags |= FT_NORMALIZE;
    }

    FourierTransform inFT (inout, flags, nThreads);
    inFT *= *this;
    inFT.inv (inout, FT_REORDER);
}
template void FourierTransform::convolveInPlace (Array<float>& inout, int flags) const;
template void FourierTransform::convolveInPlace (Array<double>& inout, int flags) const;
template void FourierTransform::convolveInPlace (Array<complex_t>& inout, int flags) const;


void FourierTransform::normalize (FourierTransform& in) {
    if (in.normalized) return;

    in *= (1.0 / (in.inputSize));
    in.inputSize = 1;
    in.normalized = true;
}


template <typename T>
void FourierTransform::reorderInto( const T* inSmall, size_t inSizeY, size_t inSizeX, T* outBig, size_t outSizeY, size_t outSizeX) {

    size_t inHalfY = inSizeY / 2;
    size_t inHalfX = inSizeX / 2;
    
    size_t chunkSize = inHalfX * sizeof (T);

    const T* southEastIn = inSmall + inHalfX;
    const T* northWestIn = inSmall + inHalfY*inSizeX;
    const T* northEastIn = inSmall + inHalfY*inSizeX + inHalfX;
    T* southEastOut = outBig + outSizeX - inHalfX;
    T* northWestOut = outBig + (outSizeY-inHalfY)*outSizeX;
    T* northEastOut = southEastOut + (outSizeY-inHalfY)*outSizeX;

    for (size_t y = 0; y < inHalfY; ++y) {
        memcpy (northEastOut, inSmall, chunkSize);
        memcpy (northWestOut, southEastIn, chunkSize);
        memcpy (outBig, northEastIn, chunkSize);
        memcpy (southEastOut, northWestIn, chunkSize);
        inSmall += inSizeX;
        northEastIn += inSizeX;
        southEastIn += inSizeX;
        northWestIn += inSizeX;
        outBig += outSizeX;
        northWestOut += outSizeX;
        southEastOut += outSizeX;
        northEastOut += outSizeX;
    }

}
template void FourierTransform::reorderInto( const double*, size_t, size_t, double*, size_t, size_t);
template void FourierTransform::reorderInto( const float*, size_t, size_t, float*, size_t, size_t);
template void FourierTransform::reorderInto( const complex_t*, size_t, size_t, complex_t*, size_t, size_t);


template <typename T>
void FourierTransform::reorder(T* in, size_t ySize, size_t xSize) {

    size_t halfY = ySize / 2;
    size_t halfX = xSize / 2;
    size_t chunkSize = halfX * sizeof (T);
    char *buf = new char[chunkSize];

    T* southEast = in + halfX;
    T* northWest = in + halfY*xSize;
    T* northEast = in + halfY*xSize + halfX;

    for (size_t y = 0; y < halfY; ++y) {
        memcpy (buf, in, chunkSize);
        memcpy (in, northEast, chunkSize);
        memcpy (northEast, buf, chunkSize);
        in += xSize;
        northEast += xSize;
        memcpy (buf, southEast, chunkSize);
        memcpy (southEast, northWest, chunkSize);
        memcpy (northWest, buf, chunkSize);
        southEast += xSize;
        northWest += xSize;
    }

    delete[] buf;

}


void FourierTransform::reorder (void) {
    if (!halfComplex) {
        reorder(get(),dimSize(0),dimSize(1));
        centered = !centered;
    } else {
        cout << "REORDER:  fix for half-complex." << endl;
    }
}


FourierTransform FourierTransform::reordered (void) const {
    FourierTransform tmp (*this);
    tmp.reorder();
    return tmp;
}


void FourierTransform::resize( size_t ySize, size_t xSize, int flags, uint8_t nT ) {
    
    if( (ySize != dimSize(0)) || (xSize != dimSize(1))) {
        Array<complex_t>::resize(ySize,xSize);
    }

    inputSize = ySize*xSize;
    normalized = false;
    centered = false;
    nThreads = nT;

    if (flags & FT_FULLCOMPLEX) {
        halfComplex = false;
        plan = getPlan( {ySize, xSize}, Plan::C2C, nThreads);
    } else {
        halfComplex = true;
        plan = getPlan( {ySize, xSize}, Plan::R2C, nThreads);
    }

    if (flags & FT_NORMALIZE) {
        normalized = true;
    }

    if (flags & FT_REORDER) {
        centered = true;
    }
    
}


const FourierTransform& FourierTransform::operator*= (const FourierTransform& rhs) {

    if (centered != rhs.centered) {
        *this *= rhs.reordered();
        return *this;
    }

    if (halfComplex == rhs.halfComplex) {
        //Array<complex_t>::operator*= (rhs);
        //std::transform (first, first+5, second, results, std::multiplies<int>());
        std::transform (ptr(), ptr()+nElements(), rhs.ptr(), ptr(), std::multiplies<complex_t>());
    } else if (halfComplex) {
        Array<complex_t> half (reinterpret_cast<const Array<complex_t>&> (rhs), dimensions());
        Array<complex_t>::operator*= (half);
    } else {
        Array<complex_t> half (reinterpret_cast<Array<complex_t>&> (*this), rhs.dimensions());

        int d0 = dimSize (0);
        int d1 = dimSize (1);
        int d2 = rhs.dimSize (1);
        int midpoint = d1 - d2 + 1; // works for both odd & even dimension sizes
        complex_t** thisPtr = redux::util::makePointers (ptr(), d0, d1);
        const complex_t** rhsPtr = redux::util::makePointers (rhs.ptr(), d0, d2);

        for (int y = 0; y < d0; ++y) {
            for (int x = 0; x < d2; ++x) {
                thisPtr[y][x] *= rhsPtr[y][x];
                if (x && x < midpoint) thisPtr[y][d1 - x] *= rhsPtr[y][x];
            }
        }

        redux::util::delPointers (thisPtr);
        redux::util::delPointers (rhsPtr);
    }

    return *this;
}


const FourierTransform& FourierTransform::operator= (const FourierTransform& rhs) {
    rhs.copy(*this);
    plan = rhs.plan;
    centered = rhs.centered;
    halfComplex = rhs.halfComplex;
    normalized = rhs.normalized;
    inputSize = rhs.inputSize;
    return *this;
}

// Explicit instantiations are needed when template implementations are not in the header.
// They should also be at the bottom to make sure all code is present when they are invoked.

template FourierTransform::FourierTransform (const Array<int16_t>&, int, uint8_t);
template FourierTransform::FourierTransform (const Array<int32_t>&, int, uint8_t);
template FourierTransform::FourierTransform (const Array<float>&, int, uint8_t);
template FourierTransform::FourierTransform (const Array<double>&, int, uint8_t);
template FourierTransform::FourierTransform (const Array<complex_t>&, int, uint8_t);

template void FourierTransform::reorder (redux::util::Array<float>& in);
template void FourierTransform::reorder (redux::util::Array<double>& in);
template void FourierTransform::reorder (redux::util::Array<complex_t>& in);


