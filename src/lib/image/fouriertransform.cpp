#include "redux/image/fouriertransform.hpp"

#include <redux/util/cache.hpp>

#include <algorithm>
#include <map>
#include <mutex>
#include <set>

using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

FourierTransform::PlansContainer FourierTransform::Plan::pc;

FourierTransform::Plan::Index::Index (const std::vector<size_t>& dims, TYPE t, uint8_t nt) : tp (t), nThreads (nt), sizes (dims) {
    sizes.erase( remove_if( sizes.begin(), sizes.end(), [](size_t i){ return i <= 1;}), sizes.end());
    if ( sizes.empty() ) {
        throw std::logic_error ("FT::Plan constructed with no non-trivial dimensions:  " + printArray (dims, "in"));
    }
}


FourierTransform::Plan::Index::Index (size_t sizeY, size_t sizeX, TYPE t, uint8_t nt) : tp (t), nThreads (nt), sizes({sizeY,sizeX}) {
    sizes.erase( remove_if( sizes.begin(), sizes.end(), [](size_t i){ return i <= 1;}), sizes.end() );
    if ( sizes.empty() ) {
        throw std::logic_error ("FT::Plan constructed with no non-trivial dimensions:  ["
        + to_string(sizeY) + ", " + to_string(sizeX) + "]");
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
    
    size_t nPix(1);
    for( auto& n: id.sizes ) nPix *= n;
    
    std::shared_ptr<complex_t> tmp1 = rdx_get_shared<complex_t>(nPix);
    std::shared_ptr<complex_t> tmp2 = rdx_get_shared<complex_t>(nPix);
    fftw_complex* ptrC1 = reinterpret_cast<fftw_complex*>(tmp1.get());
    fftw_complex* ptrC2 = reinterpret_cast<fftw_complex*>(tmp2.get());
    double* ptrD = reinterpret_cast<double*>(tmp1.get());
    
    if (id.tp == R2C) {
        if (id.sizes.size() == 2) {
            forward_plan = fftw_plan_dft_r2c_2d( id.sizes[0], id.sizes[1], ptrD, ptrC2, FFTW_MEASURE );
            backward_plan = fftw_plan_dft_c2r_2d( id.sizes[0], id.sizes[1], ptrC2, ptrD, FFTW_MEASURE );
        } else
            if (id.sizes.size() == 1) {
                forward_plan = fftw_plan_dft_r2c_1d( nPix, ptrD, ptrC2, FFTW_MEASURE);
                backward_plan = fftw_plan_dft_c2r_1d( nPix, ptrC2, ptrD, FFTW_MEASURE);
            } else {
                throw std::logic_error ("FT::Plan::init() is only implemented for 1/2 dimensions, add more when/if needed: " + printArray (id.sizes, "dims"));
            }
    } else {
        if (id.tp == C2C) {
            if (id.sizes.size() == 2) {
                forward_plan = fftw_plan_dft_2d( id.sizes[0], id.sizes[1], ptrC1, ptrC2, FFTW_FORWARD, FFTW_MEASURE );
                backward_plan = fftw_plan_dft_2d( id.sizes[0], id.sizes[1], ptrC2, ptrC1, FFTW_BACKWARD, FFTW_MEASURE );
            } else {
                if (id.sizes.size() == 1) {
                    forward_plan = fftw_plan_dft_1d( nPix, ptrC1, ptrC2, FFTW_FORWARD, FFTW_MEASURE);
                    backward_plan =  fftw_plan_dft_1d( nPix, ptrC2, ptrC1, FFTW_BACKWARD, FFTW_MEASURE);
                } else {
                    throw std::logic_error ("FT::Plan::init() is only implemented for 1/2 dimensions, add more when/if needed: " + printArray (id.sizes, "dims"));
                }
            }
        } else {
            throw std::logic_error ("FT::Plan::init() unknown tp.");
        }
    }
}


void FourierTransform::Plan::forward( double* __restrict__ in, fftw_complex* __restrict__ out ) const {
    
    if( id.tp != R2C ) throw std::logic_error ("FT::Plan::forward(double*, fftw_complex*) is only avainlable for Real-To-Complex Plans.");
    fftw_execute_dft_r2c( forward_plan, in, out );
    
}


void FourierTransform::Plan::backward( fftw_complex* __restrict__ in, double* __restrict__ out ) const {
    
    if( id.tp != R2C ) throw std::logic_error ("FT::Plan::backward(fftw_complex*, double*) is only avainlable for Real-To-Complex Plans.");
    fftw_execute_dft_c2r( backward_plan, in, out );
    
}


FourierTransform::Plan::Ptr FourierTransform::Plan::get(const std::vector<size_t>& dims, Plan::TYPE tp, uint8_t nThreads) {

    Plan::Index id(dims, tp, nThreads);
    Plan::Ptr& plan = Cache::get< Plan::Index, Plan::Ptr >( id, nullptr );

    unique_lock<mutex> lock(pc.mtx);

    if( !plan ) {  // insertion successful, i.e. new plan -> initialize it
        plan.reset( new Plan(id ) );
    }

    return plan;


}


FourierTransform::Plan::Ptr FourierTransform::Plan::get( size_t sizeY, size_t sizeX, Plan::TYPE tp, uint8_t nThreads ) {

    Plan::Index id(sizeY, sizeX, tp, nThreads);
    Plan::Ptr& plan = Cache::get< Plan::Index, Plan::Ptr >( id, nullptr );

    unique_lock<mutex> lock(pc.mtx);

    if( !plan ) {  // insertion successful, i.e. new plan -> initialize it
        plan.reset( new Plan(id ) );
    }

    return plan;

}

void FourierTransform::Plan::clear( void ) {

    Cache::clear<Plan::Index, Plan::Ptr >();

}



void FourierTransform::init( size_t ySize, size_t xSize, int flags, uint8_t nT ) {
    
    inputSize.y = ySize;
    inputSize.x = xSize;
    inPixels = ySize*xSize;
    ftSize = inputSize;
    
    currentFlags = flags;
    normalized = false;
    centered = false;
    nThreads = nT;
    
    size_t newBlockSize = ftSize.y * ftSize.x;  // always allocate for full-complex
    if( !(flags&FULLCOMPLEX) ) {
        ftSize.x = ftSize.x/2 + 1;
    }
    ftPixels = ftSize.y*ftSize.x;
    
    resize( newBlockSize );
    if( ftPixels > 0 ) {
        wrap();
        plan_full = Plan::get( {ySize, xSize}, Plan::C2C, nThreads);
        plan_half = Plan::get( {ySize, xSize}, Plan::R2C, nThreads);
    }
}


namespace redux {
    namespace image {
        
        template <typename T>
        int FourierTransform::ft( const T* __restrict__ in, complex_t* __restrict__ out, int flags ) const {
            size_t N = inputSize.y * inputSize.x;
            if( !(flags&FULLCOMPLEX) ) {
                flags &= ~REORDER_FT;     // halfComplex can't be reordered....yet.
                N = inputSize.y * (inputSize.x/2+1);
                double* __restrict__ tmpPtrD = reinterpret_cast<double*>(tmpPtr);
                if( flags&REORDER_IMG ) {
                    reorderInto( in, inputSize.y, inputSize.x, tmpPtrD, inputSize.y, inputSize.x );
                } else {
                    std::copy_n( in, inPixels, tmpPtrD );
                }
                fftw_execute_dft_r2c( plan_half->forward_plan, tmpPtrD, reinterpret_cast<fftw_complex*>(out) );
            } else {
                complex_t* __restrict__ tmpPtrC = tmpPtr;
                if( flags&REORDER_IMG ) {
                    reorderInto( in, inputSize.y, inputSize.x, tmpPtrC, inputSize.y, inputSize.x );
                } else {
                    std::copy_n( in, inPixels, tmpPtrC );
                }
                fftw_execute_dft( plan_full->forward_plan, reinterpret_cast<fftw_complex*>(tmpPtrC), reinterpret_cast<fftw_complex*>(out) );
                if( flags&REORDER_FT ) reorder( out, inputSize.y, inputSize.x );
            }
            if( flags&NORMALIZE_FT ) normalize( out, N );
            return flags;
        }
        template <> int FourierTransform::ft( const double* __restrict__ in, complex_t* __restrict__ out, int flags ) const {
            size_t N = inputSize.y * inputSize.x;
            if( !(flags&FULLCOMPLEX) ) {
                flags &= ~REORDER_FT;     // halfComplex can't be reordered....yet.
                N = inputSize.y * (inputSize.x/2+1);
                double* __restrict__ tmpPtrD = reinterpret_cast<double*>(tmpPtr);
                if( flags&REORDER_IMG ) {
                    reorderInto( in, inputSize.y, inputSize.x, tmpPtrD, inputSize.y, inputSize.x );
                } else {
                    tmpPtrD = const_cast<double*>(in);
                }
                fftw_execute_dft_r2c( plan_half->forward_plan, tmpPtrD, reinterpret_cast<fftw_complex*>(out) );
            } else {
                complex_t* tmpPtrC = tmpPtr;
                if( flags&REORDER_IMG ) {
                    reorderInto( in, inputSize.y, inputSize.x, tmpPtrC, inputSize.y, inputSize.x );
                } else {
                    std::copy_n( in, inPixels, tmpPtrC );
                }
                fftw_execute_dft( plan_full->forward_plan, reinterpret_cast<fftw_complex*>(tmpPtrC), reinterpret_cast<fftw_complex*>(out) );
                if( flags&REORDER_FT ) reorder( out, inputSize.y, inputSize.x );
            }
            if( flags&NORMALIZE_FT ) normalize( out, N );
            return flags;
        }
        template <> int FourierTransform::ft( const complex_t* __restrict__ in, complex_t* __restrict__ out, int flags ) const {
            size_t N = inputSize.y * inputSize.x;
            flags |= FULLCOMPLEX;     // calling FT with complex data forces fullComplex
            complex_t* inPtr = const_cast<complex_t*>(in);
            if( flags&REORDER_IMG ) {
                inPtr = reinterpret_cast<complex_t*>(tmpPtr);
                reorderInto( in, inputSize.y, inputSize.x, inPtr, inputSize.y, inputSize.x );
            }
            fftw_execute_dft( plan_full->forward_plan, reinterpret_cast<fftw_complex*>(inPtr), reinterpret_cast<fftw_complex*>(out) );
            if( flags&REORDER_FT ) reorder( out, inputSize.y, inputSize.x );
            if( flags&NORMALIZE_FT ) normalize( out, N );
            return flags;
        }
        template <typename T>
        int FourierTransform::ift( const complex_t* in, T* out, int flags ) const {
            complex_t* inPtr = const_cast<complex_t*>(in);        // cast away the const, since FFTW needs non-const arrays.
            if( !(flags&FULLCOMPLEX) ) {
                double* tmpPtrD = reinterpret_cast<double*>(tmpPtr);
                fftw_execute_dft_c2r( plan_half->backward_plan, reinterpret_cast<fftw_complex*>(inPtr), tmpPtrD );
                std::copy_n( tmpPtrD, inPixels, out );
            } else {
                complex_t* thisFtPtr = inPtr;
                if( flags&REORDER_FT ) {    // N.B. ift uses both tmpData & tmpData2 if FULLCOMPLEX && REORDER_FT !!
                    reorderInto( in, inputSize.y, inputSize.x, tmpPtr2 );
                    thisFtPtr = tmpPtr2;
                }
                complex_t* tmpPtrC = tmpPtr;
                fftw_execute_dft( plan_full->backward_plan, reinterpret_cast<fftw_complex*>(thisFtPtr),
                                  reinterpret_cast<fftw_complex*>(tmpPtrC) );
                std::transform( tmpPtrC, tmpPtrC+inPixels, out, [](const complex_t& c){ return c.real(); } );
            }
            if( flags&NORMALIZE_FT ) normalize( out, inPixels );
            if( flags&REORDER_IMG ) reorder( out, inputSize.y, inputSize.x );
            return flags;
        }
        template <> int FourierTransform::ift( const complex_t* in, double* out, int flags ) const {
            complex_t* inPtr = const_cast<complex_t*>(in);        // cast away the const, since FFTW needs non-const arrays.
            if( !(flags&FULLCOMPLEX) ) {
                fftw_execute_dft_c2r( plan_half->backward_plan, reinterpret_cast<fftw_complex*>(inPtr), out );
            } else {
                complex_t* thisFtPtr = inPtr;
                if( flags&REORDER_FT ) {    // N.B. ift uses both tmpData & tmpData2 if FULLCOMPLEX && REORDER_FT !!
                    reorderInto( in, inputSize.y, inputSize.x, tmpPtr2 );
                    thisFtPtr = tmpPtr2;
                }
                complex_t* tmpPtrC = tmpPtr;
                fftw_execute_dft( plan_full->backward_plan, reinterpret_cast<fftw_complex*>(thisFtPtr),
                                  reinterpret_cast<fftw_complex*>(tmpPtrC) );
                std::transform( tmpPtr, tmpPtr+inPixels, out, [](const complex_t& c){ return c.real(); } );
            }
            if( flags&NORMALIZE_FT ) normalize( out, inPixels );
            if( flags&REORDER_IMG ) reorder( out, inputSize.y, inputSize.x );
            return flags;
        }
        template <> int FourierTransform::ift( const complex_t* in, complex_t* out, int flags ) const {
            complex_t* inPtr = const_cast<complex_t*>(in);        // cast away the const, since FFTW needs non-const arrays.
            if( !(flags&FULLCOMPLEX) ) {
                double* tmpPtrD = reinterpret_cast<double*>(tmpPtr);
                fftw_execute_dft_c2r( plan_half->backward_plan, reinterpret_cast<fftw_complex*>(inPtr), tmpPtrD );
                std::copy_n( tmpPtrD, inPixels, out );
            } else {
                complex_t* thisFtPtr = inPtr;
                if( flags&REORDER_FT ) {
                    reorderInto( in, inputSize.y, inputSize.x, tmpPtr );
                    thisFtPtr = tmpPtr;
                }
                fftw_execute_dft( plan_full->backward_plan, reinterpret_cast<fftw_complex*>(thisFtPtr),
                                  reinterpret_cast<fftw_complex*>(out) );
            }
            if( flags&NORMALIZE_FT ) normalize( out, inPixels );
            if( flags&REORDER_IMG ) reorder( out, inputSize.y, inputSize.x );
            return flags;
        }
        template <typename T>
        void FourierTransform::init (const T* in, size_t ySize, size_t xSize, int flags, uint8_t nT) {
            init( ySize, xSize, flags, nT );
            currentFlags = ft<T>( in, ftPtr, flags );
            centered = currentFlags&REORDER_FT;
            normalized = currentFlags&NORMALIZE_FT;
        }
        template <>
        void FourierTransform::init( const std::complex<double>* in, size_t ySize, size_t xSize, int flags, uint8_t nT ) {
            flags |= FULLCOMPLEX;        // force full-complex for complex input.
            init( ySize, xSize, flags, nT );
            currentFlags = ft( in, ftPtr, flags );
            centered = currentFlags&REORDER_FT;
            normalized = currentFlags&NORMALIZE_FT;
        }
        
        template <>
        void FourierTransform::init( const std::complex<float>* in, size_t ySize, size_t xSize, int flags, uint8_t nT ) {
            flags |= FULLCOMPLEX;        // force full-complex for complex input.
            init( ySize, xSize, flags, nT );
            complex_t* __restrict__ ptr = tmpPtr2;
            std::copy_n( in, inPixels, ptr );
            currentFlags = ft( ptr, ftPtr, flags );
            centered = currentFlags&REORDER_FT;
            normalized = currentFlags&NORMALIZE_FT;
        }
    }
}
template int FourierTransform::ft( const int8_t*, complex_t*, int ) const;
template int FourierTransform::ft( const int16_t*, complex_t*, int ) const;
template int FourierTransform::ft( const int32_t*, complex_t*, int ) const;
template int FourierTransform::ft( const unsigned int*, complex_t*, int ) const;
template int FourierTransform::ft( const float*, complex_t*, int ) const;
template int FourierTransform::ift( const complex_t*, float* out, int ) const;
template void FourierTransform::init( const int8_t*, size_t, size_t, int, uint8_t );
template void FourierTransform::init( const uint8_t*, size_t, size_t, int, uint8_t );
template void FourierTransform::init( const int16_t*, size_t, size_t, int, uint8_t );
template void FourierTransform::init( const uint16_t*, size_t, size_t, int, uint8_t );
template void FourierTransform::init( const int32_t*, size_t, size_t, int, uint8_t );
template void FourierTransform::init( const uint32_t*, size_t, size_t, int, uint8_t );
template void FourierTransform::init( const int64_t*, size_t, size_t, int, uint8_t );
template void FourierTransform::init( const uint64_t*, size_t, size_t, int, uint8_t );
template void FourierTransform::init( const float*, size_t, size_t, int, uint8_t );
template void FourierTransform::init( const double*, size_t, size_t, int, uint8_t );


FourierTransform::FourierTransform() : centered(false), normalized(true),
    currentFlags(0), nThreads(1), inputSize(0), ftSize(0), inPixels(0), ftPixels(0), currentBlockSize(0),
    ftPtr(nullptr), tmpPtr(nullptr), tmpPtr2(nullptr) {

}


FourierTransform::FourierTransform( size_t ySize, size_t xSize, int flags, uint8_t nT ) :
    centered(false), normalized(false), currentFlags(flags), nThreads(nT),
    inputSize(0), ftSize(0), inPixels(0), ftPixels(0), currentBlockSize(0),
    ftPtr(nullptr), tmpPtr(nullptr), tmpPtr2(nullptr) {

    init( ySize, xSize, flags, nT );
    
    // since we were called without data to initialize on: set "centered" if REORDER_FT was passed.
    if( currentFlags&FULLCOMPLEX && currentFlags&REORDER_FT ) {
        centered = true;
    }

}


FourierTransform::FourierTransform( const FourierTransform& rhs ) : plan_full(rhs.plan_full), plan_half(rhs.plan_half),
    centered(rhs.centered), normalized(rhs.normalized),
    currentFlags(rhs.currentFlags), nThreads(rhs.nThreads), inputSize(rhs.inputSize), ftSize(rhs.ftSize),
    inPixels(rhs.inPixels), ftPixels(rhs.ftPixels), currentBlockSize(0),
    ftPtr(rhs.ftPtr), tmpPtr(rhs.tmpPtr), tmpPtr2(rhs.tmpPtr2) {
        
    resize( inPixels );
    std::copy_n( rhs.ftPtr, ftPixels, ftPtr );
    wrap();

}


FourierTransform::FourierTransform( FourierTransform&& rhs ) : plan_full(rhs.plan_full), plan_half(rhs.plan_half),
    centered(rhs.centered), normalized(rhs.normalized), nThreads(rhs.nThreads), inputSize(rhs.inputSize),
    ftSize(rhs.ftSize), inPixels(rhs.inPixels), ftPixels(rhs.ftPixels), currentBlockSize(0),
    ftData(rhs.ftData), ftPtr(rhs.ftPtr), tmpPtr(rhs.tmpPtr), tmpPtr2(rhs.tmpPtr2) {

    wrap();

}


void FourierTransform::getIFTx( double* out ) {

    fftw_execute_dft_c2r( plan_half->backward_plan, reinterpret_cast<fftw_complex*>(ftPtr), out);

}


void FourierTransform::getIFTx( complex_t* out ) const {

    if( centered ) {
        FourierTransform tmp(*this);
        tmp.decenter(true);
        tmp.getIFTx( out );
        return;
    }
    fftw_complex* dataPtr = reinterpret_cast<fftw_complex*>(const_cast<complex_t*>(ftPtr));
    fftw_execute_dft( plan_full->backward_plan, dataPtr, reinterpret_cast<fftw_complex*>(out) );

}


void FourierTransform::getIFTn( complex_t* out ) const {

    complex_t* thisFtPtr = const_cast<complex_t*>(ftPtr);
    if( centered ) {
        thisFtPtr = const_cast<complex_t*>(tmpPtr);
        reorderInto( ftPtr, inputSize.y, inputSize.x, thisFtPtr );
    }

    fftw_execute_dft( plan_full->backward_plan, reinterpret_cast<fftw_complex*>(thisFtPtr),
                      reinterpret_cast<fftw_complex*>(out) );

}


template <typename T>
Array<T> FourierTransform::power( bool center ) const {
    Array<T> tmp( ftSize.y, ftSize.x );
    power( tmp.get(), center );
    return tmp;
}
template Array<float> FourierTransform::power(bool) const;
template Array<double> FourierTransform::power(bool) const;


double FourierTransform::noise( int mask, double limit ) const {

    double noise_power = 0.0;

    int N = 0;
    double xx, yy;
    const complex_t* dataPtr = ftPtr;

    if( mask < 0 ) {         // specifying mask < 0 will match MvN's noise_level computation   ( mask = (np/6) )
        mask = std::max (ftSize.y, ftSize.x) / 6 + 1;
    }

    double limits = limit * limit;
    if( limit < 0 ) {        // specifying limit < 0 will match MvN's noise_level computation  ( sqrt((npY*npX)/4.0) )
        limits = inPixels / 4.0;
    }

    int endX = ftSize.x;
    int endY = ftSize.y;

    if( mask > 0 ) {
        endY = ftSize.y - mask + 1;
        if( (currentFlags&FULLCOMPLEX) ) {
            endX = ftSize.x - mask + 1;
        }
    }

    for( int x = mask; x < endX; ++x ) {
        if( !(currentFlags&FULLCOMPLEX) ) {
            xx = x;
        } else {
            xx = std::min (x, ftSize.x - x);
        }

        xx *= xx;
        int count = 1 + (!(currentFlags&FULLCOMPLEX) && x > 0 && x < ftSize.x - 1); // if !fullComplex and x is not first/last column, value is counted twice.

        for (int y = mask; y < endY; ++y) {
            yy = std::min (y, ftSize.y - y);
            if ( (xx + yy * yy) < limits) continue;         // inside region to ignore
            noise_power += count * std::norm(dataPtr[y*ftSize.x+x]);
            N += count;
        }
    }
    if( N && inPixels ) noise_power = std::sqrt (noise_power / (N * inPixels));
    return noise_power;

}


void FourierTransform::normalize( FourierTransform& in ) {
    if( in.normalized ) return;
    in *= (1.0 / (in.inPixels));
    in.normalized = true;
}


void FourierTransform::center( bool force ) {
    if( force || !centered ) {
        reorder();
        centered = true;
    }
}


void FourierTransform::decenter( bool force ) {
    if( force || centered ) {
        reorder();
        centered = false;
    }
}


FourierTransform FourierTransform::reordered (void) const {
    FourierTransform tmp (*this);
    tmp.center();
    return tmp;
}


void FourierTransform::resize( size_t blockSize ) {

    if( blockSize == 0 ) {
        ftData.reset();
        ftSize = inputSize = 0;
        inPixels = ftPixels = 0;
    } else {
        ftData = rdx_get_shared<complex_t>(3*blockSize);
        ftPtr = ftData.get();
        tmpPtr = ftPtr + blockSize;
        tmpPtr2 = tmpPtr + blockSize;
    }
    currentBlockSize = blockSize;
}


void FourierTransform::wrap( void ) {
    
    Array<complex_t>::wrap( ftPtr, ftSize.y, ftSize.x );
    
}


FourierTransform& FourierTransform::operator*= (const FourierTransform& rhs) {

    if( inputSize != rhs.inputSize ) { // size mismatch, don't do anything.
        return *this;
    }
    // FIXME: This only wokrs for non-centered FT's now !!!
    complex_t* thisPtr = ftPtr;
    const complex_t* rhsPtr = rhs.ftPtr;
    // Note: ftSize.y is always equal to rhs.ftSize.y if inputSizes matched above!
    if( ftPixels == rhs.ftPixels ) {                // same size, direct multiplication
        std::transform( thisPtr, thisPtr+ftPixels, rhsPtr, thisPtr, std::multiplies<complex_t>() ); 
    } else if( ftPixels < rhs.ftPixels ) {          // this is halfComplex (ftSize.x < rhs.ftSize.x), use half of rhs
        for( int y(0); y<ftSize.y; ++y ) {
            std::transform( thisPtr, thisPtr+ftSize.x, rhsPtr, thisPtr, std::multiplies<complex_t>() ); 
            thisPtr += ftSize.x;
            rhsPtr += rhs.ftSize.x;
        }
    } else {                                        // rhs is halfComplex (ftSize.x > rhs.ftSize.x), use half of this
        auto thisPtr2D = redux::util::makePointers(thisPtr, ftSize.y, ftSize.x );
        auto rhsPtr2D = redux::util::makePointers(rhsPtr, ftSize.y, rhs.ftSize.x );
        for (int y = 0; y < ftSize.y; ++y) {
            for (int x = 0; x < rhs.ftSize.x; ++x) {
                complex_t val = rhsPtr2D[y][x];
                thisPtr2D[y][x] *= val;
                if( x ) {
                    val = std::conj(val);
                    if( y ) thisPtr2D[ftSize.y-y][ftSize.x-x] *= val;
                    else thisPtr2D[0][ftSize.x-x] *= val;
                }
            }
        }
        redux::util::delPointers(thisPtr2D);
        redux::util::delPointers(rhsPtr2D);
    }

    return *this;
}


FourierTransform& FourierTransform::operator=(const FourierTransform& rhs) {
    
    inputSize = rhs.inputSize;
    inPixels = rhs.inPixels;
    ftSize = rhs.ftSize;
    plan_full = rhs.plan_full;
    plan_half = rhs.plan_half;
    centered = rhs.centered;
    currentFlags = rhs.currentFlags;
    normalized = rhs.normalized;
    size_t newBlockSize = inputSize.y * inputSize.x;
    resize( newBlockSize );
    std::copy_n( rhs.ftPtr, newBlockSize, ftPtr );
    wrap();
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


