#include "redux/image/fouriertransform.hpp"

#include <mutex>
#include <set>

using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

FourierTransform::Plan::Plan( const std::vector<size_t>& dims, TYPE tp ) : tp(tp) {
    if ( !dims.empty() ) {
        sizes.clear();
        for( auto it: dims ) {
            if( it > 1 ) {
                sizes.push_back(it);
            }
        }
        if( ! sizes.empty() ) init();
        else throw std::logic_error( "FT::Plan constructed with only trivial dimensions: " + printArray( dims, "dims" ) );
    } else throw std::logic_error( "FT::Plan constructed without dimensions: " + printArray( dims, "dims" ) );
}


FourierTransform::Plan::~Plan() {
    if( ! sizes.empty() ) {
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(backward_plan);
    }
}


void FourierTransform::Plan::init( void ) {
    if( tp == R2C ) {
        if( sizes.size() == 2 ) {
            auto in = sharedArray<double>( sizes[0], sizes[1] );
            auto out = sharedArray<fftw_complex>( sizes[0], sizes[1] / 2 + 1 );
            forward_plan = fftw_plan_dft_r2c_2d( sizes[0], sizes[1], *in.get(), *out.get(), FFTW_MEASURE );
            backward_plan = fftw_plan_dft_c2r_2d( sizes[0], sizes[1], *out.get(), *in.get(), FFTW_MEASURE );
        }
        else if (sizes.size() == 1) {
            auto in = sharedArray<double>( sizes[0] );
            auto out = sharedArray<fftw_complex>( sizes[0] / 2 + 1 );
            forward_plan = fftw_plan_dft_r2c_1d( sizes[0], in.get(), out.get(), FFTW_MEASURE );
            backward_plan = fftw_plan_dft_c2r_1d( sizes[0] / 2 + 1, out.get(), in.get(), FFTW_MEASURE );
        } else throw std::logic_error( "FT::Plan::init() is only implemented for 1/2 dimensions, add more when/if needed: " + printArray( sizes, "dims" ) );
    } else if( tp == C2C )  {
        if( sizes.size() == 2 ) {
            auto in = sharedArray<fftw_complex>( sizes[0], sizes[1] );
            auto out = sharedArray<fftw_complex>( sizes[0], sizes[1] );
            forward_plan = fftw_plan_dft_2d( sizes[0], sizes[1], *in.get(), *out.get(), FFTW_FORWARD, FFTW_MEASURE );
            backward_plan = fftw_plan_dft_2d( sizes[0], sizes[1], *in.get(), *out.get(), FFTW_BACKWARD, FFTW_MEASURE );
        }
        else if (sizes.size() == 1) {
            auto in = sharedArray<fftw_complex>( sizes[0] );
            auto out = sharedArray<fftw_complex>( sizes[0] );
            forward_plan = fftw_plan_dft_1d( sizes[0], in.get(), out.get(), FFTW_FORWARD, FFTW_MEASURE );
            backward_plan =  fftw_plan_dft_1d( sizes[0], in.get(), out.get(), FFTW_BACKWARD, FFTW_MEASURE );
        } else throw std::logic_error( "FT::Plan::init() is only implemented for 1/2 dimensions, add more when/if needed: " + printArray( sizes, "dims" ) );
    } else throw std::logic_error( "FT::Plan::init() unknown tp." );
}

namespace {
    struct PlanCompare {
        bool operator()( const FourierTransform::Plan::Ptr &a, const FourierTransform::Plan::Ptr &b ) const { return ( *a < *b ); }
    };
}

FourierTransform::Plan::Ptr FourierTransform::getPlan( const std::vector<size_t>& dims, Plan::TYPE tp ) {
    static std::map<std::pair<std::vector<size_t>, Plan::TYPE>,Plan::Ptr> plans;
    static mutex mtx;
    unique_lock<mutex> lock( mtx );
    auto plan = plans.emplace(make_pair(dims, tp),nullptr);
    if( plan.second ) { // insertion successful, i.e. new plan -> initialize it
        plan.first->second.reset( new Plan(dims, tp) );
    }
    return plan.first->second;
}


template <typename T>
FourierTransform::FourierTransform( const Array<T>& rhs, int flags ) : centered(false), halfComplex(true), normalized(false), reOrdered(false) {

    vector<size_t> dims = rhs.dimensions();
    if( dims.empty() ) throw logic_error( "FourierTransform called with no non-trivial dimensions: " + printArray(rhs.dimensions(), "dims") );
    else if( dims.size() > 2 ) throw logic_error( "FourierTransform only supports 1/2 dimensions at the moment: " + printArray(dimensions(), "dims") );

    Array<double> tmp;
    rhs.copy( tmp );
    
    nInputElements = rhs.nElements();
    
    if( flags & FT_REORDER ) {
        reorder(tmp);
        reOrdered = true;
    }
    
    init(tmp);

    if( flags & FT_NORMALIZE ) {
       normalize();
    }
    

    
}
namespace redux { namespace image {
template <>
    FourierTransform::FourierTransform( const Array<complex_t>& rhs, int flags ) : centered(false), halfComplex(false), normalized(false), reOrdered(false) {

    vector<size_t> dims = rhs.dimensions();
    if( dims.empty() ) throw logic_error( "FourierTransform called with no non-trivial dimensions: " + printArray(rhs.dimensions(), "dims") );
    else if( dims.size() > 2 ) throw logic_error( "FourierTransform only supports 1/2 dimensions at the moment: " + printArray(dimensions(), "dims") );
    
    Array<complex_t> tmp;
    rhs.copy( tmp );
    
    nInputElements = rhs.nElements();

    if( flags & FT_REORDER ) {
        reorder(tmp);
        reOrdered = true;
    }
    
    init(tmp);
    
    if( flags & FT_NORMALIZE ) {
       normalize();
    }

}
}}
template FourierTransform::FourierTransform( const Array<int16_t>&, int );
template FourierTransform::FourierTransform( const Array<int32_t>&, int );
template FourierTransform::FourierTransform( const Array<float>&, int );
template FourierTransform::FourierTransform( const Array<double>&, int );
template FourierTransform::FourierTransform( const Array<complex_t>&, int );


void FourierTransform::init( Array<double>& rhs ) {

    vector<size_t> dims = rhs.dimensions();
    dims.back() >>= 1;      // for r2c transforms, the last dimension has approx. half size =(n/2+1)
    dims.back() += 1; 
    
    resize( dims );
    plan = getPlan( rhs.dimensions(), Plan::R2C );
    fftw_execute_dft_r2c( plan->forward_plan, rhs.ptr(), reinterpret_cast<fftw_complex*>( ptr() ) );

}

void FourierTransform::init( Array<complex_t>& rhs ) {
    
    resize( rhs.dimensions() );
    
    plan = getPlan( rhs.dimensions(), Plan::C2C );
    fftw_execute_dft( plan->forward_plan, reinterpret_cast<fftw_complex*>( rhs.ptr() ), reinterpret_cast<fftw_complex*>( ptr() ) );

}

template <typename T>
void FourierTransform::inv( Array<T>& out ) {
    vector<size_t> dims = dimensions();
    if(halfComplex) {
        dims.back() -= 1; 
        dims.back() <<= 1;      // for r2c transforms, the last dimension has half size (n/2+1)
    }
    if( centered ) reorder();   // TBD: temporary fft-copy or reorder twice ?
    if(halfComplex) {
        Array<double> tmp(dims);
        fftw_execute_dft_c2r( plan->backward_plan, reinterpret_cast<fftw_complex*>( ptr() ), tmp.ptr());
        out.resize( dims );
        tmp.copy(out);
    } else {
        dims.push_back(2);      // add extra dim for real/imaginary
        Array<double> tmp(dims);
        out.resize( dims );
        fftw_execute_dft( plan->backward_plan, reinterpret_cast<fftw_complex*>( ptr() ), reinterpret_cast<fftw_complex*>( tmp.ptr() ));
        tmp.copy(out);   // TBD: Add extra dimension of size 2, or discard the imaginary part?
    }
    
    if( !normalized ) {
        out /= out.nElements();
    }
    
/*
    if( flags & FT_REORDER ) {
        reorder(rhs);
    }
    
    if( flags & FT_NORMALIZE ) {
        normalize(rhs);
    }
    
*/

}
template void FourierTransform::inv( Array<float>& in );
template void FourierTransform::inv( Array<double>& in );
template void FourierTransform::inv( Array<complex_t>& in );
            
template <typename T>
Array<T> FourierTransform::correlate( const Array<T>& in ) const {
    
    int f(0);
    if( ! reOrdered ) {    // only reorder inFT if *this is not already reordered
        f = ~ FT_REORDER;
    }
    FourierTransform inFT( in, f );
    //inFT *= *this;
    auto it = inFT.begin();
    for( auto& it2: *this ) *it++ *= conj(it2);
    
    Array<T> out( in.dimensions() );
    inFT.inv(out);
    reorder(out);
    
    return out;
    
}
template Array<float> FourierTransform::correlate( const Array<float>& ) const;
template Array<double> FourierTransform::correlate( const Array<double>& ) const;
template Array<complex_t> FourierTransform::correlate( const Array<complex_t>& ) const;
            


void FourierTransform::autocorrelate( void ) {
    for( auto& it: *this ) it = norm(it);
}


Array<double> FourierTransform::power( void ) const {
    
    Array<double> tmp(dimensions());
    auto it = tmp.begin();
    for( auto& it2: *this ) *it++ = norm(it2);
   
    return tmp;
    
}

void FourierTransform::convolveInPlace( Array<double>& in, int flags ) const {
    
    int f = flags;
    if( reOrdered == (f & FT_REORDER)) {    // only reorder inFT if *this is not already reordered
        f &= ~(f & FT_REORDER);
    }
    FourierTransform inFT( in, f );
    if( centered != inFT.centered ) inFT.reorder();
    inFT *= *this;
    if( inFT.centered ) inFT.reorder();
    fftw_execute_dft_c2r( plan->backward_plan, reinterpret_cast<fftw_complex*>( inFT.ptr() ), in.ptr() );
    if( (reOrdered == inFT.reOrdered) != (flags & FT_REORDER)) reorder(in);

}

void FourierTransform::convolveInPlace( Array<complex_t>& in, int flags ) const {
    int f = flags;
    if( reOrdered == (f & FT_REORDER)) {    // only reorder inFT if *this is not already reordered
        f &= ~(f & FT_REORDER);
    }
    FourierTransform inFT( in, f );
    if( centered != inFT.centered ) inFT.reorder();
    inFT *= *this;
    if( inFT.centered ) inFT.reorder();
    fftw_execute_dft( plan->backward_plan, reinterpret_cast<fftw_complex*>( inFT.ptr() ), reinterpret_cast<fftw_complex*>( in.ptr() ) );
    if( (reOrdered == inFT.reOrdered) != (flags & FT_REORDER)) reorder(in);
}


void FourierTransform::convolveInPlaceHC( Array<complex_t>& in, int flags ) const {
    int f = flags;
    if( reOrdered == (f & FT_REORDER)) {    // only reorder inFT if *this is not already reordered
        f &= ~(f & FT_REORDER);
    }
    FourierTransform inFT( in, f );
    if( centered != inFT.centered ) inFT.reorder();
cout << "FIX MULTIPLICATION" << endl;
return;
    inFT *= *this;
    if( inFT.centered ) inFT.reorder();
    fftw_execute_dft( plan->backward_plan, reinterpret_cast<fftw_complex*>( inFT.ptr() ), reinterpret_cast<fftw_complex*>( in.ptr() ) );
    if( (reOrdered == inFT.reOrdered) != (flags & FT_REORDER)) reorder(in);

}

Array<double> FourierTransform::convolve( const Array<double>& in, int flags ) const {
    
    auto dims = in.dimensions();
    Array<double> ret( dims );
    if( dims == dimensions() ) {    // matching full complex dimensions -> convolving real function with complex -> result will be complex
        Array<complex_t> tmp( dims );
        in.copy( tmp );
        convolveInPlace(tmp,flags);
        auto it = ret.begin();
        for( auto& it2: tmp ) *it++ = real(it2);
    } else {
        dims.back() >>= 1;      // for r2c transforms, the last dimension has half size =(n/2+1)
        dims.back() += 1; 
        if( dims == dimensions() ) {    // matching halfcomplex dimensions -> convolving two real functions
            in.copy( ret );
            convolveInPlace(ret,flags);
        } else cout << "Dims doesn't match real/complex fft, returning empty result !!!"<<printArray(dims,"  dims") << printArray(dimensions(),"  ldims") << endl;
    }

    return ret;
    
}

Array<complex_t> FourierTransform::convolve( const Array<complex_t>& in, int flags ) const {

cout << "cvlv 4" << endl;

    auto dims = in.dimensions();
    Array<complex_t> ret( dims );
    if( dims == dimensions() ) {    // matching full complex dimensions -> convolving complex functions
        convolveInPlace(ret,flags);
    } else {
        dims.back() >>= 1;      // for r2c transforms, the last dimension has approx. half size =(n/2+1)
        dims.back() += 1; 
        if( dims == dimensions() ) {    // matching halfcomplex dimensions -> convolving half-complex w. full-complex
            convolveInPlaceHC(ret,flags);
        } else cout << "Dims doesn't match real/complex fft, returning empty result !!!" <<printArray(dims,"  dims") << printArray(dimensions(),"  ldims") << endl;
    }
    return ret;

}

/*
template <typename T>
Array<T> FourierTransform::convolve( const Array<T>& in, int flags ) const {
cout << "cvlv 5  centered " << centered << "  reOrdered " << reOrdered << endl;
    Array<double> tmp;
    in.copy( tmp );
    convolveInPlace(tmp,flags);
    Array<T> ret( tmp.dimensions() );
    tmp.copy( ret );
    return ret;
}

template <typename T>
void FourierTransform::convolveInPlace( Array<T>& in, int flags ) const {
cout << "cvlv 6  centered " << centered << "  reOrdered " << reOrdered << endl;
    Array<double> tmp;
    in.copy( tmp );
    convolveInPlace(tmp,flags);
    in.resize( tmp.dimensions() );
    tmp.copy( in );
}
template Array<float> FourierTransform::convolve( const Array<float>& in, int ) const;
template Array<int> FourierTransform::convolve( const Array<int>& in, int ) const;
template void FourierTransform::convolveInPlace( Array<float>& in, int ) const;
template void FourierTransform::convolveInPlace( Array<int>& in, int ) const;
*/

void FourierTransform::normalize( FourierTransform& in ) {
    in *= (1.0/in.nInputElements);
    in.normalized = true;
}

template <typename T>
void FourierTransform::reorder( redux::util::Array<T>& in ) {

    size_t halfY = in.dimSize(0) / 2;
    size_t stride = in.dimSize(1);
    size_t halfX = stride / 2;
    size_t chunkSize = halfX * sizeof( T );
    char *buf = new char[chunkSize];
    
    T* southWest = in.ptr();
    T* southEast = in.ptr(0,halfX);
    T* northWest = in.ptr(halfY,0);
    T* northEast = in.ptr(halfY,halfX);
    for( size_t y=0; y<halfY; ++y ) {
        memcpy( buf, southWest, chunkSize );
        memcpy( southWest, northEast, chunkSize );
        memcpy( northEast, buf, chunkSize );
        southWest += stride;
        northEast += stride;
        memcpy( buf, southEast, chunkSize );
        memcpy( southEast, northWest, chunkSize );
        memcpy( northWest, buf, chunkSize );
        southEast += stride;
        northWest += stride;
    }
    delete[] buf;

}
template void FourierTransform::reorder( redux::util::Array<float>& in );
template void FourierTransform::reorder( redux::util::Array<double>& in );
template void FourierTransform::reorder( redux::util::Array<complex_t>& in );



