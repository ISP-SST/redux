#include "redux/image/fouriertransform.hpp"

#include <mutex>
#include <set>

using namespace redux::image;
using namespace redux::util;
using namespace std;

FourierTransform::Plan::Plan( const std::vector<size_t>& dims ) {
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
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
}


void FourierTransform::Plan::init( void ) {
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
}


const FourierTransform::Plan& FourierTransform::getPlan( const std::vector<size_t>& dims ) {
    static std::set<Plan> plans;
    static mutex mtx;
    unique_lock<mutex> lock( mtx );
    return *( plans.emplace( dims ).first );
}


template <typename T>
FourierTransform::FourierTransform( const Array<T>& rhs, int flags ) {

    Array<double> tmp;
    rhs.copy( tmp );
    vector<size_t> dims = tmp.dimensions();
    if( dims.empty() ) throw logic_error( "FourierTransform called with no non-trivial dimensions: " + printArray(rhs.dimensions(), "dims") );
    else if( dims.size() > 2 ) logic_error( "FourierTransform only supports 1/2 dimensions at the moment: " + printArray(dimensions(), "dims") );

    dims[dims.size()-1] = dims[dims.size()-1]/2 + 1; 
    
    resize( dims );

    if( flags & FT_REORDER ) {
        reorder(tmp);
    }
    
    if( flags & FT_NORMALIZE ) {
        normalize(tmp);
    }
    
    const Plan& plan = getPlan( tmp.dimensions() );
    fftw_execute_dft_r2c( plan.forward_plan, tmp.ptr(), reinterpret_cast<fftw_complex*>( ptr() ) );

}
template FourierTransform::FourierTransform( const Array<float>&, int );
template FourierTransform::FourierTransform( const Array<double>&, int );


template <typename T>
Array<T> FourierTransform::convolve( const Array<T>& in ) const {

    Array<double> tmp;
    in.copy( tmp );
    Array<T> ret( tmp.dimensions() );
    const Plan& plan = getPlan( tmp.dimensions() );
    FourierTransform ft( tmp );
    ft *= *this;
    fftw_execute_dft_c2r( plan.backward_plan, reinterpret_cast<fftw_complex*>( ft.ptr() ), tmp.ptr() );
    tmp.copy( ret );

    return ret;
}

template <>
Array<double> FourierTransform::convolve( const Array<double>& in ) const {

    Array<double> ret( in.dimensions() );
    const Plan& plan = getPlan( in.dimensions() );
    FourierTransform ft( in );
    ft *= *this;
    fftw_execute_dft_c2r( plan.backward_plan, reinterpret_cast<fftw_complex*>( ft.ptr() ), ret.ptr() );

    return ret;
}
template Array<float> FourierTransform::convolve( const Array<float>& in ) const;
template Array<double> FourierTransform::convolve( const Array<double>& in ) const;

template <typename T>
void FourierTransform::normalize( Array<T>& in ) {
    double sum = 0;
    for( auto it: in ) sum += it;
    in /= (sum*in.nElements());
}
template void FourierTransform::normalize( Array<float>& );
template void FourierTransform::normalize( Array<double>& );


template <typename T>
void FourierTransform::reorder( Array<T>& in ) {
    
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
}
template void FourierTransform::reorder( Array<float>& );
template void FourierTransform::reorder( Array<double>& );
template void FourierTransform::reorder( Array<complex_t>& );
