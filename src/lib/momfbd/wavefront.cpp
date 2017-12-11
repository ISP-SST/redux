#include "redux/momfbd/wavefront.hpp"

#include "redux/momfbd/object.hpp"

#include "redux/file/fileana.hpp"
#include "redux/util/stringutil.hpp"

using namespace redux::momfbd;
using namespace redux::file;
using namespace redux::util;
using namespace std;


WaveFront::WaveFront(uint32_t wfi) : wfIndex(wfi), alpha(nullptr), grad_alpha(nullptr), nAlpha(0) {

}


void WaveFront::setData( double* a, double* g, size_t n ) {
    
    alpha = a;
    grad_alpha=g;
    nAlpha = n;
    for( auto& im: images ) {
        im->setData(alpha);
    }
    
}


void WaveFront::clear( void ) {
    images.clear();
}


string WaveFront::print(void) {
    
    string ret = "[";
    bool first(true);
    for( auto& im: images ) {
        if(!first) ret += ", ";
        first = false;
        ret += im->idString() + (string)im->imageShift;
    }
    ret += "]";
    return std::move(ret);
    
}


void WaveFront::addImage( const shared_ptr<SubImage>& im ) {
    if( im ) {
        images.push_back(im);
    }
}


template <typename T>
void WaveFront::adjustShifts( const T* a ) const {
    for( auto& im: images ) {
        Object& o = im->object;
        o.imgShifted.fetch_or( im->adjustShifts(a) );
        ++o.progWatch;
    }
}
template void WaveFront::adjustShifts( const double* a ) const;
template void WaveFront::adjustShifts( const float* a ) const;


template <typename T>
void WaveFront::applyAlpha( const T* a ) const {
    for( auto& im: images ) {
        im->calcPhi(a);
        im->calcPFOTF();
    }
}
template void WaveFront::applyAlpha( const double* a ) const;
template void WaveFront::applyAlpha( const float* a ) const;


// void WaveFront::gradient( const bool* enabledModes, grad_tt grad ) {
// 
//     for( const shared_ptr<SubImage>& im: images ) {
//         im->calcVogelWeight();
//         grad( *im, grad_alpha, enabledModes );
//     }
// 
// }



void WaveFront::dump( std::string tag, bool dumpImages ) const {
    
    if( nAlpha ) {
        string wftag = tag + "_wf" + to_string(wfIndex);
        Ana::write( wftag + "_alpha.f0", alpha, nAlpha );
    }
    
    if( dumpImages ) {
        for( const shared_ptr<SubImage>& im: images ) {
            im->dump( tag );
        }
    }
    
}

