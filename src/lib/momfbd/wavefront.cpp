#include "redux/momfbd/wavefront.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/object.hpp"

#include "redux/file/fileana.hpp"
#include "redux/logging/logger.hpp"
#include "redux/util/stringutil.hpp"

using namespace redux::momfbd;
using namespace redux::file;
using namespace redux::logging;
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
        ret += im->idString() + (string)im->currentShift;
    }
    ret += "]";
    return ret;
    
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

WaveFronts::WaveFronts( MomfbdJob& j ) : myJob(j), logger(j.logger), nModes(0), nWaveFronts(0), coefficients() {


}


void WaveFronts::maybeInitializeStorage( void ) {
    
    lock_guard<mutex> lock(mtx);
    if( coefficients.nElements() ) return;        // already allocated/loaded
    
    size_t nPatchesX = myJob.subImagePosX.size();
    size_t nPatchesY = myJob.subImagePosY.size();
    nModes = myJob.nModes;
    nWaveFronts = myJob.waveFrontList.size();
    
    if( !(myJob.runFlags&RF_NOSWAP) ) {     // should we mmap?
        if( bfs::exists( bfs::path(cacheFile) ) ) {
            LOG_DEBUG << "opening: storage: " << cacheFile << ende;
            coefficients.openMmap( cacheFile, nPatchesY, nPatchesX, nWaveFronts, nModes );
        } else {
            LOG_DEBUG << "creating: storage: " << cacheFile << ende;
            coefficients.createMmap( cacheFile, nPatchesY, nPatchesX, nWaveFronts, nModes );
        }
    } else {
        coefficients.resize(  nPatchesY, nPatchesX, nWaveFronts, nModes );
    }
    
}


void WaveFronts::getStorage( PatchData& pData ) {

    maybeInitializeStorage();
    WavefrontData& wfData = pData.waveFronts;
    wfData.alpha.wrap(reinterpret_cast<redux::util::Array<float>&>(coefficients),
         pData.index.y, pData.index.y, pData.index.x, pData.index.x, 0, nWaveFronts-1, 0, nModes-1 );

}


void WaveFronts::loadInit( boost::asio::io_service& service, Array<PatchData::Ptr>& patches ) {

    if( !(myJob.runFlags&RF_NOSWAP) ) {    // unless swap is deactivated for this job
        cacheFile = myJob.cachePath + "wavefronts";
    }

    bool hasInit(false);
    if( hasInit ) {
        maybeInitializeStorage();
        //loadInit( service, patches );
        for( unsigned int y=0; y<patches.dimSize(0); ++y ) {
            for( unsigned int x=0; x<patches.dimSize(1); ++x ) {
                service.post( std::bind( &WaveFronts::getStorage, this, std::ref(*(patches(y,x).get())) ) );
            }
        }
    }
}
