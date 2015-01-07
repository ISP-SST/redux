#include "redux/momfbd/workspace.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/defines.hpp"
#include "redux/image/utils.hpp"

#include "redux/file/fileana.hpp"

using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

const std::string thisChannel = "workspace";

}


WorkSpace::WorkSpace( const MomfbdJob& j, PatchData::Ptr d ) : data( d ), cfg( j ) {
    //std::cout << "WorkSpace():  first = (" << d->first.x << "," << d->first.y << ")  last = (" << d->last.x << "," << d->last.y;
    //std::cout << ")   id = (" << d->index.x << "," << d->index.y << ")" << std::endl;
    int nPixels = std::min( d->nPixelsY(),d->nPixelsX() );
    window.resize( nPixels,nPixels );
    window = 1.0;
    window = redux::image::apodize( window,nPixels/8 );
    //redux::file::Ana::write( "window.f0", window );

    int md = std::min( 256, nPixels );                       // new size (maximum=256)
    md -= ( md%2 );
    //int m_offs = std::max( ( nPixels-md )/2,0 );            // centering needed?
    noiseWindow.resize( md,md );
    noiseWindow = 1.0;
    noiseWindow = redux::image::apodize( noiseWindow, md/16 );
    //redux::file::Ana::write( "noisewindow.f0", noiseWindow );

    for( auto & it: cfg.getObjects() ) {
        ObjPtr obj( new ObjectData( *this, it ) );
        objects.push_back( obj );
    }


}


WorkSpace::~WorkSpace() {
    clear();
  //  std::cout << "~WorkSpace():  1   wfsz = " << wavefronts.size() << " osz = " << objects.size() << std::endl;

}


void WorkSpace::init( boost::asio::io_service& service ) {
    for( auto & it: objects ) {
        it->prepareData( data, wavefronts, service );
    }
}


void WorkSpace::clear( void ) {
    for( auto & it: objects ) {
        it->clear();
    }
    objects.clear();
    wavefronts.clear();
}


void WorkSpace::collectResults( void ) {

    //data->images.resize();      // don't send raw data back.
}






