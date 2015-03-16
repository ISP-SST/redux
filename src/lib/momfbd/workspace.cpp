#include "redux/momfbd/workspace.hpp"

#include "redux/momfbd/momfbdjob.hpp"

#include "redux/image/utils.hpp"
#include "redux/logger.hpp"

using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

const std::string thisChannel = "workspace";

}


WorkSpace::WorkSpace( const MomfbdJob& j ) : result(new PatchResult(j)), cfg(j) {
    //std::cout << "WorkSpace():  first = (" << d->first.x << "," << d->first.y << ")  last = (" << d->last.x << "," << d->last.y;
    //std::cout << ")   id = (" << d->index.x << "," << d->index.y << ")" << std::endl;
    int nPixels = 1;//std::min( d->nPixelsY(),d->nPixelsX() );
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
/*
    for( auto & it: cfg.getObjects() ) {
        ObjPtr obj( new ObjectData( *this, it ) );
        objects.push_back( obj );
    }
*/

}


WorkSpace::~WorkSpace() {
    clear();
  //  std::cout << "~WorkSpace():  1   wfsz = " << wavefronts.size() << " osz = " << objects.size() << std::endl;

}


void WorkSpace::run( PatchData::Ptr p, boost::asio::io_service& service, uint8_t nThreads ) {
    data = p;
    result->id = data->id;
    result->index = data->index;
    result->pos = data->pos;
    LOG << "WorkSpace::run()  patch#" << data->id << "   index=" << data->index << " pos=" << data->pos;
usleep(10000);
    /*for( auto & it: objects ) {
        it->prepareData( data, wavefronts, service );
    }*/
}


void WorkSpace::clear( void ) {
    //for( auto & it: objects ) {
        //it->clear();
    //}
    objects.clear();
    wavefronts.clear();
}


PatchResult::Ptr& WorkSpace::getResult( void ) {
    return result;
    //data->images.resize();      // don't send raw data back.
}






