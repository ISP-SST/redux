#include "redux/momfbd/workspace.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/wavefront.hpp"

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


WorkSpace::WorkSpace( const MomfbdJob& j ) : objects( j.getObjects() ), result( new PatchResult(j) ), job(j), nFreeParameters(0) {
    //std::cout << "WorkSpace():  first = (" << d->first.x << "," << d->first.y << ")  last = (" << d->last.x << "," << d->last.y;
    //std::cout << ")   id = (" << d->index.x << "," << d->index.y << ")" << std::endl;

}


WorkSpace::~WorkSpace() {
    clear();
    //  std::cout << "~WorkSpace():  1   wfsz = " << wavefronts.size() << " osz = " << objects.size() << std::endl;

}


void WorkSpace::init( void ) {

    int nPixels = job.patchSize;//std::min( d->nPixelsY(),d->nPixelsX() );
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

    tilts.clear();
    wavefronts.clear();

    for( auto & it : job.getObjects() ) {
        it->initProcessing( shared_from_this() );
    }

    for( auto& it: job.modeNumbers ) {
    }


}


void WorkSpace::run( PatchData::Ptr p, boost::asio::io_service& service, uint8_t nThreads ) {

if (p->id != 1) return;


    for( auto& it: wavefronts ) {
        if(it.second) {
        } else LOG << "WorkSpace::run(1)  it.second is NULL";
    }
    p->initPatch();

    data = p;
    result->id = data->id;
    result->index = data->index;
    result->pos = data->pos;
    LOG << "WorkSpace::run()  patch#" << data->id << "   index=" << data->index << " pos=" << data->pos;

    /*

    ModeLoop:
    if shifted/new: Cut out sub-image and apply window to cutout, get stats, compute fft and add to fftSum
    Initialize PQ, get initial metric
    IterLoop
        Calculate gradient
            OTF/psf
            Line search
            Calculate metric/diff
        Check improvement
            Keep

            Discard & Restore previous
        Keep
            Calculate new alignment from tilt-coefficients
        Discard & Restore previous


    */



    coefficientMetric( service );

}


double WorkSpace::coefficientMetric(boost::asio::io_service& service) {
    double sum(0);
    for( auto& it: wavefronts ) {
        if(it.second) {
           // service.post( [it,&sum] {
                //sum.fetch_add(it.second->coefficientMetric());
                sum += it.second->coefficientMetric();
           // } );
        } else LOG << "WorkSpace::coefficientMetric()  it.second is NULL";
    }
    //runThreadsAndWait(service, job.info.maxThreads);
    cout << "WorkSpace::coefficientMetric()   sum = " << sum << endl;
    return sum;
}


double WorkSpace::objectMetric(boost::asio::io_service& service) {
    for( const shared_ptr<Object>& o: objects ) {
        if(o) {
            // async: clear and accumulate P & Q for all objects (blockwise?)
            service.post( [o] {
                o->initPQ();
                o->addAllPQ();
                //sum.fetch_add(it.second->coefficientMetric());
            } );
            // service.post( std::bind( &Object::getPQ, o.get() ) );
        } else LOG << "WorkSpace::objectMetric()  object is NULL";
    }
    runThreadsAndWait(service, job.info.maxThreads);
    double sum(0);
    for( const shared_ptr<Object>& o: objects ) {
        if(o) {
            // sync(or async?): calc/add metric for all objects (blockwise?)
        } else LOG << "WorkSpace::run()  object is NULL";
    }
    return sum;
}
           


void WorkSpace::clear( void ) {
    //for( auto & it: objects ) {
    //it->clear();
    //}

    wavefronts.clear();
}


void WorkSpace::resetAlpha( void ) {
    for( auto& it: wavefronts ) {
        if(it.second) {
            //it.second->setAlpha( alpha_init );
        } else LOG << "WorkSpace::resetAlpha()  it.second is NULL";
    }
}


PatchResult::Ptr& WorkSpace::getResult( void ) {
    return result;
    //data->images.resize();      // don't send raw data back.
}






