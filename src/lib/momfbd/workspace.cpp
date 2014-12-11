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


WorkSpace::Image::Image( const WorkSpace& ws, const ObjPtr& obj, const ChPtr& ch, const redux::util::Array<float>& stack,
                         uint32_t index, int firstY, int lastY, int firstX, int lastX )
    : Array<float>( stack, index, index, firstY, lastY, firstX, lastX ),
      index( index ), ws( ws ), object( obj ), channel( ch ) {

    // std::cout << "Image():  1   " << hexString(this) << std::endl;

}


WorkSpace::Image::~Image(void) {

  //  std::cout << "~Image():  1   " << hexString(this) << std::endl;
}


void WorkSpace::Image::init( void ) {

    auto dims = dimensions( true );
    img.resize( dims );
    for( auto& it: dims ) {
        it *= 2;
    }
    SJ.resize( dims );

    stats.getStats(*this, ST_VALUES);
    string str = "Image::init():  mean=" + to_string(stats.mean);

    img.zero();
    Array<double>::const_iterator wit = ws.window.begin();
    Array<float>::const_iterator dit = this->begin();
    for( auto& iit: img ) {                     // windowing: subtract and re-add mean afterwards
        iit = (*dit++ - stats.mean) * *wit++ + stats.mean;
    }
    
    stats.getStats(img, ST_VALUES|ST_RMS);
    str += "  mean2=" + to_string( stats.mean ) + "  std=" + to_string( stats.stddev );
    
    ft.reset( img );
//    FourierTransform::reorder( ft );
    stats.noise = channel->cfg->noiseFudge * ft.noise(-1,-1);       // mask/cutoff < 0 will revert to hardcoded values used by MvN
    str += "  noise=" + to_string( stats.noise );
    object->addToFT( ft );
/*    
    redux::file::Ana::write( "windowed_" + to_string( ws.data->index.x ) + "_" + to_string( ws.data->index.y ) +
                             "_" + to_string( index ) + ".f0", img );
    redux::file::Ana::write( "windowedft_" + to_string( ws.data->index.x ) + "_" + to_string( ws.data->index.y ) +
                             "_" + to_string( index ) + ".f0", ft );
*/
    //cout << str << endl;
    //std::cout << "Image::init():  E  " << std::endl;
}


void WorkSpace::Image::clear(void) {
    object.reset();
    channel.reset();
    wfg.reset();
}

WorkSpace::Object::Object( const WorkSpace& ws, const std::shared_ptr<ObjectCfg>& c ) : ws( ws ), cfg( c ) {
    //  std::cout << "Object():  " << hexString(this) << std::endl;

}


WorkSpace::Object::~Object() {
   // std::cout << "~Object():  1   "  << hexString(this) << printArray( ftSum.dimensions(),"   ftSum.dims" ) << std::endl;
    static int count( 0 );
    if( ftSum.nDimensions() > 1 ) {
        redux::file::Ana::write( "ftSum_" + to_string( count++ ) + ".f0", ftSum );
    }
}




void WorkSpace::Object::prepareData( const PatchData::Ptr& data, std::map<uint32_t, WfPtr>& wf, boost::asio::io_service& service  ) {
    for( auto & it: cfg->getChannels() ) {
        channels.push_back( ChPtr( new Channel( ws, shared(), it ) ) );
    }

    for( auto & it: channels ) {
        it->prepareData( data, wf, service );
    }
}


void WorkSpace::Object::clear(void) {
    for(auto &it: channels) it->clear();
    channels.clear();
    cfg.reset();
}


void WorkSpace::Object::addToFT( const redux::image::FourierTransform& ft ) {
    unique_lock<mutex> lock( mtx );
    if( !ftSum.sameSizes(ft) ) {
        ftSum.resize(ft.dimensions(true));
        ftSum.zero();
    }
    //  std::cout << "Object::addToFT():  1   " << hexString(this) << printArray(ft.dimensions(), "  ft.dims") << printArray(ftSum.dimensions(), "  ftSum.dims") << std::endl;
    ftSum += ft;
//    std::cout << "Object::addToFT():  2   " << printArray(ft.dimensions(), "ft.dims") << printArray(ftSum.dimensions(), "  ftSum.dims") << std::endl;
}


void WorkSpace::Object::addToPQ( const redux::image::FourierTransform& ft, const Array<complex_t> sj ) {
//   std::cout << "Object::addToPQ():  1   " << hexString(this) << printArray(ft.dimensions(), "  ft.dims") << printArray(sj.dimensions(), "  sj.dims") << std::endl;
    typedef Array<complex_t>::const_iterator cit;
    unique_lock<mutex> lock( mtx );
    cit ftit = ft.begin();
    cit sjit = sj.begin();
    auto qit = Q.begin();

    for( auto& pit: P ) {
        *qit++ += norm( *sjit );            // Q += sj.re^2 + sj.im^2 = norm(sj)
        pit += *ftit++ * *sjit++;           // P += ft * sj
    }
}


WorkSpace::Channel::Channel( const WorkSpace& ws, const ObjPtr& o, const std::shared_ptr<ChannelCfg>& c ) : ws( ws ), cfg( c ), obj( o ) {
    //  std::cout << "Channel():  1   " << hexString(this) << std::endl;

}


WorkSpace::Channel::~Channel(void) {
    // std::cout << "~Channel():  " << hexString(this) << std::endl;

}


void WorkSpace::Channel::prepareData( const PatchData::Ptr& data, std::map<uint32_t, WfPtr>& wf, boost::asio::io_service& service ) {
    uint32_t count=0;

    for( auto& it: cfg->imageNumbers ) {
        //   std::cout << "Channel::prepareData()  imn = " << it << std::endl;
        ImPtr im = ImPtr( new Image( ws, obj, shared(), data->images, cfg->imageOffset+count++,data->first.y,data->last.y,data->first.x,data->last.x ) );
        images.push_back( im->shared() );
//        service.post( std::bind( &WorkSpace::Image::init, im.get() ) );
im->init();

        WfPtr newWF( new WaveFront( im->shared() ) );
        im->setWaveFront(newWF);

        auto ret = wf.emplace( it, newWF );

        if( !ret.second ) {   // not inserted -> already existed -> add image to the same wavefront group
            ret.first->second->addImage( im );
        }
    }

//   std::cout << "Channel::prepareData()  done." << std::endl;
}


void WorkSpace::Channel::clear(void) {
    for(auto &it: images) it->clear();
    images.clear();
    cfg.reset();
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
        ObjPtr obj( new Object( *this, it ) );
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






