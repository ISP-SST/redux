#include "redux/momfbd/data.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/workspace.hpp"

#include "redux/image/utils.hpp"

#include "redux/file/fileana.hpp"

using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

const std::string thisChannel = "data";

}


ImageData::ImageData( const ObjPtr& obj, const ChPtr& ch, const redux::util::Array<float>& stack,
                         uint32_t index, int firstY, int lastY, int firstX, int lastX )
    : Array<float>( stack, index, index, firstY, lastY, firstX, lastX ),
      index( index ), object( obj ), channel( ch ) {

    // std::cout << "Image():  1   " << hexString(this) << std::endl;

}


ImageData::~ImageData(void) {

  //  std::cout << "~Image():  1   " << hexString(this) << std::endl;
}


void ImageData::init( void ) {

    auto dims = dimensions( true );
    img.resize( dims );
    for( auto& it: dims ) {
        it *= 2;
    }
    SJ.resize( dims );

    stats.getStats(*this, ST_VALUES);
    string str = "Image::init():  mean=" + to_string(stats.mean);

    img.zero();
/*    Array<double>::const_iterator wit = ws.window.begin();
    Array<float>::const_iterator dit = this->begin();
    for( auto& iit: img ) {                     // windowing: subtract and re-add mean afterwards
        iit = (*dit++ - stats.mean) * *wit++ + stats.mean;
    }
*/    
    stats.getStats(img, ST_VALUES|ST_RMS);
    str += "  mean2=" + to_string( stats.mean ) + "  std=" + to_string( stats.stddev );
    
    ft.reset( img );
//    FourierTransform::reorder( ft );
    stats.noise = channel->channel->noiseFudge * ft.noise(-1,-1);       // mask/cutoff < 0 will revert to hardcoded values used by MvN
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


void ImageData::clear(void) {
    object.reset();
    channel.reset();
    wfg.reset();
}
