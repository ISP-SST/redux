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


uint64_t ChannelData::size( void ) const {
    static uint64_t sz = offset.size() + residualOffset.size();
    return sz+images.size();
}


uint64_t ChannelData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = offset.pack(ptr);
    count += residualOffset.pack(ptr+count);
    count += images.pack(ptr+count);
    return count;
}


uint64_t ChannelData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = offset.unpack(ptr, swap_endian);
    count += residualOffset.unpack(ptr+count, swap_endian);
    count += images.unpack(ptr+count, swap_endian);
    return count;
}


uint64_t ObjectData::size( void ) const {
    uint64_t sz = 0;     // nChannels
    for( const ChannelData& cd: channels ) {
        sz += cd.size();
    }
    return sz;
}


uint64_t ObjectData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = 0;
    for( const ChannelData& cd: channels ) {
        count += cd.pack( ptr+count );
    }
    return count;
}


uint64_t ObjectData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = 0;
    for( ChannelData& cd: channels ) {
        count += cd.unpack( ptr+count, swap_endian );
    }
    return count;
}


PatchData::PatchData( const MomfbdJob& j, uint16_t yid, uint16_t xid) : myJob(j), index(yid,xid) {
    
    objects.resize(j.getObjects().size());
    
}


uint64_t PatchData::size( void ) const {
    uint64_t sz = Part::size();
    sz += index.size() + pos.size();
    for( const ObjectData& obj: objects ) {
        sz += obj.size();
    }
    return sz;
}


uint64_t PatchData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = Part::pack(ptr);
    count += index.pack(ptr+count);
    count += pos.pack(ptr+count);
    for( const ObjectData& obj: objects ) {
        count += obj.pack(ptr+count);
    }
    return count;
}


uint64_t PatchData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Part::unpack(ptr, swap_endian);
    count += index.unpack(ptr+count, swap_endian);
    count += pos.unpack(ptr+count, swap_endian);
    const vector<shared_ptr<Object>>& job_obj = myJob.getObjects();
    cout << "PatchData::unpack(): nObj = " << job_obj.size() << endl;
    objects.resize(job_obj.size());
    auto it = job_obj.begin();
    for( ObjectData& objData: objects ) {
        //ObjectData::Ptr objData(new ObjectData(obj));
        //objects.push_back(objData);
    cout << "PatchData::unpack(): nCh = " << (*it)->channels.size() << endl;
        objData.channels.resize((*it)->channels.size());
        count += objData.unpack( ptr+count, swap_endian );
        it++;
    }
    return count;
}


bool PatchData::operator==(const PatchData& rhs) {
    if(Part::operator==(rhs)) {
        return (index == rhs.index);
    }
    return false;
}


void GlobalData::fetch(const Cache::ModeID& id) {
    Cache& cache = Cache::getCache();
    modes.emplace(id, cache.mode(id));
}


uint64_t GlobalData::size( void ) const {
    uint64_t sz = sizeof(uint16_t); // nModes
    for( auto& it: modes) {
        sz += it.first.size();
        sz += it.second->size();
    }
    return sz;
}


uint64_t GlobalData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack(ptr,(uint16_t)modes.size());
    for(auto& it: modes) {
        count += it.first.pack(ptr+count);
        count += it.second->pack(ptr+count);
    }
    cout << "GlobalData::pack():  packed " << modes.size() << " modes." << endl;
    return count;
}


uint64_t GlobalData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint16_t tmp;
    uint64_t count = unpack(ptr,tmp,swap_endian);
    if(tmp) {
        Cache::ModeID id;
        PupilMode::Ptr mode(new PupilMode);
        while( tmp-- > 0 ) {
            count += id.unpack(ptr+count,swap_endian);
            count += mode->unpack(ptr+count,swap_endian);
            modes.emplace(id, mode);
        }
    }
    cout << "GlobalData::unpack():  got " << modes.size() << " modes." << endl;
    return count;
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
    stats.noise = channel->cfg->noiseFudge * ft.noise(-1,-1);       // mask/cutoff < 0 will revert to hardcoded values used by MvN
    str += "  noise=" + to_string( stats.noise );
    //object->addToFT( ft );
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
