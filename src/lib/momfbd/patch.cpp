#include "redux/momfbd/patch.hpp"

#include "redux/momfbd/momfbdjob.hpp"

#include "redux/util/datautil.hpp"


using namespace redux::momfbd;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

    const string thisChannel = "momfbdpatch";

    typedef std::shared_ptr<ChannelData> ChPtr;
    typedef std::shared_ptr<ObjectData> ObjPtr;
    
}


ChannelData::ChannelData( const ObjPtr& o, const std::shared_ptr<Channel>& c ) : channel( c ), object( o ) {
    //  std::cout << "Channel():  1   " << hexString(this) << std::endl;

}


ChannelData::~ChannelData(void) {
    // std::cout << "~Channel():  " << hexString(this) << std::endl;

}

uint64_t ChannelData::size( void ) const {
    uint64_t sz = images.size();
    return sz;
}


uint64_t ChannelData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = 0;
    count += images.pack(ptr+count);
    return count;
}


uint64_t ChannelData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = 0;
    count += images.unpack(ptr+count, swap_endian);
    return count;
}


void ChannelData::init(uint16_t yid, uint16_t xid) {
    channel->getPatchData(shared_from_this(),yid,xid);
}


void ChannelData::clear(void) {
    images.resize();
    channel.reset();
}


ObjectData::ObjectData( const std::shared_ptr<Object>& c ) : object( c ) {
    //  std::cout << "Object():  " << hexString(this) << std::endl;

}


ObjectData::~ObjectData() {
   // std::cout << "~Object():  1   "  << hexString(this) << printArray( ftSum.dimensions(),"   ftSum.dims" ) << std::endl;
    static int count( 0 );
    if( ftSum.nDimensions() > 1 ) {
        //redux::file::Ana::write( "ftSum_" + to_string( count++ ) + ".f0", ftSum );
    }
}


uint64_t ObjectData::size( void ) const {
    uint64_t sz = 0;     // nChannels
    for( const ChannelData::Ptr& ch: channels ) {
        sz += ch->size();
    }
    return sz;
}


uint64_t ObjectData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = 0;
    for( const ChannelData::Ptr& cd: channels ) {
        count += cd->pack( ptr+count );
    }
    return count;
}


uint64_t ObjectData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = 0;
    channels.clear();
    for( shared_ptr<Channel>& ch: object->channels ) {
        ChannelData::Ptr cd(new ChannelData(shared_from_this(),ch));
        count += cd->unpack( ptr+count, swap_endian );
        channels.push_back(cd);
    }
    return count;
}


void ObjectData::init(uint16_t yid, uint16_t xid) {
    for( const shared_ptr<Channel>& ch : object->channels ) {
        ChannelData::Ptr chData(new ChannelData(shared_from_this(),ch));
        channels.push_back(chData);
        chData->init(yid,xid);
    }
}


void ObjectData::clear(void) {
    for(auto &it: channels) it->clear();
    channels.clear();
    object.reset();
}


void ObjectData::addToFT( const redux::image::FourierTransform& ft ) {
    unique_lock<mutex> lock( mtx );
    if( !ftSum.sameSizes(ft) ) {
        ftSum.resize(ft.dimensions(true));
        ftSum.zero();
    }
    //  std::cout << "Object::addToFT():  1   " << hexString(this) << printArray(ft.dimensions(), "  ft.dims") << printArray(ftSum.dimensions(), "  ftSum.dims") << std::endl;
    ftSum += ft;
//    std::cout << "Object::addToFT():  2   " << printArray(ft.dimensions(), "ft.dims") << printArray(ftSum.dimensions(), "  ftSum.dims") << std::endl;
}


void ObjectData::addToPQ( const redux::image::FourierTransform& ft, const Array<complex_t> sj ) {
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


/*
Patch::Patch(int y, int x, uint32_t sz) : index(), first(), last() {
    if( sz > 1 ) {
        int halfSize = sz / 2;
        y = std::max(y-halfSize,0);
        x = std::max(x-halfSize,0);
    }
    first = Point(y,x);
    last = Point(y+sz-1,x+sz-1);
}
        ObjectData::Ptr objData(new ObjectData(obj));
        patch->data.push_back(objData);
        objData->init(patch->index.y,patch->index.x);


*/

void PatchData::setIndex(uint16_t yid, uint16_t xid) {
    index.x = xid;
    index.y = yid;
}


uint64_t PatchData::size( void ) const {
    uint64_t sz = Part::size();
    sz += index.size() + pos.size();
    for( const ObjectData::Ptr& obj: objects ) {
        sz += obj->size();
    }
    return sz;
}


uint64_t PatchData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = Part::pack(ptr);
    count += index.pack(ptr+count);
    count += pos.pack(ptr+count);
    for( const ObjectData::Ptr& obj: objects ) {
        count += obj->pack(ptr+count);
    }
    return count;
}


uint64_t PatchData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Part::unpack(ptr, swap_endian);
    count += index.unpack(ptr+count, swap_endian);
    count += pos.unpack(ptr+count, swap_endian);
    objects.clear();
    for( const std::shared_ptr<Object>& obj: myJob.getObjects() ) {
        ObjectData::Ptr objData(new ObjectData(obj));
        objects.push_back(objData);
        count += objData->unpack( ptr+count, swap_endian );
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
