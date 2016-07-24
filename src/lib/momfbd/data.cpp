#include "redux/momfbd/data.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/solver.hpp"

#include "redux/image/utils.hpp"

using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;


ChannelData::ChannelData( std::shared_ptr<Channel> c ) : myChannel(c) {

}


ChannelData::~ChannelData() {
    cacheRemove();        // remove cache file
    cclear();
}


void ChannelData::setPath(const std::string& path) {
    string mypath = path + "_" + to_string(myChannel->id());
    CacheItem::setPath(mypath);
}


void ChannelData::initPatch(void) {
    myChannel->initPatch(*this);
}


uint64_t ChannelData::size( void ) const {
    uint64_t sz = channelOffset.size() + offset.size() + residualOffset.size();
    sz += images.size();
    return sz;
}


uint64_t ChannelData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = channelOffset.pack(ptr);
    count += offset.pack(ptr+count);
    count += residualOffset.pack(ptr+count);
    count += images.pack(ptr+count);
    return count;
}


uint64_t ChannelData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = channelOffset.unpack(ptr, swap_endian);
    count += offset.unpack(ptr+count, swap_endian);
    count += residualOffset.unpack(ptr+count, swap_endian);
    count += images.unpack(ptr+count, swap_endian);
    if( images.size() > sizeof(uint64_t) ) isLoaded = true;
    else isLoaded = false;
    return count;
}


void ChannelData::cclear(void) {
    images.clear();
}


ObjectData::ObjectData( std::shared_ptr<Object> o ) : myObject(o), channels(o->getChannels().begin(),o->getChannels().end()) {

}


ObjectData::~ObjectData() {
    cacheRemove();        // remove cache file
    cclear();
}


void ObjectData::setPath(const std::string& path) {
    string mypath = path + "_" + to_string(myObject->id());
    CacheItem::setPath(mypath);
    for( auto& ch: channels ) {
        ch.setPath(mypath);
    }
}


void ObjectData::initPatch(void) {
    myObject->initPatch(*this);
    myObject->fitAvgPlane(*this);
    for( auto& ch: channels ) {
        ch.initPatch();
    }
    myObject->initPQ();
    myObject->addAllPQ();
}


uint64_t ObjectData::size( void ) const {
    uint64_t sz = img.size();
    sz += psf.size();
    sz += cobj.size();
    sz += res.size();
    sz += alpha.size();
    sz += div.size();
    for( auto& cd: channels ) {
        // use explicit scope to bypass compression (only compress the patch, not the individual channels)
        sz += cd.ChannelData::size();
    }
    return sz;
}


uint64_t ObjectData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = img.pack(ptr);
    count += psf.pack(ptr+count);
    count += cobj.pack(ptr+count);
    count += res.pack(ptr+count);
    count += alpha.pack(ptr+count);
    count += div.pack(ptr+count);
    for( auto& cd: channels ) {
        // use explicit scope to bypass compression (only compress the patch, not the individual channels)
        count += cd.ChannelData::pack( ptr+count );
    }
    return count;
}


uint64_t ObjectData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = img.unpack(ptr,swap_endian);
    count += psf.unpack(ptr+count,swap_endian);
    count += cobj.unpack(ptr+count,swap_endian);
    count += res.unpack(ptr+count,swap_endian);
    count += alpha.unpack(ptr+count,swap_endian);
    count += div.unpack(ptr+count,swap_endian);
    if( count > 6*sizeof(uint64_t) ) isLoaded = true;
    else isLoaded = false;
    for( auto& chan: channels ) {
        // use explicit scope to bypass compression (only compress the patch, not the individual channels)
        count += chan.ChannelData::unpack( ptr+count, swap_endian );
    }
    return count;
}


bool ObjectData::cacheLoad(bool removeAfterLoad) {
    bool loaded(false);
    loaded |= CacheItem::cacheLoad(removeAfterLoad);
    for( auto& cd: channels ) {
        loaded |= cd.cacheLoad(removeAfterLoad);
    }
    return loaded;      // true if any part was loaded, TODO better control of partial fails.
}


bool ObjectData::cacheStore(bool clearAfterStore) {
    bool stored(false);
    for( auto& cd: channels ) {
        stored |= cd.cacheStore(clearAfterStore);
    }
    stored |= CacheItem::cacheStore(clearAfterStore);      // true if any part was stored, TODO better control of partial fails.
    return stored;
}


void ObjectData::cclear(void) {
    img.clear();
    psf.clear();
    cobj.clear();
    res.clear();
    alpha.clear();
    div.clear();
}


PatchData::PatchData( const MomfbdJob& j, uint16_t yid, uint16_t xid) : myJob(j), objects( j.getObjects().begin(), j.getObjects().end() ), index(yid,xid) {
    
}


PatchData::~PatchData() {
}


void PatchData::setPath(const std::string& path) {
    string mypath = path+"/patch_"+(string)index;
    CacheItem::setPath(mypath);
    for( auto& obj: objects ) {
        obj.setPath(mypath);
    }
}


void PatchData::initPatch(void) {
    for( auto& obj: objects ) {
        obj.initPatch();
    }
}


uint64_t PatchData::size( void ) const {
    if( !isLoaded && cachedSize ) return cachedSize;
    uint64_t sz = Part::size();
    sz += index.size() + position.size() + roi.size();
    sz += sizeof(float);
    for( auto& obj: objects ) {
        // use explicit scope to bypass compression (only compress the patch, not the individual objects/channels)
        sz += obj.ObjectData::size();
    }
    return sz;
}


uint64_t PatchData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = Part::pack(ptr);
    count += index.pack(ptr+count);
    count += position.pack(ptr+count);
    count += roi.pack(ptr+count);
    count += pack( ptr+count, finalMetric );
    for( auto& obj: objects ) {
        // use explicit scope to bypass compression (only compress the patch, not the individual objects/channels)
        count += obj.ObjectData::pack(ptr+count);
    }
    return count;
}


uint64_t PatchData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Part::unpack(ptr, swap_endian);
    count += index.unpack(ptr+count, swap_endian);
    count += position.unpack(ptr+count, swap_endian);
    count += roi.unpack(ptr+count, swap_endian);
    count += unpack( ptr+count, finalMetric, swap_endian );
    for( auto& obj: objects ) {
        // use explicit scope to bypass compression (only compress the patch, not the individual objects/channels)
        count += obj.ObjectData::unpack( ptr+count, swap_endian );
    }
    return count;
}


void PatchData::cclear(void) {
    for( auto& obj: objects ) {
        obj.cclear();
    }
}

bool PatchData::cacheLoad(bool removeAfterLoad) {
    bool loaded(false);
    for( auto& obj: objects ) {
        loaded |= obj.cacheLoad(removeAfterLoad);
    }
    return loaded;
}

bool PatchData::cacheStore(bool clearAfterStore) {
    bool stored(false);
    for( auto& obj: objects ) {
        stored |= obj.cacheStore(clearAfterStore);
    }
    //stored |= CacheItem::cacheStore(clearAfterStore);
    return stored;
}

bool PatchData::operator==(const PatchData& rhs) {
    if(Part::operator==(rhs)) {
        return (index == rhs.index);
    }
    return false;
}


ModeSet& GlobalData::get(const ModeInfo& id, const ModeSet& ms) {

    unique_lock<mutex> lock(mtx);
    auto it = modes.find(id);
    if( it == modes.end() ){
        ModeSet& ret = redux::util::Cache::get<ModeInfo,ModeSet>(id, ms);
        return modes.emplace(id, ret).first->second;
    } else return it->second;

}



Pupil& GlobalData::get(const PupilInfo& id, const Pupil& ms) {

    unique_lock<mutex> lock(mtx);
    auto it = pupils.find(id);
    if( it == pupils.end() ){
        Pupil& ret = redux::util::Cache::get<PupilInfo,Pupil>(id, ms);
        return pupils.emplace(id, ret).first->second;
    } else return it->second;

}


uint64_t GlobalData::size( void ) const {
    uint64_t sz = 2*sizeof(uint16_t); // nModes & nPupils
    for( auto& mode: modes) {
        sz += mode.first.size();
        sz += mode.second.size();
    }
    for( auto& pupil: pupils) {
        sz += pupil.first.size();
        sz += pupil.second.size();
    }
    sz += constraints.size();
    return sz;
}


uint64_t GlobalData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack(ptr,(uint16_t)modes.size());
    for(auto& mode: modes) {
        count += mode.first.pack(ptr+count);
        count += mode.second.pack(ptr+count);
    }
    count += pack(ptr+count,(uint16_t)pupils.size());
    for(auto& pupil: pupils) {
        count += pupil.first.pack(ptr+count);
        count += pupil.second.pack(ptr+count);
    }
    count += constraints.pack(ptr+count);
    return count;
}


uint64_t GlobalData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint16_t tmp;
    uint64_t count = unpack(ptr,tmp,swap_endian);
    if(tmp) {
        ModeInfo id("");
        while( tmp-- > 0 ) {
            count += id.unpack(ptr+count,swap_endian);
            ModeSet& ms = redux::util::Cache::get<ModeInfo,ModeSet>(id);
            count += ms.unpack(ptr+count,swap_endian);
            modes.emplace(id, ms);
        }
    }
    count += unpack(ptr+count,tmp,swap_endian);
    if(tmp) {
        PupilInfo id("");
        while( tmp-- > 0 ) {
            count += id.unpack(ptr+count,swap_endian);
            Pupil& pup = redux::util::Cache::get<PupilInfo,Pupil>(id);
            count += pup.unpack(ptr+count,swap_endian);
            pupils.emplace(id, pup);
        }
    }
    count += constraints.unpack(ptr+count,swap_endian);
    return count;
}

