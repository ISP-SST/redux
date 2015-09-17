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


ChannelData::ChannelData( std::shared_ptr<Channel> c, const std::string& cachePath ) : myChannel(c) {
    string path = cachePath + "_" + to_string(c->id());
    setPath(path);
}


void ChannelData::getPatchData(const PatchData& patch) {
    myChannel->getPatchData(*this,patch);
    isLoaded = true;
    cacheStore(true);
}


void ChannelData::initPatch(void) {
    myChannel->initPatch(*this);
}


void ChannelData::collectResults(void) {
    myChannel->getResults(*this);
    images.clear();         // don't need input data anymore.
    isLoaded = false;
}


uint64_t ChannelData::size( void ) const {
    uint64_t sz = shift.size() + offset.size() + residualOffset.size();
    sz += images.size();
    return sz;
}


uint64_t ChannelData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = shift.pack(ptr);
    count += offset.pack(ptr+count);
    count += residualOffset.pack(ptr+count);
    count += images.pack(ptr+count);
    return count;
}


uint64_t ChannelData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = shift.unpack(ptr, swap_endian);
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


ObjectData::ObjectData( std::shared_ptr<Object> o, const std::string& cachePath ) : myObject(o) {

    string path = cachePath + "_" + to_string(o->id());
    setPath(path);
    for( auto& it: o->getChannels() ) {
        channels.push_back(ChannelData(it,path));
    }

}


void ObjectData::getPatchData(const PatchData& patch) {
    for( ChannelData& ch: channels ) {
        ch.getPatchData(patch);
    }
}


void ObjectData::initPatch(void) {
    myObject->initPatch(*this);
    for( ChannelData& ch: channels ) {
        ch.initPatch();
    }
    myObject->initPQ();
    myObject->addAllPQ();
    myObject->fitAvgPlane(*this);
}


void ObjectData::collectResults(void) {
    myObject->getResults(*this);
    if(size()>6*sizeof(uint64_t)) isLoaded = true;
    else isLoaded = false;
    for( ChannelData& ch: channels ) {
        ch.collectResults();
    }
}


uint64_t ObjectData::size( void ) const {
    uint64_t sz = img.size();
    sz += psf.size();
    sz += cobj.size();
    sz += res.size();
    sz += alpha.size();
    sz += div.size();
    for( const ChannelData& cd: channels ) {
        sz += cd.size();
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
    for( const ChannelData& cd: channels ) {
        count += cd.pack( ptr+count );
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
    channels.clear();
    for( auto& it: myObject->getChannels() ) {
        ChannelData chan(it,path());
        count += chan.unpack( ptr+count, swap_endian );
        channels.push_back(chan);
    }
    return count;
}


bool ObjectData::cacheLoad(bool removeAfterLoad) {
    bool loaded(false);
    for( ChannelData& cd: channels ) {
        loaded |= cd.cacheLoad(removeAfterLoad);
    }
    loaded |= CacheItem::cacheLoad(removeAfterLoad);
    return loaded;      // true if any part was loaded, TODO better control of partial fails.
}


bool ObjectData::cacheStore(bool clearAfterStore) {
    bool stored(false);
    for( ChannelData& cd: channels ) {
        stored |= cd.cacheStore(clearAfterStore);
    }
    stored |= CacheItem::cacheStore(clearAfterStore);      // true if any part was stored, TODO better control of partial fails.
}


void ObjectData::cclear(void) {
    for( ChannelData& cd: channels ) {
        cd.cclear();
    }
    img.clear();
    psf.clear();
    cobj.clear();
    res.clear();
    alpha.clear();
    div.clear();
}


PatchData::PatchData( const MomfbdJob& j, uint16_t yid, uint16_t xid) : myJob(j), index(yid,xid) {
    
    string path = to_string(j.info.id)+"/patch_"+(string)index;
    setPath(path); 
    for( auto& it: j.getObjects() ) {
        objects.push_back(ObjectData(it,path));
    }

    
}


void PatchData::getData(void) {
    for( ObjectData& obj: objects ) {
        obj.getPatchData(*this);
    }
    isLoaded = true;
    cacheStore(true);
}


void PatchData::initPatch(void) {
    //cout << "Initializing patch # " << index << endl;
    for( ObjectData& obj: objects ) {
        obj.initPatch();
    }
}


void PatchData::collectResults(void) {
    for( ObjectData& obj: objects ) {
        obj.collectResults();
    }
}


uint64_t PatchData::size( void ) const {
    uint64_t sz = Part::size();
    sz += index.size() + roi.size();
    for( const ObjectData& obj: objects ) {
        sz += obj.size();
    }
    return sz;
}


uint64_t PatchData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = Part::pack(ptr);
    count += index.pack(ptr+count);
    count += roi.pack(ptr+count);
    for( const ObjectData& obj: objects ) {
        count += obj.pack(ptr+count);
    }
    return count;
}


uint64_t PatchData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Part::unpack(ptr, swap_endian);
    count += index.unpack(ptr+count, swap_endian);
    count += roi.unpack(ptr+count, swap_endian);
    objects.clear();
    for( auto& it: myJob.getObjects() ) {
        ObjectData obj(it,path());
        count += obj.unpack( ptr+count, swap_endian );
        objects.push_back(obj);
    }
    return count;
}


void PatchData::cclear(void) {
    for( ObjectData& obj: objects ) {
        obj.cclear();
    }
}

bool PatchData::cacheLoad(bool removeAfterLoad) {
    bool loaded(false);
    for( ObjectData& obj: objects ) {
        loaded |= obj.cacheLoad(removeAfterLoad);
    }
    //loaded |= CacheItem::cacheLoad(removeAfterLoad);
    return loaded;
}

bool PatchData::cacheStore(bool clearAfterStore) {
    bool stored(false);
    for( ObjectData& obj: objects ) {
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


const pair<Array<double>, double>& GlobalData::fetch(uint16_t pupilPixels, double pupilRadiusInPixels ) {
    //cout << "GlobalData::fetch():   pupilPixels=" << pupilPixels << "  pupilRadiusInPixels=" << pupilRadiusInPixels << endl;
    unique_lock<mutex> lock(mtx);
    //cout << "GlobalData::fetch2():   pupilPixels=" << pupilPixels << "  pupilRadiusInPixels=" << pupilRadiusInPixels << "  pupils.sz=" << pupils.size() << endl;
    return pupils.emplace( make_pair(pupilPixels,pupilRadiusInPixels),
                           Cache::getCache().pupil(pupilPixels, pupilRadiusInPixels) ).first->second;

}


const PupilMode::Ptr GlobalData::fetch(const Cache::ModeID& id) {
//cout << " GlobalData::fetch1(mode)  id=" << id.modeNumber << endl;
    unique_lock<mutex> lock(mtx);
    Cache& cache = Cache::getCache();
    return modes.emplace( id, cache.mode(id) ).first->second;
}


uint64_t GlobalData::size( void ) const {
    uint64_t sz = sizeof(uint16_t); // nModes
    for( auto& it: modes) {
        sz += it.first.size();
        sz += it.second->size();
    }
    sz += constraints.size();
    return sz;
}


uint64_t GlobalData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack(ptr,(uint16_t)modes.size());
    for(auto& it: modes) {
        count += it.first.pack(ptr+count);
        count += it.second->pack(ptr+count);
    }
  //  cout << "GlobalData::pack():  packed " << modes.size() << " modes." << endl;
    count += constraints.pack(ptr+count);
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
 //   cout << "GlobalData::unpack():  got " << modes.size() << " modes." << endl;
    count += constraints.unpack(ptr+count,swap_endian);
    return count;
}

