#include "redux/momfbd/data.hpp"

#include "redux/logging/logger.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/solver.hpp"
#include "redux/util/cache.hpp"

#include "redux/file/fileana.hpp"
#include "redux/image/utils.hpp"

using namespace redux::file;
using namespace redux::logging;
using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;


ChannelData::ChannelData( std::shared_ptr<Channel> c ) : myChannel(c) {
    if( !c ) throw logic_error("Cannot construct ChannelData from a null ChannelPtr.");
}


ChannelData::~ChannelData() {
    
}


void ChannelData::setPath(const std::string& path) {
    string mypath = path + "_" + to_string(myChannel->id());
    CacheItem::setPath(mypath);
}


void ChannelData::initPatch(void) {
    myChannel->initPatch(*this);
}


void ChannelData::clear( void ) {
    
    try {
        cacheRemove();        // remove cache file
    } catch( exception& e ) {
        LLOG(myChannel->logger) << "Exception while removing the cache-file: \"" << getFullPath() << "\" what:" << e.what() << ende;
    }
    cclear();
    
}


uint64_t ChannelData::size( void ) const {
    static uint64_t fixed_sz = channelOffset.size() + patchStart.size() + residualOffset.size() + 1;
    return fixed_sz + images.size();
}


uint64_t ChannelData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = channelOffset.pack(ptr);
    count += patchStart.pack(ptr+count);
    count += residualOffset.pack(ptr+count);
    uint8_t hasImages(0);
    if( images.size() > sizeof(uint64_t) ) hasImages = 1;
    count += pack(ptr+count,hasImages);
    if( hasImages ) {
        count += images.pack(ptr+count);
    }
    return count;
}


uint64_t ChannelData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = channelOffset.unpack(ptr, swap_endian);
    count += patchStart.unpack(ptr+count, swap_endian);
    count += residualOffset.unpack(ptr+count, swap_endian);
    uint8_t hasImages(0);
    count += unpack(ptr+count, hasImages, swap_endian);
    if( hasImages ) {
        count += images.unpack(ptr+count, swap_endian);
    }
    if( images.size() > sizeof(uint64_t) ) isLoaded = true;
    else isLoaded = false;
    return count;
}


void ChannelData::cclear(void) {
    images.clear();
    isLoaded = false;
}


const ChannelData& ChannelData::operator=( const ChannelData& rhs ) {

    //images = rhs.images;
    cutout = rhs.cutout;
    channelOffset = rhs.channelOffset;
    patchStart = rhs.patchStart;
    residualOffset = rhs.residualOffset;

    isLoaded = images.nElements();
    
    return *this;
    
}


void ChannelData::copyResults( const ChannelData& rhs ) {
    
    isLoaded = false;

}


void ChannelData::dump( string tag ) {

    tag += "_c"+to_string(myChannel->ID);
    if( images.nElements() ) Ana::write( tag + "_images.f0", images);
 
}


ObjectData::ObjectData( std::shared_ptr<Object> o ) : myObject(o) {
    if( !o ) throw logic_error("Cannot construct ObjectData from a null ObjectPtr.");
    for( auto& c: o->getChannels() ) {
        channels.push_back( make_shared<Compressed<ChannelData,5>>(c) );
    }

}


ObjectData::~ObjectData() {

}


void ObjectData::setPath(const std::string& path) {
    string mypath = path + "_" + to_string(myObject->id());
    CacheItem::setPath(mypath);
    for( auto& cd: channels ) {
        if(cd) cd->setPath(mypath);
    }
}


void ObjectData::initPatch(void) {
    myObject->initPatch();
    for( auto& cd: channels ) {
        if(cd) cd->initPatch();
    }
    myObject->fitAvgPlane();
    myObject->initPQ();
    myObject->addAllPQ();
}


void ObjectData::clear( void ) {
    
    for( auto& cd: channels ) {
        if(cd) cd->clear();
    }
    try {
        cacheRemove();        // remove cache file
    } catch( exception& e ) {
        LLOG(myObject->logger) << "Exception while removing the cache-file: \"" << getFullPath() << "\" what:" << e.what() << ende;
    }
    cclear();
    
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
        if(cd) sz += cd->ChannelData::size();
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
        if(cd) count += cd->ChannelData::pack( ptr+count );
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
    for( auto& cd: channels ) {
        // use explicit scope to bypass compression (only compress the patch, not the individual channels)
        if(cd) count += cd->ChannelData::unpack( ptr+count, swap_endian );
    }
    return count;
}


bool ObjectData::cacheLoad(bool removeAfterLoad) {
    bool loaded(false);
    loaded |= CacheItem::cacheLoad(removeAfterLoad);
    for( auto& cd: channels ) {
        if(cd) loaded |= cd->cacheLoad(removeAfterLoad);
    }
    return loaded;      // true if any part was loaded, TODO better control of partial fails.
}


bool ObjectData::cacheStore(bool clearAfterStore) {
    bool stored(false);
    for( auto& cd: channels ) {
        if(cd) stored |= cd->cacheStore(clearAfterStore);
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
    for( auto& cd: channels ) {
        cd->cclear();
    }
    isLoaded = false;
}


const ObjectData& ObjectData::operator=( const ObjectData& rhs ) {
    
    if( channels.size() != rhs.channels.size() ) {
        throw runtime_error("Can not copy ObjectData if they contain different number of channels.");
    }

    img = rhs.img;
    psf = rhs.psf;
    cobj = rhs.cobj;
    res = rhs.res;
    alpha = rhs.alpha;
    div = rhs.div;
    
    isLoaded = (img.nElements() || psf.nElements() || cobj.nElements()
             || res.nElements() || alpha.nElements() || div.nElements() );

    for( size_t i=0; i<channels.size(); ++i ) {
        channels[i] = rhs.channels[i];
    }
    
    return *this;
    
}


void ObjectData::copyResults( const ObjectData& rhs ) {
    
    if( channels.size() != rhs.channels.size() ) {
        throw runtime_error("Can not copy ObjectResults if they contain different number of channels.");
    }

    img = rhs.img;
    psf = rhs.psf;
    cobj = rhs.cobj;
    res = rhs.res;
    alpha = rhs.alpha;
    div = rhs.div;
    
    isLoaded = (img.nElements() || psf.nElements() || cobj.nElements()
             || res.nElements() || alpha.nElements() || div.nElements() );

    for( size_t i=0; i<channels.size(); ++i ) {
        channels[i]->copyResults( *(rhs.channels[i]) );
    }
    
}


void ObjectData::dump( string tag ) {

    tag += "_o"+to_string(myObject->ID);
    for( auto& ch : channels ) {
        ch->dump(tag);
    }

    if( img.nElements() ) Ana::write( tag + "_img.f0", img);
    if( psf.nElements() ) Ana::write( tag + "_psf.f0", psf);
    if( cobj.nElements() ) Ana::write( tag + "_cobj.f0", cobj);
    if( res.nElements() ) Ana::write( tag + "_res.f0", res);
    if( alpha.nElements() ) Ana::write( tag + "_alpha.f0", alpha);
    if( div.nElements() ) Ana::write( tag + "_div.f0", div);
 
}


PatchData::PatchData( const MomfbdJob& j, uint16_t yid, uint16_t xid) : myJob(j), index(yid,xid), finalMetric(0.0) {
    vector<shared_ptr<Object>> objs = myJob.getObjects();
    for( auto& o: objs ) {
        if(o) objects.push_back( make_shared<Compressed<ObjectData,5>>(o) );
    }
}


PatchData::~PatchData() {
    
}


void PatchData::setPath(const std::string& path) {
    string mypath = path+"/patch_"+(string)index;
    CacheItem::setPath(mypath);
    for( auto& obj: objects ) {
        if(obj) obj->setPath(mypath);
    }
}


void PatchData::initPatch(void) {
    for( auto& obj: objects ) {
        if(obj) obj->initPatch();
    }
}


void PatchData::clear( void ) {
    
    for( auto& obj: objects ) {
        if(obj) obj->clear();
    }
    
}


uint64_t PatchData::size( void ) const {

    if( !isLoaded && cachedSize ) return cachedSize;
    uint64_t sz = Part::size();
    sz += index.size() + position.size() + roi.size();
    sz += sizeof(float);
    sz += metrics.size()*sizeof(float) + sizeof(uint64_t);

    for( auto& obj: objects ) {
        // use explicit scope to bypass compression (only compress the patch, not the individual objects/channels)
        if(obj) sz += obj->ObjectData::size();
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
    count += pack( ptr+count, metrics );
    for( auto& obj: objects ) {
        // use explicit scope to bypass compression (only compress the patch, not the individual objects/channels)
        if(obj) count += obj->ObjectData::pack(ptr+count);
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
    count += unpack( ptr+count, metrics, swap_endian );
    for( auto& obj: objects ) {
        // use explicit scope to bypass compression (only compress the patch, not the individual objects/channels)
        if(obj) count += obj->ObjectData::unpack( ptr+count, swap_endian );
    }
    return count;
}


void PatchData::cclear(void) {
    for( auto& obj: objects ) {
        if(obj) obj->cclear();
    }
}

bool PatchData::cacheLoad(bool removeAfterLoad) {
    bool loaded(false);
    for( auto& obj: objects ) {
        if(obj) loaded |= obj->cacheLoad(removeAfterLoad);
    }
    return loaded;
}

bool PatchData::cacheStore(bool clearAfterStore) {
    bool stored(false);
    for( auto& obj: objects ) {
        if(obj) stored |= obj->cacheStore(clearAfterStore);
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


const PatchData& PatchData::operator=( const PatchData& rhs ) {
    
    if( objects.size() != rhs.objects.size() ) {
        throw runtime_error("Can not copy PatchData if they contain different number of objects.");
    }
    
    index = rhs.index;
    roi = rhs.roi;
    finalMetric = rhs.finalMetric;
    metrics = rhs.metrics;

    for( size_t i=0; i<objects.size(); ++i ) {
        objects[i] = rhs.objects[i];
    }
    
    return *this;
    
}


void PatchData::copyResults( const PatchData& rhs ) {
    
    if( objects.size() != rhs.objects.size() ) {
        throw runtime_error("Can not copy PatchData if they contain different number of objects.");
    }

    finalMetric = rhs.finalMetric;
    metrics = rhs.metrics;
    nThreads = rhs.nThreads;
    runtime_wall = rhs.runtime_wall;
    runtime_cpu = rhs.runtime_cpu;

    for( size_t i=0; i<objects.size(); ++i ) {
        objects[i]->copyResults( *(rhs.objects[i]) );
    }
    
}


void PatchData::dump( string tag ) {

    tag += "_patch" + (string)index;
    
    for( auto & object : objects ) {
        object->dump(tag);
    }
     
}


shared_ptr<ModeSet> GlobalData::get(const ModeInfo& id, const shared_ptr<ModeSet>& ms) {

    unique_lock<mutex> lock(mtx);
    auto it = modes.find(id);
    if( it == modes.end() ){
        shared_ptr<ModeSet>& ret = redux::util::Cache::get<ModeInfo,shared_ptr<ModeSet>>(id, ms);
        if( !ret ) ret.reset( new ModeSet() );
        it = modes.emplace(id, ret).first;
    }
    return it->second;

}


shared_ptr<Pupil> GlobalData::get(const PupilInfo& id, const shared_ptr<Pupil>& ms) {

    unique_lock<mutex> lock(mtx);
    auto it = pupils.find(id);
    if( it == pupils.end() ){
        shared_ptr<Pupil>& ret = redux::util::Cache::get<PupilInfo,shared_ptr<Pupil>>(id, ms);
        if( !ret ) ret.reset( new Pupil() );
        it = pupils.emplace(id, ret).first;
    }
    return it->second;

}


bool GlobalData::verify( void ) const {

    if( !constraints.verify() ) return false;
    
    return true;
    
}


uint64_t GlobalData::size( void ) const {
    uint64_t sz = 2*sizeof(uint16_t); // nModes & nPupils
    for( auto& mode: modes ) {
        sz += mode.first.size();
        sz += mode.second->size();
    }
    for( auto& pupil: pupils ) {
        sz += pupil.first.size();
        sz += pupil.second->size();
    }
    sz += constraints.size();
    return sz;
}


uint64_t GlobalData::pack( char* ptr ) const {
    using redux::util::pack;
    unique_lock<mutex> lock(const_cast<mutex&>(mtx));
    uint64_t count = pack(ptr,(uint16_t)modes.size());
    for(auto& mode: modes) {
        count += mode.first.pack(ptr+count);
        count += mode.second->pack(ptr+count);
    }
    count += pack(ptr+count,(uint16_t)pupils.size());
    for(auto& pupil: pupils) {
        count += pupil.first.pack(ptr+count);
        count += pupil.second->pack(ptr+count);
    }
    count += constraints.pack(ptr+count);
    return count;
}


uint64_t GlobalData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    unique_lock<mutex> lock(mtx);
    uint16_t tmp;
    uint64_t count = unpack(ptr,tmp,swap_endian);
    if(tmp) {
        ModeInfo id("");
        while( tmp-- > 0 ) {
            count += id.unpack(ptr+count,swap_endian);
            shared_ptr<ModeSet>& ms = redux::util::Cache::get<ModeInfo,shared_ptr<ModeSet>>(id);
            if( !ms ) ms.reset( new ModeSet() );
            count += ms->unpack(ptr+count,swap_endian);
            modes.emplace(id, ms);
        }
    }
    count += unpack(ptr+count,tmp,swap_endian);
    if(tmp) {
        PupilInfo id("");
        while( tmp-- > 0 ) {
            count += id.unpack(ptr+count,swap_endian);
            shared_ptr<Pupil>& pup = redux::util::Cache::get<PupilInfo,shared_ptr<Pupil>>(id);
            if( !pup ) pup.reset( new Pupil() );
            count += pup->unpack(ptr+count,swap_endian);
            pupils.emplace(id, pup);
        }
    }
    count += constraints.unpack(ptr+count,swap_endian);
    constraints.makeRowsCols();
    return count;
}


void GlobalData::dump( string tag ) const {

    for( auto& m: modes ) {
        if( m.second->nElements() ) {
            Ana::write( tag + "_modes_" + (string)m.first, *(m.second) );
        }
    }

    for( auto& p: pupils ) {
        p.second->dump( tag + "_pupil_" + (string)p.first );
    }
    constraints.dump( tag + "_constraints" );
    
}

