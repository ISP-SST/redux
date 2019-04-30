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


void ChannelData::initPatch(void) {
    myChannel->initPatch(*this);
}


void ChannelData::clear( void ) {
    
}


void ChannelData::load( void ) {
    
    myChannel->getStorage(*this);
    
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
    return count;
}


const ChannelData& ChannelData::operator=( const ChannelData& rhs ) {

    cutout = rhs.cutout;
    channelOffset = rhs.channelOffset;
    patchStart = rhs.patchStart;
    residualOffset = rhs.residualOffset;
    
    return *this;
    
}


void ChannelData::copyResults( const ChannelData& rhs ) {
    images.clear();
}


void ChannelData::dump( string tag ) const {

    tag += "_c"+to_string(myChannel->ID);
    if( images.nElements() ) Ana::write( tag + "_images.f0", images);
 
}


ObjectData::ObjectData( std::shared_ptr<Object> o ) : myObject(o) {
    if( !o ) throw logic_error("Cannot construct ObjectData from a null ObjectPtr.");
    for( auto& c: o->getChannels() ) {
        channels.push_back( make_shared<ChannelData>(c) );
    }

}


ObjectData::~ObjectData() {

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
    
    img.clear();
    psf.clear();
    cobj.clear();
    res.clear();
    alpha.clear();
    div.clear();

}


void ObjectData::load( void ) {
    
    for( auto& cd: channels ) {
        if(cd) cd->load();
    }
    
}


uint64_t ObjectData::size( void ) const {
    uint64_t sz = img.size();
    sz += psf.size();
    sz += cobj.size();
    sz += res.size();
    sz += alpha.size();
    sz += div.size();
    for( auto& cd: channels ) {
        if(cd) sz += cd->size();
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
        if(cd) count += cd->pack( ptr+count );
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
    for( auto& cd: channels ) {
        if(cd) count += cd->unpack( ptr+count, swap_endian );
    }
    return count;
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
    
    for( size_t i=0; i<channels.size(); ++i ) {
        channels[i] = rhs.channels[i];
    }
    
    return *this;
    
}


void ObjectData::copyResults( const ObjectData& rhs ) {
    
    if( channels.size() != rhs.channels.size() ) {
        throw runtime_error("Can not copy ObjectResults if they contain different number of channels.");
    }
    size_t nEl = img.nElements();
    if( nEl && (nEl == rhs.img.nElements()) ) rhs.img.copyTo<float>( img.ptr() );
    nEl = psf.nElements();
    if( nEl && (nEl == rhs.psf.nElements()) )  rhs.psf.copyTo<float>( psf.ptr() );
    nEl = cobj.nElements();
    if( nEl && (nEl == rhs.cobj.nElements()) ) rhs.cobj.copyTo<float>( cobj.ptr() );
    nEl = res.nElements();
    if( nEl && (nEl == rhs.res.nElements()) ) rhs.res.copyTo<float>( res.ptr() );
    nEl = alpha.nElements();
    if( nEl && (nEl == rhs.alpha.nElements()) ) rhs.alpha.copyTo<float>( alpha.ptr() );
    nEl = div.nElements();
    if( nEl && (nEl == rhs.div.nElements()) ) rhs.div.copyTo<float>( div.ptr() );

    for( size_t i=0; i<channels.size(); ++i ) {
        channels[i]->copyResults( *(rhs.channels[i]) );
    }
    
}


void ObjectData::dump( string tag ) const {

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


uint64_t WavefrontData::size( void ) const {
    uint64_t sz = alpha.size();
    sz += ids.size()*sizeof(uint32_t) + sizeof(uint64_t);
    return sz;
}


uint64_t WavefrontData::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = alpha.pack(ptr);
    count += pack(ptr+count, ids);
    return count;
}


uint64_t WavefrontData::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = alpha.unpack(ptr,swap_endian);
    count += unpack(ptr+count, ids, swap_endian);
    return count;
}


const WavefrontData& WavefrontData::operator=( const WavefrontData& rhs ) {

    ids = rhs.ids;
    alpha = rhs.alpha;
    return *this;
    
}


void WavefrontData::copyResults( const WavefrontData& rhs ) {
    
    alpha.copyFrom<float>( rhs.alpha.ptr() );
    ids = rhs.ids;
    
}


void WavefrontData::dump( string tag ) const {

    if( !ids.empty() ) {
        vector<int32_t> tmp( ids.begin(), ids.end() );
        Ana::write( tag + "_wfi.f0", tmp );
    }
    if( alpha.nElements() ) Ana::write( tag + "_wf.f0", alpha );
 
}


PatchData::PatchData( MomfbdJob& j, uint16_t yid, uint16_t xid) : myJob(j), index(yid,xid), finalMetric(0.0) {
    vector<shared_ptr<Object>> objs = myJob.getObjects();
    for( auto& o: objs ) {
        if(o) objects.push_back( make_shared<ObjectData>(o) );
    }
}


PatchData::~PatchData() {
    
}


void PatchData::setPath(const std::string& path) {
    string mypath = path+"/patch_"+(string)index;
    CacheItem::setPath(mypath);
}


shared_ptr<ObjectData> PatchData::getObjectData( uint16_t id ) const {
    for( auto& o: objects ) {
        if( o && (o->myObject->ID == id) ) return o;
    }
    for( auto& o: trace_data ) {
        if( o && (o->myObject->ID == id) ) return o;
    }
    return shared_ptr<ObjectData>();
}


void PatchData::initPatch(void) {
    for( auto& obj: objects ) {
        if(obj) obj->initPatch();
    }
}


void PatchData::load( void ) {
    
    Part::load();
    for( auto& obj: objects ) {
        if(obj) obj->load();
    }
    
}


void PatchData::unload( void ) {
    
    Part::unload();
    packed.clear();
    
}

void PatchData::prePack( bool force ) {
    
    if( packed.packedSize && !force ) {
        return;
    }
    packed.size = size();
    if( !packed.size ) {
        return;
    }
    packed.data.reset( new char[packed.size] );
    packed.packedSize = pack( packed.data.get() );

}


void PatchData::clear( void ) {
    
    unload();
    for( auto& obj: objects ) {
        if(obj) obj->clear();
    }
    
}


uint64_t PatchData::size( void ) const {

    if( !isLoaded && cachedSize ) return cachedSize;
    uint64_t sz = Part::size();
    sz += index.size() + position.size() + roi.size();
    sz += sizeof(float);
    sz += metrics.size()*sizeof(float) + sizeof(uint64_t) + sizeof(uint16_t);
    for( auto& obj: objects ) {
        if(obj) sz += obj->size();
    }
    for( auto& td: trace_data ) {
        if(td) sz += td->size();
    }
    sz += waveFronts.size();
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
        if(obj) count += obj->pack(ptr+count);
    }
    count += waveFronts.WavefrontData::pack(ptr+count);
    uint16_t tmp = trace_data.size();
    count += pack( ptr+count, tmp );
    if( tmp ) {
        for( const auto& tobj: trace_data ) {
            if( tobj ) {
                count += tobj->pack( ptr+count );
            }
        }
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
        if(obj) count += obj->unpack( ptr+count, swap_endian );
    }
    count += waveFronts.WavefrontData::unpack(ptr+count, swap_endian);
    uint16_t tmp;
    count += unpack( ptr+count, tmp, swap_endian );
    trace_data.resize(tmp);
    for( auto& tobj: trace_data ) {
        tobj.reset( new ObjectData() );
        count += tobj->unpack( ptr+count, swap_endian );
    }

    return count;
}


void PatchData::cclear(void) {

}


bool PatchData::cacheLoad(bool removeAfterLoad) {
    return false;
}

bool PatchData::cacheStore(bool clearAfterStore) {
    return false;
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
    waveFronts = rhs.waveFronts;
    trace_data.resize( rhs.trace_data.size() );
    if( trace_data.size() != myJob.trace_objects.size() ) {
        throw runtime_error("Can not copy trace-data if they contain different number of objects.");
    }
    for( size_t i=0; i<trace_data.size(); ++i ) {
        trace_data[i] = rhs.trace_data[i];
        trace_data[i]->myObject = myJob.trace_objects[i];
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

    myJob.waveFronts.getStorage( *this );
    waveFronts.copyResults( rhs.waveFronts );
    for( size_t i=0; i<objects.size(); ++i ) {
        objects[i]->myObject->getStorage( *this, objects[i] );
        objects[i]->copyResults( *(rhs.objects[i]) );
    }
    trace_data.resize( rhs.trace_data.size() );
    if( trace_data.size() != myJob.trace_objects.size() ) {
        throw runtime_error("Can not copy trace-data if they contain different number of objects.");
    }
    for( size_t i=0; i<trace_data.size(); ++i ) {
        if( !rhs.trace_data[i] )  continue;
        trace_data[i].reset( new ObjectData() );
        trace_data[i]->myObject = myJob.trace_objects[i];
        trace_data[i]->myObject->getStorage( *this, trace_data[i] );
        trace_data[i]->copyResults( *(rhs.trace_data[i]) );
    }
    
}


void PatchData::dump( string tag ) const {

    tag += "_patch" + (string)index;
    
    for( auto & object : objects ) {
        object->dump(tag);
    }
    waveFronts.dump(tag);
     
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
    
    unique_lock<mutex> lock(mtx);
    using redux::util::pack;
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
        modes.clear();
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
        PupilInfo id( "", 0.0 );
        pupils.clear();
        while( tmp-- > 0 ) {
            count += id.unpack(ptr+count,swap_endian);
            shared_ptr<Pupil>& pup = redux::util::Cache::get<PupilInfo,shared_ptr<Pupil>>(id);
            if( !pup ) pup.reset( new Pupil() );
            count += pup->unpack(ptr+count,swap_endian);
            pupils.emplace(id, pup);
        }
    }
    count += constraints.unpack(ptr+count,swap_endian);

    return count;
}


void GlobalData::unload( void ) {
    
    Part::unload();
    
}

void GlobalData::prePack( bool force ) {
    
    if( packed.packedSize && !force ) {
        return;
    }
    packed.size = size();
    if( !packed.size ) {
        return;
    }
    packed.data.reset( new char[packed.size] );
    packed.packedSize = pack( packed.data.get() );
    
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

