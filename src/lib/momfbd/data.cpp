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


ChannelData::ChannelData( std::shared_ptr<Channel> c ) : myChannel(c) {

}


void ChannelData::getPatchData(Point16& id) {
    myChannel->getPatchData(*this,id);
}


void ChannelData::initPatch(void) {
    myChannel->initPatch(*this);
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


ObjectData::ObjectData( std::shared_ptr<Object> o ) : myObject(o) {

    for( auto& it: o->getChannels() ) {
        channels.push_back(ChannelData(it));
    }

}


void ObjectData::getPatchData(Point16& id) {
    for( ChannelData& ch: channels ) {
        ch.getPatchData(id);
    }
}


void ObjectData::initPatch(void) {
    for( ChannelData& ch: channels ) {
        ch.initPatch();
    }
    myObject->initPatch(*this);
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
    channels.clear();
    for( auto& it: myObject->getChannels() ) {
        ChannelData chan(it);
        count += chan.unpack( ptr+count, swap_endian );
        channels.push_back(chan);
    }
    return count;
}


PatchData::PatchData( const MomfbdJob& j, uint16_t yid, uint16_t xid) : myJob(j), index(yid,xid) {
    
    for( auto& it: j.getObjects() ) {
        objects.push_back(ObjectData(it));
    }

}


void PatchData::getData(void) {
    for( ObjectData& obj: objects ) {
        obj.getPatchData(index);
    }
}


void PatchData::initPatch(void) {
    //cout << "Initializing patch # " << index << endl;
    for( ObjectData& obj: objects ) {
        obj.initPatch();
    }
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
    objects.clear();
    for( auto& it: myJob.getObjects() ) {
        ObjectData obj(it);
        count += obj.unpack( ptr+count, swap_endian );
        objects.push_back(obj);
    }
    return count;
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

