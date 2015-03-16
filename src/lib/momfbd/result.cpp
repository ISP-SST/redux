#include "redux/momfbd/result.hpp"

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



void ObjectResult::addToFT( const redux::image::FourierTransform& ft ) {
    unique_lock<mutex> lock( mtx );
    if( !ftSum.sameSizes(ft) ) {
        ftSum.resize(ft.dimensions(true));
        ftSum.zero();
    }
    //  std::cout << "Object::addToFT():  1   " << hexString(this) << printArray(ft.dimensions(), "  ft.dims") << printArray(ftSum.dimensions(), "  ftSum.dims") << std::endl;
    ftSum += ft;
//    std::cout << "Object::addToFT():  2   " << printArray(ft.dimensions(), "ft.dims") << printArray(ftSum.dimensions(), "  ftSum.dims") << std::endl;
}


void ObjectResult::addToPQ( const redux::image::FourierTransform& ft, const Array<complex_t> sj ) {
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



PatchResult::PatchResult( const MomfbdJob& j ) : myJob(j) {
    
}


uint64_t PatchResult::size( void ) const {
    uint64_t sz = Part::size();
    sz += index.size() + pos.size();
    /*for( const ObjectData::Ptr& obj: objects ) {
        sz += obj->size();
    }*/
    return sz;
}


uint64_t PatchResult::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = Part::pack(ptr);
    count += index.pack(ptr+count);
    count += pos.pack(ptr+count);
    /*for( const ObjectData::Ptr& obj: objects ) {
        count += obj->pack(ptr+count);
    }*/
    return count;
}


uint64_t PatchResult::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Part::unpack(ptr, swap_endian);
    count += index.unpack(ptr+count, swap_endian);
    count += pos.unpack(ptr+count, swap_endian);
    /*objects.clear();
    for( const std::shared_ptr<Object>& obj: myJob.getObjects() ) {
        ObjectData::Ptr objData(new ObjectData(obj));
        objects.push_back(objData);
        count += objData->unpack( ptr+count, swap_endian );
    }*/
    return count;
}


