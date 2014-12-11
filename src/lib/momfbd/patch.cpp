#include "redux/momfbd/patch.hpp"

#include "redux/momfbd/momfbdjob.hpp"

#include "redux/util/datautil.hpp"

#include "redux/momfbd/defines.hpp"

using namespace redux::momfbd;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

    const string thisChannel = "momfbdpatch";

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
*/

void PatchData::setIndex(uint32_t yid, uint32_t xid) {
    index.x = xid;
    index.y = yid;
}


size_t PatchData::nPixels(void) {
    return (last.x-first.x+1)*(last.y-first.y+1);
}


size_t PatchData::nPixelsX(void) {
    return (last.x-first.x+1);
}


size_t PatchData::nPixelsY(void) {
    return (last.y-first.y+1);
}


size_t PatchData::size( void ) const {
    size_t sz = Part::size();
    sz += 3*Point::size() + PointF::size();
    sz += images.size();
    return sz;
}


uint64_t PatchData::pack( char* ptr ) const {

    using redux::util::pack;
    
    uint64_t count = Part::pack( ptr );
    count += index.pack( ptr+count );
    count += first.pack( ptr+count );
    count += last.pack( ptr+count );
    count += residualOffset.pack(ptr+count );
    count += images.pack(ptr);
    return count;
}


uint64_t PatchData::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;
    uint64_t count = Part::unpack( ptr, swap_endian );
    count += index.unpack( ptr+count, swap_endian );
    count += first.unpack( ptr+count, swap_endian );
    count += last.unpack( ptr+count, swap_endian );
    count += residualOffset.unpack( ptr+count, swap_endian );
    count += images.unpack(ptr+count, swap_endian);
    return count;

}


size_t PatchResult::size( void ) const {
    size_t sz = Part::size();
    sz += restoredObjects.size();
    sz += modeCoefficients.size();
    sz += PSFs.size();
    return sz;
}


uint64_t PatchResult::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = Part::pack( ptr );
    count += restoredObjects.pack(ptr+count);
    count += modeCoefficients.pack(ptr+count);
    count += PSFs.pack(ptr);
    return count;
}


uint64_t PatchResult::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Part::unpack( ptr, swap_endian );
    count += restoredObjects.unpack( ptr+count, swap_endian );
    count += modeCoefficients.unpack( ptr+count, swap_endian );
    count += PSFs.unpack( ptr+count, swap_endian );
    return count;
}


