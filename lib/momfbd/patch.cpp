#include "redux/momfbd/patch.hpp"

#include "redux/util/datautil.hpp"

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

void Patch::setIndex(uint32_t yid, uint32_t xid) {
    index.x = xid;
    index.y = yid;
}


size_t Patch::nPixels(void) {
    return (last.x-first.x+1)*(last.y-first.y+1);
}


size_t Patch::size( void ) const {
    size_t sz = Part::size();
    sz += 3*Point::size() + PointF::size();
    sz += sizeof(size_t) + dataSize;
    return sz;
}


char* Patch::pack( char* ptr ) const {

    using redux::util::pack;

    ptr = Part::pack( ptr );
    ptr = index.pack( ptr );
    ptr = first.pack( ptr );
    ptr = last.pack( ptr );
    ptr = residualTilts.pack( ptr );
    ptr = pack(ptr,dataSize);
    memcpy(ptr,data.get(),dataSize);
    return ptr+dataSize;
}


const char* Patch::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;

    ptr = Part::unpack( ptr, swap_endian );
    ptr = index.unpack( ptr, swap_endian );
    ptr = first.unpack( ptr, swap_endian );
    ptr = last.unpack( ptr, swap_endian );
    ptr = residualTilts.unpack( ptr, swap_endian );
    ptr = unpack(ptr,dataSize,swap_endian);
    data = sharedArray<char>(dataSize);
    memcpy(data.get(),ptr,dataSize);
    return ptr+dataSize;

}
