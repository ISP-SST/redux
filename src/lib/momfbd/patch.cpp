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


uint64_t Patch::pack( char* ptr ) const {

    using redux::util::pack;
    
    uint64_t count = Part::pack( ptr );
    count += index.pack( ptr+count );
    count += first.pack( ptr+count );
    count += last.pack( ptr+count );
    count += residualTilts.pack( ptr+count );
    count += pack(ptr+count,dataSize);
    memcpy(ptr+count,data.get(),dataSize);
    return count+dataSize;
}


uint64_t Patch::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;
    uint64_t count = Part::unpack( ptr, swap_endian );
    count += index.unpack( ptr+count, swap_endian );
    count += first.unpack( ptr+count, swap_endian );
    count += last.unpack( ptr+count, swap_endian );
    count += residualTilts.unpack( ptr+count, swap_endian );
    count += unpack(ptr+count,dataSize,swap_endian);
    data = sharedArray<char>(dataSize);
    memcpy(data.get(),ptr+count,dataSize);
    return count+dataSize;

}
