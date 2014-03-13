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

Patch::Patch(uint32_t y, uint32_t x, uint32_t sz) : index(), roi() {
    if( sz > 0 ) {
        uint32_t halfsize = sz / 2;
        if (halfsize<y) {
            y -= halfsize;
        } else {
            y = 0;
        }
        if (halfsize<x) {
            x -= halfsize;
        } else {
            x = 0;
        }
        
        roi.first = Point(y,x);
        roi.last = Point(y+sz-1,x+sz-1);
    } else {
        // TBD: throw or print error ?
    }
}


void Patch::setIndex(uint32_t yid, uint32_t xid) {
    index.x = xid;
    index.y = yid;
}


size_t Patch::size( void ) const {
    size_t sz = Part::size();
    sz += index.size() + roi.size();
    return sz;
}


char* Patch::pack( char* ptr ) const {

    using redux::util::pack;

    ptr = Part::pack( ptr );
    ptr = index.pack( ptr );
    ptr = roi.pack( ptr );

    return ptr;
}


const char* Patch::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;

    ptr = Part::unpack( ptr, swap_endian );
    ptr = index.unpack( ptr, swap_endian );
    ptr = roi.unpack( ptr, swap_endian );

    return ptr;

}
