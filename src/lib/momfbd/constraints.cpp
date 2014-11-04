#include "redux/momfbd/constraints.hpp"

#include "redux/logger.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"

#include <string>

using namespace redux::momfbd;
using namespace redux::util;
using namespace std;

#define lg Logger::mlg
namespace {

    const string thisChannel = "momfbdch";

}


Constraints::Constraints() {
    LOG_DEBUG << "Constraints::Constraints()";
}

Constraints::~Constraints() {
    LOG_DEBUG << "Constraints::~Constraints()";
}


size_t Constraints::size( void ) const {
    size_t sz = ( 4 + 2 * nGroups + n_omap + n_imap ) * sizeof( uint32_t ); // ng,nv,nh,n_omap,n_imap,omap_n,imap_n,nrv
    for( uint32_t g=0; g < nGroups; ++g ) sz += nVerticals[g] * nHorizontals[g] * sizeof( double );
    for( uint32_t o=0; o < n_omap; ++o ) sz += omap_n[o] * sizeof( uint32_t );
    for( uint32_t i=0; i < n_imap; ++i ) sz += imap_n[i] * sizeof( uint32_t );
    sz += 5 * nrv * sizeof( uint32_t ); // rv
    return sz;
}


uint64_t Constraints::pack( char* ptr ) const {

    using redux::util::pack;

    uint64_t count = pack( ptr, nGroups );
    count += pack( ptr+count, nVerticals );
    count += pack( ptr+count, nHorizontals );
    for( uint32_t g = 0; g < nGroups; ++g ) count += pack( ptr+count, *( groups[g].get() ), nVerticals[g] * nHorizontals[g] );
    count += pack( ptr+count, n_omap );
    count += pack( ptr+count, omap_n );
    for( uint32_t o = 0; o < n_omap; ++o ) count += pack( ptr+count, o_map[o] );
    count += pack( ptr+count, n_imap );
    count += pack( ptr+count, imap_n );
    for( uint32_t i = 0; i < n_imap; ++i ) count += pack( ptr+count, i_map[i] );
    count += pack( ptr+count, nrv );
    uint32_t gvoi[4];
    for( uint32_t k = 0; k < nrv; ++k ) {
        memset(gvoi,0,4*sizeof(int32_t));
        for( uint32_t ig=0; ig < nGroups; ++ig ) // search for nullspace index
            for( uint32_t iv=0; iv < nVerticals[ig]; ++iv )
                if( rv[k].k == groups[ig].get()[iv] ) {
                    gvoi[0] = ig;
                    gvoi[1] = iv;
                    iv = nVerticals[ig];
                    ig = nGroups; // break;
                }
        for( uint32_t io=0; io < n_omap; ++io ) // search for object map index
            if( rv[k].o_map == o_map[io].data() ) {
                gvoi[2] = io;
                io = n_omap;
            }
        for( uint32_t ii=0; ii < n_imap; ++ii ) // search for image map index
            if( rv[k].i_map == i_map[ii].data() ) {
                gvoi[3] = ii;
                ii = n_imap; // break;
            }
        count += pack( ptr+count, gvoi, 4 );
        count += pack( ptr+count, rv[k].m );
    }
    
    return count;

}

uint64_t Constraints::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;

    uint64_t count = unpack( ptr, nGroups, swap_endian );
    nVerticals.resize( nGroups );
    nHorizontals.resize( nGroups );
    groups.resize( nGroups );
    count += unpack( ptr+count, nVerticals, swap_endian );
    count += unpack( ptr+count, nHorizontals, swap_endian );
    for( uint32_t g = 0; g < nGroups; ++g ) {
        groups[g] = sharedArray<double>( nVerticals[g], nHorizontals[g] );
        count += unpack( ptr+count, *( groups[g].get() ), nVerticals[g] * nHorizontals[g], swap_endian );
    }
    count += unpack( ptr+count, n_omap, swap_endian );
    omap_n.resize( n_omap );
    o_map.resize( n_omap );
    count += unpack( ptr+count, omap_n, swap_endian );
    for( uint32_t o = 0; o < n_omap; ++o ) {
        o_map[o].resize( omap_n[o] );
        count += unpack( ptr+count, o_map[o], swap_endian );
    }
    count += unpack( ptr+count, n_imap, swap_endian );
    imap_n.resize( n_imap );
    i_map.resize( n_imap );
    count += unpack( ptr+count, imap_n, swap_endian );
    for( uint32_t i = 0; i < n_imap; ++i ) {
        i_map[i].resize( imap_n[i] );
        count += unpack( ptr+count, i_map[i], swap_endian );
    }
    count += unpack( ptr+count, nrv, swap_endian );
    rv.resize( nrv );
    for( uint32_t k = 0; k < nrv; ++k ) {
        uint32_t g, v, o, i;
        count += unpack( ptr+count, g, swap_endian );
        count += unpack( ptr+count, v, swap_endian );
        count += unpack( ptr+count, o, swap_endian );
        count += unpack( ptr+count, i, swap_endian );
        count += unpack( ptr+count, rv[k].m, swap_endian );
        rv[k].o_map = o_map[o].data();
        rv[k].i_map = i_map[i].data();
        rv[k].k = groups[g].get()[v];
        rv[k].nk = nHorizontals[g];
    }
    return count;

}
