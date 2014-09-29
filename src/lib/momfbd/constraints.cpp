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


char* Constraints::pack( char* ptr ) const {

    using redux::util::pack;

    ptr = pack( ptr, nGroups );
    ptr = pack( ptr, nVerticals );
    ptr = pack( ptr, nHorizontals );
    for( uint32_t g = 0; g < nGroups; ++g ) ptr = pack( ptr, *( groups[g].get() ), nVerticals[g] * nHorizontals[g] );
    ptr = pack( ptr, n_omap );
    ptr = pack( ptr, omap_n );
    for( uint32_t o = 0; o < n_omap; ++o ) ptr = pack( ptr, o_map[o] );
    ptr = pack( ptr, n_imap );
    ptr = pack( ptr, imap_n );
    for( uint32_t i = 0; i < n_imap; ++i ) ptr = pack( ptr, i_map[i] );
    ptr = pack( ptr, nrv );
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
        ptr = pack( ptr, gvoi, 4 );
        ptr = pack( ptr, rv[k].m );
    }
    
    return ptr;

}

const char* Constraints::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;

    ptr = unpack( ptr, nGroups, swap_endian );
    nVerticals.resize( nGroups );
    nHorizontals.resize( nGroups );
    groups.resize( nGroups );
    ptr = unpack( ptr, nVerticals, swap_endian );
    ptr = unpack( ptr, nHorizontals, swap_endian );
    for( uint32_t g = 0; g < nGroups; ++g ) {
        groups[g] = sharedArray<double>( nVerticals[g], nHorizontals[g] );
        ptr = unpack( ptr, *( groups[g].get() ), nVerticals[g] * nHorizontals[g], swap_endian );
    }
    ptr = unpack( ptr, n_omap, swap_endian );
    omap_n.resize( n_omap );
    o_map.resize( n_omap );
    ptr = unpack( ptr, omap_n, swap_endian );
    for( uint32_t o = 0; o < n_omap; ++o ) {
        o_map[o].resize( omap_n[o] );
        ptr = unpack( ptr, o_map[o], swap_endian );
    }
    ptr = unpack( ptr, n_imap, swap_endian );
    imap_n.resize( n_imap );
    i_map.resize( n_imap );
    ptr = unpack( ptr, imap_n, swap_endian );
    for( uint32_t i = 0; i < n_imap; ++i ) {
        i_map[i].resize( imap_n[i] );
        ptr = unpack( ptr, i_map[i], swap_endian );
    }
    ptr = unpack( ptr, nrv, swap_endian );
    rv.resize( nrv );
    for( uint32_t k = 0; k < nrv; ++k ) {
        uint32_t g, v, o, i;
        ptr = unpack( ptr, g, swap_endian );
        ptr = unpack( ptr, v, swap_endian );
        ptr = unpack( ptr, o, swap_endian );
        ptr = unpack( ptr, i, swap_endian );
        ptr = unpack( ptr, rv[k].m, swap_endian );
        rv[k].o_map = o_map[o].data();
        rv[k].i_map = i_map[i].data();
        rv[k].k = groups[g].get()[v];
        rv[k].nk = nHorizontals[g];
    }
    return ptr;

}
