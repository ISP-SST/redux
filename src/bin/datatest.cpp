#include "redux/file/filefits.hpp"
#include "redux/util/array.hpp"
#include "redux/util/ricecompress.hpp"

#include <iostream>
#include <string>

using namespace redux;
using namespace redux::file;
using namespace redux::util;
using namespace std;

const string bdir = "/scratch/tomas/debug/";
//const string bdir = "/scratch/tomas_local/rice/";

const string files[] = { bdir + "/bad/file.fits",
                         bdir + "/good/file.fits",
                         bdir + "/dark/file.fits",
                         bdir + "/flat/file.fits",
                         bdir + "/ph/file.fits"
};

string use_file = files[1];


int main( int argc, char *argv[] ) {

    //static void read( const std::string& filename, redux::util::Array<T>& data, std::shared_ptr<redux::file::Fits>& hdr );
    // implement this !!!
    // static void write( const std::string& filename, const char* data, std::shared_ptr<redux::file::Fits> hdr, bool compress = false, int slice=5 );
    int dbg_it = -1;
    if (argc > 1)
        dbg_it = atoi(argv[1]);
    
    Array<int16_t> data;
    shared_ptr<Fits> hdr;

    Fits::read( use_file, data, hdr );
    
    size_t nY = data.dimSize(1);
    size_t nX = data.dimSize(2);
    size_t frameSize = nX*nY;
    size_t dataSize = frameSize*2;
    size_t blockSize = dataSize*10; // should be more than enough
    cout << "nX=" << nX << " nY=" << nY << " frameSize=" << dataSize << " blockSize=" << blockSize << endl;
    if( !blockSize ) return EXIT_FAILURE;
    
    
    unique_ptr<uint8_t[]> tmpData( new uint8_t[blockSize] );
    
    uint8_t* ptr_good = tmpData.get();
    uint8_t* ptr_bad = ptr_good + 5*frameSize;
    
    int16_t* const inPtr = data.get();
    int16_t* out_ref = inPtr + frameSize;
    int16_t* out_good = out_ref + frameSize;      // just mangle frames > 1
    int16_t* out_bad = out_good + frameSize;
    int16_t* out_fix = out_bad + frameSize;
    
    size_t rice_bsize = 32;
    //int rice_comp16( const int16_t* in, size_t inSize, uint8_t* out, size_t outSize, size_t blockSize );
    
    int csize_good = rice_comp16( inPtr, frameSize, ptr_good, 5*frameSize, rice_bsize );
    int csize_bad = rice_comp16_bad( inPtr, frameSize, ptr_bad, 5*frameSize, rice_bsize );
    
    cout << __LINE__ << "  csize_good = " << csize_good << "  csize_bad = " << csize_bad << endl;
    if( csize_bad > dataSize) {
        cout << __LINE__ << " WARN:   csize_bad = " << csize_bad << " bigger might have overflown buffer !!! " << endl;
    }
    //cout << printBits( ptr_good, 10 ) << endl;
    //cout << printBits( ptr_bad, 10 ) << endl;
    
    //int rice_decomp16( const uint8_t* in, size_t inSize, int16_t* out, size_t outSize, size_t blockSize );
    size_t decomp_size = frameSize;

    if( decomp_size > 2000 ) dbg_it = -1;
cout << "----------------------------------------------------------------------------" << endl;
    int dsize_ref = rice_decomp16( ptr_good, 5*frameSize, out_ref, decomp_size, rice_bsize );
cout << "----------------------------------------------------------------------------" << endl;
    int dsize_good = rice_decomp16_dbg( ptr_good, 5*frameSize, out_good, decomp_size, rice_bsize );
cout << "----------------------------------------------------------------------------" << endl;
    int dsize_bad = rice_decomp16( ptr_bad, 5*frameSize, out_bad, decomp_size, rice_bsize );
cout << "----------------------------------------------------------------------------" << endl;
    int dsize_fix = rice_decomp16_fix( ptr_bad, 5*frameSize, out_fix, decomp_size, rice_bsize );
cout << "----------------------------------------------------------------------------" << endl;

// decomp_size = 2*rice_bsize; //frameSize;
// dsize_good = rice_decomp16_dbg( ptr_good, frameSize, out_good, decomp_size, rice_bsize);
// dsize_bad = rice_decomp16( ptr_bad, frameSize, out_bad, decomp_size, rice_bsize );
// dsize_fix = rice_decomp16_fix( ptr_bad, frameSize, out_fix, decomp_size, rice_bsize );
    
    size_t offset=0;
    size_t pos=0;
    size_t cmp_bad(0), cmp_fix(0), cmp_ref(0);
    int val_in_bad(0), val_in_fix(0), val_bad(0), val_fix(0), val_in_ref(0), val_ref(0);
    for( size_t i(0); i<decomp_size; ++i ) {
        if( !cmp_ref && (inPtr[i] != out_ref[i]) ) {
            pos = cmp_ref = i;
            val_in_ref = inPtr[i];
            val_ref = out_bad[i];
        }
        if( false && (inPtr[i] != out_bad[i]) ) {
            pos = cmp_bad = i;
            val_in_bad = inPtr[i];
            val_bad = out_bad[i];
        }
        if( !cmp_fix && (inPtr[i] != out_fix[i]) ) {
            pos = cmp_fix = i;
            val_in_fix = inPtr[i];
            val_fix = out_fix[i];
        }
        if( pos ) {
            break;
        }
    }
    //cmp_good = memcmp( inPtr, out_good, frameSize );
    //cmp_bad = memcmp( inPtr, out_bad, frameSize );
    //cmp_fix = memcmp( inPtr, out_fix, frameSize );
//     while( abs(cmp_good)+abs(cmp_bad)+abs(cmp_fix) < 1 ) {
//         cmp_good = memcmp( inPtr+offset, out_good, rice_bsize );
//         cmp_bad = memcmp( inPtr+offset, out_bad, rice_bsize );
//         cmp_fix = memcmp( inPtr+offset, out_fix, rice_bsize );
//         offset += rice_bsize;
//     }

    cout << "  in_vs_ref @" << cmp_ref << " (" << val_in_ref << " vs. " << val_ref << ")"
         << "  in_vs_bad @" << cmp_bad << " (" << val_in_bad << " vs. " << val_bad << ")"
         << "  in_vs_fix @" << cmp_fix << " (" << val_in_fix << " vs. " << val_fix << ")" << endl;
    cout << __LINE__ << "  dsize_ref = " << dsize_ref << "  dsize_good = " << dsize_good << "  dsize_bad = " << dsize_bad << "  dsize_fix = " << dsize_fix << endl;
    if( dbg_it >= 0 ) {
        cout << __LINE__ << " DSFGSDFGSF" << endl;
        offset = (dbg_it/rice_bsize)*rice_bsize;
    } else {
        offset = (pos/rice_bsize)*rice_bsize;
    }
    decomp_size = 2*rice_bsize;
    cout << __LINE__ << "  offset = " << offset <<  "  decomp_size = " << decomp_size << endl;
    cout << endl << printArray( inPtr+offset, decomp_size, "in" ) << endl;
    //cout << endl << printArray( out_good+offset, decomp_size, "good" ) << endl;
    cout << endl << printArray( out_bad+offset, decomp_size, "bad" ) << endl;
    cout << endl << printArray( out_fix+offset, decomp_size, "fix" ) << endl;
    
    
    return EXIT_SUCCESS;

}

