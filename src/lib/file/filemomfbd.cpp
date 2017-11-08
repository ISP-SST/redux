#include "redux/file/filemomfbd.hpp"

#include <iostream>
#include <cmath>     // NAN
#include <cstdlib>     // atof


using namespace redux::file;
using namespace redux::util;
using namespace std;

#if REDUX_BYTE_ORDER == REDUX_LITTLE_ENDIAN
static const int system_is_big_endian = 0;
#elif REDUX_BYTE_ORDER == REDUX_BIG_ENDIAN
static const int system_is_big_endian = 1;
#else
#error REDUX_BYTE_ORDER not set
#endif


uint8_t FileMomfbd::PatchInfo::parse ( ifstream& file, const bool& swapNeeded, const float& version ) {

    int64_t tmpSize, patchSize ( 0 );
    uint8_t dataMask( 0 );

    offset = file.tellg();

    patchSize += readOrThrow ( file, region, 4, "PatchInfo:region" );
    
    if( version >= 20110714.0 ){
        patchSize += readOrThrow ( file, &offx, 1, "PatchInfo:offx" );
        patchSize += readOrThrow ( file, &offy, 1, "PatchInfo:offy" );
    }
    
    patchSize += readOrThrow ( file, &nChannels, 1, "PatchInfo:nChannels" );
    if ( swapNeeded ) {
        swapEndian ( &region, 4 );
        swapEndian ( nChannels );
    }

    if ( region[0] > region[1] ) swap ( region[0], region[1] );
    if ( region[2] > region[3] ) swap ( region[2], region[3] );

    nPixelsX = region[1] - region[0] + 1;
    nPixelsY = region[3] - region[2] + 1;

    nim.reset ( new int32_t [nChannels], [] ( int32_t * p ) { delete[] p; } );
    dx.reset ( new int32_t [nChannels], [] ( int32_t * p ) { delete[] p; } );
    dy.reset ( new int32_t [nChannels], [] ( int32_t * p ) { delete[] p; } );

    patchSize += readOrThrow ( file, nim.get(), nChannels, "PatchInfo:nim" );
    patchSize += readOrThrow ( file, dx.get(), nChannels, "PatchInfo:dx" );
    patchSize += readOrThrow ( file, dy.get(), nChannels, "PatchInfo:dy" );

    if ( swapNeeded ) {
        swapEndian ( nim.get(), nChannels );
        swapEndian ( dx.get(), nChannels );
        swapEndian ( dy.get(), nChannels );
    }

    char hasImage;
    patchSize += readOrThrow ( file, &hasImage, 1, "PatchInfo:hasImage" );

    if ( hasImage ) {
        dataMask |= MOMFBD_IMG;
        imgPos = file.tellg();
        tmpSize = nPixelsX * nPixelsY * sizeof ( float );
        patchSize += tmpSize;
        file.seekg ( tmpSize, ios_base::cur );
    }

    npsf = 0;
    patchSize += readOrThrow ( file, &npsf, 1, "PatchInfo:npsf" );
    if ( swapNeeded ) {
        swapEndian ( npsf );
    }

    if ( npsf ) {
        dataMask |= MOMFBD_PSF;
        psfPos = file.tellg();
        tmpSize = npsf * nPixelsX * nPixelsY * sizeof ( float );
        patchSize += tmpSize;
        file.seekg ( tmpSize, ios_base::cur );
    }

    nobj = 0;
    patchSize += readOrThrow ( file, &nobj, 1, "PatchInfo:nobj" );
    if ( swapNeeded ) swapEndian ( nobj );
    if ( nobj ) {
        dataMask |= MOMFBD_OBJ;
        objPos = file.tellg();
        tmpSize = nobj * nPixelsX * nPixelsY * sizeof ( float );
        patchSize += tmpSize;
        file.seekg ( tmpSize, ios_base::cur );
    }

    nres = 0;
    patchSize += readOrThrow ( file, &nres, 1, "PatchInfo:nres" );
    if ( swapNeeded ) swapEndian ( nres );
    if ( nres ) {
        dataMask |= MOMFBD_RES;
        resPos = file.tellg();
        tmpSize = nres * nPixelsX * nPixelsY * sizeof ( float );
        patchSize += tmpSize;
        file.seekg ( tmpSize, ios_base::cur );
    }

    nalpha = nm = 0;
    patchSize += readOrThrow ( file, &nalpha, 1, "PatchInfo:nalpha" );
    if ( swapNeeded ) swapEndian ( nalpha );
    if ( nalpha ) {
        dataMask |= MOMFBD_ALPHA;
        patchSize += readOrThrow ( file, &nm, 1, "PatchInfo:nm" );
        if ( swapNeeded ) swapEndian ( nm );
        alphaPos = file.tellg();
        tmpSize = nalpha * nm * sizeof ( float );
        patchSize += tmpSize;
        file.seekg ( tmpSize, ios_base::cur );
    }

    ndiv = 0;
    if ( version >= 20100726.0 ) {  // check if this is a version that can have div info...
        patchSize += readOrThrow ( file, &ndiv, 1, "PatchInfo:ndiv" );
        if ( swapNeeded ) swapEndian ( ndiv );
    }

    nphx = nPixelsX / 2;
    nphy = nPixelsY / 2;
    if ( ndiv ) {
        dataMask |= MOMFBD_DIV;
        if ( version >= 20110708.0 ) {
            patchSize += readOrThrow ( file, &nphx, 1, "PatchInfo:nphx" );
            patchSize += readOrThrow ( file, &nphy, 1, "PatchInfo:nphy" );
            if( swapNeeded) {
                swapEndian ( nphx );
                swapEndian ( nphy );
            }

        } else cout << "WARNING: diversity of versions < 20110708.0 may have wrong dimensions!  (version = " << version << ")" << endl;

        diversityPos = file.tellg();
        tmpSize = ndiv * nphy *  nphx * sizeof ( float );
        if( version >= 20110916.0 ) tmpSize += ndiv;        // + byte with diversity-type.
        patchSize += tmpSize;
        file.seekg ( tmpSize, ios_base::cur );
    }

    if ( patchSize != ( file.tellg() - offset ) ) {
        throw ios_base::failure ( "Failed to read Momfbd Patch. Size mismatch: "+to_string ( patchSize ) +" <-> " + to_string ( file.tellg() - offset ) );
    }

    return dataMask;
}


size_t FileMomfbd::PatchInfo::load ( ifstream& file, char* ptr, const bool& swapNeeded, const float& version, uint8_t loadMask, int verbosity, uint8_t alignTo ) const {

    size_t count(0);
    float* fPtr;
    short int* tmpInt;

    size_t tmpSize, nxny;

    // region
    tmpInt = ( short int* ) ptr;
    for ( int i = 0; i < 4; ++i ) {
        tmpInt[i] = region[i];
    }
    if ( version >= 20110714.0 ) {
        tmpInt[4] = offx;
        tmpInt[5] = offy;
        count += 2 * sizeof ( short int );
    }
    count += 4 * sizeof ( short int );

    tmpInt = ( short int* ) (ptr+count);
    for ( int i = 0; i < nChannels; ++i ) {
        tmpInt[ i + 0 * nChannels ] = dx.get() [ i ];
        tmpInt[ i + 1 * nChannels ] = dy.get() [ i ];
        tmpInt[ i + 2 * nChannels ] = nim.get() [ i ];
    }
    count += 3 * nChannels * sizeof ( short int );

    while ( count % alignTo ) count++;

    nxny = nPixelsX * nPixelsY;
    if ( ( loadMask & MOMFBD_IMG ) && imgPos ) {
        file.seekg ( imgPos );
        fPtr = reinterpret_cast<float*> ( ptr+count );
        count += readOrThrow ( file, fPtr, nxny, "MomfbdPatch:img" );

        if ( swapNeeded ) {
            swapEndian ( fPtr, nxny );
        }
    }

    if ( ( loadMask & MOMFBD_PSF ) && npsf ) {
        file.seekg ( psfPos );
        fPtr = reinterpret_cast<float*> ( ptr+count );
        count += readOrThrow ( file, fPtr, npsf * nxny, "MomfbdPatch:psf" );
        if ( swapNeeded ) {
            swapEndian ( fPtr, npsf * nxny );
        }
    }

    if ( ( loadMask & MOMFBD_OBJ ) && nobj ) {
        file.seekg ( objPos );
        fPtr = reinterpret_cast<float*> ( ptr+count );
        count += readOrThrow ( file, fPtr, nobj * nxny, "MomfbdPatch:obj" );
        if ( swapNeeded ) {
            swapEndian ( fPtr, nobj * nxny );
        }
    }

    if ( ( loadMask & MOMFBD_RES ) && nres ) {
        file.seekg ( resPos );
        fPtr = reinterpret_cast<float*> ( ptr+count );
        count += readOrThrow ( file, fPtr, nres * nxny, "MomfbdPatch:res" );
        if ( swapNeeded ) {
            swapEndian ( fPtr, nres * nxny );
        }
    }
    if ( ( loadMask & MOMFBD_ALPHA ) && nalpha ) {
        tmpSize = nalpha * nm;
        file.seekg ( alphaPos );
        fPtr = reinterpret_cast<float*> ( ptr+count );
        size_t cnt = readOrThrow ( file, fPtr, tmpSize, "MomfbdPatch:alpha" );
        count += cnt;
        if ( swapNeeded ) {
            swapEndian ( fPtr, tmpSize );
        }
    }

    if ( ( loadMask & MOMFBD_DIV ) && ndiv ) {
        tmpSize = ndiv * nphy * nphx;
        file.seekg ( diversityPos );
        fPtr = reinterpret_cast<float*> ( ptr+count );
        if( version >= 20110916.0 ) {           // + byte with diversity-type.
            for( int d=0; d<ndiv; ++d ) {
                char typ;
                readOrThrow ( file, &typ, 1, "MomfbdPatch:div-typ" );   // just discard div-type for now.
                count += readOrThrow ( file, ptr+count, nphy*nphx*sizeof(float), "MomfbdPatch:div" );
            }
        } else {
            count += readOrThrow ( file, ptr+count, tmpSize*sizeof(float), "MomfbdPatch:div" );
        }
        if ( swapNeeded ) {
            swapEndian ( fPtr, tmpSize );
        }

    }

    return count;
}

void FileMomfbd::PatchInfo::write ( ofstream& file, const char* data, const float& version, uint8_t writeMask ) {

    char tmp8;
    const float* fPtr;
    
    if ( region[0] > region[1] ) swap ( region[0], region[1] );
    if ( region[2] > region[3] ) swap ( region[2], region[3] );

    nPixelsX = region[1] - region[0] + 1;
    nPixelsY = region[3] - region[2] + 1;
    
    writeOrThrow ( file, region, 4, "PatchInfo:region" );
    
    if( version >= 20110714.0 ){
        writeOrThrow ( file, &offx, 1, "PatchInfo:offx" );
        writeOrThrow ( file, &offy, 1, "PatchInfo:offy" );
    }
    
    writeOrThrow ( file, &nChannels, 1, "PatchInfo:nChannels" );
    
    writeOrThrow ( file, nim.get(), nChannels, "PatchInfo:nim" );
    writeOrThrow ( file, dx.get(), nChannels, "PatchInfo:dx" );
    writeOrThrow ( file, dy.get(), nChannels, "PatchInfo:dy" );
    
    tmp8 = ( writeMask&MOMFBD_IMG );
    writeOrThrow ( file, &tmp8, 1, "FileMomfbd:withImage" );
    if ( tmp8 ) {
        fPtr = reinterpret_cast<const float*> ( data + imgPos );
        writeOrThrow ( file, fPtr, nPixelsX*nPixelsY, "FileMomfbd:IMG" );
    }
    
    if ( !(writeMask&MOMFBD_PSF) ) npsf = 0;
    writeOrThrow ( file, &npsf, 1, "FileMomfbd:npsf" );
    if ( npsf ) {
        fPtr = reinterpret_cast<const float*> ( data + psfPos );
        writeOrThrow ( file, fPtr, npsf*nPixelsX*nPixelsY, "FileMomfbd:PSF" );
    }
    
    if ( !(writeMask&MOMFBD_OBJ) ) nobj = 0;
    writeOrThrow ( file, &nobj, 1, "FileMomfbd:nobj" );
    if ( nobj ) {
        fPtr = reinterpret_cast<const float*> ( data + objPos );
        writeOrThrow ( file, fPtr, nobj*nPixelsX*nPixelsY, "FileMomfbd:OBJ" );
    }
    
    if ( !(writeMask&MOMFBD_RES) ) nres = 0;
    writeOrThrow ( file, &nres, 1, "FileMomfbd:nres" );
    if ( nres ) {
        fPtr = reinterpret_cast<const float*> ( data + resPos );
        writeOrThrow ( file, fPtr, nres*nPixelsX*nPixelsY, "FileMomfbd:RES" );
    }
    
    if( !(writeMask&MOMFBD_ALPHA) ) nalpha = 0;
    writeOrThrow ( file, &nalpha, 1, "FileMomfbd:nalpha" );
    if ( nalpha ) {
        writeOrThrow ( file, &nm, 1, "FileMomfbd:nm" );
        if(nm) {
            fPtr = reinterpret_cast<const float*> ( data + alphaPos );
            writeOrThrow ( file, fPtr, nalpha*nm, "FileMomfbd:ALPHA" );
        }
    }
    
    if( !(writeMask&MOMFBD_DIV) ) ndiv = 0;
    if ( version >= 20100726.0 ) {  // check if this is a version that can have div info...
        writeOrThrow ( file, &ndiv, 1, "FileMomfbd:ndiv" );
        if ( ndiv ) {
            if ( version >= 20110708.0 ) {
                writeOrThrow ( file, &nphx, 1, "PatchInfo:nphx" );
                writeOrThrow ( file, &nphy, 1, "PatchInfo:nphy" );
            }
            fPtr = reinterpret_cast<const float*> ( data + diversityPos );
            if( version >= 20110916.0 ) {           // + byte with diversity-type.
                for( int d=0; d<ndiv; ++d ) {
                    char typ=0;             // TODO: "real" type
                    writeOrThrow ( file, &typ, 1, "MomfbdPatch:div-typ" );   // just discard div-type for now.
                    writeOrThrow ( file, fPtr, nphy*nphx, "MomfbdPatch:div" );
                    fPtr += nphy*nphx;
                }
            } else {
                writeOrThrow ( file, fPtr, ndiv * nphy *  nphx, "MomfbdPatch:div" );
            }
        }
    }

}


FileMomfbd::FileMomfbd ( void ) : version ( 0 ), pix2cf(NAN), cf2pix(NAN), 
    nChannels ( 0 ), nFileNames ( 0 ), nPH ( 0 ), nModes ( 0 ), nPatchesX( 0 ), nPatchesY( 0 ), nPoints(0), phOffset ( 0 ),
    modesOffset ( 0 ), filenameOffset ( 0 ), patchDataSize ( 0 ), headerSize ( 0 ), dataMask(0), swapNeeded ( false )  {

}


FileMomfbd::FileMomfbd ( const std::string& filename ) : version ( 0 ), pix2cf(NAN), cf2pix(NAN),
    nChannels ( 0 ), nFileNames ( 0 ), nPH ( 0 ), nModes ( 0 ), nPatchesX( 0 ), nPatchesY( 0 ), nPoints(0), phOffset ( -1 ),
    modesOffset ( -1 ), filenameOffset ( -1 ), patchDataSize ( 0 ), headerSize ( 0 ), dataMask(0), swapNeeded ( false ) {

    read ( filename );
    
}


size_t FileMomfbd::getPatchSize( const FileMomfbd* const info, uint8_t loadMask, const float& version, size_t alignTo ) {
    
    if ( info->patches.nElements() == 0 ) return 0;                            // no patches
    
    size_t patchSize  = 4 * sizeof(short int);                                 // xl, yl, xh, yl
    if ( version >= 20110714.0 ) {
        patchSize += 2*sizeof(short int);                                      // off, offy
    }
    patchSize += 3*info->nChannels * sizeof( short int );                      // dx, dy, nim
    while( patchSize % alignTo ) patchSize++;                                  // pad if not on boundary (needed if nChannels is odd)

    if ( !(loadMask & MOMFBD_PATCH) ) return patchSize;
    
    const FileMomfbd::PatchInfo& tmpPatch = info->patches(0);

    size_t nFloats = 0;
    if ( loadMask & MOMFBD_IMG ) {
        nFloats += tmpPatch.nPixelsX * tmpPatch.nPixelsY;
    }
    if ( loadMask & MOMFBD_PSF ) {
        nFloats += tmpPatch.npsf * tmpPatch.nPixelsX * tmpPatch.nPixelsY;
    }
    if ( loadMask & MOMFBD_OBJ ) {
        nFloats += tmpPatch.nobj * tmpPatch.nPixelsX * tmpPatch.nPixelsY;
    }
    if ( loadMask & MOMFBD_RES ) {
        nFloats += tmpPatch.nPixelsX * tmpPatch.nPixelsY * tmpPatch.nres;
    }
    if ( loadMask & MOMFBD_ALPHA ) {
        nFloats += tmpPatch.nm * tmpPatch.nalpha;
    }
    if ( loadMask & MOMFBD_DIV ) {
        nFloats += tmpPatch.nphx * tmpPatch.nphy * tmpPatch.ndiv;
    }
    patchSize += nFloats * sizeof(float);

    return patchSize;

}



void FileMomfbd::clear(void) {
    
    version = 0;
    memset( region, 0, 4*sizeof(int32_t) );
    nChannels = 0;
    nFileNames = 0;
    nPH = 0;
    nModes = 0;
    nPatchesX = nPatchesY = nPoints = 0;
    phOffset = modesOffset = filenameOffset = -1;
    patchDataSize = 0;
    headerSize = 0;
    swapNeeded = false;
    
    dateString.clear();
    timeString.clear();
    versionString.clear();
    
    fileNames.clear();

    clipStartX = clipEndX = clipStartY = clipEndY = 0;

    patches.resize();
    
}


void FileMomfbd::read ( std::ifstream& file ) {

    headerSize = 0;
    file.seekg ( 0 );

    char tmp8;
    int32_t tmp32;
    vector<char> tmpStr ( 1024, 0 );

    // first byte = endian-flag
    readOrThrow ( file, &tmp8, 1, "FileMomfbd:endian" );
    swapNeeded = system_is_big_endian ^ tmp8;

    // version string
    readOrThrow ( file, &tmp32, 1, "FileMomfbd:version-length" );
    if ( swapNeeded ) swapEndian ( &tmp32 );
    tmpStr.reserve ( tmp32 );   // strings are stored including \0 termination, so no additional char needed
    readOrThrow ( file, & ( tmpStr[0] ), tmp32, "FileMomfbd:version-string" );
    versionString.assign ( tmpStr.begin(), tmpStr.begin() + tmp32 );
    version = atof ( versionString.c_str() );

    // time string
    readOrThrow ( file, &tmp32, 1, "FileMomfbd:time-length" );
    if ( swapNeeded ) swapEndian ( &tmp32 );
    tmpStr.reserve ( tmp32 );   // strings are stored including \0 termination, so no additional char needed
    readOrThrow ( file, & ( tmpStr[0] ), tmp32, "FileMomfbd:time-string" );
    timeString.assign ( tmpStr.begin(), tmpStr.begin() + tmp32 );

    // date string
    readOrThrow ( file, &tmp32, 1, "FileMomfbd:date-length" );
    if ( swapNeeded ) swapEndian ( &tmp32 );
    tmpStr.reserve ( tmp32 );   // strings are stored including \0 termination, so no additional char needed
    readOrThrow ( file, & ( tmpStr[0] ), tmp32, "FileMomfbd:date-string" );
    dateString.assign ( tmpStr.begin(), tmpStr.begin() + tmp32 );

    readOrThrow ( file, &tmp8, 1, "FileMomfbd:hasModes" );
    if ( tmp8 ) {
        dataMask |= MOMFBD_MODES;
        if ( version >= 20110714.0 ) {
            readOrThrow ( file, &pix2cf, 1, "FileMomfbd:pix2cf" );
            readOrThrow ( file, &cf2pix, 1, "FileMomfbd:cf2pix" );
            if ( swapNeeded ) {
                swapEndian ( pix2cf );
                swapEndian ( cf2pix );
            }
        }
        
        readOrThrow ( file, &nPH, 1, "FileMomfbd:nPH" );
        readOrThrow ( file, &nModes, 1, "FileMomfbd:nModes" );
        if ( swapNeeded ) {
            swapEndian ( nPH );
            swapEndian ( nModes );
        }
        size_t tmpSz = nPH*nPH*sizeof(float);
        phOffset = file.tellg();
        file.seekg ( tmpSz, ios_base::cur );
        modesOffset = file.tellg();
        if(version >= 20120201.0) tmpSz++;              // Add a byte to specify mode-type
        file.seekg ( nModes*tmpSz, ios_base::cur );
    } else nModes = nPH = 0;
    
    readOrThrow ( file, &nChannels, 1, "FileMomfbd:nChannels" );
    if ( swapNeeded ) {
        swapEndian ( nChannels );
    }

    // global clip
    clipStartX.reset ( new int16_t[ nChannels ], [] ( int16_t * p ) { delete[] p; } );
    clipEndX.reset ( new int16_t[ nChannels ], [] ( int16_t * p ) { delete[] p; } );
    clipStartY.reset ( new int16_t[ nChannels ], [] ( int16_t * p ) { delete[] p; } );
    clipEndY.reset ( new int16_t[ nChannels ], [] ( int16_t * p ) { delete[] p; } );

    readOrThrow ( file, clipStartX.get(), nChannels, "FileMomfbd:clipStartX" );
    readOrThrow ( file, clipEndX.get(), nChannels, "FileMomfbd:clipEndX" );
    readOrThrow ( file, clipStartY.get(), nChannels, "FileMomfbd:clipStartY" );
    readOrThrow ( file, clipEndY.get(), nChannels, "FileMomfbd:clipEndY" );
    if ( swapNeeded ) {
        swapEndian ( clipStartX.get(), nChannels );
        swapEndian ( clipEndX.get(), nChannels );
        swapEndian ( clipStartY.get(), nChannels );
        swapEndian ( clipEndY.get(), nChannels );
    }

    readOrThrow ( file, &nPatchesX, 1, "FileMomfbd:nPatchesX" );    // slow index of patch-array stored first
    readOrThrow ( file, &nPatchesY, 1, "FileMomfbd:nPatchesY" );
    readOrThrow ( file, &nPoints, 1, "FileMomfbd:nPoints" );

    if ( swapNeeded ) {
        swapEndian ( nPatchesX );
        swapEndian ( nPatchesY );
        swapEndian ( nPoints );
    }

    region[0] = region[2] = numeric_limits<int32_t>::max();
    region[1] = region[3] = numeric_limits<int32_t>::min();
    patches.resize ( nPatchesX, nPatchesY );
    for ( int x = 0; x < nPatchesX; ++x ) {
        for ( int y = 0; y < nPatchesY; ++y ) {
            FileMomfbd::PatchInfo* patch = patches.ptr ( x, y );
            dataMask |= patch->parse ( file, swapNeeded, version );
            region[0] = std::min( region[0], patch->region[0] );
            region[1] = std::max( region[1], patch->region[1] );
            region[2] = std::min( region[2], patch->region[2] );
            region[3] = std::max( region[3], patch->region[3] );
        }
    }
    
    try {
        readOrThrow ( file, &nFileNames, 1, "FileMomfbd:nFileNames" );
        dataMask |= MOMFBD_NAMES;
        if ( swapNeeded ) {
            swapEndian ( nFileNames );
        }
        filenameOffset = file.tellg();
    } catch ( const std::ios_base::failure& e ) {
        nFileNames = 0;
    }


}

void FileMomfbd::read ( const std::string& filename ) {
    ifstream file ( filename );
    read ( file );
}
#include "redux/file/fileana.hpp"
void FileMomfbd::write ( std::ofstream& file, const char* data, uint8_t writeMask, int verbosity ) {

    headerSize = 0;
    file.seekp ( 0 );

    char tmp8;
    int32_t tmp32;
    vector<char> tmpStr ( 1024, 0 );
    int32_t dummySize(0);
    shared_ptr<float> dummy;
    if((phOffset<0) && writeMask&MOMFBD_MODES) dummySize = nPH*nPH;
    if((modesOffset<0) && writeMask&MOMFBD_MODES) dummySize = std::max(nModes*(nPH*nPH+1), dummySize);

    if(dummySize) {
        dummy = sharedArray<float>(dummySize);
        memset(dummy.get(),0,dummySize*sizeof(float));
    }


    // first byte = endian-flag
    tmp8 = system_is_big_endian;
    writeOrThrow ( file, &tmp8, 1, "FileMomfbd:endian" );

    // version string
    tmp32 = versionString.length() + 1;
    writeOrThrow ( file, &tmp32, 1, "FileMomfbd:version-length" );
    writeOrThrow ( file, versionString.c_str(), tmp32, "FileMomfbd:version-string" );
    
    // time string
    tmp32 = timeString.length() + 1;
    writeOrThrow ( file, &tmp32, 1, "FileMomfbd:time-length" );
    writeOrThrow ( file, timeString.c_str(), tmp32, "FileMomfbd:time-string" );

    // date string
    tmp32 = dateString.length() + 1;
    writeOrThrow ( file, &tmp32, 1, "FileMomfbd:date-length" );
    writeOrThrow ( file, dateString.c_str(), tmp32, "FileMomfbd:date-string" );

    tmp8 = (writeMask&MOMFBD_MODES);
    writeOrThrow ( file, &tmp8, 1, "FileMomfbd:hasModes" );
    if ( tmp8 ) {
        if ( version >= 20110714.0 ) {
            writeOrThrow ( file, &pix2cf, 1, "FileMomfbd:pix2cf" );
            writeOrThrow ( file, &cf2pix, 1, "FileMomfbd:cf2pix" );
        }
        writeOrThrow ( file, &nPH, 1, "FileMomfbd:nPH" );
        writeOrThrow ( file, &nModes, 1, "FileMomfbd:nModes" );
        const float* fPtr;
        if(phOffset>=0) {
            fPtr = reinterpret_cast<const float*> ( data + phOffset );
        } else fPtr = dummy.get();
        writeOrThrow ( file, fPtr, nPH * nPH, "FileMomfbd:PH-data" );
        if ( nModes ) {
            if( modesOffset>=0 ) {
                fPtr = reinterpret_cast<const float*> ( data + modesOffset );
            } else fPtr = dummy.get();
            uint8_t tmp8(1);
            for( int i=0; i<nModes; ++i ) {
                if(version >= 20120201.0) writeOrThrow ( file, &tmp8, 1, "FileMomfbd:Mode-type" );   // set mode 1 (= Pupil-mode)
                writeOrThrow ( file, fPtr, nPH * nPH, "FileMomfbd:Mode-data" );
                fPtr += nPH * nPH;
            }
        }
    }


    writeOrThrow ( file, &nChannels, 1, "FileMomfbd:nChannels" );
    if ( nChannels ) {
        // global clip
        writeOrThrow ( file, clipStartX.get(), nChannels, "FileMomfbd:clipStartX" );
        writeOrThrow ( file, clipEndX.get(), nChannels, "FileMomfbd:clipEndX" );
        writeOrThrow ( file, clipStartY.get(), nChannels, "FileMomfbd:clipStartY" );
        writeOrThrow ( file, clipEndY.get(), nChannels, "FileMomfbd:clipEndY" );
    } else cout << "nChannels = 0!!" << endl;

    writeOrThrow ( file, &nPatchesX, 1, "FileMomfbd:nPatchesX" );
    writeOrThrow ( file, &nPatchesY, 1, "FileMomfbd:nPatchesY" );
    writeOrThrow ( file, &nPoints, 1, "FileMomfbd:nPoints" );
    
    for ( int x = 0; x < nPatchesX; ++x ) {
        for ( int y = 0; y < nPatchesY; ++y ) {
            FileMomfbd::PatchInfo* patch = patches.ptr ( x, y );
            patch->write ( file, data, version, writeMask );
        }
    }

    if( (writeMask & MOMFBD_NAMES) && fileNames.size() ) {
        nFileNames = fileNames.size();
        writeOrThrow ( file, &nFileNames, 1, "FileMomfbd:nFileNames" );
        for( auto & fn:fileNames ) {
            int32_t nameLength = fn.length() + 1;
            writeOrThrow ( file, &nameLength, 1, "FileMomfbd:nameLength" );
            writeOrThrow ( file, fn.c_str(), nameLength, "FileMomfbd:FileName" );
        }
    }
}


void FileMomfbd::write ( const std::string& filename, const char* data, uint8_t writeMask, int verbosity ) {
    ofstream file ( filename, ofstream::binary );
    write ( file, data, writeMask, verbosity );
}



size_t FileMomfbd::load ( ifstream& file, char* ptr, uint8_t loadMask, int verbosity, uint8_t alignTo ) {

    size_t count(0);
    file.clear();

    // pupil & modes
    if ( ( loadMask & MOMFBD_MODES ) && nPH ) {

        // pupil
        if ( phOffset > 0 ) {
            file.seekg ( phOffset );
            float* fPtr = reinterpret_cast<float*> ( ptr+count );
            if( verbosity > 2 ) {
                cout << "FileMomfbd::loading pupil.    size = (" << nPH << "," << nPH << ")" << endl;
            }
            size_t n = readOrThrow ( file, fPtr, nPH * nPH, "MomfbdData:pupil" );
            if( n != (nPH * nPH * sizeof(float))) {
                cout << "FileMomfbd:pupil:  size mismatch: " << n << " != " << (nPH * nPH * sizeof(float)) << endl;
            }
            if ( swapNeeded ) {
                swapEndian ( fPtr, nPH * nPH );
            }
            count += nPH * nPH * sizeof(float);
        }

        // modes
        if ( nModes && (modesOffset>0) ) {
            // Load data
            file.seekg ( modesOffset );
            float* fPtr = reinterpret_cast<float*> ( ptr+count );
            if( verbosity > 2 ) {
                cout << "FileMomfbd::loading modes.    size = (" << nModes << "," << nPH << "," << nPH << ")  modesOffset = " << modesOffset << endl;
            }
            for (int i=0; i< nModes; ++i ) {
                file.seekg ( 1, ios_base::cur ); // FIXME: skip first byte for now.
                size_t n = readOrThrow ( file, fPtr, nPH * nPH, "MomfbdData:modes" );
                if( n != (nPH * nPH * sizeof(float))) {
                    cout << "FileMomfbd:modes:  size mismatch: " << n << " != " << (nPH * nPH * sizeof(float)) << endl;
                }
                fPtr += nPH * nPH;
            }
            if ( swapNeeded ) {
                swapEndian ( fPtr, nModes * nPH * nPH );
            }
            count += nModes * nPH * nPH * sizeof(float);
        }

    }

    if ( loadMask & MOMFBD_PATCH && nPatchesX > 0 && nPatchesY > 0 ) {
        const FileMomfbd::PatchInfo* patch;

        if ( verbosity > 1 ) {
            cout << "Total patches: " << nPatchesX << " x " << nPatchesY << endl;
        }

        for ( int x = 0; x < nPatchesX; ++x ) {
            for ( int y = 0; y < nPatchesY; ++y ) {
                patch = patches.ptr ( x, y );
                if ( patch ) {
                    if ( verbosity > 1 ) {
                        cout << "Loading patch (" << x << "," << y << ")   \r" << flush;
                    }
                    count += patch->load ( file, ptr+count, swapNeeded, version, loadMask, verbosity, alignTo );
                    //cout << "ptr = " << hexString(ptr) << "   ptr2 = " << hexString(ptr2) << "   diff = " << (ptrdiff_t)(ptr2-ptr0)<< endl;
                }   //  if( patch )

            }   // x-loop

        }   // y-loop

        if ( verbosity > 1 ) {
            cout << endl;
        }

    }   // if(MOMFBD_PATCH)


    if ( ( loadMask & MOMFBD_NAMES ) && (filenameOffset>0) ) {
        file.seekg ( filenameOffset );
        int32_t nameLength;
        vector<char> tmpStr ( 1024, 0 );
        int i = 0;
        while ( i < nFileNames ) {
            readOrThrow ( file, &nameLength, 1, "MomfbdData:namelength #" + to_string ( i ) );
            if ( swapNeeded ) swapEndian ( nameLength );
            tmpStr.reserve ( nameLength );
            readOrThrow ( file, & ( tmpStr[0] ), nameLength,  "MomfbdData:name #" + to_string ( i ) );
            fileNames.push_back ( string ( tmpStr.begin(), tmpStr.begin() + nameLength ) );
            i++;
        }
    } else fileNames.clear();

    return count;
}




std::shared_ptr<FileMomfbd> redux::file::readMomfbdInfo ( const std::string& filename ) {

    std::shared_ptr<FileMomfbd> hdr ( new FileMomfbd ( filename ) );
    return hdr;

}


