#include "redux/file/filemomfbd.hpp"

#include <iostream>
#include <cstdlib>     // atof
//#include <stdlib.h>     // atof

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


namespace {

    size_t parsePatch ( ifstream& file, FileMomfbd::PatchInfo* p, const bool& swapNeeded, const float& version ) {

        int64_t tmpSize, patchSize(0);

        p->offset = file.tellg();

        patchSize += readOrThrow ( file, p->region, 4, "PatchInfo:region" );
        patchSize += readOrThrow ( file, &p->nChannels, 1, "PatchInfo:nChannels" );
        if ( swapNeeded ) {
            swapEndian ( &p->region, 4 );
            swapEndian ( p->nChannels );
        }

        if ( p->region[0] > p->region[1] ) swap ( p->region[0], p->region[1] );
        if ( p->region[2] > p->region[3] ) swap ( p->region[2], p->region[3] );

        p->nPixelsX = p->region[1] - p->region[0] + 1;
        p->nPixelsY = p->region[3] - p->region[2] + 1;

        p->nim.reset ( new int32_t [p->nChannels], [] ( int32_t * p ) {
            delete[] p;
        } );
        p->dx.reset ( new int32_t [p->nChannels], [] ( int32_t * p ) {
            delete[] p;
        } );
        p->dy.reset ( new int32_t [p->nChannels], [] ( int32_t * p ) {
            delete[] p;
        } );

        patchSize += readOrThrow ( file, p->nim.get(), p->nChannels, "PatchInfo:nim" );
        patchSize += readOrThrow ( file, p->dx.get(), p->nChannels, "PatchInfo:dx" );
        patchSize += readOrThrow ( file, p->dy.get(), p->nChannels, "PatchInfo:dy" );

        if ( swapNeeded ) {
            swapEndian ( p->nim.get(), p->nChannels );
            swapEndian ( p->dx.get(), p->nChannels );
            swapEndian ( p->dy.get(), p->nChannels );
        }

        char hasImage;
        patchSize += readOrThrow ( file, &hasImage, 1, "PatchInfo:hasImage" );

        if ( hasImage ) {
            p->imgPos = file.tellg();
            tmpSize = p->nPixelsX * p->nPixelsY * sizeof ( float );
            patchSize += tmpSize;
            file.seekg ( tmpSize, ios_base::cur );
        }

        p->npsf = 0;
        patchSize += readOrThrow ( file, &p->npsf, 1, "PatchInfo:npsf" );
        if ( swapNeeded ) {
            swapEndian ( p->npsf );
        }

        if ( p->npsf ) {
            p->psfPos = file.tellg();
            tmpSize = p->npsf * p->nPixelsX * p->nPixelsY * sizeof ( float );
            patchSize += tmpSize;
            file.seekg ( tmpSize, ios_base::cur );
        }

        p->nobj = 0;
        patchSize += readOrThrow ( file, &p->nobj, 1, "PatchInfo:nobj" );
        if ( swapNeeded ) swapEndian ( p->nobj );
        if ( p->nobj ) {
            p->objPos = file.tellg();
            tmpSize = p->nobj * p->nPixelsX * p->nPixelsY * sizeof ( float );
            patchSize += tmpSize;
            file.seekg ( tmpSize, ios_base::cur );
        }

        p->nres = 0;
        patchSize += readOrThrow ( file, &p->nres, 1, "PatchInfo:nres" );
        if ( swapNeeded ) swapEndian ( p->nres );
        if ( p->nres ) {
            p->resPos = file.tellg();
            tmpSize = p->nres * p->nPixelsX * p->nPixelsY * sizeof ( float );
            patchSize += tmpSize;
            file.seekg ( tmpSize );
        }

        p->nalpha = p->nm = 0;
        patchSize += readOrThrow ( file, &p->nalpha, 1, "PatchInfo:nalpha" );
        if ( swapNeeded ) swapEndian ( p->nalpha );
        if ( p->nalpha ) {
            patchSize += readOrThrow ( file, &p->nm, 1, "PatchInfo:nm" );
            if ( swapNeeded ) swapEndian ( p->nm );
            p->alphaPos = file.tellg();
            tmpSize = p->nalpha * p->nm * sizeof ( float );
            patchSize += tmpSize;
            file.seekg ( tmpSize, ios_base::cur );
        }

        p->ndiv = 0;
        if ( version >= 20100726.0 ) {  // check if this is a version that can have div info...
            patchSize += readOrThrow ( file, &p->ndiv, 1, "PatchInfo:ndiv" );
            if ( swapNeeded ) swapEndian ( p->ndiv );
        }

        p->nphx = p->nPixelsX / 2;
        p->nphy = p->nPixelsY / 2;
        if ( p->ndiv ) {
            p->diversityPos = file.tellg();
            tmpSize = p->ndiv * p->nphy *  p->nphx * sizeof ( float );
            patchSize += tmpSize;
            file.seekg ( tmpSize, ios_base::cur );
        }

        if ( patchSize != ( file.tellg() - p->offset ) ) {
            throw ios_base::failure ( "Failed to read Momfbd Patch. Size mismatch: "+to_string(patchSize)+" <-> " + to_string(file.tellg() - p->offset) );
        }

        return patchSize;
    }

}   // anon namespace


FileMomfbd::FileMomfbd ( void ) : version ( 0 ), startX ( 0 ), endX ( 0 ), startY ( 0 ), endY ( 0 ),
    nChannels ( 0 ), nFileNames ( 0 ), nPH ( 0 ), nModes ( 0 ), phOffset ( 0 ),
    modesOffset ( 0 ), filenameOffset ( 0 ), patchDataSize ( 0 ), headerSize ( 0 ) {

}

FileMomfbd::FileMomfbd ( const std::string& filename ) : version ( 0 ), startX ( 0 ), endX ( 0 ), startY ( 0 ), endY ( 0 ),
    nChannels ( 0 ), nFileNames ( 0 ), nPH ( 0 ), nModes ( 0 ), phOffset ( 0 ),
    modesOffset ( 0 ), filenameOffset ( 0 ), patchDataSize ( 0 ), headerSize ( 0 ), swapNeeded ( false ) {

    read ( filename );
}

void FileMomfbd::read ( std::ifstream& file ) {

    headerSize = 0;
    file.seekg ( 0 );

    /*size_t rawSize = sizeof( struct raw_header );  // = 512 by construction
    memset( &m_Header, 0, rawSize );

    file.read( reinterpret_cast<char*>( &m_Header ), rawSize );
    if( !file.good() ) {
        throw ios_base::failure( "Failed to read ANA header" );
    }

    hdrSize = file.tellg();
    bool hasMagic = ( m_Header.synch_pattern == MAGIC_ANA );
    bool hasReversedMagic = ( m_Header.synch_pattern == MAGIC_ANAR );

    */


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
        readOrThrow ( file, &nPH, 1, "FileMomfbd:nPH" );
        readOrThrow ( file, &nModes, 1, "FileMomfbd:nModes" );

        if ( swapNeeded ) {
            swapEndian ( nPH );
            swapEndian ( nModes );
        }

        phOffset = file.tellg();
        file.seekg ( nPH * nPH * sizeof ( float ), ios_base::cur );
        modesOffset = file.tellg();
        file.seekg ( nModes * nPH * nPH * sizeof ( float ), ios_base::cur );

    } else nModes = nPH = 0;

    readOrThrow ( file, &nChannels, 1, "FileMomfbd:nChannels" );
    if ( swapNeeded ) {
        swapEndian ( nChannels );
    }

    // global clip
    clipStartX.reset ( new int16_t[ nChannels ] );
    clipEndX.reset ( new int16_t[ nChannels ] );
    clipStartY.reset ( new int16_t[ nChannels ] );
    clipEndY.reset ( new int16_t[ nChannels ] );

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

    readOrThrow ( file, &nPatchesX, 1, "FileMomfbd:nPatchesX" );
    readOrThrow ( file, &nPatchesY, 1, "FileMomfbd:nPatchesY" );
    readOrThrow ( file, &nPoints, 1, "FileMomfbd:nPoints" );

    if ( swapNeeded ) {
        swapEndian ( nPatchesX );
        swapEndian ( nPatchesY );
        swapEndian ( nPoints );
    }

    patches.resize ( nPatchesX, nPatchesY );
    for ( int x = 0; x < nPatchesX; ++x ) {
        for ( int y = 0; y < nPatchesY; ++y ) {
            FileMomfbd::PatchInfo* ptr = patches.ptr ( x, y );
            //cout << "ptr = " << hexString( ptr ) << endl;
            parsePatch ( file, ptr, swapNeeded, version );
        }
    }

    try {
        readOrThrow ( file, &nFileNames, 1, "FileMomfbd:nFileNames" );
        if ( swapNeeded ) {
            swapEndian ( nFileNames );
        }
        filenameOffset = file.tellg();
    } catch ( ... ) {
        nFileNames = 0;
    }


}

void FileMomfbd::read ( const std::string& filename ) {
    ifstream file ( filename );
    read ( file );
}

void FileMomfbd::write ( ofstream& file ) {

    throw invalid_argument( "FileMomfbd::write not implemented yet" );

}



std::shared_ptr<FileMomfbd> redux::file::readMomfbdInfo ( const std::string& filename ) {

    std::shared_ptr<FileMomfbd> hdr ( new FileMomfbd ( filename ) );
    return hdr;

}

