#include "redux/momfbd/channel.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/object.hpp"
#include "redux/momfbd/subimage.hpp"
#include "redux/momfbd/util.hpp"

#ifdef DEBUG_
#   define TRACE_THREADS
#endif

#include "redux/constants.hpp"
#include "redux/file/fileana.hpp"
#include "redux/math/functions.hpp"
#include "redux/image/cachedfile.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/descatter.hpp"
#include "redux/image/utils.hpp"
#include "redux/image/zernike.hpp"
#include "redux/logging/logger.hpp"
#include "redux/translators.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/util/projective.hpp"
#include "redux/util/trace.hpp"

#include <functional>
#include <math.h>
#include <numeric>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace bfs = boost::filesystem;
using namespace redux::file;
using namespace redux::image;
using namespace redux::logging;
using namespace redux::momfbd;
using namespace redux::util;
using namespace redux;
using namespace std;


Channel::Channel (Object& o, MomfbdJob& j, uint16_t id) : otfNormalization(0), ID (id),
    flipX(false), flipY(false), nTotalFrames(0), myObject(o), myJob(j), logger(j.logger) {

}


Channel::~Channel() {
    THREAD_MARK
    cleanup();
    THREAD_MARK
}


void Channel::parsePropertyTree (bpt::ptree& tree, redux::logging::Logger& logger) {

    ChannelCfg::parseProperties( tree, logger, myObject );       // parse using our parent-object as default-values.
    checkFlip();    // just to set flipX/Y

}



bpt::ptree Channel::getPropertyTree( bpt::ptree& tree, bool showAll ) {

    bpt::ptree node;
    ChannelCfg::getProperties( node, myObject, showAll );         // generate using our parent-object as default-values.
    tree.push_back (bpt::ptree::value_type ("channel", node));
    return node;

}


void Channel::cleanup(void) {
    
    THREAD_MARK
    imageStats.clear();
    images.clear();
    dark.clear();
    gain.clear();
    ccdResponse.clear();
    ccdScattering.clear();
    psf.clear();
    modulationMatrix.clear();
    xOffset.clear();
    yOffset.clear();
    subImages.clear();
    phi_fixed.clear();
    phi_channel.clear();
    
    THREAD_MARK
    if( !cacheFile.empty() ) {
        bfs::path tmpP(cacheFile);
        if( bfs::exists(tmpP) ) {
            try {
                bfs::remove(tmpP);
            } catch( exception& e ) {
                LOG_ERR << "Channel::cleanup: failed to remove storage: " << cacheFile << "  reason: " << e.what() << endl;
            }
        }
    }
    THREAD_MARK

}


size_t Channel::size (void) const {

    size_t sz = ChannelCfg::size();
    sz += sizeof (uint16_t) + sizeof(uint32_t) + 2;          // ID + nTotalFrames + flipX/Y
    sz += nFrames.size()*sizeof( size_t ) + sizeof( uint64_t );
    sz += sizeof( size_t );                               // frameNumbersPerFile.size()
    for( const auto& fn: frameNumbersPerFile )  sz += fn.size()*sizeof( size_t )+ sizeof( size_t );
    sz += imgSize.size();
    sz += imageStats.size() * ArrayStats().size() + sizeof (uint16_t);
    return sz;
    
}


uint64_t Channel::pack (char* ptr) const {
    using redux::util::pack;
    uint64_t count = ChannelCfg::pack (ptr);
    count += pack (ptr + count, ID);
    count += pack (ptr + count, nTotalFrames);
    count += pack (ptr + count, nFrames );
    count += pack (ptr + count, frameNumbersPerFile.size() );
    for( const auto& fn: frameNumbersPerFile )  count += pack( ptr + count, fn );
    count += pack (ptr + count, flipX );
    count += pack (ptr + count, flipY );
    count += imgSize.pack (ptr + count);
    uint16_t statSize = imageStats.size();
    count += pack (ptr + count, statSize);
    for (auto & stat : imageStats)
        count += stat->pack (ptr + count);
    if (count != size()) {
        LOG_ERR << "(" << hexString (this) << "): Packing failed, there is a size mismatch:  count = " << count << "  sz = " << size() << ende;
    }
    return count;
}


uint64_t Channel::unpack (const char* ptr, bool swap_endian) {
    using redux::util::unpack;
    uint64_t count = ChannelCfg::unpack (ptr, swap_endian);
    count += unpack (ptr + count, ID, swap_endian);
    count += unpack (ptr + count, nTotalFrames, swap_endian);
    count += unpack (ptr + count, nFrames, swap_endian );
    size_t nFpF;
    count += unpack (ptr + count, nFpF, swap_endian );
    frameNumbersPerFile.resize( nFpF );
    for( auto& fn: frameNumbersPerFile )  count += unpack( ptr + count, fn,swap_endian  );
    count += unpack (ptr + count, flipX );
    count += unpack (ptr + count, flipY );
    count += imgSize.unpack (ptr + count, swap_endian);
    uint16_t statSize;
    count += unpack (ptr + count, statSize, swap_endian);
    imageStats.resize (statSize);
    for (auto & stat : imageStats) {
        stat.reset (new ArrayStats());
        count += stat->unpack (ptr + count, swap_endian);
    }
    return count;
}


void Channel::checkFlip( void ) {
    
    flipX = flipY = false;
    if( alignMapX.size() == 16 && alignMapY.size() == 16 ) {
        if( alignMapX[4] < 0 ) flipX = true;
        if( alignMapY[1] < 0 ) flipY = true;
    } else if( alignMap.size() == 9 ) {
        if( alignMap[0] < 0 ) flipX = true;
        if( alignMap[4] < 0 ) flipY = true;
    } else if( alignClip.size() == 4) {
        if( alignClip[0] > alignClip[1] ) flipX = true;
        if( alignClip[2] > alignClip[3] ) flipY = true;
    }
    
}


bool Channel::checkCfg (void) {

    // Do we have a correct filename template ?
    if (imageTemplate.empty()) {
        LOG_ERR << "No filename template specified." << ende;
        return false;
    }
    size_t nWild = std::count (imageTemplate.begin(), imageTemplate.end(), '%');
    if (nWild > 2) {
        LOG_ERR << "Filename template contains too many wildcards: \"" << imageTemplate << "\"" << ende;
        return false;
    } else if (nWild == 1 && fileNumbers.empty()) {
        LOG_ERR << "Filename template contains wildcard and no image-numbers given (with IMAGE_NUM)" << ende;
        return false;
    } /*else if( nWild == 2 && sequenceNumber == 0 ) {
        LOG_ERR << "Filename template contains 2 wildcards and no sequence-number given (with SEQUENCE_NUM)" << ende;
        return false;
    }*/


    if( darkTemplate.empty() != gainFile.empty() ) {
        LOG_ERR << "Either BOTH or NONE of DARK_TEMPLATE/GAIN_FILE has to be specified!!!" << ende;
        return false;
    } else if (darkTemplate.empty()) {
        LOG_WARN << "No dark/gain files specified, assuming the data is pre-processed already." << ende;
        if ( !responseFile.empty() ) {
            LOG_ERR << "Detector response correction only possible when flatfielding, It will NOT be performed!!!" << ende;
        }
    } else {
        nWild = std::count (darkTemplate.begin(), darkTemplate.end(), '%');
        if (nWild > 1) {
            LOG_ERR << "Dark template contains too many wildcards: \"" << darkTemplate << "\"" << ende;
            return false;
        } else if (nWild == 1 && darkNumbers.empty()) {
            LOG_ERR << "Dark template contains wildcard and no dark-numbers given (with DARK_NUM)" << ende;
            return false;
        } else if (nWild == 0 && darkNumbers.size()) {
            //LOG_WARN << "Dark template contains no wildcard AND dark-numbers specified. Numbers will be ignored and the dark-template used as a single filename." << ende;
            darkNumbers.clear();    // TODO: fix this properly, numbers might reappear after transfer (because of inheritance)
        }
    }
    
    if( !alignClip.empty() && alignClip.size() != 4 ) {
        LOG_WARN << "ALIGN_CLIP does not seem to contain 4 integers. Whole image area will be used." << ende;
        alignClip.clear();
    }
    
    if( alignClip.size() == 4 ) {
        if( !alignClip[0] || !alignClip[1] || !alignClip[2] || !alignClip[3] ) {
            LOG_ERR << "ALIGN_CLIP values should be 1-based, i.e. first pixel is (1,1)." << ende;
            return false;
        }
        for( auto& c: alignClip ) c -= 1;           // keep the clips 0-based internally
    }
    
    checkFlip();    // just to set flipX/Y
    
    if( discard.size() > 2 ) {
        LOG_WARN << "DISCARD only uses 2 numbers, how many frames to drop at beginning/end of a cube. It will be truncated." << ende;
    }
    discard.resize(2,0);

    if( diversityValues.size() != diversityModes.size() ) {
        LOG_WARN << "Number of diversity orders (" << diversityModes.size()
                 << ") does not match number of diversity coefficients (" << diversityValues.size() << ")!" << ende;
        return false;
    }

    // Check if data directory exists (doesn't have to, it might only be visible by the manager, not necessarily on the machine where the job was submitted.
    bfs::path fn = bfs::path( boost::str(boost::format (imageTemplate) % (imageNumberOffset)) ).parent_path();
    if( !bfs::is_directory(fn) ) {    // try alternative path for data directory
        fn = bfs::path(imageDataDir) / bfs::path (boost::str(boost::format (imageTemplate) % (imageNumberOffset))).parent_path();
    }
    if( incomplete && bfs::is_directory(fn) ) {   // if directory is present, and INCOMPLETE is specified, check which image numbers exists
        for( size_t i(0); i < fileNumbers.size(); ) {
            fn = bfs::path( boost::str( boost::format(imageTemplate) % (imageNumberOffset + fileNumbers[i]) ) );
            if( !bfs::is_regular_file(fn) ) {
                fn = bfs::path(imageDataDir) / bfs::path( boost::str( boost::format(imageTemplate) % (imageNumberOffset + fileNumbers[i]) ) );
                if (!bfs::is_regular_file(fn)) {
                    //LOG_TRACE << "File not found: \"" << fn.string() << "\", removing from list of image numbers." << ende;
                    fileNumbers.erase(fileNumbers.begin() + i);
                    continue;
                }
            }
            ++i;
        }
        if( fileNumbers.empty() ) {
            LOG_FATAL << boost::format("No files found for incomplete channel with filename template \"%s\" in directory \"%s\"") % imageTemplate % imageDataDir << ende;
            return false;
        }
    }
    
    return true;

}


bool Channel::checkData( bool verbose ) {

    // Images
    if( incomplete ) {   // check if files are present
        for( size_t i(0); i < fileNumbers.size(); ) {
            bfs::path fn = bfs::path( boost::str( boost::format(imageTemplate) % (imageNumberOffset + fileNumbers[i]) ) );
            if( !bfs::is_regular_file (fn)) {
                fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % (imageNumberOffset + fileNumbers[i])));
                if( !bfs::is_regular_file(fn) ) {
                    //LOG_TRACE << "File not found: \"" << fn.string() << "\", removing from list of image numbers." << ende;
                    fileNumbers.erase(fileNumbers.begin() + i);
                    continue;
                }
            }
            ++i;
        }
        if( fileNumbers.empty() ) {
            LOG_FATAL << boost::format("No files found for incomplete channel with filename template \"%s\" in directory \"%s\"") % imageTemplate % imageDataDir << ende;
            return false;
        }
    }
    
    size_t nFiles = std::max<size_t>(1,fileNumbers.size());         // If no numbers, load template as single file
    nFrames.resize( nFiles, 1 );                                    // Assume 1 frame per file by default
    frameNumbersPerFile.resize( nFiles );                           // for storing the frameNumbers obtained from each file.

    if (fileNumbers.empty()) {         // single file
        bfs::path fn = bfs::path(imageDataDir) / bfs::path(imageTemplate);
        if (!bfs::is_regular_file(fn)) {
            LOG_ERR << boost::format ("Input-file %s not found!") % fn << ende;
            return false;
        }
        try {
            Image<float> tmp;
            //CachedFile::load( tmp, fn.string(), true );       // Only read metadata
            redux::file::readFile( fn.string(), tmp, true );       // Only read metadata
            uint8_t nDims = tmp.meta->nDims();
            vector<size_t> fNums = tmp.meta->getFrameNumbers();
            if( discard[0] ) {
                fNums.erase( fNums.begin(), fNums.begin()+discard[0] );
            }
            if( discard[1] ) {
                fNums.erase( fNums.end()-discard[1], fNums.end() );
            }
            frameNumbersPerFile[0] = fNums;
            //CachedFile::unload<float>( fn.string() );
            if( nDims == 3 ) {
                nFrames[0] = tmp.meta->dimSize(0) - (discard[0]+discard[1]);
                //imgSize = Point16( tmp.meta->dimSize(1), tmp.meta->dimSize(2) );
            } else if( nDims != 2 ) {
                LOG_ERR << boost::format ("Input-file %s not 2D or 3D!") % imageTemplate << ende;
                return false;
            }
        } catch( const std::exception& e ) {
            LOG_ERR << boost::format ("Failed to get nFrames from file %s what(): %s") % imageTemplate % e.what() << ende;
            return false;
        } catch( ... ) {
            LOG_ERR << boost::format ("Failed to get nFrames from file %s") % imageTemplate << ende;
            return false;
        }
    } else {                            // template + numbers
        for( size_t i=0; i<nFiles; ++i) {
            uint32_t number = fileNumbers[i] + imageNumberOffset;
            string thisFile = boost::str( boost::format(imageTemplate) % (number) );
            bfs::path fn = bfs::path (imageDataDir) / bfs::path( thisFile );
            if (!bfs::is_regular_file(fn)) {
                LOG_ERR << boost::format ("Input-file %s not found!") % fn << ende;
                return false;
            }
            try {
                Image<float> tmp;
                //CachedFile::load( tmp, fn.string(), true );       // Only read metadata
                redux::file::readFile( fn.string(), tmp, true );       // Only read metadata
                uint8_t nDims = tmp.meta->nDims();
                vector<size_t> fNums = tmp.meta->getFrameNumbers(number);
                if( discard[0] ) {
                    fNums.erase( fNums.begin(), fNums.begin()+discard[0] );
                }
                if( discard[1] ) {
                    fNums.erase( fNums.end()-discard[1], fNums.end() );
                }
                frameNumbersPerFile[i] = fNums;
                //CachedFile::unload<float>( fn.string() );
                if( nDims == 3 ) {
                    nFrames[i] = tmp.meta->dimSize(0) - (discard[0]+discard[1]);
                    //imgSize = Point16( tmp.meta->dimSize(1), tmp.meta->dimSize(2) );
                } else if( nDims != 2 ) {
                    LOG_ERR << boost::format ("Input-file %s not 2D or 3D!") % thisFile << ende;
                    return false;
                }
            } catch( const std::exception& e ) {
                LOG_ERR << boost::format ("Failed to read nFrames from file %s what(): %s") % thisFile % e.what() << ende;
                return false;
            } catch( ... ) {
                LOG_ERR << boost::format ("Failed to read nFrames from file %s") % thisFile << ende;
                return false;
            }
            
        }
    }
    
    nTotalFrames = std::accumulate( nFrames.begin(), nFrames.end(), 0 );
    if( waveFrontList.size() != nFiles ) {
        waveFrontList.resize( nFiles, imageNumberOffset );
    }

    if( waveFrontList.size() == nFiles ) {
        std::transform( waveFrontList.begin(), waveFrontList.end(), fileNumbers.begin(), waveFrontList.begin(), std::plus<uint32_t>() );
    }

    if( nFiles != nTotalFrames ) {
        vector<uint32_t> tmpV;
        for( size_t i=0; i<nFiles; ++i ) {
            for( const auto& fn: frameNumbersPerFile[i] ) {
                tmpV.push_back( fn );
            }
        }
        std::swap( waveFrontList, tmpV );
    }
    string wfStr = redux::util::uIntsToString( waveFrontList );
    if( verbose ) LOG_DETAIL << "Channel " << myObject.ID << ":" << ID << " waveFronts: " << wfStr << ende;

    // Dark(s)
    if (!darkTemplate.empty()) {
        size_t nWild = std::count(darkTemplate.begin(), darkTemplate.end(), '%');
        if (nWild == 0 || darkNumbers.empty()) {          // single file, DARK_NUM will be ignored if no wildcard in the template
            if (!bfs::is_regular_file(darkTemplate)) {
                bfs::path fn = bfs::path (imageDataDir) / bfs::path (darkTemplate);
                if (!bfs::is_regular_file(fn)) {
                    logAndThrow("Dark-file " + darkTemplate + " not found!");
                } else darkTemplate = fn.string();
            }
        } else {                            // template
            for (auto & number : darkNumbers) {
                bfs::path fn = bfs::path(boost::str(boost::format(darkTemplate) % number));
                if (!bfs::is_regular_file(fn)) {
                    fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (darkTemplate) % number));
                    if (!bfs::is_regular_file(fn)) {
                        logAndThrow("Dark-file " + (boost::format(darkTemplate) % number).str() + " not found!");
                    } else darkTemplate = (bfs::path(imageDataDir) / bfs::path(darkTemplate)).string();      // found in imageDataDir, prepend it to darkTemplate
                }
            }
        }
    }

    // Gain
    if (!gainFile.empty()) {
        if (!bfs::is_regular_file(gainFile)) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (gainFile);
            if (!bfs::is_regular_file(fn)) {
                logAndThrow("Gain-file " + gainFile + " not found!");
            } else gainFile = fn.string();
        }
    }

    if (!responseFile.empty()) {
        if (!bfs::is_regular_file(responseFile)) {
            bfs::path fn = bfs::path(imageDataDir) / bfs::path(responseFile);
            if (!bfs::is_regular_file(fn)) {
                logAndThrow("Response-file " + responseFile + " not found!");
            } else responseFile = fn.string();
        }
    }

    if (!backgainFile.empty()) {
        if (!bfs::is_regular_file(backgainFile)) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (backgainFile);
            if (!bfs::is_regular_file(fn)) {
                logAndThrow("Backgain-file " + backgainFile + " not found!");
            } else backgainFile = fn.string();
        }
    }

    if (!psfFile.empty()) {
        if (!bfs::is_regular_file(psfFile)) {
            bfs::path fn = bfs::path(imageDataDir) / bfs::path(psfFile);
            if (!bfs::is_regular_file(fn)) {
                logAndThrow("PSF-file " + psfFile + " not found!");
            } else psfFile = fn.string();
        }
    }

    if (!mmFile.empty()) {
        if (!bfs::is_regular_file(mmFile)) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (mmFile);
            if (!bfs::is_regular_file(fn)) {
                logAndThrow("Modulation-matrix file " + mmFile + " not found!");
            } else mmFile = fn.string();
        }
    }

    if (!xOffsetFile.empty()) {
        if (!bfs::is_regular_file(xOffsetFile)) {
            bfs::path fn = bfs::path(imageDataDir) / bfs::path(xOffsetFile);
            if (!bfs::is_regular_file(fn)) {
                logAndThrow("Offset-file " + xOffsetFile + " not found!");
            } else xOffsetFile = fn.string();
        }
    }

    if (!yOffsetFile.empty()) {
        if (!bfs::is_regular_file(yOffsetFile)) {
            bfs::path fn = bfs::path(imageDataDir) / bfs::path(yOffsetFile);
            if (!bfs::is_regular_file(fn)) {
                logAndThrow("Offset-file " + yOffsetFile + " not found!");
            } else yOffsetFile = fn.string();
        }
    }

    return true;
}



void Channel::initChannel (void) {

    ModeInfo info( myJob.klMinMode, myJob.klMaxMode, 0, myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, myJob.klCutoff );
    
    if( !diversityModeFile.empty() ) {      // A file with modes was supplied
        
        info = ModeInfo( diversityModeFile, myObject.pupilPixels, divModeFileNormalize );

        const shared_ptr<ModeSet>& modes = myJob.globalData->get( info );
        lock_guard<mutex> lock( modes->mtx );
        if( modes->empty() && bfs::is_regular_file(diversityModeFile) ) {   // load file
            LOG_DEBUG << "initChannel(" << idString() << "):   Empty diversity ModeSet: loading \"" << diversityModeFile << "\"" << ende;
            modes->load( diversityModeFile, myObject.pupilPixels );
            if( !modes->empty() ) {
                if( myObject.pupil ) modes->getNorms( *(myObject.pupil) );
                else modes->getNorms();
                if( divModeFileNormalize ) modes->normalize();
                LOG_DEBUG << "initChannel(" << idString() << "): normalize: " << (divModeFileNormalize?"Yes":"No") << ende;
                LOG_DEBUG << "initChannel(" << idString() << "): " << printArray( modes->norms,"measured norms" ) << ende;
            }
        }
        if( modes->empty() ) {
            LOG_ERR << "initChannel(" << idString() << "):  Could not load diversityModeFile=\"" << diversityModeFile << "\"" << ende;
        }
        return;
    }

    
    for( size_t i(0); i < diversityModes.size(); ++i ) {
        
        uint16_t modeNumber = diversityModes[i].mode;
        int modeType = diversityModes[i].type;
        
        ModeInfo mi2 = info;
        if( modeNumber == 2 || modeNumber == 3 || modeType == ZERNIKE ) {
            mi2.firstMode = mi2.lastMode = 0;
        }
        mi2.modeNumber = modeNumber;
        const shared_ptr<ModeSet>& dmode = myJob.globalData->get(mi2);
        unique_lock<mutex> lock( dmode->mtx);
        if( dmode->empty() ) {    // this (single-mode-) set was inserted, so it is not generated yet.
            if( modeType == ZERNIKE) {
                LOG_DEBUG << "initChannel(" << idString() << "):   Empty diversity mode: generating Zernike mode." << ende;
                dmode->generate( myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, diversityModes, Zernike::NORMALIZE );
            } else {
            LOG_DEBUG << "initChannel(" << idString() << "):   Empty diversity mode: generating KL mode." << ende;
                dmode->generate( myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, myJob.klMinMode, myJob.klMaxMode, diversityModes, myJob.klCutoff, Zernike::NORMALIZE );
            }
            if( dmode->nDimensions() != 3 || dmode->dimSize(1) != myObject.pupilPixels || dmode->dimSize(2) != myObject.pupilPixels ) {    // mismatch
                LOG_ERR << "Generated diversity-mode does not match. This should NOT happen!!" << ende;
            } else {
                if( myObject.pupil && myObject.pupil->valid() ) dmode->getNorms( *(myObject.pupil) );
                else dmode->getNorms();
                LOG_DETAIL << "Generated diversity-mode: " << mi2 << printArray( dmode->norms,"\n\t\t  norm" ) << ende;
                dmode->normalize();
            }
        }  
    }
    
}


void Channel::loadCalib( boost::asio::io_service& service ) {     // load through cache-functionality, so the same data is recycled for other objects

    //LOG_TRACE << "Channel::loadCalib() << ende";
    // TODO: absolute/relative paths
    // TODO: cache files and just fetch shared_ptr

    if( !darkTemplate.empty() ) {        // needs to be read synchronously because of adding/normalization
        Image<float> tmp;
        size_t nDarkFrames(0);
        size_t nWild = std::count (darkTemplate.begin(), darkTemplate.end(), '%');
        if (nWild == 0 || darkNumbers.empty()) {
            CachedFile::load( tmp, darkTemplate );
            nDarkFrames = tmp.meta->getNumberOfFrames();
        } else {
            Image<float> tmp2;
            for (size_t di = 0; di < darkNumbers.size(); ++di) {
                bfs::path fn = bfs::path (boost::str (boost::format (darkTemplate) % darkNumbers[di]));
                if (!di) {
                    CachedFile::load( tmp, fn.string() );
                    nDarkFrames = tmp.meta->getNumberOfFrames();
                } else {
                    CachedFile::load( tmp2, fn.string() );
                    tmp += tmp2;
                    nDarkFrames += tmp2.meta->getNumberOfFrames();
                }
            }
            
            //   NOTE: no normalization here, modify metadata instead. TODO: += operator for meta to add frames??

        }
        dark = tmp.copy();
        if (nDarkFrames) dark *= 1.0/nDarkFrames;
    }


    if( !gainFile.empty() ) {
        CachedFile::load<float>( gain, gainFile );
        gainMask = rdx_get_shared<uint8_t>(imgSize.y*imgSize.x);
        make_mask( gain.get(), gainMask.get(), imgSize.y, imgSize.x, 0, 8, true, true ); // filter away larger features than ~8 pixels
    }


    if( !responseFile.empty() ) {
        CachedFile::load<float>( ccdResponse, responseFile );
    }

    if (!backgainFile.empty()) {
        CachedFile::load<float>( ccdScattering, backgainFile );
    }

    if( !psfFile.empty() ) {
        CachedFile::load<float>( psf, psfFile );
    }

    if( !mmFile.empty() ) {
        CachedFile::load<float>( modulationMatrix, mmFile );
    }

    if( !xOffsetFile.empty() ) {
        CachedFile::load<int16_t>( xOffset, xOffsetFile );
        if( alignClip.size() == 4 ) {
            Point clipSize( abs(alignClip[3]-alignClip[2])+1, abs(alignClip[1]-alignClip[0])+1 );
            if( clipSize.x != xOffset.dimSize(1) ||
                clipSize.y != xOffset.dimSize(0) ) {
                LOG_ERR << "Size of offset file: " << xOffsetFile << " does not match the align_clip." << ende;
                xOffset.clear();
            }
        }
    }

    if( !yOffsetFile.empty() ) {
        CachedFile::load<int16_t>( yOffset, yOffsetFile );
        if( alignClip.size() == 4 ) {
            Point clipSize( abs(alignClip[3]-alignClip[2])+1, abs(alignClip[1]-alignClip[0])+1 );
            if( clipSize.x != yOffset.dimSize(1) ||
                clipSize.y != yOffset.dimSize(0) ) {
                LOG_ERR << "Size of offset file: " << yOffsetFile << " does not match the align_clip." << ende;
                yOffset.clear();
            }
        }
    }

    
}


void Channel::loadData( boost::asio::io_service& service, redux::util::Array<PatchData::Ptr>& patches ) {

    size_t nFiles = std::max<size_t>( 1, fileNumbers.size() );       // If no numbers, load template as single file

    startT = bpx::pos_infin;
    endT = bpx::neg_infin;
    
    if( nTotalFrames == 0 ) {
        throw logic_error("No input images for channel " + to_string(myObject.ID) + ":" + to_string(ID));
    }
    
    // Prepare needed storage
    if( !(myJob.runFlags&RF_NOSWAP) ) {    // unless swap is deactivated for this job
        cacheFile = myJob.cachePath + "images_" + to_string(myObject.ID) + "_" + to_string(ID);
    }
    
    maybeLoadImages();
    imageStats.resize( nTotalFrames );
    
    bool saveFFData = (myObject.saveMask&SF_SAVE_FFDATA) || (myJob.runFlags&RF_FLATFIELD);
    if( saveFFData ) {
        myObject.progWatch.increaseTarget( nFiles );    // add saving of calibrated datafiles.
    }
    
    loadCalib(service);
    
    progWatch.set(patches.dimSize(0)*patches.dimSize(1)+nFiles);
    progWatch.setHandler( [this](){
        if( !(myJob.runFlags&RF_NOSWAP) ) {
            images.clear();
        }
    });
    for( unsigned int y=0; y<patches.dimSize(0); ++y ) {
        for( unsigned int x=0; x<patches.dimSize(1); ++x ) {
            service.post( [this,&patches,y,x]() {
                auto oData = patches(y,x)->getObjectData(myObject.ID);
                if( !oData ) throw runtime_error("patches(y,x)->getObjectData() returned a null pointer !");
                auto chData = oData->channels[ID];
                adjustCutout( *chData, patches(y,x) );
                ++progWatch;
            } );
        }
    }
    
    size_t nPreviousFrames(0);
    for( size_t i=0; i<nFiles; ++i ) {
        service.post( [&,i,nPreviousFrames,saveFFData](){
            try {
                loadFile(i,nPreviousFrames);
                if( imgSize.y < 1 || imgSize.x < 1 ) throw logic_error("Image size is zero.");
                size_t nF = nFrames[i];
                for( size_t j=0; j<nF; ++j ) {
                    preprocessImage( nPreviousFrames+j );
                    //service.post(std::bind( &Channel::preprocessImage, this, nPreviousFrames+j) );
                }
                if( saveFFData ) {
                    bfs::path fn;
                    try {
                        fn = bfs::path (myJob.info.outputDir) / bfs::path (boost::str (boost::format (imageTemplate) % fileNumbers[i])).filename();
                        bfs::path ext = fn.extension();
                        if(ext.string().empty() || ext.string().length() > 5 ) {     // we assume the filename does not have a proper extenstion, add a temporary dummy
                            fn = bfs::path( fn.string() + ".ext" );
                        }
                        fn.replace_extension(".cor.f0");
                        LOG_DEBUG << boost::format("Saving dark/flat corrected data (%d:%d:%d) as %s.") % myObject.ID % ID % i % fn.string() << ende;
                        Image<float> view( images, nPreviousFrames, nPreviousFrames+nF-1, 0, imgSize.y-1, 0, imgSize.x-1 );
                        redux::file::Ana::write( fn.string(), view.copy() );   // TODO: other formats
                    } catch ( const std::exception& e ) {
                        LOG_ERR << "Failed to save corrected file: " << fn << "  reason: " << e.what() << ende;
                    }
                    ++myObject.progWatch;
                }
                ++progWatch;
                return;
            } catch ( const std::exception& e ) {
                LOG_ERR << "Failed to load/preprocess file. reason: " << e.what() << ende;
            } catch ( ... ) {
                LOG_ERR << "Failed to load/preprocess file for unknown reason." << ende;
            }
            Job::moveTo( &myJob, Job::JSTATE_ERR );
            myJob.progWatch.clear();
            progWatch.clear();
            myJob.updateProgressString();
        });
        nPreviousFrames += nFrames[i];
    }

}


void Channel::unloadCalib(void) {               // unload what was accessed through the cache, this should be called when all objects are done pre-processing.

    dark.clear();
    gain.clear();
    ccdResponse.clear();
    ccdScattering.clear();
    psf.clear();
    modulationMatrix.clear();
    xOffset.clear();
    yOffset.clear();

    if (!darkTemplate.empty()) {
        size_t nWild = std::count (darkTemplate.begin(), darkTemplate.end(), '%');
        if (nWild == 0 || darkNumbers.empty()) {
            CachedFile::unload<float>(darkTemplate);
        } else {
            for (size_t di = 0; di < darkNumbers.size(); ++di) {
                bfs::path fn = bfs::path (boost::str (boost::format (darkTemplate) % darkNumbers[di]));
                CachedFile::unload<float>(fn.string());
            }
        }
    }
    if( !gainFile.empty() ) CachedFile::unload<float>(gainFile);
    if( !responseFile.empty() ) CachedFile::unload<float>(responseFile);
    if( !backgainFile.empty() ) CachedFile::unload<float>(backgainFile);
    if( !psfFile.empty() ) CachedFile::unload<float>(psfFile);
    if( !mmFile.empty() ) CachedFile::unload<float>(mmFile);
    if( !xOffsetFile.empty() ) CachedFile::unload<int16_t>(xOffsetFile);
    if( !yOffsetFile.empty() ) CachedFile::unload<int16_t>(yOffsetFile);

}


double Channel::getMaxMean(void) {
    double mm = std::numeric_limits<double>::lowest();
    for ( const ArrayStats::Ptr &stat : imageStats ) {
        if( stat ) mm = std::max( mm, stat->mean );
    }
    return mm;
}


double Channel::getMaxMedian(void) {
    double mm = 0.0;
    for ( const ArrayStats::Ptr &stat : imageStats ) {
        if( stat ) mm = std::max( mm, stat->median );
    }
    return mm;
}


double Channel::getMedianMedian(void) {
    vector<double> medians;
    for ( const ArrayStats::Ptr &stat : imageStats ) {
        if( stat ) medians.push_back( stat->median );
    }
    double mm = 0.0;
    size_t nM = medians.size();
    if( nM ) {
        size_t nHalf = (nM>>1) - !(nM&1);
        auto mid = medians.begin() + nHalf;
        std::nth_element( medians.begin(), mid, medians.end());
        mm = *mid;
        if( !(nM&1) ) {
            std::nth_element( mid, mid+1, medians.end());
            mm += *(mid+1);
            mm /= 2;
        }
    }
    return mm;
}


void Channel::getFileNames( std::vector<std::string>& files, const std::vector<uint32_t>& wf ) const {

    set<uint32_t> wfSet( wf.begin(), wf.end() );
    size_t nFiles = std::min( nFrames.size(), fileNumbers.size() );
    for( size_t i=0; i<nFiles; ++i ) {
        bool inList(false);
        for( const auto& fn: frameNumbersPerFile[i] ) {
            if( wfSet.count( fn ) ) {
                inList = true;
                break;
            }
        }
        if( inList ) {
            bfs::path fn = bfs::path(imageDataDir)/bfs::path(boost::str(boost::format(imageTemplate) % fileNumbers[i]));
            files.push_back(fn.string());
        }
    }
    
}


void Channel::initProcessing( Solver& solver ) {
    
    initPhiFixed();

    subImages.resize(nTotalFrames);
    for (uint16_t i=0; i < nTotalFrames; ++i) {
        if( !subImages[i] ) subImages[i].reset( new SubImage(myObject, *this, solver.window, solver.noiseWindow) );
        auto ret = solver.wavefronts.emplace( waveFrontList[i], nullptr );
        if( !ret.first->second ) ret.first->second.reset( new WaveFront(waveFrontList[i]) );
        ret.first->second->addImage( subImages[i] );
    }
    
}


void Channel::initPatch (ChannelData& cd) {

    size_t nImages = cd.images.dimSize();
    if (nTotalFrames != nImages) {
        string msg = "Number of images in stack does not match the number of total frames.  " +to_string(nTotalFrames)
        + "!=" + printArray(cd.images.dimensions(),"dims");
        LOG_ERR << msg << ende;
        throw logic_error(msg);
    }

    uint16_t patchSize = myObject.patchSize;
    Point16 dataSize( cd.images.dimSize(1), cd.images.dimSize(2) );
    Region16 dataRegion( dataSize-1 );
    Region16 patchRegion( patchSize-1, patchSize-1 );
    patchRegion += cd.patchStart;
    
    if( cd.patchStart+patchSize > dataSize ) {
        LOG_ERR << "initPatch ch=" << myObject.ID << ":" << ID << ": patchStart(" << cd.patchStart
                 << ")+patchSize(" << patchSize << ") > dataSize(" << dataSize << "). Something is very strange here..." << ende;
    }
    
    size_t dataPixels = dataSize.x * dataSize.y;
    if( myObject.normalizeTo > 0.0 ) {
        if( (imageStats.size() == nImages) && (cd.images.nDimensions() == 3) ) {
            float* dataPtr = cd.images.get();
            for( uint16_t i=0; i<nImages; ++i ) {
                double scale = myObject.normalizeTo;
                if( myJob.normType == NORM_OBJ_MAX_MEAN ) {
                    scale /= imageStats[i]->mean;
                } else if( (myJob.normType == NORM_OBJ_MAX_MEDIAN) || (myJob.normType == NORM_OBJ_MEDIAN_MEDIAN) ) {
                    scale /= imageStats[i]->median;
                } else if( myJob.normType == NORM_NONE ) scale = 1.0;
                transform( dataPtr, dataPtr+ dataPixels, dataPtr, [scale]( const float &a ){ return a*scale; } );
                dataPtr += dataPixels;
            }
        }
    }

    std::shared_ptr<float*> arrayPtr = cd.images.reshape( nTotalFrames*dataSize.y, dataSize.x );
    float** imgPtr = arrayPtr.get();
    for( size_t i=0; i<nTotalFrames; ++i) {
        if( flipX ) redux::util::reverseX(imgPtr, dataSize.y, dataSize.x);
        if( flipY ) redux::util::reverseY(imgPtr, dataSize.y, dataSize.x);
        imgPtr += dataSize.y;
    }
    
    if( flipX ) {
        cd.residualOffset.x *= -1;
    }
    if( flipY ) {
        cd.residualOffset.y *= -1;
    }
    
    patchRegion = dataRegion.getAsGlobal( patchRegion );
    

    for (uint16_t i=0; i < nImages; ++i) {
        subImages[i]->wrap( cd.images, i, i, patchRegion.first.y, patchRegion.last.y, patchRegion.first.x, patchRegion.last.x );
        subImages[i]->setPatchInfo( i, cd.patchStart, cd.residualOffset, patchSize, dataSize.x, myObject.pupilPixels, myJob.nModes );
        subImages[i]->stats.getStats( cd.images.ptr(i,0,0), dataPixels, ST_VALUES|ST_RMS );
    }

    phi_fixed.copy( phi_channel );
    
    double* phiPtr = phi_channel.get();
    size_t pupilSize2 = myObject.pupilPixels*myObject.pupilPixels;
    
    PointD residualAlpha = myObject.modes->shiftToAlpha (cd.residualOffset);

    LOG_DEBUG << "initPatch ch=" << myObject.ID << ":" << ID << ":  patchRegion=" << patchRegion << "  dataRegion=" << dataRegion
             << "\n    patchStart=" << cd.patchStart << "  exactPatchPosition=" << cd.exactPatchPosition
             << "  cutoutRegion=" << cd.cutoutRegion << "  residualOffset=" << cd.residualOffset
             << "  residualAlpha=" << residualAlpha << ende;


    int32_t mIndex = myObject.modes->tiltMode.x;
    if( mIndex >= 0 && fabs(cd.residualOffset.x) > 0 ) {
        const double* modePtr = myObject.modes->modePointers[mIndex];
        double res = -residualAlpha.x;   // positive coefficient shifts image to the left
        transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
            [res](const double& p, const double& m) {
                return p + res*m;
            });
    }
    
    mIndex = myObject.modes->tiltMode.y;
    if( mIndex >= 0 && fabs(cd.residualOffset.y) > 0 ) {
        const double* modePtr = myObject.modes->modePointers[mIndex];
        double res = -residualAlpha.y;   // positive coefficient shifts image downwards
        transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
            [res](const double& p, const double& m) {
                return p + res*m;
            });
    }

    otfNormalization = 1;       // measure the normalization (because calculating it as sqrt(1.0 / object.pupil.area / otfSize2) is a bit inaccurate)
    if( !subImages.empty() && subImages[0] ) {
        subImages[0]->zeroPhi();
        subImages[0]->calcPFOTF();
        complex_t* otfPtr = subImages[0]->OTF.get();
        double otfMax(0);
        pupilSize2 *= 4;    // OTF-size is 2 x pupilSize
        for( size_t i(0); i<pupilSize2; ++i ){
            otfMax = max(otfMax,abs(otfPtr[i]));
        }
        otfNormalization = 1.0/sqrt(otfMax);
    }
}


void Channel::initPhiFixed(void) {

    phi_fixed.resize( myObject.pupilPixels, myObject.pupilPixels );
    phi_fixed.zero();
    double* phiPtr = phi_fixed.get();
    size_t pupilSize2 = myObject.pupilPixels*myObject.pupilPixels;
    
    if( diversityModes.empty() ) {
        return;
    }
    
    ModeInfo minfo( myJob.klMinMode, myJob.klMaxMode, 0, myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, myJob.klCutoff );
    double scale4 = util::def2cf( myJob.telescopeD / myObject.telescopeF );     // conversion factor from m to radians, only used for focus-term.
    
    if( ! diversityModeFile.empty() ) {     // using mode-file
        minfo = ModeInfo( diversityModeFile, myObject.pupilPixels, divModeFileNormalize );
    
        const shared_ptr<ModeSet>& modes = myJob.globalData->get( minfo );
        lock_guard<mutex> lock( modes->mtx );
        if( !modes->empty() ) {
            for( size_t i(0); i < diversityModes.size(); ++i ) {
                uint16_t modeNumber = diversityModes[i].mode;
                if( modeNumber >= modes->modePointers.size() ) {
                    LOG_WARN << "initPhiFixed(" << idString() << "): modeNumber = " << modeNumber
                            << " >= nModes in modefile \"" << diversityModeFile << "\"\n\t\tThis Diversity mode will be ignored!!!" << ende;
                    continue;
                }
                double alpha = diversityValues[i].coefficient;
                if( diversityValues[i].physical ) {
                    alpha *= scale4/myObject.wavelength;
                }
                const double* modePtr = modes->modePointers[static_cast<size_t>(modeNumber)];
                transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
                    [alpha](const double& p, const double& m) {
                        return p + alpha*m;
                    });
            }
        } else {
            LOG_ERR << "initPhiFixed: DIV_MODE_FILE=\"" << diversityModeFile
                    << "\", but ModeSet is empty!! This should NOT happen!!" << ende;
        }

        return;
    }
    
    for( size_t i(0); i < diversityModes.size(); ++i ) {
        
        ModeInfo mi2 = minfo;
        mi2.modeNumber = diversityModes[i].mode;
        int modeType = diversityModes[i].type;
        double alpha = diversityValues[i].coefficient;
        if( diversityValues[i].physical ) {
            alpha *= scale4/myObject.wavelength;
        }

        if( mi2.modeNumber == 2 || mi2.modeNumber == 3 || modeType == ZERNIKE ) {   // generated tilts always Zernike
            mi2.firstMode = mi2.lastMode = 0;
        }
        
        const shared_ptr<ModeSet>& ms = myJob.globalData->get(mi2);

        if( ms->empty() ) {    // generate
            if( modeType == ZERNIKE ) {
                ms->generate( myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, myJob.modeList, Zernike::NORMALIZE );
            } else {
                ms->generate( myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, myJob.klMinMode, myJob.klMaxMode, myJob.modeList, myJob.klCutoff, Zernike::NORMALIZE );
            }
            if( ms->nDimensions() != 3 || ms->dimSize(1) != myObject.pupilPixels || ms->dimSize(2) != myObject.pupilPixels ) {    // mismatch
                LOG_ERR << "Generated diversity-mode does not match. This should NOT happen!!" << ende;
            } else {
                LOG_DETAIL << "Generated diversity-mode: " << mi2 << printArray( ms->norms,"\n\t\t  norm" ) << ende;
            }
        }        
        auto it = std::find( ms->modeList.begin(), ms->modeList.end(), mi2.modeNumber );
        if( it != ms->modeList.end() ) {
            const double* modePtr = ms->modePointers[static_cast<size_t>(it-ms->modeList.begin())];
            transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
                [alpha](const double& p, const double& m) {
                    return p + alpha*m;
                });
        }

    }

}


void Channel::addAllFT (redux::util::Array<double>& ftsum) {
    for (auto& subimage : subImages) {
        subimage->addFT (ftsum);
    }
}


void Channel::addTimeStamps( const bpx::ptime& newStart, const bpx::ptime& newEnd ) {

    unique_lock<mutex> lock(mtx);
    
    if(startT.is_special()) startT = newStart;
    else startT = std::min( startT, newStart );
    if(endT.is_special()) endT = newEnd;
    else endT = std::max( endT, newEnd );

}


void Channel::loadFile( size_t fileIndex, size_t offset ) {

    bfs::path fn = bfs::path(imageDataDir) / bfs::path( );
    Image<float> tmpImg;
    
    try {
        
        if( fileIndex >= fileNumbers.size() ) {
            string msg = "File-index is out of range. (" + to_string(fileIndex) + "/" + to_string(fileNumbers.size()) + ").";
            throw logic_error(msg);
        }
        if( !fileNumbers.empty() ) {
            fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( imageTemplate ) % fileNumbers[fileIndex] ) );
        }
    
        string imStr = to_string(fileIndex);
        redux::file::readFile( fn.string(), tmpImg, false );
        size_t nF = nFrames[fileIndex];
        Image<float> view( images, offset, offset+nF-1, 0, imgSize.y-1, 0, imgSize.x-1 );
        if( nF > 1 ) {
            Image<float> tmpView( tmpImg, discard[0], discard[0]+nF-1, 0, imgSize.y-1, 0, imgSize.x-1 );
            view.assign( reinterpret_cast<redux::util::Array<float>&>(tmpView) );
            imStr = to_string(offset)+"-"+to_string(offset+nF-1);
        } else {
            view.assign( reinterpret_cast<redux::util::Array<float>&>(tmpImg) );           
        }
        if( tmpImg.meta ) {
            addTimeStamps( tmpImg.meta->getStartTime(), tmpImg.meta->getEndTime() );
        }
        LOG_DETAIL << boost::format ("Loaded file "+imageTemplate+"  (%d:%d:%s)") % fileNumbers[fileIndex] % myObject.ID % ID % imStr << ende;
    } catch ( const std::exception& ) {
        throw;
    } catch ( ... ) {
        string msg = "Failed to load file " + fn.string() + " for unknown reason.";
        throw runtime_error(msg);
    }
    

}


void Channel::preprocessImage( size_t i ) {

    try {
        
        Image<float> view( images, i, i, 0, imgSize.y-1, 0, imgSize.x-1 );
        Array<double> tmpImg( imgSize.y, imgSize.x );       // local temporary with double-precision
        tmpImg = view.copy<double>();
        
        LOG_DEBUG << boost::format ("Preprocessing image (%d:%d:%d)  %s") % myObject.ID % ID % i % printArray(tmpImg.dimensions(),"dims") << ende;
        
    //         bfs::path fn = bfs::path (boost::str (boost::format (imageTemplate) % i));
    //         fn = bfs::path(fn.leaf().string() + ".inp");
    //         redux::file::Ana::write( fn.string(), tmpImg );   // TODO: other formats

        /*
        // Michiel's method for detecting bitshifted Sarnoff images.    TODO make this an SST-specific "filter" to be applied on raw data
        if (imgMean > 4 * avgMean) {
            LOG_WARN << boost::format ("Image bit shift detected for image %s (mean > 4*avgMean). adjust factor=0.625 (keep your fingers crossed)!") % fn << ende;
            tmpImg *= 0.625;
            modified = true;
        } else if (imgMean < 0.25 * avgMean) {
            LOG_WARN << boost::format ("Image bit shift detected for image %s (mean < 0.25*avgMean). adjust factor=16 (keep your fingers crossed)!") % fn << ende;
            tmpImg *= 16;
            modified = true;
        }*/

        if (dark.valid() && gain.valid()) {
            if (! tmpImg.sameSize (dark)) {
                LOG_ERR << boost::format ("Dimensions of dark (%s) does not match this image (%s), skipping flatfielding !!")
                        % printArray (dark.dimensions(), "") % printArray (tmpImg.dimensions(), "") << ende;
                return;
            }
            if (! tmpImg.sameSize (gain)) {
                LOG_ERR << boost::format ("Dimensions of gain (%s) does not match this image (%s), skipping flatfielding !!")
                        % printArray (gain.dimensions(), "") % printArray (tmpImg.dimensions(), "") << ende;
                return;
            }
            if (ccdResponse.valid() && !tmpImg.sameSize (ccdResponse)) {
                LOG_WARN << boost::format ("Dimensions of ccd-response (%s) does not match this image (%s), will not be used !!")
                        % printArray (ccdResponse.dimensions(), "") % printArray (tmpImg.dimensions(), "") << ende;
                ccdResponse.resize();
            }
            
            double n;
            if(dark.meta && ((n=dark.meta->getNumberOfFrames()) > 1)) {
                tmpImg.subtract(dark,1.0/n);
            } else {
                tmpImg -= dark;
            }
            
            if (ccdResponse.valid()) {   // correct for the detector response (this should not contain the gain correction and must be done before descattering)
                tmpImg *= ccdResponse;
            }

            if (ccdScattering.valid() && psf.valid()) {           // apply backscatter correction
                if (tmpImg.sameSize (ccdScattering) && tmpImg.sameSize (psf)) {
                    LOG_DEBUG << boost::format("Applying correction for CCD transparency for image (%d:%d:%d)") % myObject.ID % ID % i << ende;
                    redux::image::descatter( tmpImg, ccdScattering, psf );
                } else {
                    LOG_ERR << boost::format ("Dimensions of ccdScattering (%s) or psf (%s) does not match this image (%s), skipping flatfielding !!")
                            % printArray (ccdScattering.dimensions(), "") % printArray (psf.dimensions(), "") % printArray (tmpImg.dimensions(), "") << ende;
                }
            }

            tmpImg *= gain;
    //         fn = bfs::path (boost::str (boost::format (imageTemplate) % i));
    //         fn = bfs::path(fn.leaf().string() + ".dg");
    //         redux::file::Ana::write( fn.string(), tmpImg );   // TODO: other formats

            namespace sp = std::placeholders;
            size_t sy = tmpImg.dimSize(0);
            size_t sx = tmpImg.dimSize(1);

            shared_ptr<double*> array = tmpImg.reshape(sy,sx);
            shared_ptr<uint8_t*> mask2D;
            if( gainMask ) {
                mask2D = reshapeArray( gainMask.get(), sy, sx );
            }
            double** arrayPtr = array.get();
            switch (myJob.fillpixMethod) {
                case FPM_HORINT: {
                    //LOG_DETAIL << "Filling bad pixels using horizontal interpolation." << ende;
                    function<double (size_t, size_t) > func = bind( horizontalInterpolation<double>, arrayPtr, sy, sx, sp::_1, sp::_2);
                    fillPixels (arrayPtr, sy, sx, func, std::bind( std::less_equal<double>(), sp::_1, myJob.badPixelThreshold), mask2D.get());
                    break;
                }
                case FPM_MEDIAN: {
                    // TODO: median method
                    break;
                }
                case FPM_INVDISTWEIGHT:       // inverse distance weighting is the default method, so fall through
                default: {
                    //LOG_DETAIL << "Filling bad pixels using inverse distance weighted average." << ende;
                    function<double (size_t, size_t) > func = bind( inverseDistanceWeight<double>, arrayPtr, sy, sx, sp::_1, sp::_2);
                    fillPixels (arrayPtr, sy, sx, func, std::bind( std::less_equal<double>(), sp::_1, myJob.badPixelThreshold), mask2D.get());
                }
            }

            // Fill larger features that the mask will exclude. This will fill e.g. black borders.
            // TBD: Should this be skipped and force the user to be stricter with the clip/ROI instead?
            function<double (size_t, size_t) > func = bind( inverseDistanceWeight<double>, arrayPtr, sy, sx, sp::_1, sp::_2);
            fillPixels (arrayPtr, sy, sx, func, std::bind( std::less_equal<double>(), sp::_1, myJob.badPixelThreshold));
            
            // FIXME: This is a hack to create truncated values as the old code!!
            //for(size_t i=0; i<sy*sx; ++i) arrayPtr[0][i] = (int)arrayPtr[0][i];
            
        }

        view.assign(tmpImg);                            // copy back to image
        
        
        if( borderClip ) {
            std::vector<int64_t> first, last;
            const std::vector<size_t>& dims = tmpImg.dimensions();
            size_t nDims = dims.size();
            for( size_t i = 0; i < nDims; ++i ) {
                if( dims[i] == 1 ) {            // if dimSize == 1 we don't clip it.
                    first.push_back( tmpImg.first()[i] );
                    last.push_back( tmpImg.last()[i] );
                }
                else {
                    if( 2 * borderClip > dims[i] ) {    // all data clipped, nothing to do
                        return;
                    }
                    first.push_back( tmpImg.first()[i] + borderClip );
                    last.push_back( tmpImg.last()[i] - borderClip );
                }
            }
            tmpImg = Array<double>( tmpImg, first, last ).copy();
        }
        
        imageStats[i].reset( new ArrayStats() );
        imageStats[i]->getStats( tmpImg );    // get stats for corrected data
        
        // compute median value
        double* tmpPtr = tmpImg.ptr();
        size_t nEl = tmpImg.nElements();
        size_t halfEl = (nEl>>1) - !(nEl&1);
        double* mid = tmpPtr + halfEl;
        std::nth_element( tmpPtr, mid, tmpPtr+nEl);
        imageStats[i]->median = *mid;
        if( !(nEl&1) ) {
            std::nth_element( mid, mid+1, tmpPtr+nEl );
            imageStats[i]->median += *(mid+1);
            imageStats[i]->median /= 2;
        }
    } catch ( const std::exception& e ) {
        LOG_ERR << boost::format("Failed to preprocess image (%d:%d:%d)  %s") % myObject.ID % ID % i % printArray(images.dimensions(),"dims")
                << ".  reason: " << e.what() << ende;
    } catch ( ... ) {
        LOG_ERR << boost::format("Failed to preprocess image (%d:%d:%d)  %s") % myObject.ID % ID % i % printArray(images.dimensions(),"dims")
                << " for unknown reason."  << ende;
    }

    ++myJob.progWatch;          // this will trigger unloading the calibration data after all objects/images are done.
    ++myObject.progWatch;       // this will trigger calculating max_mean and storing patch-data to disk
        

}


void Channel::maybeLoadImages( void ) {
    
    lock_guard<mutex> lock(mtx);
    if( images.nElements() ) return;        // already allocated/loaded
    
    if( !(myJob.runFlags&RF_NOSWAP) ) {     // should we mmap?
        if( bfs::exists( bfs::path(cacheFile) ) ) {
            LOG_TRACE << "opening: Images: " << cacheFile << ende;
            images.openMmap( cacheFile, nTotalFrames, imgSize.y, imgSize.x );
        } else {
            LOG_TRACE << "creating: Images: " << cacheFile << ende;
            images.createMmap( cacheFile, nTotalFrames, imgSize.y, imgSize.x );
        }
    } else {
        images.resize( nTotalFrames, imgSize.y, imgSize.x );
    }
    
}


void Channel::getStorage( ChannelData& chData ) {

    maybeLoadImages();
    chData.images.wrap(reinterpret_cast<redux::util::Array<float>&>(images), 0, nTotalFrames-1, chData.cutoutRegion.first.y, chData.cutoutRegion.last.y, chData.cutoutRegion.first.x, chData.cutoutRegion.last.x);
 
}


uint32_t Channel::nImages( const std::vector<uint32_t>& wf ) {
    set<uint32_t> wfSet( wf.begin(), wf.end() );
    uint32_t nIms(0);
    for( uint32_t& w: waveFrontList ) {
        if( wfSet.count( w) ) {
            nIms++;
        }
    }
    return nIms;
}


uint32_t Channel::nImages(void) {
    
    if( nTotalFrames == 0 ) {
        // just to allow a rough estimate of needed diskspace before the preprocessing starts.
        if( nFrames.empty() ) return fileNumbers.size();
        nTotalFrames = std::accumulate( nFrames.begin(), nFrames.end(), 0 );
    }
    
    return nTotalFrames;
    
}


PointF Channel::getOffsetAt( const Point16& pos, size_t sz ) const {

    PointF ret(0,0);
    if( !xOffset.valid() && !yOffset.valid() ) {     //  No offset files
        return ret;
    }
    
    sz /= 2;    // We will only use half-size below
    RegionI roi( pos.y-sz, pos.x-sz, pos.y+sz, pos.x+sz );
    roi.first.max( PointI() );      // restrict to non-negative positions
    
    Image<int16_t> patch;
    ArrayStats stats;
    if( xOffset.valid() ) {
        size_t ySz = xOffset.dimSize(0);
        size_t xSz = xOffset.dimSize(1);
        if( false && sz ) {     // FIXME: stats with stride seems buggy !!
            roi.last.min( PointI(ySz-1, xSz-1) );    // restrict to array size
            if( roi.span() != 0 ) {
                patch.wrap( xOffset, roi.first.y, roi.first.x, roi.last.y, roi.last.x );
                stats.getMinMaxMean( patch );
                ret.x = stats.mean/100.0;
            }
        } else if( (pos.y < ySz) && (pos.x < xSz) ) {
            ret.x = xOffset( pos.y, pos.x )/100.0;
        }
    }
    if( yOffset.valid() ) {
        size_t ySz = yOffset.dimSize(0);
        size_t xSz = yOffset.dimSize(1);
        if( false && sz ) {     // FIXME: stats with stride seems buggy !!
            roi.last.min( PointI(ySz-1, xSz-1) );    // restrict to array size
            if( roi.span() != 0 ) {
                patch.wrap( yOffset, roi.first.y, roi.first.x, roi.last.y, roi.last.x );
                stats.getMinMaxMean( patch );
                ret.x = stats.mean/100.0;
            }
        } else if( (pos.y < ySz) && (pos.x < xSz) ) {
            ret.y = yOffset( pos.y, pos.x )/100.0;
        }
    }
    
    return ret;
    
}


void Channel::adjustCutout( ChannelData& chData, const PatchData::Ptr& patch ) const {
    
    if( !patch ) throw runtime_error("Channel::adjustCutout() called with null patch!");

    if( imgSize == 0 ) {
        throw runtime_error("No valid imgSize when adjusting cutout, this should not happen!!!");
    }
    
    // Patch (midpoint) position. This is specified in the clipped/reference frame, if align_clip is specified.
    Region16 clipRegion( 0, 0, imgSize.y-1, imgSize.x-1 );     // defines the allowed placement for the patch
    
    if( alignClip.size() == 4 ) {   // if there is align_clip, the coordinates are given relative to that area
        clipRegion.first = Point16( alignClip[2], alignClip[0] );
        clipRegion.last  = Point16( alignClip[3], alignClip[1] );
    }   // TODO use job::roi if it exists !
    
    // Local copy of patch-position (which is in "clipped&flipped" coordinates, if applicable)
    Point16 patchPos = patch->position;
    
    // Doubles for processing.  Patch position in local/camera coordinates. 
    PointD localPos = clipRegion.getAsGlobal( patch->position );
    PointD localOffset;
    clipRegion.normalize();    // make sure first < last
    
    const uint16_t patchSize = myObject.patchSize;
    const uint16_t halfPatch = patchSize/2;

    chData.channelOffset = 0;
    
    if( alignMapX.size() == 16 && alignMapY.size() == 16 ) {
        localPos = pointWarp( alignMapX, alignMapY, localPos );
    } else if( alignMap.size() == 9 ) {
        ProjectiveMap map( alignMap );
        // the align-map is always in global coordinates
        localPos = map * (localPos-0.5) + 0.5;
        // localPos is now in global (camera) coordinates.
    } else {                                // old style alignment with clips & offsetfiles
        // getOffsetAt() must be called with "clipped&flipped" coordinates, and returns the offset, also in "clipped&flipped"
        localOffset = getOffsetAt( localPos-clipRegion.first, myObject.patchSize );
        chData.channelOffset = localOffset.round();
        localPos += localOffset + PointI(flipY,flipX);
        // localPos is now global (camera) coordinates, channelOffset holds the integer "offset" from files
    }
    chData.exactPatchPosition = localPos;
    chData.cutoutPosition = localPos.round(); // + PointI(flipY,flipX);
    chData.residualOffset = localPos.remainder();
    
    RegionI desiredPatch( myObject.patchSize-1, myObject.patchSize-1 );         // patch region
    desiredPatch += chData.cutoutPosition - halfPatch;                          // ... now centered on cutoutPosition, this is our patch
    desiredPatch.expand( myObject.maxLocalShift );

    clipRegion.first += halfPatch;
    clipRegion.last -= (halfPatch-1);
    RegionI actualPatch = desiredPatch;

    bool alreadyWarned(false);
    if( !clipRegion.isInside(chData.cutoutPosition) ) {
        LOG_WARN << "Patch " << patchPos << " (mapped to: " << chData.cutoutPosition
                << "), does not lie within allowed patch placement (" << clipRegion << ") in channel " << myObject.ID << ":" << ID
                << ". Difference=" << clipRegion.diff(localPos) << "  This *will* cause severe artifacts!!!" << ende;
        alreadyWarned = true;
    }
    
    // N.B. align-clip should be ignored for the restriction of max_local_shift
    clipRegion.first = halfPatch + myObject.maxLocalShift;
    clipRegion.last = imgSize - halfPatch - myObject.maxLocalShift;
    if( !clipRegion.isInside(chData.cutoutPosition) ) {
        PointD diff = clipRegion.diff(localPos);
        if( !alreadyWarned ) {
            LOG_NOTICE << "Patch " << patchPos << " (mapped to: " << chData.cutoutPosition
                      << "), lies closer to the border than MAX_LOCAL_SHIFT in channel " << myObject.ID
                      << ":" << ID << ", this *might* cause artifacts if image motion is large."
                      << "\n\tPlace patch within (" << clipRegion << ") to avoid this potential issue."
                      << "\n\tDifference=" << diff << "." << ende;
        }
        actualPatch.shrinkSigned( diff.round() );
    }
    chData.cutoutRegion = actualPatch;
    chData.patchStart = (chData.cutoutPosition - halfPatch) - actualPatch.first;
    string flipStr = (flipX||flipY)?"  flip: ":"";
    if( flipX ) {
        flipStr += "X";
        chData.patchStart.x = actualPatch.last.x - (chData.cutoutPosition.x + halfPatch - 1);
    }
    if( flipY ) {
        flipStr += flipX?"/Y":"Y";
        chData.patchStart.y = actualPatch.last.y - (chData.cutoutPosition.y + halfPatch - 1);
    }
    string actStr = "";
    if( chData.cutoutRegion != desiredPatch ) actStr = "  actual="+(string)chData.cutoutRegion;
    if( chData.cutoutRegion != desiredPatch ) actStr = "  actual="+(string)chData.cutoutRegion;
    LOG_DETAIL<< "AdjustCutout ch=" << myObject.ID << ":" << ID << ": specifiedPosition" << patchPos
              << "  mappedPosition=" << chData.exactPatchPosition
              << "  cutoutPosition=" << chData.cutoutPosition
              << "\n    residualPosition=" << chData.residualOffset
              << "  desiredPatch=" << desiredPatch << actStr
              << "  localOffset=" << chData.channelOffset <<"  actualPatch=" << actualPatch
              << "   start=" << chData.patchStart << flipStr << ende;
         
}


void Channel::adjustCutouts( Array<PatchData::Ptr>& patches ) {
    
    if( patches.nDimensions() == 2 ) {
        size_t nPatchesY = patches.dimSize(0);
        size_t nPatchesX = patches.dimSize(1);
        for( unsigned int py=0; py<nPatchesY; ++py ) {
            for( unsigned int px=0; px<nPatchesX; ++px ) {
                PatchData::Ptr& patch( patches(py,px) );
                if( !patch ) {
                    throw runtime_error("Channel::adjustCutout() patch pointer is NULL! (" + to_string(py) + "," + to_string(px) + ")");
                }
                ChannelData::Ptr chData = patch->getChannelData( myObject.ID, ID );
                if( !chData ) {
                    string msg = "patch->getChannelData() returned a null pointer!";
                    msg += "\n  objID = " + to_string(myObject.ID) + "  chanID = " + to_string(ID)
                        + " patch.position = " + (string)patch->position;
                    throw runtime_error( msg );
                }
                adjustCutout( *chData, patch );
            }
        }
    }
    
}


Point16 Channel::getImageSize( bool force ) {

    if( force || (imgSize == 0) ) {
//         if( alignClip.size() == 4 ) {   // we have align-clip, but no mapping => reference channel.
//             imgSize = Point16(abs(alignClip[3]-alignClip[2])+1, abs(alignClip[1]-alignClip[0])+1);
//         } else {                        //  No align-map or align-clip, get full image size.
            checkCfg();
            string thisFile;
            if( !fileNumbers.empty() ) {
                thisFile = boost::str(boost::format (imageTemplate) % fileNumbers[0]);
            } else {
                thisFile = imageTemplate;
            }
            bfs::path fn = bfs::path(imageDataDir) / bfs::path(thisFile);
            
            if( bfs::is_regular_file(fn) ) {
                try {
                    Image<float> tmp;
                    redux::file::readFile( fn.string(), tmp, true );       // Only read metadata
                    uint8_t nDims = tmp.meta->nDims();
                    if( nDims == 3 ) {
                        imgSize = Point16( tmp.meta->dimSize(1), tmp.meta->dimSize(2) );
                    } else if( nDims == 2 ) {
                        imgSize = Point16( tmp.meta->dimSize(0), tmp.meta->dimSize(1) );
                    } else LOG_ERR << "Image " << fn << " is not 2D or 3D." << ende;
                } catch( const std::exception& e ) {
                    LOG_ERR << boost::format ("Failed to get imgSize from file %s what(): %s") % thisFile % e.what() << ende;
                } catch( ... ) {
                    LOG_ERR << boost::format ("Failed to get imgSize from file %s") % thisFile << ende;
                }
            }
//        }
            if( imgSize.x && imgSize.y ) { 
                uint16_t old_clip = borderClip;
                borderClip = min<uint16_t>( borderClip, imgSize.x/8 );     // ensure we are not clipping too much.
                borderClip = min<uint16_t>( borderClip, imgSize.y/8 );
                if( borderClip != old_clip ) {
                    LOG_WARN << "BORDER_CLIP has been reduced from " << old_clip << " to " << borderClip << " to avoid clipping too much of the images." << ende;
                }
            }

    }

    return imgSize;
 
}
            

void Channel::logAndThrow( string msg ) {
    msg = to_string(myObject.ID)+":"+to_string(ID)+": "+msg;
    LOG_ERR << msg << ende;
    throw job_error(msg);
    
}


string Channel::idString( void ) const {
    return to_string(myObject.ID) + ":" + to_string(ID);
}


void Channel::dump (std::string tag) {
    
    tag += ":"+to_string(ID);
    
    Ana::write (tag + "_phi_fixed.f0", phi_fixed);
    Ana::write (tag + "_phi_channel.f0", phi_channel);
    
    uint16_t patchSize = myObject.patchSize;
    uint16_t otfSize = myObject.pupilPixels<<1;
    size_t blockSize = patchSize*patchSize;
    Array<float> tmpF( subImages.size(), patchSize, patchSize );
    Array<double> tmpD( patchSize, patchSize );
    Array<complex_t> tmpC( subImages.size(), patchSize, patchSize );
    FourierTransform tmpFT( patchSize, patchSize );
    Array<float> statArr( subImages.size(), 4 );
    Array<int16_t> shiftArr( subImages.size(), 2 );
    ArrayStats s;

    if( subImages.size() && subImages[0] ) {
        Array<float> wrap(reinterpret_cast<Array<float>&>(*subImages[0]));
        wrap.resetLimits();     // all subimages share a datablock, reset the limits to the "full" data and write.
        Ana::write( tag + "_data.f0", wrap );
    }
    
    float* fPtr = tmpF.get();
    for( shared_ptr<SubImage>& im: subImages ) {
        im->copyTo<float>(fPtr);
        fPtr += blockSize;
    }
    Ana::write( tag + "_imgs.f0", tmpF );
    
    fPtr = tmpF.get();
    complex_t* cPtr = tmpC.get();
    int idx(0);
    for( shared_ptr<SubImage>& im: subImages ) {
        im->getWindowedImg( tmpD, s, true );
        tmpFT.init( tmpD.get(), patchSize, patchSize, FULLCOMPLEX );
        FourierTransform::reorder(tmpFT);
        tmpD.copyTo<float>(fPtr);
        tmpFT.copyTo<complex_t>(cPtr);
        const vector<int64_t>& first = im->first();
        shiftArr(idx,0) = first[1];
        shiftArr(idx,1) = first[2];
        statArr(idx,0) = s.min;
        statArr(idx,1) = s.max;
        statArr(idx,2) = s.mean;
        statArr(idx++,3) = s.stddev;
        fPtr += blockSize;
        cPtr += blockSize;
    }
    Ana::write( tag + "_wimgs.f0", tmpF );
    Ana::write( tag + "_fts.f0", tmpC );
    Ana::write( tag + "_stat.f0", statArr );
    Ana::write( tag + "_shift.f0", shiftArr );
    tmpC.resize();      // free some memory
    
    tmpD.resize( otfSize, otfSize );
    tmpF.resize( subImages.size(), otfSize, otfSize );
    blockSize = otfSize*otfSize;
    fPtr = tmpF.get();
    for( shared_ptr<SubImage>& im: subImages ) {
        im->getPSF( tmpD.get() );
        tmpD.copyTo<float>(fPtr);
        fPtr += blockSize;
    }
    Ana::write( tag + "_psfs.f0", tmpF );
    tmpF.resize();      // free some memory

    
}

