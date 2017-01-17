#include "redux/momfbd/channel.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/object.hpp"
#include "redux/momfbd/subimage.hpp"
#include "redux/momfbd/util.hpp"

#include "redux/constants.hpp"
#include "redux/file/fileana.hpp"
#include "redux/math/functions.hpp"
#include "redux/image/cachedfile.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/descatter.hpp"
#include "redux/image/utils.hpp"
#include "redux/logging/logger.hpp"
#include "redux/translators.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/util/projective.hpp"

#include <functional>
#include <math.h>
#include <numeric>
#include <string>

#include <boost/algorithm/string.hpp>
//#include <boost/range/algorithm.hpp>
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


namespace {
    
    double def2cf( double pd_defocus, double telescope_r ) { // defocus distance in meters
        static const double tmp = M_PI/(8.0 * sqrt(3.0));
        return pd_defocus * telescope_r * telescope_r * tmp;
    }
   
}


Channel::Channel (Object& o, MomfbdJob& j, uint16_t id) : ID (id), nTotalFrames(0), myObject(o),
    myJob(j), logger(j.logger) {

    setLogChannel(myJob.getLogChannel());
    
}


Channel::~Channel() {
    cleanup();
}


void Channel::parsePropertyTree (bpt::ptree& tree) {

    ChannelCfg::parseProperties (tree, myObject);       // parse using our parent-object as default-values.

}



bpt::ptree Channel::getPropertyTree (bpt::ptree& tree) {

    bpt::ptree node;
    ChannelCfg::getProperties (node, myObject);         // generate using our parent-object as default-values.
    tree.push_back (bpt::ptree::value_type ("channel", node));
    return node;

}


void Channel::cleanup(void) {
    
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

}


size_t Channel::size (void) const {

    size_t sz = ChannelCfg::size();
    sz += sizeof (uint16_t) + sizeof(uint32_t);          // ID + nTotalFrames
    sz += imgSize.size();
    sz += imageStats.size() * ArrayStats::size() + sizeof (uint16_t);
    return sz;
}


uint64_t Channel::pack (char* ptr) const {
    using redux::util::pack;
    uint64_t count = ChannelCfg::pack (ptr);
    count += pack (ptr + count, ID);
    count += pack (ptr + count, nTotalFrames);
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
            LOG_WARN << "Dark template contains no wildcard AND dark-numbers specified. Numbers will be ignored and the dark-template used as a single filename." << ende;
            darkNumbers.clear();    // TODO: fix this properly, numbers might reappear after transfer (because of inheritance)
        }
    }
    
    if( !alignClip.empty() && alignClip.size() != 4 ) {
        LOG_WARN << "ALIGN_CLIP does not seem to contain 4 integers. Whole image area will be used." << ende;
        alignClip.clear();
    }
    
    if( alignClip.size() == 4 ) {
        if( !alignClip[0] || !alignClip[1] || !alignClip[2] || !alignClip[3] ) {
            LOG_WARN << "ALIGN_CLIP values should be 1-based, i.e. first pixel is (1,1)." << ende;
            alignClip.clear();
        }
    }
    
    if( discard.size() > 2 ) {
        LOG_WARN << "DISCARD only uses 2 numbers, how many frames to drop at beginning/end of a cube. It will be truncated." << ende;
    }
    discard.resize(2,0);

    if( diversity.size() == diversityModes.size() ) {
        for( unsigned int i=0; i<diversity.size(); ++i ) {
            if( diversityModes[i] == 4 && physicalDefocusDistance ) {   // focus term, convert from physical length (including mm/cm) to coefficient
                diversity[i] = def2cf( diversity[i], myJob.telescopeD / myObject.telescopeF );
            }
        }
    } else {
        LOG_WARN << "Number of diversity orders does not match number of diversity coefficients!" << ende;
    }

    
    return true;

}


bool Channel::checkData (void) {

    // Images
    if ( incomplete ) {   // check if files are present
        for (size_t i (0); i < fileNumbers.size();) {
            bfs::path fn = bfs::path (boost::str (boost::format (imageTemplate) % (imageNumberOffset + fileNumbers[i])));
            if (!bfs::is_regular_file (fn)) {
                fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % (imageNumberOffset + fileNumbers[i])));
                if (!bfs::is_regular_file(fn)) {
                    //LOG_TRACE << "File not found: \"" << fn.string() << "\", removing from list of image numbers." << ende;
                    fileNumbers.erase(fileNumbers.begin() + i);
                    continue;
                }
            }
            ++i;
        }
        if (fileNumbers.empty()) {
            LOG_FATAL << boost::format ("No files found for incomplete object with filename template \"%s\" in directory \"%s\"") % imageTemplate % imageDataDir << ende;
            return false;
        }
    }
    
    size_t nFiles = std::max<size_t>(1,fileNumbers.size());         // If no numbers, load template as single file
    nFrames.resize( nFiles, 1 );                                    // Assume 1 frame per file by default

    if (fileNumbers.empty()) {         // single file
        bfs::path fn = bfs::path(imageDataDir) / bfs::path(imageTemplate);
        if (!bfs::is_regular_file(fn)) {
            LOG_ERR << boost::format ("Input-file %s not found!") % fn << ende;
            return false;
        }
        Image<float> tmp;
        CachedFile::load( tmp, fn.string(), true );       // Only read metadata
        uint8_t nDims = tmp.meta->nDims();
        //CachedFile::unload<float>( fn.string() );
        if( nDims == 3 ) {
            nFrames[0] = tmp.meta->dimSize(0) - (discard[0]+discard[1]);
            //imgSize = Point16( tmp.meta->dimSize(1), tmp.meta->dimSize(2) );
        } else if( nDims != 2 ) {
            LOG_ERR << boost::format ("Input-file %s not 2D or 3D!") % imageTemplate << ende;
            return false;
        }
    } else {                            // template + numbers
        for( size_t i=0; i<fileNumbers.size(); ++i) {
            uint32_t number = fileNumbers[i];
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % (imageNumberOffset + number)));
            if (!bfs::is_regular_file(fn)) {
                LOG_ERR << boost::format ("Input-file %s not found!") % boost::str (boost::format (imageTemplate) % (imageNumberOffset + number)) << ende;
                return false;
            }
            Image<float> tmp;
            CachedFile::load( tmp, fn.string(), true );       // Only read metadata
            uint8_t nDims = tmp.meta->nDims();
            //CachedFile::unload<float>( fn.string() );
            if( nDims == 3 ) {
                nFrames[i] = tmp.meta->dimSize(0) - (discard[0]+discard[1]);
                //imgSize = Point16( tmp.meta->dimSize(1), tmp.meta->dimSize(2) );
            } else if( nDims != 2 ) {
                LOG_ERR << boost::format ("Input-file %s not 2D or 3D!") % boost::str (boost::format (imageTemplate) % (imageNumberOffset + number)) << ende;
                return false;
            }
            
        }
    }
    
    nTotalFrames = std::accumulate( nFrames.begin(), nFrames.end(), 0 );
    if( waveFronts.size() != nFiles ) {
        waveFronts.resize( nFiles, imageNumberOffset );
    }

    if( waveFronts.size() == fileNumbers.size() ) {
        std::transform( waveFronts.begin(), waveFronts.end(), fileNumbers.begin(), waveFronts.begin(), std::plus<uint32_t>() );
    }

    if( nFiles != nTotalFrames ) {
        vector<uint32_t> tmpV;
        for( size_t i=0; i<nFiles; ++i ) {
            for( size_t j=0; j<nFrames[i]; ++j ) {
                tmpV.push_back( waveFronts[i]+j+discard[0] );
            }
        }
        std::swap( waveFronts, tmpV );
    }
    string wfStr = redux::util::uIntsToString( waveFronts );
    LOG_DEBUG << "Channel " << myObject.ID << ":" << ID << " contains waveFronts: " << wfStr << ende;
    myJob.info.progress[1] += std::max(fileNumbers.size(),1UL)*3;  // load, process, split & store (TODO more accurate progress reporting)

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


//#define INDEX_THRESHOLD  0
#define INDEX_THRESHOLD  1E-12      // cutoff to get rid of some fft noise
// TBD: this should be a parameter somewhere...
void Channel::initCache (void) {


    ModeInfo mi(myJob.klMinMode, myJob.klMaxMode, 0, myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, myJob.klCutoff);
    for (unsigned int i=0; i < diversityModes.size(); ++i) {
        uint16_t modeNumber = diversityModes[i];
        ModeInfo mi2 = mi;
        if (modeNumber == 2 || modeNumber == 3 || diversityTypes[i] == ZERNIKE) {
            mi2.firstMode = mi2.lastMode = 0;
        }
        ModeSet& ret = myJob.globalData->get(mi2);
        unique_lock<mutex> lock(ret.mtx);
        if( ret.empty() ) {    // this set was inserted, so it is not generated yet.
            if(diversityTypes[i] == ZERNIKE) {
                ret.generate( myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, diversityModes );
            } else {
                ret.generate( myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, myJob.klMinMode, myJob.klMaxMode, diversityModes, myJob.klCutoff );
            }
            if( ret.nDimensions() != 3 || ret.dimSize(1) != myObject.pupilPixels || ret.dimSize(2) != myObject.pupilPixels ) {    // mismatch
                LOG_ERR << "Generated ModeSet does not match. This should NOT happen!!" << ende;
            } else {
                LOG_DEBUG << "Generated Modeset with " << ret.dimSize(0) << " modes. (" << myObject.pupilPixels << "x" << myObject.pupilPixels
                << "  radius=" << myObject.pupilRadiusInPixels << ")" << ende;
                ret.getNorms( myObject.pupil );
            }
        }  
    }
    
}


void Channel::loadCalib( boost::asio::io_service& service ) {     // load through cache-functionality, so the same data is recycled for other objects

    //LOG_TRACE << "Channel::loadCalib() << ende";
    // TODO: absolute/relative paths
    // TODO: cache files and just fetch shared_ptr

    if( !darkTemplate.empty() ) {        // needs to be read synchronously because of adding/normalization
        progWatch.increaseTarget( std::max<int>(darkNumbers.size(),1) );
        service.post( [this](){
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
            ++progWatch;
        });
    }


    if( !gainFile.empty() ) {
        progWatch.increaseTarget( 1 );
        service.post([this](){
            CachedFile::load<float>( gain, gainFile );
            gainMask.reset( new uint8_t[imgSize.y*imgSize.x], []( uint8_t* p ){ delete[] p; } );
            make_mask( gain.get(), gainMask.get(), imgSize.y, imgSize.x, 0, 8, true, false ); // filter away larger features than ~8 pixels
            ++progWatch;
        });
    }


    if( !responseFile.empty() ) {
        progWatch.increaseTarget( 1 );
        service.post([this](){
            CachedFile::load<float>( ccdResponse, responseFile );
            ++progWatch;
        });
    }

    if (!backgainFile.empty()) {
        progWatch.increaseTarget( 1 );
        service.post([this](){
            CachedFile::load<float>( ccdScattering, backgainFile );
            ++progWatch;
        });
    }

    if( !psfFile.empty() ) {
        progWatch.increaseTarget( 1 );
        service.post([this](){
            CachedFile::load<float>( psf, psfFile );
            ++progWatch;
        });
    }

    if( !mmFile.empty() ) {
        progWatch.increaseTarget( 1 );
        service.post([this](){
            CachedFile::load<float>( modulationMatrix, mmFile );
            ++progWatch;
        });
    }

    if( !xOffsetFile.empty() ) {
        progWatch.increaseTarget( 1 );
        service.post([this](){
            CachedFile::load<int16_t>( xOffset, xOffsetFile );
            if( alignClip.size() == 4 ) {
                Point clipSize( abs(alignClip[3]-alignClip[2])+1, abs(alignClip[1]-alignClip[0])+1 );
                if( clipSize.x != xOffset.dimSize(1) ||
                    clipSize.y != xOffset.dimSize(0) ) {
                    LOG_ERR << "Size of offset file: " << xOffsetFile << " does not match the align_clip." << ende;
                    xOffset.clear();
                }
                //cout << "Loaded offset file: " << xOffsetFile << printArray(xOffset.dimensions(),"  dims") << endl;
                //clipImage( xOffset, alignClip );
            }
            ++progWatch;
        });
    }

    if( !yOffsetFile.empty() ) {
        progWatch.increaseTarget( 1 );
        service.post([this](){
            CachedFile::load<int16_t>( yOffset, yOffsetFile );
            if( alignClip.size() == 4 ) {
                Point clipSize( abs(alignClip[3]-alignClip[2])+1, abs(alignClip[1]-alignClip[0])+1 );
                if( clipSize.x != yOffset.dimSize(1) ||
                    clipSize.y != yOffset.dimSize(0) ) {
                    LOG_ERR << "Size of offset file: " << yOffsetFile << " does not match the align_clip." << ende;
                    yOffset.clear();
                }
                //cout << "Loaded offset file: " << yOffsetFile << printArray(yOffset.dimensions(),"  dims") << endl;
                //clipImage( yOffset, alignClip );
            }
            ++progWatch;
        });
    }
    
    ++progWatch;
    
}


void Channel::loadData( boost::asio::io_service& service, redux::util::Array<PatchData::Ptr>& patches ) {

   //LOG_TRACE << "Channel::loadData()" << ende;
    size_t nFiles = std::max<size_t>( 1, fileNumbers.size() );       // If no numbers, load template as single file

    startT = bpx::pos_infin;
    endT = bpx::neg_infin;
    
    if( nTotalFrames == 0 ) {
        throw logic_error("No input images for channel " + to_string(myObject.ID) + ":" + to_string(ID));
    }
    
    //myJob.progWatch.increaseTarget(nTotalFrames);
    //myObject.progWatch.increaseTarget(nTotalFrames);
    
    // Prepare needed storage
    images.resize( nTotalFrames, imgSize.y, imgSize.x );
    imageStats.resize( nTotalFrames );
    
    progWatch.set(1);
    progWatch.setHandler([this,nFiles,&service,&patches](){
        //cout << "Calib loaded.   " << myObject.ID << ":" << ID << endl;
        for( unsigned int y=0; y<patches.dimSize(0); ++y ) {
            for( unsigned int x=0; x<patches.dimSize(1); ++x ) {
                service.post( [this,&patches,y,x](){
                    auto chData = patches(y,x)->objects[myObject.ID]->channels[ID];
                    adjustCutout( *chData, patches(y,x) );
                    ++myJob.progWatch;
                } );
            }
        }
        size_t nPreviousFrames(0);
        for( size_t i=0; i<nFiles; ++i ) {

            service.post( [&,i,nPreviousFrames](){
                try {
                    loadFile(i,nPreviousFrames);
                    
                    if( imgSize.y < 1 || imgSize.x < 1 ) throw logic_error("Image size is zero.");
                    
                    size_t nF = nFrames[i];

                    for( size_t j=0; j<nF; ++j ) {
                        preprocessImage( nPreviousFrames+j );
                        //service.post(std::bind( &Channel::preprocessImage, this, nPreviousFrames+j) );
                    }
                    if( myObject.saveMask & SF_SAVE_FFDATA ) {
                        myJob.progWatch.increaseTarget(1);
                        progWatch.increaseTarget(1);
                        service.post( [&,i,nF,nPreviousFrames](){
                            bfs::path fn;
                            try {
                                fn = bfs::path (myJob.info.outputDir) / bfs::path (boost::str (boost::format (imageTemplate) % fileNumbers[i])).leaf();
                                bfs::path ext = fn.extension();
                                if(ext.string().empty() || ext.string().length() > 5 ) {     // we assume the filename does not have a proper extenstion, add a temporary dummy
                                    fn = bfs::path( fn.string() + ".ext" );
                                }
                                fn.replace_extension(".cor.f0");
                                LOG_DEBUG << boost::format ("Saving dark/flat corrected file %s.") % fn.string() << ende;
                                LOG_DEBUG << "Saving dark/flat corrected file " << printArray(images.dimensions(),"  imgdims") << ende;
                                Image<float> view( images, nPreviousFrames, nPreviousFrames+nF-1, 0, imgSize.y-1, 0, imgSize.x-1 );
                                redux::file::Ana::write( fn.string(), view.copy() );   // TODO: other formats
                            } catch ( const std::exception& e ) {
                                LOG_ERR << "Failed to save corrected file: " << fn << "  reason: " << e.what() << ende;
                            }
                            ++progWatch;
                            ++myJob.progWatch;
                        });
                    } 
                } catch ( const std::exception& e ) {
                    LOG_ERR << "Failed to load/preprocess file. reason: " << e.what() << ende;
                } catch ( ... ) {
                    LOG_ERR << "Failed to load/preprocess file for unknown reason." << ende;
                }
            });
            nPreviousFrames += nFrames[i];
        }
        
    });
    loadCalib(service);

}


void Channel::storePatches(boost::asio::io_service& service, Array<PatchData::Ptr>& patches) {

    size_t nPatchesY = patches.dimSize(0);
    size_t nPatchesX = patches.dimSize(1);
    progWatch.set(nPatchesX*nPatchesY);
    progWatch.setHandler([this](){
        images.clear();                                                                 // release resources after storing the patch-data.
    });

    for(unsigned int py=0; py<nPatchesY; ++py) {
        for(unsigned int px=0; px<nPatchesX; ++px) {
            service.post( [this,&patches,py,px](){
                auto chData = patches(py,px)->objects[myObject.ID]->channels[ID];
                copyImagesToPatch(*chData);
                chData->cacheStore(true);    // store to disk and clear array
                ++progWatch;
            });
        }
    }
    
}


void Channel::unloadData(void) {               // unload what was accessed through the cache, this should be called when all objects are done pre-processing.
    
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

    if( fileNumbers.empty() ) {         // single file
        bfs::path fn = bfs::path(imageDataDir) / bfs::path(imageTemplate);
        CachedFile::unload<float>( fn.string() );
    } else {
        size_t nFiles = std::max<size_t>( 1, fileNumbers.size() );       // If no numbers, load template as single file
        for( size_t i=0; i<nFiles; ++i ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( imageTemplate ) % fileNumbers[i] ) );
            CachedFile::unload<float>( fn.string() );
        }
    }

}


double Channel::getMaxMean(void) {
    double maxMean = std::numeric_limits<double>::lowest();
    for ( const ArrayStats::Ptr &stat : imageStats ) {
        if( stat ) maxMean = std::max(maxMean,stat->mean);
    }
    return maxMean;
}


void Channel::getFileNames(std::vector<std::string>& files) const {

    for (auto &num: fileNumbers) {
        bfs::path fn = bfs::path(imageDataDir)/bfs::path(boost::str(boost::format(imageTemplate) % num));
        files.push_back(fn.string());
    }
    
}


void Channel::initProcessing( const Solver& solver ) {
    
    initPhiFixed();

    subImages.resize(nTotalFrames);
    for (uint16_t i=0; i < nTotalFrames; ++i) {
        if( !subImages[i] ) subImages[i].reset( new SubImage(myObject, *this, solver.window, solver.noiseWindow) );
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

    size_t blockSizeY = cd.images.dimSize(1);
    size_t blockSizeX = cd.images.dimSize(2);
    if( (imageStats.size() == nImages) && (cd.images.nDimensions() == 3) ) {
        float* dataPtr = cd.images.get();
        size_t imgPixels = blockSizeY*blockSizeX;
        for( uint16_t i=0; i<nImages; ++i ) {
            double scale = myObject.objMaxMean/imageStats[i]->mean;
            //double scale = 1.0/imageStats[i]->mean;
            transform(dataPtr, dataPtr+imgPixels, dataPtr, bind1st(multiplies<double>(),scale) );
            dataPtr += imgPixels;
        }
    }
    
    uint16_t patchSize = myObject.patchSize;
    
    PointF localShift; 
    int firstY = max( cd.offset.y, 0 );
    int lastY = min<int>( firstY+patchSize-1, blockSizeY-1);
    firstY = max( lastY-patchSize+1, 0 );   // this should actually never go out-of-bounds unless patchSize > blockSize
    localShift.y = firstY-cd.offset.y;
    int firstX = max( cd.offset.x, 0 );
    int lastX = min<int>( firstX+patchSize-1, blockSizeX-1);
    firstX = max( lastX-patchSize+1, 0 );
    localShift.x = firstX-cd.offset.x;

    for (uint16_t i=0; i < nImages; ++i) {
        subImages[i]->setPatchInfo( i, cd.channelOffset, cd.residualOffset, patchSize, myObject.pupilPixels, myJob.modeNumbers.size() );
        subImages[i]->wrap( cd.images, i, i, firstY, lastY, firstX, lastX );
        subImages[i]->stats.getStats( cd.images.ptr(i,0,0), blockSizeY*blockSizeX, ST_VALUES|ST_RMS );
    }

    phi_fixed.copy( phi_channel );
    
    double* phiPtr = phi_channel.get();
    size_t pupilSize2 = myObject.pupilPixels*myObject.pupilPixels;
    
    PointF totalshift = localShift + cd.residualOffset;
    
    int32_t mIndex = myObject.modes.tiltMode.x;
    if( mIndex >= 0 && fabs(totalshift.x) > 0 ) {
        const double* modePtr = myObject.modes.modePointers[mIndex];
        float res = -totalshift.x*myObject.shiftToAlpha.x;   // positive coefficient shifts image to the left
        transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
            [res](const double& p, const double& m) {
                return p + res*m;
            });
    }
    
    mIndex = myObject.modes.tiltMode.y;
    if( mIndex >= 0 && fabs(totalshift.y) > 0 ) {
        const double* modePtr = myObject.modes.modePointers[mIndex];
        float res = -totalshift.y*myObject.shiftToAlpha.y;   // positive coefficient shifts image downwards
        transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
            [res](const double& p, const double& m) {
                return p + res*m;
            });
    }

    
}


void Channel::initPhiFixed(void) {

    phi_fixed.resize( myObject.pupilPixels, myObject.pupilPixels );
    phi_fixed.zero();
    double* phiPtr = phi_fixed.get();
    size_t pupilSize2 = myObject.pupilPixels*myObject.pupilPixels;
    
    ModeInfo mi (myJob.klMinMode, myJob.klMaxMode, 0, myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, myJob.klCutoff);
    for (unsigned int i = 0; i < diversityModes.size(); ++i) {
        uint16_t modeNumber = diversityModes[i];
        ModeInfo mi2 = mi;
        if (modeNumber == 2 || modeNumber == 3 || diversityTypes[i] == ZERNIKE) {
            mi2.firstMode = mi2.lastMode = 0;
        }
        mi2.modeNumber = modeNumber;
        ModeSet& ms = myJob.globalData->get(mi2);

        if( ms.empty() ) {    // generate
            if( diversityTypes[i] == ZERNIKE ) {
                ms.generate( myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, myJob.modeNumbers );
            } else {
                ms.generate( myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, myJob.klMinMode, myJob.klMaxMode, myJob.modeNumbers, myJob.klCutoff );
            }
            if( ms.nDimensions() != 3 || ms.dimSize(1) != myObject.pupilPixels || ms.dimSize(2) != myObject.pupilPixels ) {    // mismatch
                LOG_ERR << "Generated ModeSet does not match. This should NOT happen!!" << ende;
            } else {
                LOG_DEBUG << "Generated Modeset with " << ms.dimSize(0) << " modes. (" << myObject.pupilPixels << "x" << myObject.pupilPixels << "  radius=" << myObject.pupilRadiusInPixels << ")" << ende;

//                ms.getNorms( myObject.pupil );
//cout << "initPhiFixed " << __LINE__ << endl;
//                ms.normalize( 1.0/myObject.wavelength );
//cout << "initPhiFixed " << __LINE__ << endl;
            }
        }        
        auto it = std::find(ms.modeNumbers.begin(), ms.modeNumbers.end(), modeNumber);
        if( it != ms.modeNumbers.end() ) {
            double div = diversity[i]/myObject.wavelength;     // minus/sign is just to keep old cfg-files usable.
            const double* modePtr = ms.modePointers[static_cast<size_t>(it-ms.modeNumbers.begin())];
            transform( phiPtr, phiPtr+pupilSize2, modePtr, phiPtr,
                [div](const double& p, const double& m) {
                    return p + div*m;
                });
        }

    }

}


void Channel::addAllFT (redux::util::Array<double>& ftsum) {
    for (auto& subimage : subImages) {
        subimage->addFT (ftsum);
    }
}

// void Channel::addAllPQ(void) const {
//     for (const auto& subimage : subImages) {
//         myObject.addToPQ( subimage->imgFT.get(), subimage->OTF.get() );
//     }
// }


void Channel::addTimeStamps( const bpx::ptime& newStart, const bpx::ptime& newEnd ) {

    unique_lock<mutex> lock(mtx);
    
    if(startT.is_special()) startT = newStart;
    else startT = std::min( startT, newStart );
    if(endT.is_special()) endT = newEnd;
    else endT = std::max( endT, newEnd );

}


void Channel::loadFile( size_t fileIndex, size_t offset ) {
    static mutex lmtx;
    unique_lock<mutex> lock(lmtx);
    //cout << ("imageDataDir="+imageDataDir+" imageTemplate="+imageTemplate) << flush;
    bfs::path fn = bfs::path(imageDataDir) / bfs::path( );
    //cout << "...ok"  << endl;
    
    try {
        
        if( fileIndex >= fileNumbers.size() ) {
            string msg = "File-index is out of range. (" + to_string(fileIndex) + "/" + to_string(fileNumbers.size()) + ").";
            throw logic_error(msg);
        }
        if( !fileNumbers.empty() ) {
            fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( imageTemplate ) % fileNumbers[fileIndex] ) );
        }
    
        string imStr = to_string(fileIndex);
        Image<float> tmpImg;
        CachedFile::load( tmpImg, fn.string() );
        size_t nF = nFrames[fileIndex];
        Image<float> view( images, offset, offset+nF-1, 0, imgSize.y-1, 0, imgSize.x-1 );
        if( nF > 1 ) {
            Image<float> tmpView( tmpImg, discard[0], discard[0]+nF-1, 0, imgSize.y-1, 0, imgSize.x-1 );
            view.assign( reinterpret_cast<redux::util::Array<float>&>(tmpView) );
            imStr = to_string(offset)+"-"+to_string(offset+nF-1);
        } else {
            view.assign( reinterpret_cast<redux::util::Array<float>&>(tmpImg) );           
        }
        LOG << boost::format ("Loaded file "+imageTemplate+"  (%d:%d:%s)") % fileNumbers[fileIndex] % myObject.ID % ID % imStr << ende;
        myJob.info.progress[0]++;
        if( tmpImg.meta ) {
            addTimeStamps( tmpImg.meta->getStartTime(), tmpImg.meta->getEndTime() );
        }

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
        
        LOG_TRACE << boost::format ("Preprocessing image (%d:%d:%d)  %s") % myObject.ID % ID % i % printArray(tmpImg.dimensions(),"dims") << ende;
        
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
                    LOG_DETAIL << "Applying correction for CCD transparency." << ende;
                    redux::image::descatter (tmpImg, ccdScattering, psf);
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
                    function<double (size_t, size_t) > func = bind (horizontalInterpolation<double>, arrayPtr, sy, sx, sp::_1, sp::_2);
                    fillPixels (arrayPtr, sy, sx, func, std::bind2nd (std::less_equal<double>(), myJob.badPixelThreshold), mask2D.get());
                    break;
                }
                case FPM_MEDIAN: {
                    // TODO: median method
                    break;
                }
                case FPM_INVDISTWEIGHT:       // inverse distance weighting is the default method, so fall through
                default: {
                    //LOG_DETAIL << "Filling bad pixels using inverse distance weighted average." << ende;
                    function<double (size_t, size_t) > func = bind (inverseDistanceWeight<double>, arrayPtr, sy, sx, sp::_1, sp::_2);
                    fillPixels (arrayPtr, sy, sx, func, std::bind2nd (std::less_equal<double>(), myJob.badPixelThreshold), mask2D.get());
                }
            }

            // FIXME: This is a hack to create truncated values as the old code!!
            //for(size_t i=0; i<sy*sx; ++i) arrayPtr[0][i] = (int)arrayPtr[0][i];
            
        }

        view.assign(tmpImg);                            // copy back to image
        imageStats[i].reset( new ArrayStats() );
        imageStats[i]->getStats( borderClip, tmpImg );    // get stats for corrected data
        myJob.info.progress[0]++;
    
        ++myJob.progWatch;          // this will trigger calculating max_mean and storing patches in Object::loadData()
        ++myObject.progWatch;       // this will trigger unloading the calibration data after all objects/images are done.
        
    } catch ( const std::exception& e ) {
        LOG_ERR << "Failed to preprocess image #" << i << ". reason: " << e.what() << ende;
    } catch ( ... ) {
        LOG_ERR << "Failed to preprocess image #" << i << " for unknown reason."  << ende;
    }

}


void Channel::copyImagesToPatch(ChannelData& chData) {

    Array<float> block(reinterpret_cast<redux::util::Array<float>&>(images), 0, nTotalFrames-1, chData.cutout.first.y, chData.cutout.last.y, chData.cutout.first.x, chData.cutout.last.x);
    Array<float> tmp = block.copy(true);
    
    bool flipX(false);
    bool flipY(false);
    if( alignMap.size() == 9 ) {
        if( alignMap[0] < 0 ) flipX = true;
        if( alignMap[0] < 4 ) flipY = true;
    } else if( alignClip.size() == 4) {
        if( alignClip[0] > alignClip[1] ) flipX = true;
        if( alignClip[2] > alignClip[3] ) flipY = true;
    }
    if( flipX || flipY ) {
        size_t sy = chData.cutout.last.y - chData.cutout.first.y + 1;
        size_t sx = chData.cutout.last.x - chData.cutout.first.x + 1;
        std::shared_ptr<float*> arrayPtr = tmp.reshape(nTotalFrames*sy, sx);
        float** imgPtr = arrayPtr.get();
        for( size_t i=0; i<nTotalFrames; ++i) {
            if( flipX ) redux::util::reverseX(imgPtr, sy, sx);
            if( flipY ) redux::util::reverseY(imgPtr, sy, sx);
            imgPtr += sy;
        }
    }
    
    chData.images = std::move(tmp);       // chData.images will share datablock with the "images" stack, so minimal RAM usage.
    chData.setLoaded();
    
}


uint32_t Channel::nImages(void) {
    
    if( nTotalFrames == 0 ) {
        nTotalFrames = std::accumulate( nFrames.begin(), nFrames.end(), 0 );
    }
    
    return nTotalFrames;
    
}


void Channel::adjustCutout( ChannelData& chData, const PatchData::Ptr& patch ) const {

    Point16 refPos = patch->position;       // position in clipped reference-channel coordinates
    
    PointF localPos;                        // position in this channel in non-clipped coordinates
    uint16_t halfPatch = myObject.patchSize/2;
    RegionI desiredCutout( myObject.patchSize-1, myObject.patchSize-1 );                    // local cut-out
    desiredCutout -= PointI( halfPatch, halfPatch );                                        // ...now centered on 0,0

    chData.residualOffset = 0;
    chData.offset = myObject.maxLocalShift;
    
    
    bool flipX(false);
    bool flipY(false);
    
    if( alignMap.size() == 9 ) {
        if( alignMap[0] < 0 ) flipX = true;
        if( alignMap[0] < 4 ) flipY = true;
        ProjectiveMap map( alignMap );
        localPos = map*(refPos + myJob.roi.first);  // position in reference-channel, global coordinates
    } else {                                // old style alignment with clips & offsetfiles
        localPos = refPos;
        if( xOffset.valid() ) {
            float shift = xOffset(refPos.y, refPos.x)/100.0;
            chData.residualOffset.x = shift;
            localPos.x += shift;
        }
        if( yOffset.valid() ) {
            float shift = yOffset(refPos.y, refPos.x)/100.0;
            chData.residualOffset.y = shift;
            localPos.y += shift;
        }
        if( alignClip.size() == 4 ) {
            if( alignClip[0] > alignClip[1] ) {
                flipX = true;
                localPos.x = alignClip[0] - localPos.x - 1 + myObject.patchSize%2;
            } else {
                localPos.x += alignClip[0] - 1;
            }
            if( alignClip[2] > alignClip[3] ) {
                flipY = true;
                localPos.y = alignClip[2] - localPos.y - 1 + myObject.patchSize%2;
            } else {
                localPos.y += alignClip[2] - 1;
            }
        }
    }
    
    if( imgSize == 0 ) {
        LOG_ERR << "No valid imgSize when adjusting cutout, that should not happen..." << ende;
        return;
    }
    
    if( flipX ) localPos.x += 1;
    if( flipY ) localPos.y += 1;
    
    PointI finalPos = localPos+0.5;
    desiredCutout += finalPos;                                                  // ...now centered on localPos, this is the cutout we want
    RegionI imgBoundary(0,0,imgSize.y-1, imgSize.x-1);                          // restrict to lie inside image
    
    RegionI tmpCutout = desiredCutout;
    tmpCutout.restrict(imgBoundary);
    bool alreadyWarned(false);
    if( tmpCutout != desiredCutout ) {
        LOG_WARN << "Patch " << patch->index << " does not lie completely within image in channel " << myObject.ID << ":" << ID
        << ": desired: " << desiredCutout << "  This will likely cause severe artifacts!!!" << ende;
        alreadyWarned = true;
        if( !flipY && (tmpCutout.first.y != desiredCutout.first.y) ) {
            chData.offset.y -= (tmpCutout.first.y - desiredCutout.first.y);
            finalPos.y += (tmpCutout.first.y - desiredCutout.first.y);
        }
        if( !flipX && (tmpCutout.first.x != desiredCutout.first.x) ) {
            chData.offset.x -= (tmpCutout.first.x - desiredCutout.first.x);
            finalPos.x += (tmpCutout.first.x - desiredCutout.first.x);
        }
        if( flipY && (tmpCutout.last.y != desiredCutout.last.y) ) {
            chData.offset.y += (tmpCutout.last.y - desiredCutout.last.y);
            finalPos.y -= (tmpCutout.last.y - desiredCutout.last.y);
        }
        if( flipX && (tmpCutout.last.x != desiredCutout.last.x) ) {
            chData.offset.x += (tmpCutout.last.x - desiredCutout.last.x);
            finalPos.x -= (tmpCutout.last.x - desiredCutout.last.x);
        }
    }
    
    chData.channelOffset = PointI( lround(chData.residualOffset.y), lround(chData.residualOffset.x) );  // possibly apply shift from offsetfiles.
    chData.channelOffset -= imgBoundary.outside( tmpCutout+chData.channelOffset );                      // restrict the shift inside the image, leave the rest in "residualOffset" to be dealt with using Zernike tilts.
    chData.residualOffset = localPos - finalPos;
    
    tmpCutout.grow( myObject.maxLocalShift );                               // add maxLocalshift
    desiredCutout = tmpCutout;
    tmpCutout.restrict(imgBoundary);
    if( tmpCutout != desiredCutout ) {
        if( !alreadyWarned )
            LOG_DETAIL << "Patch " << patch->index << " + maxLocalShift does not lie completely within image in channel " <<
            myObject.ID << ":" << ID << ":  desired: " << desiredCutout << "  actual: " << tmpCutout << ende;
        if( !flipY && (tmpCutout.first.y != desiredCutout.first.y) ) {
            chData.offset.y -= (tmpCutout.first.y - desiredCutout.first.y);
        }
        if( !flipX && (tmpCutout.first.x != desiredCutout.first.x) ) {
            chData.offset.x -= (tmpCutout.first.x - desiredCutout.first.x);
        }
        if( flipY && (tmpCutout.last.y != desiredCutout.last.y) ) {
            chData.offset.y += (tmpCutout.last.y - desiredCutout.last.y);
        }
        if( flipX && (tmpCutout.last.x != desiredCutout.last.x) ) {
            chData.offset.x += (tmpCutout.last.x - desiredCutout.last.x);
        }
    }
    
    chData.cutout = tmpCutout;
    string actStr = "";
    if( chData.cutout != desiredCutout ) actStr = "  actual="+(string)chData.cutout;
    LOG_TRACE << "AdjustCutout ch=" <<myObject.ID << ":" << ID << ": patch=" << patch->index << refPos << "  mapped=" << localPos
              << "  desired=" << desiredCutout << actStr
              << "  localOffset=" << chData.channelOffset << "   start=" << chData.offset
              << "   residual=" << chData.residualOffset << ende;
         
}


void Channel::adjustCutouts( Array<PatchData::Ptr>& patches ) {
    
    if( patches.nDimensions() == 2 ) {
        size_t nPatchesY = patches.dimSize(0);
        size_t nPatchesX = patches.dimSize(1);
        for( unsigned int py=0; py<nPatchesY; ++py ) {
            for( unsigned int px=0; px<nPatchesX; ++px ) {
                PatchData::Ptr& patch( patches(py,px) );
                auto chData = patch->objects[myObject.ID]->channels[ID];
                adjustCutout( *chData, patch );
            }
        }
    }
    
}


Point16 Channel::getImageSize(void) {

    if( imgSize == 0 ) {
//         if( alignClip.size() == 4 ) {   // we have align-clip, but no mapping => reference channel.
//             imgSize = Point16(abs(alignClip[3]-alignClip[2])+1, abs(alignClip[1]-alignClip[0])+1);
//         } else {                        //  No align-map or align-clip, get full image size.
            bfs::path fn;
            if( nImages() ) {
                fn = bfs::path(imageDataDir) / bfs::path(boost::str(boost::format (imageTemplate) % fileNumbers[0]));
            } else {
                fn = bfs::path(imageDataDir) / bfs::path(imageTemplate);
            }
            Image<float> tmp;
            CachedFile::load( tmp, fn.string(), true );       // Image will be cached, so this is not a "wasted" load.

            uint8_t nDims = tmp.meta->nDims();

            if( nDims == 3 ) {
                imgSize = Point16( tmp.meta->dimSize(1), tmp.meta->dimSize(2) );
            } else if( nDims == 2 ) {
                imgSize = Point16( tmp.meta->dimSize(0), tmp.meta->dimSize(1) );
            } else LOG_ERR << "Image " << fn << " is not 2D or 3D." << ende;
//        }
    }

    return imgSize;
 
}


void Channel::logAndThrow( string msg ) {
    msg = to_string(myObject.ID)+":"+to_string(ID)+": "+msg;
    LOG_ERR << msg << ende;
    throw job_error(msg);
    
}


void Channel::dump (std::string tag) {
    
    tag += "_c"+to_string(ID);
    
    Ana::write (tag + "_phi_fixed.f0", phi_fixed);
    Ana::write (tag + "_phi_channel.f0", phi_channel);
    
    uint16_t patchSize = myObject.patchSize;
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
    
    float* tmpPtr = tmpF.get();
    for( shared_ptr<SubImage>& im: subImages ) {
        im->copyTo<float>(tmpPtr);
        tmpPtr += blockSize;
    }
    Ana::write( tag + "_img.f0", tmpF );
    
    tmpPtr = tmpF.get();
    complex_t* ftPtr = tmpC.get();
    int idx(0);
    for( shared_ptr<SubImage>& im: subImages ) {
        im->getWindowedImg( tmpD, s );
        tmpFT.reset( tmpD.get(), patchSize, patchSize, FT_FULLCOMPLEX );
        FourierTransform::reorder(tmpFT);
        tmpD.copyTo<float>(tmpPtr);
        tmpFT.copyTo<complex_t>(ftPtr);
        const vector<int64_t>& first = im->first();
        shiftArr(idx,0) = first[1];
        shiftArr(idx,1) = first[2];
        statArr(idx,0) = s.min;
        statArr(idx,1) = s.max;
        statArr(idx,2) = s.mean;
        statArr(idx++,3) = s.stddev;
        tmpPtr += blockSize;
        ftPtr += blockSize;
    }
    Ana::write( tag + "_wimg.f0", tmpF );
    Ana::write( tag + "_ft.f0", tmpC );
    Ana::write( tag + "_stat.f0", statArr );
    Ana::write( tag + "_shift.f0", shiftArr );
    tmpC.resize();      // free some memory
    
    tmpPtr = tmpF.get();
    for( shared_ptr<SubImage>& im: subImages ) {
        im->getPSF( tmpD.get() );
        tmpD.copyTo<float>(tmpPtr);
        tmpPtr += blockSize;
    }
    Ana::write( tag + "_psf.f0", tmpF );
    tmpF.resize();      // free some memory

    
}

