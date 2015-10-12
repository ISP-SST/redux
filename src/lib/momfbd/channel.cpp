#include "redux/momfbd/channel.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/object.hpp"
#include "redux/momfbd/subimage.hpp"
#include "redux/momfbd/util.hpp"

#include "redux/constants.hpp"
#include "redux/file/fileana.hpp"
#include "redux/file/fileio.hpp"
#include "redux/math/functions.hpp"
#include "redux/image/utils.hpp"
#include "redux/logger.hpp"
#include "redux/translators.hpp"
#include "redux/util/stringutil.hpp"

#include <functional>
#include <math.h>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "redux/image/fouriertransform.hpp"

namespace bfs = boost::filesystem;
using namespace redux::file;
using namespace redux::image;
using namespace redux::momfbd;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

    const string thisChannel = "channel";
    
    struct ClippedFile {
        string filename;
        const vector<int16_t> clip;
        bool symmetricClip;

        ClippedFile(const string& fn, const vector<int16_t>& cl, bool sym=false) :
            filename(fn), clip(cl), symmetricClip(sym) { }

        bool operator<(const ClippedFile& rhs) const {
            if(filename == rhs.filename) {
                if(symmetricClip == rhs.symmetricClip) {
                    return (clip < rhs.clip);
                }
                return (symmetricClip < rhs.symmetricClip);
            }
            return (filename < rhs.filename);
        }

        template <typename T>
        void loadImage(Image<T>& img, bool normalize=false) const {

            redux::file::readFile(filename, img);
            
            if( clip.size() == 4 ) {
                bool flipX = false, flipY = false;
                vector<int16_t> alignClip = clip;
                if (alignClip[0] > alignClip[1]) {     // we have the y (row/slow) dimension first, momfbd cfg-files (and thus alignClip) has x first. TBD: how should this be ?
                    std::swap (alignClip[0], alignClip[1]);
                    flipX = true;
                }
                if (alignClip[2] > alignClip[3]) {
                    std::swap (alignClip[2], alignClip[3]);
                    flipY = true;
                }
                for (auto & index : alignClip)
                    --index;       // NOTE: momfbd cfg files uses 1-based indexes, internally we start with 0.
                size_t sy = alignClip[3] - alignClip[2] + 1;
                size_t sx = alignClip[1] - alignClip[0] + 1;
                if(symmetricClip) {
                    const std::vector<size_t>& dims = img.dimensions();
                    int skewY = (dims[0] - sy) / 2  - alignClip[2];
                    int skewX = (dims[1] - sx) / 2  - alignClip[0];
                    alignClip[0] += skewX;
                    alignClip[1] += skewX;
                    alignClip[2] += skewY;
                    alignClip[3] += skewY;
                }
                img.setLimits (alignClip[2], alignClip[3], alignClip[0], alignClip[1]);
                img.trim();

                if (flipX || flipY) {
                    shared_ptr<T*> arrayPtr = img.reshape(sy, sx);
                    T** imgPtr = arrayPtr.get();
                    if (flipX) reverseX(imgPtr, sy, sx);
                    if (flipY) reverseY(imgPtr, sy, sx);
                }
            }

        }

        template <typename T>
        static void load(Image<T>& img, const string& fn, const vector<int16_t>& cl, bool norm=false, bool sym=false) {
            ClippedFile cf(fn,cl,sym);
            Image<T>& cimg = redux::util::Cache::get<ClippedFile, Image<T> >(cf);
            {
                unique_lock<mutex> lock(cimg.imgMutex);
                if(cimg.nElements() == 0) cf.loadImage(cimg,norm);
            }
            img = cimg;
        }

        template <typename T>
        static void unload(const string& fn, const vector<int16_t>& cl, bool sym=false) {
            ClippedFile cf(fn,cl,sym);
            redux::util::Cache::erase<ClippedFile, Image<T>>(cf);
        }
    };

}


Channel::Channel (Object& o, MomfbdJob& j, uint16_t id) : ID (id), myObject (o), myJob (j) {

}

Channel::~Channel() {

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


size_t Channel::size (void) const {

    size_t sz = ChannelCfg::size();
    sz += sizeof (uint16_t);          // ID
    sz += imgSize.size();
    sz += imageStats.size() * ArrayStats::size() + sizeof (uint16_t);
    return sz;
}


uint64_t Channel::pack (char* ptr) const {
    using redux::util::pack;
    uint64_t count = ChannelCfg::pack (ptr);
    count += pack (ptr + count, ID);
    count += imgSize.pack (ptr + count);
    uint16_t statSize = imageStats.size();
    count += pack (ptr + count, statSize);
    for (auto & stat : imageStats)
        count += stat->pack (ptr + count);
    if (count != size()) {
        LOG_ERR << "(" << hexString (this) << "): Packing failed, there is a size mismatch:  count = " << count << "  sz = " << size();
    }
    return count;
}


uint64_t Channel::unpack (const char* ptr, bool swap_endian) {
    using redux::util::unpack;

    uint64_t count = ChannelCfg::unpack (ptr, swap_endian);
    count += unpack (ptr + count, ID, swap_endian);
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
        LOG_ERR << "No filename template specified.";
        return false;
    }
    size_t nWild = std::count (imageTemplate.begin(), imageTemplate.end(), '%');
    if (nWild > 2) {
        LOG_ERR << "Filename template contains too many wildcards: \"" << imageTemplate << "\"";
        return false;
    } else if (nWild == 1 && imageNumbers.empty()) {
        LOG_ERR << "Filename template contains wildcard and no image-numbers given (with IMAGE_NUM)";
        return false;
    } /*else if( nWild == 2 && sequenceNumber == 0 ) {
        LOG_ERR << "Filename template contains 2 wildcards and no sequence-number given (with SEQUENCE_NUM)";
        return false;
    }*/


    if( darkTemplate.empty() != gainFile.empty() ) {
        LOG_ERR << "Either BOTH or NONE of DARK_TEMPLATE/GAIN_FILE has to be specified!!!";
        return false;
    } else if (darkTemplate.empty()) {
        LOG_WARN << "No dark/gain files specified, assuming the data is pre-processed already.";
        if ( !responseFile.empty() ) {
            LOG_ERR << "Detector response correction only possible when flatfielding, It will NOT be performed!!!";
        }
    } else {
        nWild = std::count (darkTemplate.begin(), darkTemplate.end(), '%');
        if (nWild > 1) {
            LOG_ERR << "Dark template contains too many wildcards: \"" << darkTemplate << "\"";
            return false;
        } else if (nWild == 1 && darkNumbers.empty()) {
            LOG_ERR << "Dark template contains wildcard and no dark-numbers given (with DARK_NUM)";
            return false;
        } else if (nWild == 0 && darkNumbers.size()) {
            LOG_WARN << "Dark template contains no wildcard AND dark-numbers specified. Numbers will be ignored and the dark-template used as a single filename.";
            darkNumbers.clear();    // TODO: fix this properly, numbers might reappear after transfer (because of inheritance)
        }
    }
    
    if( !alignClip.empty() && alignClip.size() != 4 ) {
        LOG_WARN << "ALIGN_CLIP does not seem to contain 4 integers. Whole image area will be used.";
        alignClip.clear();
    }
    
    return true;

}


bool Channel::checkData (void) {

    // Images
    if ( incomplete ) {   // check if files are present
        for (size_t i (0); i < imageNumbers.size();) {
            bfs::path fn = bfs::path (boost::str (boost::format (imageTemplate) % (imageNumberOffset + imageNumbers[i])));
            if (!bfs::exists (fn)) {
                fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % (imageNumberOffset + imageNumbers[i])));
                if (!bfs::exists (fn)) {
                    //LOG_TRACE << "File not found: \"" << fn.string() << "\", removing from list of image numbers.";
                    imageNumbers.erase(imageNumbers.begin() + i);
                    continue;
                }
            }
            ++i;
        }
        if (imageNumbers.empty()) {
            LOG_CRITICAL << boost::format ("No files found for incomplete object with filename template \"%s\" in directory \"%s\"") % imageTemplate % imageDataDir;
            return false;
        }
    }
    if (imageNumbers.empty()) {         // single file
        bfs::path fn = bfs::path (imageDataDir) / bfs::path (imageTemplate);
        if (! bfs::exists (fn)) {
            LOG_ERR << boost::format ("Image-file %s not found!") % fn;
            return false;
        }
    } else {                            // template + numbers
        for (auto & number : imageNumbers) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % (imageNumberOffset + number)));
            if (!bfs::exists (fn)) {
                LOG_ERR << boost::format ("Image-file %s not found!") % boost::str (boost::format (imageTemplate) % (imageNumberOffset + number));
                return false;
            }
        }
    }


    // Dark(s)
    size_t nWild = std::count (darkTemplate.begin(), darkTemplate.end(), '%');
    if (nWild == 0 || darkNumbers.empty()) {          // single file, DARK_NUM will be ignored if no wildcard in the template
        if (! bfs::exists (bfs::path (darkTemplate))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (darkTemplate);
            if (! bfs::exists (fn)) {
                logAndThrow("Dark-file " + darkTemplate + " not found!");
            } else darkTemplate = fn.c_str();
        }
    } else {                            // template
        for (auto & number : darkNumbers) {
            bfs::path fn = bfs::path (boost::str (boost::format (darkTemplate) % number));
            if (!bfs::exists (fn)) {
                fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (darkTemplate) % number));
                if (!bfs::exists (fn)) {
                    logAndThrow("Dark-file " + (boost::format(darkTemplate) % number).str() + " not found!");
                } else darkTemplate = fn.c_str();
            }
        }
    }


    // Gain
    if (!gainFile.empty()) {
        if (! bfs::exists (bfs::path (gainFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (gainFile);
            if (! bfs::exists (fn)) {
                logAndThrow("Gain-file " + gainFile + " not found!");
            } else gainFile = fn.c_str();
        }
    }

    if (!responseFile.empty()) {
        if (! bfs::exists (bfs::path (responseFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (responseFile);
            if (! bfs::exists (fn)) {
                logAndThrow("Response-file " + responseFile + " not found!");
            } else responseFile = fn.c_str();
        }
    }

    if (!backgainFile.empty()) {
        if (! bfs::exists (bfs::path (backgainFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (backgainFile);
            if (! bfs::exists (fn)) {
                logAndThrow("Backgain-file " + backgainFile + " not found!");
            } else backgainFile = fn.c_str();
        }
    }

    if (!psfFile.empty()) {
        if (! bfs::exists (bfs::path (psfFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (psfFile);
            if (! bfs::exists (fn)) {
                logAndThrow("PSF-file " + psfFile + " not found!");
            } else psfFile = fn.c_str();
        }
    }

    if (!mmFile.empty()) {
        if (! bfs::exists (bfs::path (mmFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (mmFile);
            if (! bfs::exists (fn)) {
                logAndThrow("Modulation-matrix file " + mmFile + " not found!");
            } else mmFile = fn.c_str();
        }
    }

    if (!xOffsetFile.empty()) {
        if (! bfs::exists (bfs::path (xOffsetFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (xOffsetFile);
            if (! bfs::exists (fn)) {
                logAndThrow("Offset-file " + xOffsetFile + " not found!");
            } else xOffsetFile = fn.c_str();
        }
    }

    if (!yOffsetFile.empty()) {
        if (! bfs::exists (bfs::path (yOffsetFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (yOffsetFile);
            if (! bfs::exists (fn)) {
                logAndThrow("Offset-file " + yOffsetFile + " not found!");
            } else yOffsetFile = fn.c_str();
        }
    }

    return true;
}


//#define INDEX_THRESHOLD  0
#define INDEX_THRESHOLD  1E-12      // cutoff to get rid of some fft noise
// TBD: this should be a parameter somewhere...
void Channel::initCache (void) {

    //LOG_DETAIL << "wavelength = " << myObject.wavelength << "   patchSize = " << patchSize << "  telescopeD = " << myJob.telescopeD << "  arcSecsPerPixel = " << arcSecsPerPixel;
//     Pupil::calculatePupilSize (frequencyCutoff, pupilRadiusInPixels, pupilPixels, myObject.wavelength, patchSize, myJob.telescopeD, arcSecsPerPixel);
//     
//     myJob.patchSize = myObject.patchSize = patchSize;     // TODO: fulhack until per-channel sizes is implemented
//     myJob.pupilPixels = myObject.pupilPixels = pupilPixels;
    //size_t otfPixels = 2 * pupilPixels;
    //LOG_DETAIL << "frequencyCutoff = " << frequencyCutoff << "  pupilSize = " << pupilPixels << "  pupilRadiusInPixels = " << pupilRadiusInPixels;
    //pupil = myJob.globalData->fetch (pupilPixels, pupilRadiusInPixels);
    // Create a temporary OTF and store the indices where the OTF/pupil are non-zero. This will be used in loops to skip irrelevant evaluations.
    /*Array<double> tmpImg (2 * pupilPixels, 2 * pupilPixels);
    tmpImg.zero();
    Array<double> tmpSubImg (tmpImg, 0, pupilPixels - 1, 0, pupilPixels - 1);
    pupil.first.copy(tmpSubImg);
    size_t cnt(0);
    for (auto & it: tmpSubImg) {      // map indices where the pupil-mask is non-zero.
        if( it > INDEX_THRESHOLD ) {
            pupilIndices.insert (cnt);
            size_t otfOffset = (cnt / pupilPixels) * otfPixels + (cnt % pupilPixels) + (pupilPixels / 2) + (pupilPixels / 2) * otfPixels;
            pupilInOTF.insert(make_pair (cnt, otfOffset));
            it = 1;
        }
        cnt++;
    }
    FourierTransform::autocorrelate (tmpImg);

    LOG_DEBUG << "Generated pupilIndices with " << pupilIndices.size() << " elements. (full size = " << pupil.first.nElements() << ", area = " << pupil.second << " )";
    double* tmpPtr = tmpImg.get();

    for (size_t index = 0; index < tmpImg.nElements(); ++index) {           // map indices where the OTF-mask (auto-correlated pupil-mask) is non-zero.
        if (fabs(tmpPtr[index]) > INDEX_THRESHOLD) {
            otfIndices.insert (index);
            //tmpPtr[index] = 1;
        } else tmpPtr[index] = 0;
    }
    //Ana::write ("otfsupport.f0", tmpImg);

    LOG_DEBUG << "Generated otfIndices with " << otfIndices.size() << " elements. (full size = " << tmpImg.nElements() << ")";
    */
  
    
    
    
/*
    Cache::ModeID id (myJob.klMinMode, myJob.klMaxMode, 0, pupilPixels, pupilRadiusInPixels, rotationAngle, myJob.klCutoff);

    bool needTiltCoeffs(false);
    for (uint i = 0; i < diversityModes.size(); ++i) {
        uint16_t modeNumber = diversityModes[i];
        Cache::ModeID id2 = id;
        if (diversityTypes[i] == ZERNIKE) {
            id2.firstMode = id2.lastMode = 0;
        }
        id2.modeNumber = modeNumber;
        myJob.globalData->fetch (id2);
    }

    bfs::path fn = bfs::path("mvn_modes.f0");
    Array<double> debugModes;
    if ( bfs::exists (fn)) {
        redux::file::readFile(fn.string(), debugModes);
    }
    
    
    if ( !modeFile.empty() ) {
        modes.setPupilSize( pupilPixels, 0.0, 0.0 );
        modes.loadFile( modeFile );
    }
    
    if( modes.nPupilPixels != pupilPixels ) {
        modes.setPupilSize( pupilPixels, pupilRadiusInPixels, rotationAngle );
        modes.generate( myJob );
        // TODO generate
    }
    
    modeNumbers.clear();
    cnt = 0;

    if( needTiltCoeffs ) {
        id.firstMode = id.lastMode = 0;
        id.modeNumber = 3;
        const auto mode = myJob.globalData->fetch (id);
        double ddx = (*mode)(pupilPixels/2,pupilPixels/2+1) - (*mode)(pupilPixels/2,pupilPixels/2);
        pixelsToAlpha = util::pix2cf(arcSecsPerPixel,myJob.telescopeD)/(0.5*frequencyCutoff*ddx);
        alphaToPixels = 1.0/pixelsToAlpha;
        
    }
    */
    //defocusToAlpha = util::def2cf(myJob.telescopeD/2.0);
    //alphaToDefocus = 1.0/pixelsToAlpha;
    
}


void Channel::loadCalib(boost::asio::io_service& service) {     // load through cache-functionality, so the same data is recycled for other objects

    //LOG_TRACE << "Channel::loadCalib()";
    // TODO: absolute/relative paths
    // TODO: cache files and just fetch shared_ptr

    if (!darkTemplate.empty()) {        // needs to be read synchronously because of adding/normalization
        size_t nWild = std::count (darkTemplate.begin(), darkTemplate.end(), '%');
        if (nWild == 0 || darkNumbers.empty()) {
            ClippedFile::load(dark,darkTemplate,alignClip,true);
        } else {
            size_t nFrames = 1;
            Image<float> tmp;
            for (size_t di = 0; di < darkNumbers.size(); ++di) {
                bfs::path fn = bfs::path (boost::str (boost::format (darkTemplate) % darkNumbers[di]));
                if (!di) {
                    ClippedFile::load(dark,fn.string(),alignClip,true);
                } else {
                    ClippedFile::load(tmp,fn.string(),alignClip,true);
                    dark += tmp;
                    nFrames++;
                }
            }
            
            //   NOTE: no normalization here, modify metadata instead. TODO: += operator for meta to add frames??

        }

    }


    if (!gainFile.empty()) {
        service.post(std::bind(ClippedFile::load<float>, std::ref(gain), gainFile, alignClip, false/*norm*/, false/*symclip*/));
    }


    if (!responseFile.empty()) {
        service.post(std::bind(ClippedFile::load<float>, std::ref(ccdResponse), responseFile, alignClip, false/*norm*/, false/*symclip*/));
    }

    if (!backgainFile.empty()) {
        service.post(std::bind(ClippedFile::load<float>, std::ref(ccdScattering), backgainFile, alignClip, false/*norm*/, false/*symclip*/));
    }

    if (!psfFile.empty()) {
        service.post(std::bind(ClippedFile::load<float>, std::ref(psf), psfFile, alignClip, false/*norm*/, true/*symclip*/));
    }

    if (!mmFile.empty()) {
        service.post(std::bind(ClippedFile::load<float>, std::ref(modulationMatrix), mmFile, alignClip, false/*norm*/, false/*symclip*/));
    }

    if (!xOffsetFile.empty()) {
        service.post(std::bind(ClippedFile::load<int16_t>, std::ref(xOffset), xOffsetFile, alignClip, false/*norm*/, false/*symclip*/));
    }

    if (!yOffsetFile.empty()) {
        service.post(std::bind(ClippedFile::load<int16_t>, std::ref(yOffset), yOffsetFile, alignClip, false/*norm*/, false/*symclip*/));
    }
    
}


void Channel::loadData(boost::asio::io_service& service, Array<PatchData::Ptr>& patches) {

   //LOG_TRACE << "Channel::loadData()";
    size_t nImages = std::max<size_t>(1,imageNumbers.size());       // If no numbers, load template as single file

    startT = bpx::pos_infin;
    endT = bpx::neg_infin;
    
    // Prepare needed storage
    images.resize(nImages,imgSize.y,imgSize.x);
    imageStats.resize(nImages);
    service.post(std::bind(&Channel::adjustCutouts, this, std::ref(patches)));
    
    // load data (and preprocess)
    imageStats[0].reset(new ArrayStats());
    service.post(std::bind(&Channel::loadImage, this, 0));           // will *not* be loaded through the cache, so only saved in "images"
    for (size_t i = 1; i < nImages; ++i) {
        imageStats[i].reset (new ArrayStats());
        service.post(std::bind(&Channel::loadImage, this, i));       // will *not* be loaded through the cache, so only saved in "images"
    }
    runThreadsAndWait(service, myJob.info.maxThreads);
    
    patchWriteFail = std::async( launch::async, [this,&patches](){                      // launch as async to do writing in the background whil loading/pre-processing the rest.
        size_t nPatchesY = patches.dimSize(0);                                          // we will synchronize all the "patchWriteFail" futures at the end of MomfbdJob::preProcess
        size_t nPatchesX = patches.dimSize(1);
        for(uint py=0; py<nPatchesY; ++py) {
            for(uint px=0; px<nPatchesX; ++px) {
                ChannelData& chData(patches(py,px)->objects[myObject.ID].channels[ID]);
                copyImagesToPatch(chData);
                chData.cacheStore(true);    // store to disk and clear array
            }
        }
        images.clear();                                                                 // release resources after storing the patch-data.
        return false;   // TODO return true if any error
    });

    
    runThreadsAndWait(service, myJob.info.maxThreads);

}


void Channel::unloadCalib(void) {               // unload what was accessed through the cache, this should be called when all objects are done pre-processing.
    
    if (!darkTemplate.empty()) {
        size_t nWild = std::count (darkTemplate.begin(), darkTemplate.end(), '%');
        if (nWild == 0 || darkNumbers.empty()) {
            ClippedFile::unload<float>(darkTemplate,alignClip);
        } else {
            for (size_t di = 0; di < darkNumbers.size(); ++di) {
                bfs::path fn = bfs::path (boost::str (boost::format (darkTemplate) % darkNumbers[di]));
                ClippedFile::unload<float>(fn.string(),alignClip);
            }
        }
    }
    if( !gainFile.empty() ) ClippedFile::unload<float>(gainFile, alignClip);
    if( !responseFile.empty() ) ClippedFile::unload<float>(responseFile, alignClip);
    if( !backgainFile.empty() ) ClippedFile::unload<float>(backgainFile, alignClip);
    if( !psfFile.empty() ) ClippedFile::unload<float>(psfFile, alignClip, true);
    if( !mmFile.empty() ) ClippedFile::unload<float>(mmFile, alignClip);
    if( !xOffsetFile.empty() ) ClippedFile::unload<int16_t>(xOffsetFile, alignClip);
    if( !yOffsetFile.empty() ) ClippedFile::unload<int16_t>(yOffsetFile, alignClip);
    
    dark.clear();
    gain.clear();
    ccdResponse.clear();
    ccdScattering.clear();
    psf.clear();
    modulationMatrix.clear();
    xOffset.clear();
    yOffset.clear();

}


double Channel::getMaxMean(void) {
    double maxMean = std::numeric_limits<double>::lowest();
    for (auto &stat : imageStats) {
        maxMean = std::max(maxMean,stat->mean);
    }
    return maxMean;
}


void Channel::getFileNames(std::vector<std::string>& files) const {

    for (auto &num: imageNumbers) {
        bfs::path fn = bfs::path(imageDataDir)/bfs::path(boost::str(boost::format(imageTemplate) % num));
        files.push_back(fn.string());
    }
    
}


void Channel::initProcessing( const Solver& solver ) {
    
    initPhiFixed();

    size_t nImages = imageNumbers.size();
    subImages.resize(nImages);
    for (uint16_t i=0; i < nImages; ++i) {
        if( !subImages[i] ) subImages[i].reset( new SubImage(myObject, *this, solver.window, solver.noiseWindow) );
    }
    
}


void Channel::initPatch (ChannelData& cd) {

    size_t nImages = cd.images.dimSize();
    if (imageNumbers.empty() || imageNumbers.size() != nImages) {
        LOG_ERR << "Number of images in stack does not match the imageNumbers.  " << imageNumbers.size() << "!=" << printArray(cd.images.dimensions(),"dims");
        return;
    }


    if( (imageStats.size() == nImages) && (cd.images.nDimensions() == 3) ) {
        Array<float> view( cd.images, 0, 0, 0, cd.images.dimSize(1)-1, 0, cd.images.dimSize(2)-1 );
        for( uint16_t i=0; i<nImages; ++i ) {
            view *= (myObject.objMaxMean/imageStats[i]->mean);
            view.shift(0,1);
        }
    }
    uint16_t patchSize = myObject.patchSize;
    for (uint16_t i=0; i < nImages; ++i) {
        subImages[i]->setPatchInfo( i, cd.offset, cd.shift, patchSize, myObject.pupilPixels, myJob.modeNumbers.size() );
        subImages[i]->wrap( cd.images, i, i, cd.offset.y, cd.offset.y+patchSize-1, cd.offset.x, cd.offset.x+patchSize-1 );
        subImages[i]->stats.getStats( cd.images.ptr(i,0,0), cd.images.dimSize(1)*cd.images.dimSize(2), ST_VALUES|ST_RMS );
        subImages[i]->init();
        //subImages[i]->dump("o"+to_string(myObject.ID)+"_c"+to_string(ID)+"_im"+to_string(i));
    }
//cout << "Ch::initP   modes @ " << hexString(myObject.modes.get()) << "   uc=" << myObject.modes.use_count() << endl;
}


void Channel::initPhiFixed (void) {
    phi_fixed.resize (myObject.pupilPixels, myObject.pupilPixels);
    phi_fixed.zero();
    ModeInfo id (myJob.klMinMode, myJob.klMaxMode, 0, myObject.pupilPixels, myObject.pupilRadiusInPixels, rotationAngle, myJob.klCutoff);
    uint16_t modeNumber;
    for (uint i = 0; i < diversityModes.size(); ++i) {
        ModeInfo id2 = id;
        modeNumber = diversityModes[i];
        cout << "Channel::initPhiFixed()  i=" << i << endl;
        if (modeNumber == 2 || modeNumber == 3 || diversityTypes[i] == ZERNIKE) {
            id2.firstMode = id2.lastMode = 0;
        }
        id2.modeNumber = modeNumber;
        //const auto mode = myJob.globalData->fetch (id2);
        //redux::file::Ana::write ("mode_" + to_string (modeNumber) + "_" + to_string (i) + ".f0", *mode);
        //phi_fixed.add (*mode, diversity[i]);
        //redux::file::Ana::write ("phi-mode_" + to_string (modeNumber) + "_" + to_string (i) + ".f0", phi_fixed);
    }
    computePhi();   // no tilts for now, just initialize once
}


void Channel::computePhi (void) {
    //cout << "Channel::computePhi()" << endl;
    phi_fixed.copy(phi_channel);
    static int bla (0);
    if (diversityModes.size()) {
        redux::file::Ana::write ("phi_" + to_string (bla++) + ".f0", phi_channel);
    }
    // TODO: add tilt corrections
}

/*
void Channel::addMode (redux::util::Array<double>& phi, uint16_t modenumber, double weight) const {
    const auto mode = modes.at (modenumber);
    // cout << "Channel::addMode()  mode = " << modenumber << "  weight = " << weight << endl;
    //redux::file::Ana::write ("mode_" + to_string (modenumber) + ".f0", *mode);
    //redux::file::Ana::write ("pupil.f0", pupil.first);
    if (mode) {
        phi.add (*mode, weight);
    }
}


void Channel::getPhi (redux::util::Array<double>& phi, const WaveFront& wf) const {
    phi_channel.copy(phi);
    //return;
    //cout << "Channel::getPhi()" << endl;
    for (auto & it : wf.modes) {
        //cout << "Channel::getPhi()  it.first = " << it.first << endl;
        const auto mode = modes.at (it.first);
        if (mode) {  //&& it.second.second ) { // TODO: possibility to enable/disable modes
            //cout << "Channel::getPhi()  mode = " << hexString(mode.get()) << endl;
            phi.add (*mode, *it.second.value);
        }
    }
}
*/

void Channel::addAllFT (redux::util::Array<double>& ftsum) {
    for (auto& subimage : subImages) {
        subimage->addFT (ftsum);
    }
}

void Channel::addAllPQ(void) const {
    for (const auto& subimage : subImages) {
        myObject.addToPQ( subimage->imgFT, subimage->OTF );
    }
}

double Channel::metric (void) {

    double sum = 0.0;
//   for(auto &im: subImages) {
//       for( auto& a: im->wf->alpha) {
//           double coeff = a.second.first;
//           sum += coeff*coeff * modes.at(a.first)->inv_atm_rms;
//       }
//   }
    return sum;

}


void Channel::addTimeStamps( const bpx::ptime& newStart, const bpx::ptime& newEnd ) {

    unique_lock<mutex> lock(mtx);
    
    if(startT.is_special()) startT = newStart;
    else startT = std::min( startT, newStart );
    if(endT.is_special()) endT = newEnd;
    else endT = std::max( endT, newEnd );

}


void Channel::loadImage(size_t i) {
    
    bfs::path fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % imageNumbers[i]));
    ClippedFile cf(fn.string(),alignClip,false);
    
    Image<float> tmpImg;
    Image<float> view(images,i,i,0,imgSize.y-1,0,imgSize.x-1);
    
    LOG_DEBUG << boost::format ("Loading image %s  (%d:%d:%d  %s)") % fn % myObject.ID % ID % i % printArray(alignClip,"clip");
    cf.loadImage(tmpImg,false);
    if(tmpImg.meta) {
        addTimeStamps( tmpImg.meta->getStartTime(), tmpImg.meta->getEndTime() );
    }
    preprocessImage(i, tmpImg);
    imageStats[i]->getStats(borderClip, tmpImg);     // get stats for corrected data
    view.assign(reinterpret_cast<redux::util::Array<float>&>(tmpImg));
 
}


void Channel::preprocessImage( size_t index, Image<float>& img ) {

    Array<double> tmpImg;       // local temporary with double-precision
    tmpImg = img.copy<double>();
    
    /*
    // Michiel's method for detecting bitshifted Sarnoff images.    TODO make this an SST-specific "filter" to be applied on raw data
    if (imgMean > 4 * avgMean) {
        LOG_WARN << boost::format ("Image bit shift detected for image %s (mean > 4*avgMean). adjust factor=0.625 (keep your fingers crossed)!") % fn;
        tmpImg *= 0.625;
        modified = true;
    } else if (imgMean < 0.25 * avgMean) {
        LOG_WARN << boost::format ("Image bit shift detected for image %s (mean < 0.25*avgMean). adjust factor=16 (keep your fingers crossed)!") % fn;
        tmpImg *= 16;
        modified = true;
    }*/

    if (dark.valid() && gain.valid()) {
        if (! tmpImg.sameSize (dark)) {
            LOG_ERR << boost::format ("Dimensions of dark (%s) does not match this image (%s), skipping flatfielding !!") % printArray (dark.dimensions(), "") % printArray (tmpImg.dimensions(), "");
            return;
        }
        if (! tmpImg.sameSize (gain)) {
            LOG_ERR << boost::format ("Dimensions of gain (%s) does not match this image (%s), skipping flatfielding !!") % printArray (gain.dimensions(), "") % printArray (tmpImg.dimensions(), "");
            return;
        }
        if (ccdResponse.valid() && !tmpImg.sameSize (ccdResponse)) {
            LOG_WARN << boost::format ("Dimensions of ccd-response (%s) does not match this image (%s), will not be used !!") % printArray (ccdResponse.dimensions(), "") % printArray (tmpImg.dimensions(), "");
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
                LOG_DETAIL << "Applying correction for CCD transparency.";
                redux::image::descatter (tmpImg, ccdScattering, psf);
            } else {
                LOG_ERR << boost::format ("Dimensions of ccdScattering (%s) or psf (%s) does not match this image (%s), skipping flatfielding !!")
                        % printArray (ccdScattering.dimensions(), "") % printArray (psf.dimensions(), "") % printArray (tmpImg.dimensions(), "");
            }
        }

        tmpImg *= gain;

        namespace sp = std::placeholders;
        size_t sy = tmpImg.dimSize(0);
        size_t sx = tmpImg.dimSize(1);

        shared_ptr<double*> array = tmpImg.reshape(sy,sx);
        double** arrayPtr = array.get();
        switch (myJob.fillpixMethod) {
            case FPM_HORINT: {
                //LOG_DETAIL << "Filling bad pixels using horizontal interpolation.";
                function<double (size_t, size_t) > func = bind (horizontalInterpolation<double>, arrayPtr, sy, sx, sp::_1, sp::_2);
                fillPixels (arrayPtr, sy, sx, func, std::bind2nd (std::less_equal<double>(), myJob.badPixelThreshold));
                break;
            }
            case FPM_MEDIAN: {
                // TODO: median method
                break;
            }
            case FPM_INVDISTWEIGHT:       // inverse distance weighting is the default method, so fall through
            default: {
                //LOG_DETAIL << "Filling bad pixels using inverse distance weighted average.";
                function<double (size_t, size_t) > func = bind (inverseDistanceWeight<double>, arrayPtr, sy, sx, sp::_1, sp::_2);
                fillPixels (arrayPtr, sy, sx, func, std::bind2nd (std::less_equal<double>(), myJob.badPixelThreshold));
            }
        }

    }

    img.assign(tmpImg);         // copy back to image
  
    // save corrected data if specified.
    if((myObject.saveMask & SF_SAVE_FFDATA)) {
        bfs::path fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % imageNumbers[index]));
        fn = bfs::path(fn.leaf().string() + ".cor");
        LOG_DETAIL << boost::format ("Saving flat/dark corrected file %s") % fn.string();
        redux::file::Ana::write(fn.string(), tmpImg.copy<float>());   // TODO: other formats
    }


}


void Channel::copyImagesToPatch(ChannelData& chData) {
    
    size_t nImages = std::max<size_t>(1,imageNumbers.size());
    Array<float> block(reinterpret_cast<redux::util::Array<float>&>(images), 0, nImages-1, chData.cutout.first.y, chData.cutout.last.y, chData.cutout.first.x, chData.cutout.last.x);
    chData.images = std::move(block);       // chData.images will share datablock with the "images" stack, so minimal RAM usage.
    chData.setLoaded();
    
}


void Channel::adjustCutout(ChannelData& chData, const Region16& origCutout) const {

    RegionI desiredCutout = origCutout;
    desiredCutout.grow(myObject.maxLocalShift);

    chData.cutout = desiredCutout;
    chData.offset = 0;
    chData.residualOffset = 0;
    
    if( imgSize == 0 ) {
        LOG_ERR << "No valid imgSize when adjusting cutout, that should not happen...";
        return;
    }
    
    RegionI imgBoundary(0,0,imgSize.y-1, imgSize.x-1);
    chData.cutout.restrict(imgBoundary);

    if (xOffset.valid()) {
        ArrayStats stats;
        stats.getStats (Image<int16_t> (xOffset, chData.cutout.first.y, chData.cutout.last.y, chData.cutout.first.x, chData.cutout.last.x), ST_VALUES);
        chData.residualOffset.x = stats.mean/100.0;
    }
    if (yOffset.valid()) {
        ArrayStats stats;
        stats.getStats (Image<int16_t> (yOffset, chData.cutout.first.y, chData.cutout.last.y, chData.cutout.first.x, chData.cutout.last.x), ST_VALUES);
        chData.residualOffset.y = stats.mean/100.0;
    }

    chData.shift = PointI(lround(chData.residualOffset.y),lround(chData.residualOffset.x));     // possibly apply shift from offsetfiles.
    chData.shift -= imgBoundary.outside(chData.cutout+chData.shift);                            // restrict the shift inside the image, leave the rest in "residualOffset" to be dealt with using Zernike tilts.
    chData.cutout += chData.shift;
    chData.residualOffset -= chData.shift;
    chData.offset = myObject.maxLocalShift;
    if (chData.cutout != desiredCutout) {
        chData.offset -= (chData.cutout.first - desiredCutout.first - chData.shift);
    }
    
}


void Channel::adjustCutouts(Array<PatchData::Ptr>& patches) {
    
    if( patches.nDimensions() == 2 ) {
        size_t nPatchesY = patches.dimSize(0);
        size_t nPatchesX = patches.dimSize(1);
        for(uint py=0; py<nPatchesY; ++py) {
            for(uint px=0; px<nPatchesX; ++px) {
                PatchData& patch(*patches(py,px));
                ChannelData& chData = patch.objects[myObject.ID].channels[ID];
                adjustCutout(chData,patch.roi);
            }
        }
    }
    
}


void Channel::storePatchData(boost::asio::io_service& service, Array<PatchData::Ptr>& patches) {
    
    if( patches.nDimensions() == 2 ) {
        size_t nPatchesY = patches.dimSize(0);
        size_t nPatchesX = patches.dimSize(1);
        cout << "storing (" << nPatchesY << "x" << nPatchesX << ") patches..." << endl;
        for(uint py=0; py<nPatchesY; ++py) {
            for(uint px=0; px<nPatchesX; ++px) {
                PatchData& patch(*patches(py,px));
                ChannelData& chData = patch.objects[myObject.ID].channels[ID];
                service.post(std::bind(&ChannelData::cacheStore,std::ref(chData),true));
            }
        }
    }
    
}

Point16 Channel::getImageSize(void) {

    if( imgSize == 0 ) {
        if( alignClip.size() == 4 ) {
            imgSize = Point16(abs(alignClip[3]-alignClip[2])+1, abs(alignClip[1]-alignClip[0])+1);
            return imgSize;
        }
        size_t nImages = imageNumbers.size();
        bfs::path fn;
        if(nImages) {
            fn = bfs::path(imageDataDir) / bfs::path(boost::str (boost::format (imageTemplate) % imageNumbers[0]));
        } else {
            nImages = 1;
            fn = bfs::path(imageDataDir) / bfs::path(imageTemplate);
        }
        Image<float> tmp;
        ClippedFile::load(tmp,fn.string(),alignClip);       // Image will be cached, so this is not a "wasted" load.
        if( tmp.nDimensions() == 2 ) {
            imgSize = Point16(tmp.dimSize(0),tmp.dimSize(1));
        }
    }
    return imgSize;
 
}


void Channel::logAndThrow( string msg ) {
    msg = to_string(myObject.ID)+":"+to_string(ID)+": "+msg;
    LOG_ERR << msg;
    throw job_check_failed(msg);
    
}


void Channel::dump (std::string tag) {
    
    tag += "_c"+to_string(ID);
    
    Ana::write (tag + "_phi_fixed.f0", phi_fixed);
    Ana::write (tag + "_phi_channel.f0", phi_channel);
    int cnt(0);
    for( auto& im: subImages ) {
        im->dump(tag+"_im"+to_string(cnt++));
    }

}

