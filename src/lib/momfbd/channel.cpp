#include "redux/momfbd/channel.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/object.hpp"
#include "redux/momfbd/subimage.hpp"
#include "redux/momfbd/util.hpp"
#include "redux/momfbd/wavefront.hpp"

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

    bool checkImageScale (double& F, double& A, double& P) {

        double rad2asec = 180.0 * 3600.0 / redux::PI;
        size_t count = F > 0 ? 1 : 0;
        count += A > 0 ? 1 : 0;
        count += P > 0 ? 1 : 0;
        if (count > 2) {
            LOG_WARN << "Too many parameters specified: replacing telescope focal length (" << F
                     << ") with computed value (" << (P * rad2asec / A) << ")";
            F = P * rad2asec / A;
            return true;
        } else if (count < 2) {
            LOG_ERR << "At least two of the parameters \"TELESCOPE_F\", \"ARCSECPERPIX\" and \"PIXELSIZE\" has to be provided.";
        } else {    // count == 2
            if (F <= 0) {
                F = P * rad2asec / A;
            } else if (A <= 0) {
                A = P * rad2asec / F;
            } else if (P <= 0) {
                P = F * A / rad2asec;
            }
            return true;
        }
        return false;
    }

    void calculatePupilSize (double &frequencyCutoff, double &pupilRadiusInPixels, uint16_t &nPupilPixels, double wavelength, uint32_t nPixels, double telescopeDiameter, double arcSecsPerPixel) {
        static double radians_per_arcsec = redux::PI / (180.0 * 3600.0);         // (2.0*redux::PI)/(360.0*3600.0)
        double radians_per_pixel = arcSecsPerPixel * radians_per_arcsec;
        double q_number = wavelength / (radians_per_pixel * telescopeDiameter);
        frequencyCutoff = (double) nPixels / q_number;
        nPupilPixels = nPixels >> 2;
        pupilRadiusInPixels = frequencyCutoff / 2.0;                   // telescope radius in pupil pixels...
        if (nPupilPixels < pupilRadiusInPixels) {            // this should only be needed for oversampled images
            uint16_t goodsizes[] = { 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45, 48, 50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128, 135, 144 };
            for (int i = 0; (nPupilPixels = max (goodsizes[i], nPupilPixels)) < pupilRadiusInPixels; ++i);     // find right size
        }
        nPupilPixels <<= 1;
    }

    struct ClippedFile {
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
            
            LOG_DEBUG << boost::format ("Loaded file \"%s\"  (%s)") % filename % printArray(clip,"clip");
 
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
                for (auto & it : alignClip) --it;       // NOTE: momfbd cfg files uses 1-based indexes, internally we start with 0.
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
            img.resize();
            Image<T>& cimg = redux::util::Cache::get<ClippedFile, Image<T> >(cf,img);
            {
                unique_lock<mutex> lock(cimg.imgMutex);
                if(cimg.nElements() == 0) cf.loadImage(cimg,norm);
            }
            img = cimg;
        }

        template <typename T>
        static void unload(const string& fn, const vector<int16_t>& cl, bool sym=false) {
            ClippedFile cf(fn,cl,sym);
            redux::util::Cache::erase<ClippedFile, Image<T> >(cf);
        }

        string filename;
        const vector<int16_t> clip;
        bool symmetricClip;
    };

}


Channel::Channel (Object& o, MomfbdJob& j, uint16_t id) : ID (id), myObject (o), myJob (j) {

}

Channel::~Channel() {

}

void Channel::parsePropertyTree (bpt::ptree& tree) {

    ChannelCfg::parseProperties (tree, myObject);

    /*if( !alignClip.empty() && alignClip.size() != 4 ) {
        LOG_ERR << "argument to ALIGN_CLIP could not be translated to 4 integers. Whole image area will be used !!";
        alignClip.clear();
    }*/

    /*if( imageDataDir.length() == 0 ) {
        imageDataDir = cleanPath( "./" );      // Nothing specified in cfg file, use current directory.
    }*/

    if (imageTemplate.length() == 0) {
        LOG_ERR << "no filename template specified.";
    }

    /*    if( darkTemplate.length() > 0 ) {
            if( darkNumbers.size() == 0 ) {
                LOG_ERR << "darkfield template specified but no dark numbers.";
            }
        }
        else if( darkNumbers.size() > 0 ) {
            LOG_ERR << "darkfield dark numbers specified but no darkfield template.";
        }
    */
    if (gainFile.length() > 0) {
        if (darkTemplate.length() == 0) {
            LOG_ERR << "a gain file name but no dark field was specified.";
        }
    } else if (darkTemplate.length() > 0) {
        LOG_ERR << "a dark field name but no gain file was specified.";
    }

    if ( (responseFile.length() > 0) && (gainFile.length() == 0)) {
        LOG_ERR << "detector response correction only possible when flatfielding.";
    }



    //MomfbdJob::maybeOverride( tree.get<bool>( "NO_RESTORE", flags & MFBD_NO_RESTORE ), flags, MFBD_NO_RESTORE );

    size_t p;
    if ( (p = imageTemplate.find_first_of ('%')) != string::npos) {
        /*     if( sequenceNumber > 0 ) {
                 size_t q;
                 if( ( q = imageTemplate.find_first_of( '%', p + 1 ) ) != string::npos ) {
                     string tmpString = boost::str( boost::format( imageTemplate.substr( 0, q ) ) % sequenceNumber );
                     imageTemplate = tmpString + imageTemplate.substr( q );
                 }
                 else  LOG_WARN << boost::format( "file name template %s does not contain a 2nd format specifier (needs 2)" ) % imageTemplate;
             }*/
    } else {
        //LOG_WARN << boost::format( "file name template %s does not contain a format specifier (needs %d)" ) % imageTemplate % ( 1 + ( sequenceNumber >= 0 ) );
    }


    //LOG_DEBUG << "Channel::parseProperties() done.";

}



bpt::ptree Channel::getPropertyTree (bpt::ptree& tree) {

    bpt::ptree node;

    /*
        if( sequenceNumber != myObject.sequenceNumber ) node.put( "SEQUENCE_NUM", sequenceNumber );
        if( mmRow ) node.put( "MMROW", mmRow );
        if( mmWidth ) node.put( "MMWIDTH", mmWidth );
        if( stokesWeights != myObject.stokesWeights ) node.put( "VECTOR", stokesWeights );
        if( !diversity.empty() ) node.put( "DIVERSITY", diversity );
        if( !diversityOrders.empty() ) node.put( "DIV_ORDERS", diversityOrders ); // TODO types missing
        if( image_num_offs != 0 ) node.put( "DT", image_num_offs );
        uint32_t dflags = flags ^ myObject.flags;
        if( dflags & MFBD_NO_RESTORE ) node.put( "NO_RESTORE", ( bool )( flags & MFBD_NO_RESTORE ) );
        if( dflags & MFBD_SAVE_FFDATA ) node.put( "SAVE_FFDATA", ( bool )( flags & MFBD_SAVE_FFDATA ) );
        // TODO "INCOMPLETE"
    */
    ChannelCfg::getProperties (node, myObject);

    tree.push_back (bpt::ptree::value_type ("channel", node));

    return node;

}


size_t Channel::size (void) const {

    size_t sz = ChannelCfg::size();
    sz += sizeof (uint16_t);        // ID;
    //sz += sizeof (uint32_t);        // dataOffset;
    sz += dark.size();
    sz += imageStats.size() * ArrayStats::size() + sizeof (uint16_t);
    return sz;
}


uint64_t Channel::pack (char* ptr) const {
    using redux::util::pack;
    uint64_t count = ChannelCfg::pack (ptr);
    count += pack (ptr + count, ID);
    //count += pack (ptr + count, dataOffset);
    count += dark.pack (ptr + count);
    uint16_t statSize = imageStats.size();
    count += pack (ptr + count, statSize);
    for (auto & it : imageStats) count += it->pack (ptr + count);
    if (count != size()) {
        LOG_ERR << "(" << hexString (this) << "): Packing failed, there is a size mismatch:  count = " << count << "  sz = " << size();
    }
    return count;
}


uint64_t Channel::unpack (const char* ptr, bool swap_endian) {
    using redux::util::unpack;

    uint64_t count = ChannelCfg::unpack (ptr, swap_endian);
    count += unpack (ptr + count, ID, swap_endian);
    //count += unpack (ptr + count, dataOffset, swap_endian);
    count += dark.unpack (ptr + count, swap_endian);
    uint16_t statSize;
    count += unpack (ptr + count, statSize, swap_endian);
    imageStats.resize (statSize);
    for (auto & it : imageStats) {
        it.reset (new ArrayStats());
        count += it->unpack (ptr + count, swap_endian);
    }
    return count;
}


bool Channel::checkCfg (void) {

    if (!checkImageScale (telescopeF, arcSecsPerPixel, pixelSize)) {
        return false;
    }

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

    // Do we have a correct dark template ?
    if (darkTemplate.empty()) {
        LOG_ERR << "No filename template specified.";
        return false;
    }
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

    return true;

}


bool Channel::checkData (void) {

    // Images
    if (incomplete) {   // check if files are present
        for (size_t i (0); i < imageNumbers.size();) {
            bfs::path fn = bfs::path (boost::str (boost::format (imageTemplate) % (imageNumberOffset + imageNumbers[i])));
            if (!bfs::exists (fn)) {
                fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % (imageNumberOffset + imageNumbers[i])));
                if (!bfs::exists (fn)) {
                    LOG_TRACE << "File not found: \"" << fn.string() << "\", removing from list of image numbers.";
                    imageNumbers.erase (imageNumbers.begin() + i);
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
        for (auto & it : imageNumbers) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % (imageNumberOffset + it)));
            if (!bfs::exists (fn)) {
                LOG_ERR << boost::format ("Image-file %s not found!") % boost::str (boost::format (imageTemplate) % (imageNumberOffset + it));
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
                LOG_ERR << boost::format ("Dark-file %s not found!") % darkTemplate;
                return false;
            } else darkTemplate = fn.c_str();
        }
    } else {                            // template
        for (auto & it : darkNumbers) {
            bfs::path fn = bfs::path (boost::str (boost::format (darkTemplate) % it));
            if (!bfs::exists (fn)) {
                fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (darkTemplate) % it));
                if (!bfs::exists (fn)) {
                    LOG_ERR << boost::format ("Dark-file %s not found!") % boost::str (boost::format (darkTemplate) % it);
                    return false;
                } else darkTemplate = fn.c_str();
            }
        }
    }


    // Gain
    if (!gainFile.empty()) {
        if (! bfs::exists (bfs::path (gainFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (gainFile);
            if (! bfs::exists (fn)) {
                LOG_ERR << boost::format ("Gain-file %s not found!") % gainFile;
                return false;
            } else gainFile = fn.c_str();
        }
    }

    if (!responseFile.empty()) {
        if (! bfs::exists (bfs::path (responseFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (responseFile);
            if (! bfs::exists (fn)) {
                LOG_ERR << boost::format ("Response-file %s not found!") % responseFile;
                return false;
            } else responseFile = fn.c_str();
        }
    }

    if (!backgainFile.empty()) {
        if (! bfs::exists (bfs::path (backgainFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (backgainFile);
            if (! bfs::exists (fn)) {
                LOG_ERR << boost::format ("Backgain-file %s not found!") % backgainFile;
                return false;
            } else backgainFile = fn.c_str();
        }
    }

    if (!psfFile.empty()) {
        if (! bfs::exists (bfs::path (psfFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (psfFile);
            if (! bfs::exists (fn)) {
                LOG_ERR << boost::format ("PSF-file %s not found!") % psfFile;
                return false;
            } else psfFile = fn.c_str();
        }
    }

    if (!mmFile.empty()) {
        if (! bfs::exists (bfs::path (mmFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (mmFile);
            if (! bfs::exists (fn)) {
                LOG_ERR << boost::format ("Modulation-matrix file %s not found!") % mmFile;
                return false;
            } else mmFile = fn.c_str();
        }
    }

    if (!xOffsetFile.empty()) {
        if (! bfs::exists (bfs::path (xOffsetFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (xOffsetFile);
            if (! bfs::exists (fn)) {
                LOG_ERR << boost::format ("Offset-file %s not found!") % xOffsetFile;
                return false;
            } else xOffsetFile = fn.c_str();
        }
    }

    if (!yOffsetFile.empty()) {
        if (! bfs::exists (bfs::path (yOffsetFile))) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (yOffsetFile);
            if (! bfs::exists (fn)) {
                LOG_ERR << boost::format ("Offset-file %s not found!") % yOffsetFile;
                return false;
            } else yOffsetFile = fn.c_str();
        }
    }

    return true;
}

void Channel::init (void) {

}

//#define INDEX_THRESHOLD  0
#define INDEX_THRESHOLD  1E-12      // cutoff to get rid of some fft noise
// TBD: this should be a parameter somewhere...
void Channel::initCache (void) {

    //LOG_DETAIL << "wavelength = " << myObject.wavelength << "   patchSize = " << patchSize << "  telescopeD = " << myJob.telescopeD << "  arcSecsPerPixel = " << arcSecsPerPixel;
    calculatePupilSize (frequencyCutoff, pupilRadiusInPixels, pupilPixels, myObject.wavelength, patchSize, myJob.telescopeD, arcSecsPerPixel);
    
    myJob.patchSize = myObject.patchSize = patchSize;     // TODO: fulhack until per-channel sizes is implemented
    myJob.pupilPixels = myObject.pupilPixels = pupilPixels;
    size_t otfPixels = 2 * pupilPixels;
    //LOG_DETAIL << "frequencyCutoff = " << frequencyCutoff << "  pupilSize = " << pupilPixels << "  pupilRadiusInPixels = " << pupilRadiusInPixels;
    pupil = myJob.globalData->fetch (pupilPixels, pupilRadiusInPixels);
    // Create a temporary OTF and store the indices where the OTF/pupil are non-zero. This will be used in loops to skip irrelevant evaluations.
    Array<double> tmpImg (2 * pupilPixels, 2 * pupilPixels);
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

    for (uint16_t & it : myJob.modeNumbers) {
        Cache::ModeID id2 = id;
        if(it == 2 || it == 3) needTiltCoeffs = true;
        if (myJob.modeBasis == ZERNIKE || it == 2 || it == 3) {     // force use of Zernike modes for all tilts
            id2.firstMode = id2.lastMode = 0;
        }
        id2.modeNumber = it;
        const PupilMode::Ptr mode = myJob.globalData->fetch (id2);
        modes.emplace (it, myJob.globalData->fetch (id2));
    }

    if( needTiltCoeffs ) {
        id.firstMode = id.lastMode = 0;
        id.modeNumber = 2;
        const PupilMode::Ptr mode = myJob.globalData->fetch (id);
        double dxdp = (*mode)(pupilPixels/2,pupilPixels/2+1) - (*mode)(pupilPixels/2,pupilPixels/2);
        pixelsToAlpha = util::pix2cf(arcSecsPerPixel,myJob.telescopeD)/(0.5*frequencyCutoff*dxdp);
        alphaToPixels = 1.0/pixelsToAlpha;
    }
    
    defocusToAlpha = util::def2cf(myJob.telescopeD/2.0);
    alphaToDefocus = 1.0/pixelsToAlpha;
    
}


void Channel::cleanup (void) {

}



namespace {
    template <typename T>
    void loadWrapper (const string& fn, T& img) {
        redux::file::readFile (fn, img);
        LOG_DETAIL << boost::format ("Loaded file \"%s\"") % fn;
    }
}


void Channel::loadData (boost::asio::io_service& service) {

    LOG_TRACE << "Channel::loadData()";
    // TODO: absolute/relative paths
    // TODO: cache files and just fetch shared_ptr

    //ClippedFile tmpFile(alignClip,"",false);
    
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

    if (!pupilFile.empty()) {
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

    size_t nImages = std::max<size_t>(1,imageNumbers.size());       // If no numbers, load template as single file
    
    images.resize(nImages);
    imageStats.resize(nImages);

    if (imageNumbers.empty()) {
        imageStats[0].reset (new ArrayStats());
        service.post( [this](){
                    bfs::path fn = bfs::path(imageDataDir) / bfs::path(imageTemplate);
                    ClippedFile::load(images[0],fn.string(),alignClip);
                    imageStats[0]->getStats(myJob.borderClip, images[0], ST_VALUES);  // only get min/max/mean
                   });
    } else {
        for (size_t i = 0; i < nImages; ++i) {
            imageStats[i].reset (new ArrayStats());
            service.post( [this,i](){
                        bfs::path fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % imageNumbers[i]));
                        ClippedFile::load(images[i],fn.string(),alignClip);
                        imageStats[i]->getStats(myJob.borderClip, images[i], ST_VALUES);  // only get min/max/mean
                       });
        }
    }


}


void Channel::unloadData(void) {
    
    size_t nImages = imageNumbers.size();
    if (nImages) {
        for (size_t i = 0; i < nImages; ++i) {
            bfs::path fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % imageNumbers[i]));
            ClippedFile::unload<float>(fn.string(),alignClip);
        }
    } else  {
        bfs::path fn = bfs::path(imageDataDir) / bfs::path(imageTemplate);
        ClippedFile::unload<float>(fn.string(),alignClip);
    }
    images.clear();
  
}


void Channel::unloadCalib(void) {
    
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
    ClippedFile::unload<float>(gainFile, alignClip);
    ClippedFile::unload<float>(responseFile, alignClip);
    ClippedFile::unload<float>(backgainFile, alignClip);
    ClippedFile::unload<float>(psfFile, alignClip, true);
    ClippedFile::unload<float>(mmFile, alignClip);
    ClippedFile::unload<int16_t>(xOffsetFile, alignClip);
    ClippedFile::unload<int16_t>(yOffsetFile, alignClip);
    dark.clear();
    gain.clear();
    ccdResponse.clear();
    ccdScattering.clear();
    psf.clear();
    modulationMatrix.clear();
    xOffset.clear();
    yOffset.clear();

}

void Channel::preprocessData (boost::asio::io_service& service) {

    size_t nImages = images.size();
    if( nImages ) {
        double avgMean = 0.0;
        startT = bpx::pos_infin;
        endT = bpx::neg_infin;
        for (size_t i = 0; i < nImages; ++i) {
            avgMean += imageStats[i]->mean;
            if(images[i].meta) {
                if(startT.is_special()) startT = images[i].meta->getStartTime();
                else startT = std::min(startT,images[i].meta->getStartTime());
                if(endT.is_special()) endT = images[i].meta->getEndTime();
                else endT = std::max(endT,images[i].meta->getEndTime());
            }
        }
        avgMean /= static_cast<double> (nImages);
        for (size_t i = 0; i < nImages; ++i) {
            service.post(std::bind (&Channel::preprocessImage, this, i, avgMean));
        }
    }

}


double Channel::getMaxMean (void) const {
    double maxMean = std::numeric_limits<double>::lowest();
    for (ArrayStats::Ptr imStat : imageStats) {
        if (imStat->mean > maxMean) maxMean = imStat->mean;
    }
    return maxMean;
}


void Channel::getFileNames(std::vector<std::string>& files) const {

    for (auto &num: imageNumbers) {
        bfs::path fn = bfs::path(imageDataDir)/bfs::path(boost::str(boost::format(imageTemplate) % num));
        files.push_back(fn.string());
    }
    
}


void Channel::initProcessing (WorkSpace::Ptr ws) {
    workspace = ws;
    initCache();        // this will initialize modes & pupil for this channel
    initPhiFixed();
    for (auto & m : modes) {           // check if the tilts are present
        if (m.first == 2 || m.first == 3) {
            shared_ptr<Tilts> tilts (new Tilts (*this, m.first));
            auto ret = workspace->tilts.emplace (m.first, tilts);
            if (ret.second) {       // new tiltmode, i.e. this will be the reference channel
                tilts->nFreeAlpha = imageNumbers.size();
                tilts->anchorChannel = true;
            } else {                // existing tiltmode, this channel is added as a "relative tilt"
                ret.first->second->addRelativeTilt (tilts);
            }
            tilts->init();
        }
    }

    for (uint16_t i = 0; i < imageNumbers.size(); ++i) {
        uint32_t imageNumber = imageNumbers[i];
        std::shared_ptr<WaveFront>& wf = workspace->wavefronts[imageNumber];
        if (!wf) wf.reset (new WaveFront());
        for (auto & m : modes) {
            if (m.first > 3) {      // tilts are treated separately
                wf->addWeight (m.first, 1.0 / (m.second->atm_rms * myObject.wavelength * myObject.wavelength));
            }
        }
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
    
    subImages.clear();
    for (uint16_t i=0; i < nImages; ++i) {
        uint32_t imageNumber = imageNumbers[i];
        std::shared_ptr<SubImage> simg (new SubImage (myObject, *this, workspace->window, workspace->noiseWindow, cd.images, i, cd.offset,
                                        patchSize, pupilPixels));   // TODO: fix offsets
        subImages.push_back (simg);
        simg->init();
        std::shared_ptr<WaveFront>& wf = workspace->wavefronts[imageNumber];
        if (!wf) cout << "Channel::initPatch(): wf = NULL for imageNumber = " << imageNumber << endl;
        wf->addImage (simg);
    }
}


void Channel::getResults (ChannelData& cd) {
    
}


void Channel::writeAna (const redux::util::Array<PatchData::Ptr>& patches) {

     
}



void Channel::initPhiFixed (void) {
    phi_fixed.resize (pupilPixels, pupilPixels);
    phi_fixed.zero();
    Cache::ModeID id (myJob.klMinMode, myJob.klMaxMode, 0, pupilPixels, pupilRadiusInPixels, rotationAngle, myJob.klCutoff);
    uint16_t modeNumber;
    for (uint i = 0; i < diversityModes.size(); ++i) {
        Cache::ModeID id2 = id;
        modeNumber = diversityModes[i];
        cout << "Channel::initPhiFixed()  i=" << i << endl;
        if (modeNumber == 2 || modeNumber == 3 || diversityTypes[i] == ZERNIKE) {
            id2.firstMode = id2.lastMode = 0;
        }
        id2.modeNumber = modeNumber;
        const PupilMode::Ptr mode = myJob.globalData->fetch (id2);
        //redux::file::Ana::write ("mode_" + to_string (modeNumber) + "_" + to_string (i) + ".f0", *mode);
        phi_fixed.add (*mode, diversity[i]);
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


void Channel::addMode (redux::util::Array<double>& phi, uint16_t modenumber, double weight) const {
    const PupilMode::Ptr mode = modes.at (modenumber);
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
        const PupilMode::Ptr mode = modes.at (it.first);
        if (mode) {  //&& it.second.second ) { // TODO: possibility to enable/disable modes
            //cout << "Channel::getPhi()  mode = " << hexString(mode.get()) << endl;
            phi.add (*mode, *it.second.value);
        }
    }
}


void Channel::addAllFT (redux::util::Array<double>& ftsum) {
    for (shared_ptr<SubImage>& it : subImages) {
        it->addFT (ftsum);
    }
}

double Channel::metric (void) {

    double sum = 0.0;
//   for(shared_ptr<SubImage> &im: subImages) {
//       for( auto& a: im->wf->alpha) {
//           double coeff = a.second.first;
//           sum += coeff*coeff * modes.at(a.first)->inv_atm_rms;
//       }
//   }
    return sum;

}


void Channel::normalizeData (boost::asio::io_service& service, double value) {
    size_t nImages = imageNumbers.size();
    for (size_t i = 0; i < nImages; ++i) {
        service.post (std::bind (&Channel::normalizeImage, this, i, value));
    }
}


void Channel::loadImage (size_t index) {
/*    Image<double> subimg (images, index, index, 0, images.dimSize (1) - 1, 0, images.dimSize (2) - 1);
    bfs::path fn = bfs::path (imageDataDir) / bfs::path (boost::str (boost::format (imageTemplate) % imageNumbers[index]));
    redux::file::readFile (fn.string(), subimg);
    LOG_DETAIL << boost::format ("Loaded file %s") % fn;
    imageStats[index]->getStats (myJob.borderClip, subimg, ST_VALUES);                       // only get min/max/mean
*/}


void Channel::preprocessImage (size_t index, double avgMean) {

    double imgMean = imageStats[index]->mean;
    bool modified = false;

    Array<double> tmpImg;
    tmpImg = images[index].copy<double>();
    //images[index].copy(tmpImg);
    bfs::path fn = bfs::path (boost::str (boost::format (imageTemplate) % imageNumbers[index]));
    LOG_TRACE << boost::format ("Pre-processing image %s") % fn;
    redux::file::Ana::write(fn.string()+"_raw.f0", images[index]);
    redux::file::Ana::write(fn.string()+"_copy.f0", tmpImg);
    // Michiel's method for detecting bitshifted Sarnoff images.
    if (imgMean > 4 * avgMean) {
        LOG_WARN << boost::format ("Image bit shift detected for image %s (mean > 4*avgMean). adjust factor=0.625 (keep your fingers crossed)!") % fn;
        tmpImg *= 0.625;
        modified = true;
    } else if (imgMean < 0.25 * avgMean) {
        LOG_WARN << boost::format ("Image bit shift detected for image %s (mean < 0.25*avgMean). adjust factor=16 (keep your fingers crossed)!") % fn;
        tmpImg *= 16;
        modified = true;
    }

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
        //tmpImg -= reinterpret_cast<const redux::util::Array<float>&>(dark);
        modified = true;
    redux::file::Ana::write(fn.string()+"_dark.f0", dark);
    redux::file::Ana::write(fn.string()+"_dark2.f0", reinterpret_cast<redux::util::Array<float>&>(dark));
    redux::file::Ana::write(fn.string()+"_darked.f0", tmpImg);
        
        if (ccdResponse.valid()) {   // correct for the detector response (this should not contain the gain correction and must be done before descattering)
            //tmpImg *= ccdResponse;
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
        //tmpImg *= gain;
        //tmpImg *= reinterpret_cast<const redux::util::Array<float>&>(gain);
    redux::file::Ana::write(fn.string()+"_gain.f0", gain);
    redux::file::Ana::write(fn.string()+"_gained.f0", tmpImg);

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
    redux::file::Ana::write(fn.string()+"_filled.f0", tmpImg);

    imageStats[index]->getStats(myJob.borderClip, tmpImg);     // get stats for corrected data
    tmpImg.copy(images[index]);
    
    if (modified && (myObject.saveMask & SF_SAVE_FFDATA)) {
        fn = bfs::path (fn.leaf().string() + ".cor");
        LOG_DETAIL << boost::format ("Saving flat/dark corrected file %s") % fn.string();
        redux::file::Ana::write (fn.string(), images[index]);   // TODO: other formats
    }

}


void Channel::normalizeImage (size_t index, double value) {
    images[index] *= (value / imageStats[index]->mean);
    double noise1 = imageStats[index]->noise;
    imageStats[index]->getStats(images[index]);
    LOG_TRACE << boost::format("Normalizing image (%d,%d,%d):  noise1 = %f   noise2 = %f") % myObject.ID % ID % (index) % noise1 % imageStats[index]->noise;
}


size_t Channel::sizeOfPatch (uint32_t npixels) const {
    size_t sz = sizeof (size_t) + imageStats.size() * sizeof (float);
    sz += npixels * images.size() * sizeof (float);
    return sz;
}


void Channel::getPatchData (ChannelData& chData, const PatchData& patch) const {

    size_t nImages = images.size();
    chData.offset = 0;
    chData.residualOffset = 0;
    
    if (nImages && (images[0].nDimensions() == 2)) {
        RegionI imgBoundary(0,0,images[0].dimSize(0)-1,images[0].dimSize(1)-1);
        RegionI desiredCutout = patch.roi;
        desiredCutout.grow(maxLocalShift);
        RegionI actualCutout = desiredCutout;
        actualCutout.restrict(imgBoundary);
        if (xOffset.valid()) {
            ArrayStats stats;
            stats.getStats (Image<int16_t> (xOffset, actualCutout.first.y, actualCutout.last.y, actualCutout.first.x, actualCutout.last.x), ST_VALUES);
            chData.residualOffset.x = stats.mean/100.0;
        }
        if (yOffset.valid()) {
            ArrayStats stats;
            stats.getStats (Image<int16_t> (yOffset, actualCutout.first.y, actualCutout.last.y, actualCutout.first.x, actualCutout.last.x), ST_VALUES);
            chData.residualOffset.y = stats.mean/100.0;
        }

        chData.shift = PointI(lround(chData.residualOffset.y),lround(chData.residualOffset.x));
        chData.shift -= imgBoundary.outside(actualCutout+chData.shift);             // restrict the shift inside the image, leave the rest in "residualOffset" to be dealt with using Zernike tilts.
        actualCutout += chData.shift;
        chData.residualOffset -= chData.shift;
        chData.offset = maxLocalShift;
        if (actualCutout != desiredCutout) {
            chData.offset -= (actualCutout.first - desiredCutout.first - chData.shift);
        }

        PointI pSize = actualCutout.last - actualCutout.first + 1;
        chData.images.resize(nImages, pSize.y, pSize.x);
        Array<float> patchImg(chData.images, 0, 0, 0, pSize.y-1, 0, pSize.x-1);
        for( uint i=0; i<nImages; ++i) {
            Image<float> tmpImage(images[i], actualCutout.first.y, actualCutout.last.y, actualCutout.first.x, actualCutout.last.x);
            patchImg.assign(reinterpret_cast<const redux::util::Array<float>&>(tmpImage));
            patchImg.shift(0,1);
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


void Channel::dump (std::string tag) {

    Ana::write (tag + "_pupil.f0", pupil.first);
    Ana::write (tag + "_fittedplane.f0", fittedPlane);
    Ana::write (tag + "_phi_fixed.f0", phi_fixed);
    Ana::write (tag + "_phi_channel.f0", phi_channel);

}

