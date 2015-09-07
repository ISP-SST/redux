#include "redux/momfbd/object.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/data.hpp"
#include "redux/momfbd/util.hpp"

#include "redux/file/filemomfbd.hpp"
#include "redux/math/functions.hpp"
#include "redux/translators.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/constants.hpp"
#include "redux/logger.hpp"
#include "redux/revision.hpp"

#include "redux/file/fileana.hpp"
#include <cstdio>
#include <limits>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/info_parser.hpp>

namespace bfs = boost::filesystem;
using namespace redux::momfbd;
using namespace redux::file;
using namespace redux::math;
using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;
using boost::algorithm::iequals;

#define lg Logger::mlg
namespace {

    const string thisChannel = "object";

}


Object::Object (MomfbdJob& j, uint16_t id) : ObjectCfg (j), myJob (j), ID (id), nObjectImages(0) {


}


Object::~Object() {

}


void Object::parsePropertyTree( bpt::ptree& tree ) {

    ObjectCfg::parseProperties(tree, myJob);

    uint16_t nCh(0);
    for( auto & it : tree ) {
        if( iequals( it.first, "CHANNEL" ) ) {
            Channel* tmpCh = new Channel( *this, myJob, nCh++ );
            tmpCh->parsePropertyTree( it.second );
            channels.push_back( shared_ptr<Channel>( tmpCh ) );
        }
    }

    //LOG_DEBUG << "Object::parseProperties() done.";

}


bpt::ptree Object::getPropertyTree( bpt::ptree& tree ) {

    bpt::ptree node;

    for( shared_ptr<Channel>& ch : channels ) {
        ch->getPropertyTree( node );
    }

    ObjectCfg::getProperties(node,myJob);

    tree.push_back( bpt::ptree::value_type( "object", node ) );

    return node;

}


size_t Object::size(void) const {
    size_t sz = ObjectCfg::size();
    sz += 2*sizeof(uint16_t);                   // channels.size() + ID
    for( const shared_ptr<Channel>& ch : channels ) {
        sz += ch->size();
    }
    sz += sizeof(nObjectImages);
    return sz;
}


uint64_t Object::pack(char* ptr) const {
    using redux::util::pack;
    uint64_t count = ObjectCfg::pack(ptr);
    count += pack(ptr+count, ID);
    count += pack(ptr+count, (uint16_t)channels.size());
    for( const shared_ptr<Channel>& ch : channels ) {
        count += ch->pack(ptr+count);
    }
    count += pack (ptr + count, nObjectImages);
    if (count != size()) {
        LOG_ERR << "(" << hexString (this) << "): Packing failed, there is a size mismatch:  count = " << count << "  sz = " << size();
    }
    return count;
}


uint64_t Object::unpack(const char* ptr, bool swap_endian) {
    using redux::util::unpack;

    uint64_t count = ObjectCfg::unpack(ptr, swap_endian);
    count += unpack(ptr+count, ID, swap_endian);
    uint16_t tmp;
    count += unpack(ptr+count, tmp, swap_endian);
    channels.resize(tmp);
    for( shared_ptr<Channel>& ch : channels ) {
        ch.reset(new Channel(*this,myJob));
        count += ch->unpack(ptr+count, swap_endian);
    }
    count += unpack (ptr + count, nObjectImages, swap_endian);
    return count;
}


size_t Object::nImages(void) const {
    if(nObjectImages) return nObjectImages;
    size_t objImages = 0;        // reset image-count
    for (const shared_ptr<Channel>& ch : channels) objImages += ch->nImages();
    return objImages;
}


void Object::calcPatchPositions(const std::vector<uint16_t>& y, const std::vector<uint16_t>& x) {
    for( shared_ptr<Channel>& ch : channels ) ch->calcPatchPositions(y,x);
}


void Object::initProcessing( WorkSpace::Ptr ws ) {
    //cout << "Object::initProcessing(" << hexString(this) << ")" << endl;
    if( patchSize && pupilPixels ) {
        P.resize(2*pupilPixels,2*pupilPixels);
        Q.resize(2*pupilPixels,2*pupilPixels);
        ftSum.resize(patchSize,patchSize);          // full-complex for now
        for( auto& ch : channels ) {
            if( patchSize == ch->patchSize && pupilPixels == ch->pupilPixels) {
                ch->initProcessing(ws);
            } else {
                LOG_ERR << "Different sub-field sizes not supported....yet.  obj.patchSize=" << patchSize
                << "  ch.patchSize=" << ch->patchSize << "  obj.pupilPixels=" << pupilPixels << "  ch.pupilPixels=" << ch->pupilPixels;
            }
        }
    } else {
        LOG_ERR << "Object patchSize is 0 !!!";
    }
}


void Object::initPatch( ObjectData& od ) {
    unique_lock<mutex> lock (mtx);
//    cout << "Object::initPatch(" << hexString (this) << ")  reg_gamma = " << myJob.reg_gamma << endl;
    reg_gamma = 0;
    ftSum.zero();
}


void Object::getResults(ObjectData& od) {
    
    unique_lock<mutex> lock (mtx);
    
    FourierTransform avgObjFT(patchSize, patchSize, FT_FULLCOMPLEX|FT_REORDER );
    Array<complex_t> tmpC(patchSize, patchSize);
    Array<double> tmpD(patchSize, patchSize);
    avgObjFT.zero();
    tmpD.zero();
    complex_t* aoPtr = avgObjFT.get();
    double* dPtr = tmpD.get();
    double avgNoiseVariance = 0.0;
    double frequencyCutoff = sqrt(patchSize*patchSize);
    for (shared_ptr<Channel>& ch : channels) {
        for (shared_ptr<SubImage>& im : ch->subImages) {
            im->restore(aoPtr, dPtr);
            avgNoiseVariance += sqr(im->stats.noise);
        }
        if (ch->frequencyCutoff < frequencyCutoff ) frequencyCutoff = ch->frequencyCutoff;
    }
    avgNoiseVariance /= nObjectImages;
    for (auto & ind : otfIndices) {
        if( fabs(dPtr[ind]) > 0 ) {
            aoPtr[ind] /= dPtr[ind];
        } else aoPtr[ind] = dPtr[ind] = 0;
    }
    if (!(myJob.runFlags&RF_NO_FILTER)) {
  //      cout << "Applying filter:  noise = " << avgNoiseVariance << " lf = " << frequencyCutoff << endl;
        ScharmerFilter(aoPtr, dPtr, patchSize, patchSize, avgNoiseVariance, 0.90 * frequencyCutoff);
    }
    
    avgObjFT.directInverse(tmpC.get());

    // TODO: implement add/subtract plane properly

    od.results.resize(patchSize, patchSize);
    od.results.assign(tmpC);


}


void Object::initPQ (void) {
//    cout << "i" << flush;
 //   cout << "Object::initPQ(" << hexString(this) << ")  reg_gamma = " << myJob.reg_gamma << endl;
    P.zero();
    Q = (reg_gamma); // * 0.1 / nImages);
}

/*
void Object::addAllFT( void ) {
    unique_lock<mutex> lock( mtx );
    ftSum.zero();
    for( shared_ptr<Channel>& ch : channels ) {
        for( shared_ptr<SubImage>& im : ch->subImages ) {
            im->addFT(ftSum);
        }
    }
}
*/

void Object::addToFT( const redux::image::FourierTransform& ft, double rg ) {
    unique_lock<mutex> lock( mtx );
    reg_gamma += 0.10*rg;///nImages;
    const complex_t* ftPtr = ft.get();
    double* ftsPtr = ftSum.get();
    for (size_t ind = 0; ind < ftSum.nElements(); ++ind) {
        ftsPtr[ind] += norm (ftPtr[ind]);
    }
}


void Object::addDiffToPQ(const redux::image::FourierTransform& ft, const Array<complex_t>& otf, const Array<complex_t>& oldotf) {

    unique_lock<mutex> lock (mtx);
//    cout << "Object::addDiffToPQ()" << endl;
    double *qPtr = Q.get();
    complex_t *pPtr = P.get();
    const complex_t *ftPtr = ft.get();
    const complex_t *otfPtr = otf.get();
    const complex_t *ootfPtr = oldotf.get();

    for( auto ind: otfIndices ) {
        qPtr[ind] += norm(otfPtr[ind]) - norm(ootfPtr[ind]);
        pPtr[ind] += ftPtr[ind] * conj(otfPtr[ind] - ootfPtr[ind]);
    }
}


void Object::addAllPQ(void) {
    for( shared_ptr<Channel>& ch : channels ) {
        for( shared_ptr<SubImage>& im : ch->subImages ) {
            unique_lock<mutex> lock( mtx );
            im->addPQ(P.get(),Q.get());
        }
    }
}


void Object::slask(void) {
//    cout << "Object::slask(void)" << endl;
//     static int bla(0);
//     unique_lock<mutex> lock( mtx );
//     redux::file::Ana::write( "ftsum_" + to_string( bla++ ) + ".f0", ftSum );
//     Array<double> img;
//     ftSum.inv(img);
//     redux::file::Ana::write( "ftsuminv_" + to_string( bla ) + ".f0", img );
}


void Object::calcMetric (void) {
    
   // cout << "Object::calcMetric(" << hexString(this) << ")  " << __LINE__ << "   otfsz = " << otfIndices.size() << endl;
//cout << "c" << flush;
    currentMetric = 0.0;
    const double* ftsPtr = ftSum.get();
    const complex_t* pPtr = P.get();
    const double* qPtr = Q.get();
    
    unique_lock<mutex> lock (mtx);
 //   Ana::write ("Qm.f0", Q);
 //   Ana::write ("Pm.f0", P);
 //   Ana::write ("FTm.f0", ftSum);
    size_t N = 4*pupilPixels*pupilPixels;
    for (size_t ind=0; ind<N; ++ind) {
    //for (auto & ind : otfIndices) {
        currentMetric += (ftsPtr[ind] - norm (pPtr[ind]) / qPtr[ind]);
    }
   // cout << "Object::calcMetric(" << hexString(this) << ")   m1 = " << currentMetric << endl;
    //currentMetric *= (weight / otfIndices.size());
    currentMetric *= (weight / N);

   // cout << "Object::calcMetric(" << hexString(this) << ")   m = " << setprecision(12) << currentMetric << endl;

}


bool Object::checkCfg (void) {

    if ( (saveMask & SF_SAVE_PSF) && (saveMask & SF_SAVE_PSF_AVG)) {
        LOG_WARN << "Both GET_PSF and GET_PSF_AVG mode specified.";
    }
    if (channels.empty()) {
        LOG_CRITICAL << "Each object must have at least 1 channel specified.";
    }

    for (shared_ptr<Channel>& ch : channels) {
        if (!ch->checkCfg()) return false;
    }

    if (outputFileName.empty()) {   // TODO: clean this up
        string tpl = channels[0]->imageTemplate;
        size_t p = tpl.find_first_of ('%');
        if (p != string::npos) {
            string tmpString = boost::str (boost::format (tpl) % 1);
            auto it = tmpString.begin();
            auto it2 = tpl.begin();
            p = 0;
            size_t i = std::min (tmpString.length(), tpl.length());
            while (p < i && tmpString[p] == tpl[p]) p++;
            it = tmpString.end();
            it2 = tpl.end();
            size_t ii = tmpString.length() - 1;
            i = tpl.length() - 1;
            while (ii && i && tmpString[ii] == tpl[i]) {
                ii--;
                i--;
            }
            tmpString.replace (p, ii - p + 1, "%d..%d");
            if (count (tmpString.begin(), tmpString.end(), '%') == 2) {
                outputFileName = boost::str (boost::format (tmpString) % *channels[0]->imageNumbers.begin() % *channels[0]->imageNumbers.rbegin());
            } else {
                LOG_CRITICAL << boost::format ("failed to generate output filename from \"%s\"  (->\"%s\").") % tpl % tmpString;
                return false;
            }
        } else LOG_CRITICAL << boost::format ("first filename template \"%s\" does not contain valid format specifier.") % tpl;
    }

    return true;
}


bool Object::checkData (void) {

    bfs::path tmpOF(outputFileName + ".ext");
    bfs::path tmpPath = tmpOF.parent_path();
    if( !bfs::exists (tmpPath) ) {
        if( !bfs::create_directories(tmpPath) ) {
            LOG_CRITICAL << boost::format ("failed to create directory for output: %s") % tmpPath;
            return false;
        } else LOG_TRACE << boost::format ("create output directory %s") % tmpPath;
    }
    try {
        bfs::path slask(tmpPath);
        slask += "_test_writability_";
        bfs::create_directory(slask);
        bfs::remove_all(slask);
    } catch( ... ) {
        LOG_CRITICAL << boost::format ("output directory %s not writable: %s") % tmpPath;
        return false;
    }
    for (int i = 1; i & FT_MASK; i <<= 1) {
        if (i & myJob.outputFileType) {  // this filetype is specified.
            tmpOF.replace_extension (FileTypeExtensions.at ( (FileType) i));
            if( bfs::exists(tmpOF) && !(myJob.runFlags & RF_FORCE_WRITE) ) {
                LOG_CRITICAL << boost::format ("output file %s already exists! Use -f (or OVERWRITE) to replace file.") % tmpOF;
                return false;
            } else {
                LOG << "Output filename: " << tmpOF;
            }
        }
    }

    for (shared_ptr<Channel>& ch : channels) {
        if (!ch->checkData()) return false;
    }

    return true;
}


void Object::init( void ) {

    for( shared_ptr<Channel>& ch : channels ) {
        ch->init();
    }

//   init( KL_cfg* kl_cfg, double lambda, double r_c, int nph_in, int basis, int nm, int *mode_num,
//              int nch, int *ndo, int **dorder, int **dtype, int kl_min_mode, int kl_max_mode, double svd_reg, double angle, double **pupil_in )

//    modes.init(coeff,lambda,r_c,nph,myJob.basis,myJob.modes);

    /*    for( int o = 1; o <= nObjects; ++o ) {
            mode[o] = new modes( kl_cfs, cfg->lambda[o], cfg->lim_freq[o] / 2.0, cfg->nph[o], cfg->basis, cfg->nModes, cfg->mode_num, nChannels[o], cfg->nDiversityOrders[o], cfg->dorder[o], cfg->dtype[o], cfg->kl_min_mode, cfg->kl_max_mode, cfg->svd_reg, cfg->angle[o], cfg->pupil[o], io );
    //              cfg->pix2cf[o]/=0.5*cfg->lambda[o]*cfg->lim_freq[o]*(mode[o]->mode[0][2][cfg->nph[o]/2+1][cfg->nph[o]/2]-mode[o]->mode[0][2][cfg->nph[o]/2][cfg->nph[o]/2]);
    //              cfg->cf2pix[o]*=0.5*cfg->lambda[o]*cfg->lim_freq[o]*(mode[o]->mode[0][2][cfg->nph[o]/2+1][cfg->nph[o]/2]-mode[o]->mode[0][2][cfg->nph[o]/2][cfg->nph[o]/2]);
            cfg->pix2cf[o] /= 0.5 * cfg->lambda[o] * cfg->lim_freq[o] * mode[o]->mode[0][2]->ddx();
            cfg->cf2pix[o] *= 0.5 * cfg->lambda[o] * cfg->lim_freq[o] * mode[o]->mode[0][2]->ddx();
        }
    */


}


void Object::initCache (void) {
    for (shared_ptr<Channel>& ch : channels) {
        ch->initCache();
        //size_t sz = pupilIndices.size();
        pupilIndices.insert (ch->pupilIndices.begin(), ch->pupilIndices.end());
        //if (sz != pupilIndices.size()) cout << "Added to pupilIndices: " << sz << " -> " <<  pupilIndices.size() << endl;
        //sz = otfIndices.size();
        otfIndices.insert (ch->otfIndices.begin(), ch->otfIndices.end());
        //if (sz != otfIndices.size()) cout << "Added to otfIndices: " << sz << " -> " <<  otfIndices.size() << endl;
    }
}


void Object::cleanup( void ) {

}


void Object::loadData( boost::asio::io_service& service ) {
    for( shared_ptr<Channel>& ch : channels ) {
        ch->loadData( service );
    }
}


void Object::preprocessData (boost::asio::io_service& service) {
    nObjectImages = nImages();
    for (shared_ptr<Channel>& ch : channels) {
        ch->preprocessData (service);
    }
}


void Object::normalize(boost::asio::io_service& service ) {

    double maxMean = std::numeric_limits<double>::lowest();
    for( shared_ptr<Channel>& ch : channels ) {
        double mM = ch->getMaxMean();
        if( mM > maxMean ) maxMean = mM;
    }
    for( shared_ptr<Channel>& ch : channels ) {
        ch->normalizeData(service, maxMean);
    }
}


void Object::prepareStorage (void) {

    bfs::path fn = bfs::path (outputFileName + ".momfbd");      // TODO: fix storage properly

    LOG << "Preparing file " << fn << " for temporary, and possibly final, storage.";

    std::shared_ptr<FileMomfbd> info (new FileMomfbd());

    // Extract date/time from the git commit.
    int day, month, year, hour;
    char buffer [15];
    sscanf (reduxCommitTime, "%4d-%2d-%2d %2d", &year, &month, &day, &hour);
    sprintf (buffer, "%4d%02d%02d.%02d", year, month, day, hour);
    info->versionString = buffer;
    info->version = atof (info->versionString.c_str());

    info->dateString = "FIXME";
    info->timeString = "FIXME";
//     if(false) {
//         for( auto& it: channels ) {
//             info->fileNames.push_back ( "FIXME" );
//         }
//         info->dataMask |= MOMFBD_NAMES;
//     }
    info->nFileNames = info->fileNames.size();

    int32_t n_img = nImages();
    int32_t nChannels = info->nChannels = channels.size();
    info->clipStartX = sharedArray<int16_t> (nChannels);
    info->clipEndX = sharedArray<int16_t> (nChannels);
    info->clipStartY = sharedArray<int16_t> (nChannels);
    info->clipEndY = sharedArray<int16_t> (nChannels);
    for (int i = 0; i < nChannels; ++i) {
        info->clipStartX.get() [ i ] = channels[i]->alignClip[0];
        info->clipEndX.get() [ i ] = channels[i]->alignClip[1];
        info->clipStartY.get() [ i ] = channels[i]->alignClip[2];
        info->clipEndY.get() [ i ] = channels[i]->alignClip[3];
    }

    info->nPH = pupilPixels;

    Array<float> tmp;

    if (saveMask & SF_SAVE_MODES && (info->nPH > 0)) {
        double pupilRadiusInPixels = pupilPixels / 2.0;
        if (channels.size()) pupilRadiusInPixels = channels[0]->pupilRadiusInPixels;
        tmp.resize (myJob.modeNumbers.size() + 1, info->nPH, info->nPH);            // +1 to also fit pupil in the array
        tmp.zero();
        Array<float> tmp_slice (tmp, 0, 0, 0, info->nPH - 1, 0, info->nPH - 1);     // subarray
        tmp_slice.assign(myJob.globalData->fetch (pupilPixels, pupilRadiusInPixels).first);     // store pupil at index 0
        info->phOffset = 0;
        if (myJob.modeNumbers.size()) {
            Cache::ModeID id (myJob.klMinMode, myJob.klMaxMode, 0, pupilPixels, pupilRadiusInPixels, rotationAngle, myJob.klCutoff);
            info->nModes = myJob.modeNumbers.size();
            info->modesOffset = pupilPixels * pupilPixels * sizeof (float);
            for (uint16_t & it : myJob.modeNumbers) {   // Note: globalData might also contain modes we don't want to save here, e.g. PhaseDiversity modes.
                if (it > 3 && myJob.modeBasis != ZERNIKE) {       // use Zernike modes for the tilts
                    id.firstMode = myJob.klMinMode;
                    id.lastMode = myJob.klMaxMode;
                } else id.firstMode = id.lastMode = 0;
                tmp_slice.shift (0, 1);     // shift subarray 1 step
                id.modeNumber = it;
                tmp_slice.assign(reinterpret_cast<const Array<double>&> (*myJob.globalData->fetch(id)));

            }
        }
    }

    /*
    info->pix2cf = pix2cf;
    info->cf2pix = cf2pix;
    */
    info->nPatchesX = 0;//nPatchesX;
    info->nPatchesY = 0;//nPatchesY;
    info->patches.resize (info->nPatchesX, info->nPatchesY);

    auto dummy = sharedArray<int32_t> (nChannels);
    for (int x = 0; x < info->nPatchesX; ++x) {
        for (int y = 0; y < info->nPatchesY; ++y) {
            info->patches(x,y).region[0] = info->patches(x,y).region[2] = 1;
            info->patches(x,y).region[1] = info->patches(x,y).region[3] = patchSize;
            info->patches(x,y).nChannels = nChannels;
            info->patches(x,y).nim = sharedArray<int32_t> (nChannels); //dummy;
            info->patches(x,y).dx = sharedArray<int32_t> (nChannels); //dummy;
            info->patches(x,y).dy = sharedArray<int32_t> (nChannels);; //dummy;
            for (int i = 0; i < nChannels; ++i) {
                info->patches(x,y).nim.get() [i] = 1000 + x * 100 + y * 10 + i;
                info->patches(x,y).dx.get() [i] = 2000 + x * 100 + y * 10 + i;
                info->patches(x,y).dy.get() [i] = 3000 + x * 100 + y * 10 + i;
            }
            info->patches(x,y).npsf = n_img;
            info->patches(x,y).nobj = n_img;
            info->patches(x,y).nres = n_img;
            info->patches(x,y).nalpha = n_img;
            info->patches(x,y).ndiv = n_img;
            info->patches(x,y).nm = info->nModes;
            info->patches(x,y).nphx = info->nPH;
            info->patches(x,y).nphy = info->nPH;

        }   // y-loop
    }   // x-loop



    uint8_t writeMask = MOMFBD_IMG;                                                 // always output image
    if (saveMask & SF_SAVE_PSF || saveMask & SF_SAVE_PSF_AVG)    writeMask |= MOMFBD_PSF;
    if (saveMask & SF_SAVE_COBJ)    writeMask |= MOMFBD_OBJ;
    if (saveMask & SF_SAVE_RESIDUAL)    writeMask |= MOMFBD_RES;
    if (saveMask & SF_SAVE_ALPHA)    writeMask |= MOMFBD_ALPHA;
    if (saveMask & SF_SAVE_DIVERSITY)    writeMask |= MOMFBD_DIV;
    if (saveMask & SF_SAVE_MODES)    writeMask |= MOMFBD_MODES;

    //cout << "prepareStorage: " << bitString(writeMask) << endl;
    info->write (fn.string(), reinterpret_cast<char*> (tmp.ptr()), writeMask);
    //cout << "prepareStorage done."  << endl;

}


void Object::writeAna (const redux::util::Array<PatchData::Ptr>& patches) {

    LOG << "Writing output to ANA.   baseName=\"" << outputFileName << "\"";
    
    LOG_WARN << "Writing to ANA still not properly implemented...";

    for (uint y = 0; y < patches.dimSize(0); ++y) {
        for (uint x = 0; x < patches.dimSize(1); ++x) {
            patches(y,x)->cacheLoad(true);
            bfs::path fn = bfs::path (outputFileName + "_img_"+to_string(x)+"_"+to_string(y)+".f0");
            Ana::write(fn.string(), patches(y,x)->objects[ID].results);
        }
    }
    
    for (shared_ptr<Channel>& ch : channels) {
        ch->writeAna(patches);
    }
    
    
}


void Object::writeFits (const redux::util::Array<PatchData::Ptr>& patches) {
    bfs::path fn = bfs::path (outputFileName + ".fits");
    LOG << "Writing output to file: " << fn;
    LOG_ERR << "Writing to FITS still not implemented...";
}


void Object::writeMomfbd (const redux::util::Array<PatchData::Ptr>& patches) {

    bfs::path fn = bfs::path (outputFileName + "_thi.momfbd");      // TODO: fix storage properly

    LOG << "Writing output to file: " << fn;

    std::shared_ptr<FileMomfbd> info (new FileMomfbd());

    // Extract date/time from the git commit.
    int day, month, year, hour;
    char buffer [15];
    sscanf (reduxCommitTime, "%4d-%2d-%2d %2d", &year, &month, &day, &hour);
    sprintf (buffer, "%4d%02d%02d.%01d", year, month, day, hour);
    info->versionString = buffer;
    info->version = atof (info->versionString.c_str());

    info->dateString = myJob.observationDate;
    info->timeString = "FIXME";
//     if(false) {
//         for( auto& it: channels ) {
//             info->fileNames.push_back ( "FIXME" );
//         }
//         info->dataMask |= MOMFBD_NAMES;
//     }
    info->nFileNames = info->fileNames.size();

    int32_t nChannels = info->nChannels = channels.size();
    info->clipStartX = sharedArray<int16_t> (nChannels);
    info->clipEndX = sharedArray<int16_t> (nChannels);
    info->clipStartY = sharedArray<int16_t> (nChannels);
    info->clipEndY = sharedArray<int16_t> (nChannels);
    for (int i = 0; i < nChannels; ++i) {
        info->clipStartX.get() [ i ] = channels[i]->alignClip[0];
        info->clipEndX.get() [ i ] = channels[i]->alignClip[1];
        info->clipStartY.get() [ i ] = channels[i]->alignClip[2];
        info->clipEndY.get() [ i ] = channels[i]->alignClip[3];
    }

    info->nPH = pupilPixels;

    Array<float> modes;

    if (saveMask & SF_SAVE_MODES && (info->nPH > 0)) {
        double pupilRadiusInPixels = pupilPixels / 2.0;
        if (channels.size()) pupilRadiusInPixels = channels[0]->pupilRadiusInPixels;
        modes.resize (myJob.modeNumbers.size() + 1, info->nPH, info->nPH);            // +1 to also fit pupil in the array
        modes.zero();
        Array<float> tmp_slice (modes, 0, 0, 0, info->nPH - 1, 0, info->nPH - 1);     // subarray
        tmp_slice.assign(myJob.globalData->fetch (pupilPixels, pupilRadiusInPixels).first);
        info->phOffset = 0;
        if (myJob.modeNumbers.size()) {
            Cache::ModeID id (myJob.klMinMode, myJob.klMaxMode, 0, pupilPixels, pupilRadiusInPixels, rotationAngle, myJob.klCutoff);
            info->nModes = myJob.modeNumbers.size();
            info->modesOffset = pupilPixels * pupilPixels * sizeof (float);
            for (uint16_t & it : myJob.modeNumbers) {   // Note: globalData might also contain modes we don't want to save here, e.g. PhaseDiversity modes.
                if (it > 3 && myJob.modeBasis != ZERNIKE) {       // use Zernike modes for the tilts
                    id.firstMode = myJob.klMinMode;
                    id.lastMode = myJob.klMaxMode;
                } else id.firstMode = id.lastMode = 0;
                tmp_slice.shift (0, 1);     // shift subarray 1 step
                id.modeNumber = it;
                tmp_slice.assign(reinterpret_cast<const Array<double>&>(*myJob.globalData->fetch(id)));

            }
        }
    }

    /*
    info->pix2cf = pix2cf;
    info->cf2pix = cf2pix;
    */
    info->nPatchesY = patches.dimSize (0);
    info->nPatchesX = patches.dimSize (1);
    info->patches.resize(info->nPatchesX, info->nPatchesY);

    auto dummy = sharedArray<int32_t> (nChannels);
    for (int x = 0; x < info->nPatchesX; ++x) {
        for (int y = 0; y < info->nPatchesY; ++y) {
            patches(y,x)->cacheLoad(true);
            info->patches(x,y).region[0] = patches(y,x)->roi.first.x;
            info->patches(x,y).region[1] = patches(y,x)->roi.last.x;
            info->patches(x,y).region[2] = patches(y,x)->roi.first.y;
            info->patches(x,y).region[3] = patches(y,x)->roi.last.y;
            info->patches(x,y).nChannels = nChannels;
            info->patches(x,y).nim = sharedArray<int32_t> (nChannels); //dummy;
            info->patches(x,y).dx = sharedArray<int32_t> (nChannels); //dummy;
            info->patches(x,y).dy = sharedArray<int32_t> (nChannels);; //dummy;
            for (int i = 0; i < nChannels; ++i) {
                info->patches(x,y).nim.get()[i] = channels[i]->nImages();
                info->patches(x,y).dx.get()[i] = patches(y,x)->objects[ID].channels[i].shift.x;
                info->patches(x,y).dy.get()[i] = patches(y,x)->objects[ID].channels[i].shift.y;
            }
            info->patches(x,y).npsf = nObjectImages;
            info->patches(x,y).nobj = nObjectImages;
            info->patches(x,y).nres = nObjectImages;
            info->patches(x,y).nalpha = nObjectImages;
            info->patches(x,y).ndiv = nObjectImages;
            info->patches(x,y).nm = info->nModes;
            info->patches(x,y).nphx = info->nPH;
            info->patches(x,y).nphy = info->nPH;
        }   // y-loop
    }   // x-loop



    uint8_t writeMask = MOMFBD_IMG;                                                 // always output image
    int64_t imgSize = patchSize*patchSize*sizeof(float);
    size_t patchDataSize = imgSize;
    if (false) { // not implemented yet
        if (saveMask & SF_SAVE_PSF || saveMask & SF_SAVE_PSF_AVG) {
            writeMask |= MOMFBD_PSF;
            patchDataSize += patchSize * patchSize * sizeof (float);
        }
        if (saveMask & SF_SAVE_COBJ)                               writeMask |= MOMFBD_OBJ;
        if (saveMask & SF_SAVE_RESIDUAL)                           writeMask |= MOMFBD_RES;
        if (saveMask & SF_SAVE_ALPHA)                              writeMask |= MOMFBD_ALPHA;
        if (saveMask & SF_SAVE_DIVERSITY)                          writeMask |= MOMFBD_DIV;
        if (saveMask & SF_SAVE_MODES)                              writeMask |= MOMFBD_MODES;
    }

    size_t modeSize = modes.nElements() * sizeof (float);
    size_t totalSize = modeSize + patches.nElements() * patchDataSize;

    auto tmp = sharedArray<char> (totalSize);
    memcpy(tmp.get(), modes.get(), modeSize);
    char* tmpPtr = tmp.get();
    int64_t offset = modeSize;
    for (int x = 0; x < info->nPatchesX; ++x) {
        for (int y = 0; y < info->nPatchesY; ++y) {
            memcpy(tmpPtr+offset, patches(y,x)->objects[ID].results.get(), imgSize);
            info->patches(x,y).imgPos = offset;
            offset += imgSize;
        }
    }

    //cout << "prepareStorage: " << bitString(writeMask) << endl;
    info->write (fn.string(), reinterpret_cast<char*> (tmp.get()), writeMask);
    //cout << "prepareStorage done."  << endl;

}

void Object::writeResults (const redux::util::Array<PatchData::Ptr>& patches) {
    if (myJob.outputFileType & FT_ANA) writeAna (patches);
    if (myJob.outputFileType & FT_FITS) writeFits (patches);
    if (myJob.outputFileType & FT_MOMFBD) writeMomfbd (patches);
}

void Object::storePatches (WorkInProgress& wip, boost::asio::io_service& service, uint8_t nThreads) {

    bfs::path fn = bfs::path (outputFileName);
    fn.replace_extension ("momfbd");
    std::shared_ptr<FileMomfbd> info (new FileMomfbd (fn.string()));

    LOG_DEBUG << "storePatches()";

    for (auto & it : wip.parts) {
        PatchData::Ptr patch = static_pointer_cast<PatchData> (it);
        LOG_DEBUG << "storePatches() index: (" << patch->index.x << "," << patch->index.y << ")  offset = "
                  << info->patches (patch->index.x , patch->index.y).offset;
        patch->step = MomfbdJob::JSTEP_COMPLETED;
    }

}


size_t Object::sizeOfPatch (uint32_t npixels) const {
    size_t sz (0);
    for (const shared_ptr<Channel>& ch : channels) {
        sz += ch->sizeOfPatch (npixels);
    }
    return sz;
}


Point16 Object::getImageSize (void) {
    Point16 sizes;
    for (shared_ptr<Channel>& ch : channels) {
        Point16 tmp = ch->getImageSize();
        if (sizes.x == 0) {
            sizes = tmp;
        } else if (tmp != sizes) {
            throw std::logic_error ("The images have different sizes for the different channels, please verify the ALIGN_CLIP values.");
        }
    }
    return sizes;
}


void Object::dump (std::string tag) {
//    cout << "Dumping object  #" << ID << "  this=" << hexString (this) << " with tag=" << tag << endl;
    tag += "_obj_" + to_string (ID);
    Ana::write (tag + "_ftsum.f0", ftSum);
    Ana::write (tag + "_q.f0", Q);
    Ana::write (tag + "_p.f0", P);
    for (shared_ptr<Channel>& ch : channels) {
        ch->dump (tag);
    }

}
