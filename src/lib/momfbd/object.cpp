#include "redux/momfbd/object.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/data.hpp"
#include "redux/momfbd/util.hpp"

#include "redux/file/filemomfbd.hpp"
#include "redux/translators.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/constants.hpp"
#include "redux/logger.hpp"
#include "redux/revision.hpp"

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
using namespace redux::util;
using namespace redux;
using namespace std;
using boost::algorithm::iequals;

#define lg Logger::mlg
namespace {

const string thisChannel = "object";

}


Object::Object( MomfbdJob& j ) : ObjectCfg(j), myJob(j) {


}


Object::~Object() {

}


void Object::parsePropertyTree( bpt::ptree& tree ) {

    ObjectCfg::parseProperties(tree, myJob);

    for( auto & it : tree ) {
        if( iequals( it.first, "CHANNEL" ) ) {
            Channel* tmpCh = new Channel( *this, myJob );
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
    sz += sizeof(uint16_t);                   // channels.size()
    for( const shared_ptr<Channel>& ch : channels ) {
        sz += ch->size();
    }
    return sz;
}


uint64_t Object::pack(char* ptr) const {
    using redux::util::pack;
    uint64_t count = ObjectCfg::pack(ptr);
    count += pack(ptr+count, (uint16_t)channels.size());
    for( const shared_ptr<Channel>& ch : channels ) {
        count += ch->pack(ptr+count);
    }
    if(count != size()) {
        LOG_ERR << "(" << hexString(this) << "): Packing failed, there is a size mismatch:  count = " << count << "  sz = " << size();
    }
    return count;
}


uint64_t Object::unpack(const char* ptr, bool swap_endian) {
    using redux::util::unpack;

    uint64_t count = ObjectCfg::unpack(ptr, swap_endian);
    uint16_t tmp;
    count += unpack(ptr+count, tmp, swap_endian);
    channels.resize(tmp);
    for( shared_ptr<Channel>& ch : channels ) {
        ch.reset(new Channel(*this,myJob));
        count += ch->unpack(ptr+count, swap_endian);
    }
    return count;
}


size_t Object::nImages(size_t offset) {
    size_t n(0);
    for( shared_ptr<Channel>& ch : channels ) n += ch->nImages(offset+n);
    return n;
}


void Object::collectImages(redux::util::Array<float>& stack) const {
    for( const shared_ptr<Channel>& ch : channels ) ch->collectImages(stack);
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
    cout << "Object::initPatch(" << hexString(this) << ")  reg_gamma = " << myJob.reg_gamma << endl;
    addAllFT();
}


void Object::initPQ( void ) {
    unique_lock<mutex> lock( mtx );
    P.zero();
    Q = myJob.reg_gamma;
}


void Object::addAllFT( void ) {
    unique_lock<mutex> lock( mtx );
    ftSum.zero();
    for( shared_ptr<Channel>& ch : channels ) ch->addAllFT(ftSum);
}


void Object::addToFT( const redux::image::FourierTransform& ft ) {
    unique_lock<mutex> lock( mtx );
    if( !ftSum.sameSizes(ft) ) {
        ftSum.resize(ft.dimensions(true));
        ftSum.zero();
    }
    //  std::cout << "Object::addToFT():  1   " << hexString(this) << printArray(ft.dimensions(), "  ft.dims") << printArray(ftSum.dimensions(), "  ftSum.dims") << std::endl;
   // ftSum += ft;
//    std::cout << "Object::addToFT():  2   " << printArray(ft.dimensions(), "ft.dims") << printArray(ftSum.dimensions(), "  ftSum.dims") << std::endl;
}


void Object::addToPQ( const redux::image::FourierTransform& ft, const Array<complex_t> sj ) {
//   std::cout << "Object::addToPQ():  1   " << hexString(this) << printArray(ft.dimensions(), "  ft.dims") << printArray(sj.dimensions(), "  sj.dims") << std::endl;
    typedef Array<complex_t>::const_iterator cit;
    unique_lock<mutex> lock( mtx );
    cit ftit = ft.begin();
    cit sjit = sj.begin();
    auto qit = Q.begin();

    for( auto& pit: P ) {
        *qit++ += norm( *sjit );            // Q += sj.re^2 + sj.im^2 = norm(sj)
        pit += *ftit++ * *sjit++;           // P += ft * sj
    }
}


void Object::addAllPQ(void) {
    for( shared_ptr<Channel>& ch : channels ) {
        for( shared_ptr<SubImage>& im : ch->subImages ) {
            unique_lock<mutex> lock( mtx );
            im->addPQ(P,Q);
        }
    }
}


#include "redux/file/fileana.hpp"

void Object::slask(void) {
//    cout << "Object::slask(void)" << endl;
//     static int bla(0);
//     unique_lock<mutex> lock( mtx );
//     redux::file::Ana::write( "ftsum_" + to_string( bla++ ) + ".f0", ftSum );
//     Array<double> img;
//     ftSum.inv(img);
//     redux::file::Ana::write( "ftsuminv_" + to_string( bla ) + ".f0", img );
}


double Object::metric(void) {

    double m = 0.0;
    uint nY = ftSum.dimSize(0);
    uint nX = ftSum.dimSize(1);
    if(nY != P.dimSize(0) || nX != P.dimSize(1) ||
       nY != Q.dimSize(0) || nX != Q.dimSize(1) ) {
        cout << "Object::metric()   dim mismatch: " << printArray(ftSum.dimensions(),"ftdims")
                                                    << printArray(P.dimensions()," pdims")
                                                    << printArray(Q.dimensions()," qdims") << endl;
    }
    double** fs = makePointers(ftSum.ptr(), nY, nX);
    complex_t** p = makePointers(P.ptr(), nY, nX);
    double** q = makePointers(Q.ptr(), nY, nX);
    
    for(uint y=0; y<nY; ++y) {
        for(uint x=0; x<nX; ++x) {
            m += fs[y][x]-norm(p[y][x])/q[y][x]; // warning: |q| may be < 1E-20 (i.e. not "safe")?
        }
    }
    m/=(double)(nY*nX);

    cout << "Object::metric()   m = " << m << endl;
    delPointers(fs);
    delPointers(p);
    delPointers(q);
//   for(int x=1;x<=nPixels;++x)
//     for(int y=1;y<=nPixels;++y)
//       l+=ftSum[x][y]-(p[x][y].re*p[x][y].re+p[x][y].im*p[x][y].im)/q[x][y]; // warning: |q| may be < 1E-20 (i.e. not "safe")?
//   l/=(double)(nPixels*nPixels);
//   if(reg_alpha){
//     double sum=0.0;
//     for(int i=1;i<=nImages;++i) sum+=subImages[i]->metric(alpha[i]); // alpha is passed down
//     l+=0.5*reg_alpha*sum;
//   }
  return m;

}


bool Object::checkCfg(void) {

    if( ( saveMask & SF_SAVE_PSF ) && ( saveMask & SF_SAVE_PSF_AVG ) ) {
        LOG_WARN << "Both GET_PSF and GET_PSF_AVG mode specified.";
    }
    if( channels.empty() ) {
        LOG_CRITICAL << "Each object must have at least 1 channel specified.";
    }

    for( shared_ptr<Channel>& ch: channels ) {
        if( !ch->checkCfg() ) return false;
    }

    if( outputFileName.empty() ) {  // TODO: clean this up
        string tpl = channels[0]->imageTemplate;
        size_t p = tpl.find_first_of( '%' );
        if( p != string::npos ) {
            string tmpString = boost::str( boost::format( tpl ) % 1 );
            auto it = tmpString.begin();
            auto it2 = tpl.begin();
            p = 0;
            size_t i = std::min( tmpString.length(), tpl.length() );
            while( p < i && tmpString[p] == tpl[p] ) p++;
            it = tmpString.end();
            it2 = tpl.end();
            size_t ii = tmpString.length() - 1;
            i = tpl.length() - 1;
            while( ii && i && tmpString[ii] == tpl[i] ) {
                ii--;
                i--;
            }
            tmpString.replace( p, ii - p + 1, "%d..%d" );
            if( count( tmpString.begin(), tmpString.end(), '%' ) == 2 ) {
                outputFileName = boost::str( boost::format( tmpString ) % *channels[0]->imageNumbers.begin() % *channels[0]->imageNumbers.rbegin() );
            } else {
                LOG_CRITICAL << boost::format( "failed to generate output file from \"%s\"  (->\"%s\")." ) % tpl % tmpString;
                return false;
            }
        }
        else LOG_CRITICAL << boost::format( "first filename template \"%s\" does not contain valid format specifier." ) % tpl;
    }

    bfs::path tmpOF( outputFileName+".ext" );
    for( int i = 1; i&FT_MASK; i <<= 1 ) {
        if( i & myJob.outputFileType ) { // this filetype is specified.
            tmpOF.replace_extension( FileTypeExtensions.at( ( FileType )i ) );
            if( !(myJob.runFlags&RF_FORCE_WRITE) && bfs::exists( tmpOF ) ) {
                LOG_CRITICAL << boost::format( "output file %s already exists! Use -f (or OVERWRITE) to replace file." ) % tmpOF;
                return false;
            } else {
                LOG << "Output filename: " << tmpOF;
            }
        }
    }

    return true;
}


bool Object::checkData(void) {

//   fn = bfs::path( outputFileName );
    /*    if( !bfs::exists( fn ) ) {
            LOG_CRITICAL << "Not found !!! \"" << fn.string() << "\"";
            imageNumbers.erase( imageNumbers.begin() + i );
            continue;
        }
      */

    for( shared_ptr<Channel>& ch : channels ) {
        if( !ch->checkData() ) return false;
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


void Object::initCache( void ) {
    for( shared_ptr<Channel>& ch : channels ) {
        ch->initCache();
    }
}


void Object::cleanup( void ) {

}


void Object::loadData( boost::asio::io_service& service ) {
    for( shared_ptr<Channel>& ch : channels ) {
        ch->loadData( service );
    }
}


void Object::preprocessData(boost::asio::io_service& service ) {
    for( shared_ptr<Channel>& ch : channels ) {
        ch->preprocessData(service);
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


void Object::prepareStorage(void) {

    bfs::path fn = bfs::path( outputFileName + ".momfbd" );     // TODO: fix storage properly

    LOG << "Preparing file " << fn << " for temporary, and possibly final, storage.";

    std::shared_ptr<FileMomfbd> info ( new FileMomfbd() );

    // Extract date/time from the git commit.
    int day, month, year, hour;
    char buffer [15];
    sscanf(reduxCommitTime, "%4d-%2d-%2d %2d", &year, &month, &day, &hour);
    sprintf (buffer, "%4d%02d%02d.%02d", year, month, day, hour);
    info->versionString = buffer;
    info->version = atof ( info->versionString.c_str() );

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
    info->clipStartX = sharedArray<int16_t>(nChannels);
    info->clipEndX = sharedArray<int16_t>(nChannels);
    info->clipStartY = sharedArray<int16_t>(nChannels);
    info->clipEndY = sharedArray<int16_t>(nChannels);
    for(int i=0; i<nChannels; ++i ) {
        info->clipStartX.get() [ i ] = channels[i]->alignClip[0];
        info->clipEndX.get() [ i ] = channels[i]->alignClip[1];
        info->clipStartY.get() [ i ] = channels[i]->alignClip[2];
        info->clipEndY.get() [ i ] = channels[i]->alignClip[3];
    }

    info->nPH = pupilPixels;

    Array<float> tmp;

    if(saveMask&SF_SAVE_MODES && (info->nPH>0)) {
        double pupilRadiusInPixels = pupilPixels/2.0;
        if( channels.size() ) pupilRadiusInPixels = channels[0]->pupilRadiusInPixels;
        tmp.resize(myJob.modeNumbers.size()+1,info->nPH,info->nPH);                 // +1 to also fit pupil in the array
        tmp.zero();
        Array<float> tmp_slice(tmp, 0, 0, 0, info->nPH-1, 0, info->nPH-1);          // subarray
        tmp_slice = myJob.globalData->fetch(pupilPixels,pupilRadiusInPixels).first;
        info->phOffset = 0;
        if(myJob.modeNumbers.size()) {
            Cache::ModeID id(myJob.klMinMode, myJob.klMaxMode, 0, pupilPixels, pupilRadiusInPixels, wavelength, rotationAngle);
            if(myJob.modeBasis == ZERNIKE) {    // Needed because klMinMode/klMaxMode might be non-zero even if we are using Zerikes
                id.firstMode = id.lastMode = 0;
            }
            info->nModes = myJob.modeNumbers.size();
            info->modesOffset = pupilPixels*pupilPixels*sizeof(float);
            for( uint16_t& it: myJob.modeNumbers ) {    // Note: globalData might also contain modes we don't want to save here, e.g. PhaseDiversity modes.
                tmp_slice.shift(0,1);       // shift subarray 1 step (in the first dimension)
                id.modeNumber = it;
                tmp_slice = *(myJob.globalData->fetch(id));

            }
        }
    }

    /*
    info->pix2cf = pix2cf;
    info->cf2pix = cf2pix;
    */
    info->nPatchesX = 0;//nPatchesX;
    info->nPatchesY = 0;//nPatchesY;
    info->patches.resize ( info->nPatchesY, info->nPatchesX );

    auto dummy = sharedArray<int32_t>(nChannels);
    for ( int x = 0; x < info->nPatchesX; ++x ) {
        for ( int y = 0; y < info->nPatchesY; ++y ) {
            info->patches(x,y).region[0] = info->patches(x,y).region[2] = 1;
            info->patches(x,y).region[1] = info->patches(x,y).region[3] = patchSize;
            info->patches(x,y).nChannels = nChannels;
            info->patches(x,y).nim = sharedArray<int32_t>(nChannels); //dummy;
            info->patches(x,y).dx = sharedArray<int32_t>(nChannels); //dummy;
            info->patches(x,y).dy = sharedArray<int32_t>(nChannels);; //dummy;
            for(int i=0; i<nChannels; ++i) {
                info->patches(x,y).nim.get()[i] = 1000+x*100+y*10+i;
                info->patches(x,y).dx.get()[i] = 2000+x*100+y*10+i;
                info->patches(x,y).dy.get()[i] = 3000+x*100+y*10+i;
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
    if(saveMask&SF_SAVE_PSF || saveMask&SF_SAVE_PSF_AVG)    writeMask |= MOMFBD_PSF;
    if(saveMask&SF_SAVE_COBJ)    writeMask |= MOMFBD_OBJ;
    if(saveMask&SF_SAVE_RESIDUAL)    writeMask |= MOMFBD_RES;
    if(saveMask&SF_SAVE_ALPHA)    writeMask |= MOMFBD_ALPHA;
    if(saveMask&SF_SAVE_DIVERSITY)    writeMask |= MOMFBD_DIV;
    if(saveMask&SF_SAVE_MODES)    writeMask |= MOMFBD_MODES;

    //cout << "prepareStorage: " << bitString(writeMask) << endl;
    info->write ( fn.string(), reinterpret_cast<char*>(tmp.ptr()), writeMask );
    //cout << "prepareStorage done."  << endl;

}


void Object::storePatches( WorkInProgress& wip, boost::asio::io_service& service, uint8_t nThreads) {

    bfs::path fn = bfs::path( outputFileName );
    fn.replace_extension( "momfbd" );
    std::shared_ptr<FileMomfbd> info ( new FileMomfbd(fn.string()) );

    LOG_DEBUG << "storePatches()";

    for( auto& it: wip.parts ) {
        PatchData::Ptr patch = static_pointer_cast<PatchData>(it);
        LOG_DEBUG << "storePatches() index: (" << patch->index.x << "," << patch->index.y << ")  offest = "
                  << info->patches( patch->index.x ,patch->index.y).offset;
        patch->step = MomfbdJob::JSTEP_COMPLETED;
    }

}


size_t Object::sizeOfPatch(uint32_t npixels) const {
    size_t sz(0);
    for( const shared_ptr<Channel>& ch : channels ) {
        sz += ch->sizeOfPatch(npixels);
    }
    return sz;
}


Point16 Object::clipImages(void) {
    Point16 sizes;
    for( shared_ptr<Channel>& ch : channels ) {
        Point16 tmp = ch->clipImages();
        if(sizes.x == 0) {
            sizes = tmp;
        } else if( tmp != sizes ) {
            throw std::logic_error("The clipped images have different sizes for the different channels, please verify the ALIGN_CLIP values.");
        }
    }
    return sizes;
}

