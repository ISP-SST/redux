#include "redux/momfbd/objectcfg.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/defines.hpp"
#include "redux/momfbd/channelcfg.hpp"
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

    const string thisChannel = "momfbdobj";

    void getObjectPupilSize( double &lim_freq, double &r_c, uint32_t &nPupilPixels, double wavelength, uint32_t nPixels, double telescopeDiameter, double arcSecsPerPixel ) {
        double radians_per_arcsec = redux::PI/(180.0*3600.0);             // (2.0*redux::PI)/(360.0*3600.0)
        double radians_per_pixel = arcSecsPerPixel * radians_per_arcsec;
        double q_number = wavelength / ( radians_per_pixel * telescopeDiameter );
        lim_freq = ( double )nPixels / q_number;
        nPupilPixels = nPixels>>2;
        r_c = lim_freq / 2.0;                   // telescope radius in pupil pixels...
        if( nPupilPixels < r_c ) {           // this should only be needed for oversampled images
            uint32_t goodsizes[] = { 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45, 48, 50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128, 135, 144 };
            for( int i = 0; ( nPupilPixels = max( goodsizes[i], nPupilPixels ) ) < r_c; ++i ); // find right size
        }
        nPupilPixels <<= 1;
    }

}


ObjectCfg::ObjectCfg( const MomfbdJob& j ) : fillpix_method(j.fillpix_method), output_data_type(j.output_data_type), objectSize(j.patchSize),
                                sequenceNumber(j.sequenceNumber), objectPupilSize(0), flags(j.flags), reg_gamma(j.reg_gamma),
                                weight(1), angle(0), wavelength(0), lim_freq(), r_c(), cf2pix(), pix2cf(), cf2def(), 
                                imageNumbers(j.imageNumbers), darkNumbers(j.darkNumbers), stokesWeights(j.stokesWeights),
                                imageDataDir(j.imageDataDir), pupil(j.pupil), myJob(j) {
    
    
}

ObjectCfg::~ObjectCfg() {
}

void ObjectCfg::parseProperties( bpt::ptree& tree, const string& fn ) {


    reg_gamma = tree.get<double>( "REG_GAMMA", myJob.reg_gamma );
    weight = tree.get<double>( "WEIGHT", DEF_WEIGHT );
    angle = tree.get<double>( "ANGLE", DEF_ANGLE );
    wavelength = tree.get<double>( "WAVELENGTH" );
    imageNumbers = tree.get<vector<uint32_t>>( "IMAGE_NUM", myJob.imageNumbers );
    sequenceNumber = tree.get<uint32_t>( "SEQUENCE_NUM", myJob.sequenceNumber );
    darkNumbers = tree.get<vector<uint32_t>>( "DARK_NUM", myJob.darkNumbers );
    wf_num = tree.get<vector<uint32_t>>( "WFINDEX", vector<uint32_t>() );

    LOG_ERR << "myJob." << printArray(myJob.imageNumbers,"imageNumbers" );
    imageDataDir = cleanPath( tree.get<string>( "IMAGE_DATA_DIR", myJob.imageDataDir ) );
    LOG_ERR << "obj." << printArray(imageNumbers,"imageNumbers" );

    string tmpString;
    fillpix_method = MomfbdJob::getFromMap( tmpString = tree.get<string>( "FPMETHOD", "" ), fillpixMap );
    if( fillpix_method == 0 ) {
        if( tmpString.length() > 0 ) {
            LOG_ERR << "unknown fillpix method \"" << tmpString << "\"\n  Valid entries currently are: "
                    "\"median\", \"invdistweight\" or \"horint\"";
        }
        fillpix_method = myJob.fillpix_method;
    }

    getObjectPupilSize( lim_freq, r_c, objectPupilSize, wavelength, objectSize, myJob.telescopeDiameter, myJob.arcSecsPerPixel );

    pupilFile = cleanPath( tree.get<string>( "PUPIL", myJob.pupilFile ), imageDataDir );
    if( pupilFile != "" ) {
        util::loadPupil( pupilFile, pupil, objectPupilSize );
        pupil.resize();
        if(pupil.nDimensions()>0) pupilFile.clear();    // pupil already loaded, master doesn't need to look for it.
    }

    
    auto filetypes = tree.get<vector<FileType>>( "FILE_TYPE", vector<FileType>( 1, FT_NONE ) );
    if( filetypes.empty() ) LOG_ERR << "\"FILE_TYPE\" has to be one of ANA/FITS/MOMFBD.";
    else {
        uint32_t tmp = 0;
        for( auto & it : filetypes ) tmp |= it;
        if( tmp ) flags = ( flags&~MFBD_FT_MASK ) | tmp; // override global setting
    }

    tmpString = tree.get<string>( "DATA_TYPE", "" );
    if( iequals( tmpString, "FLOAT" ) ) {
        output_data_type = MFBD_F32T;
    }
    else if( iequals( tmpString, "SHORT" ) ) {
        output_data_type = MFBD_I16T;
    }
    else if( tmpString == "" ) {
        output_data_type = myJob.output_data_type;
    }
    else {
        LOG_WARN << "\"DATA_TYPE\" unrecognised data type \"" << tmpString << "\", using INT16";
        output_data_type = MFBD_I16T;
    }

    flags = myJob.flags;
    MomfbdJob::maybeOverride( tree.get<bool>( "FIT_PLANE", flags & MFBD_FIT_PLANE ), flags, MFBD_FIT_PLANE );
    MomfbdJob::maybeOverride( tree.get<bool>( "GET_ALPHA", flags & MFBD_GET_ALPHA ), flags, MFBD_GET_ALPHA );
    MomfbdJob::maybeOverride( tree.get<bool>( "GET_COBJ", flags & MFBD_GET_COBJ ), flags, MFBD_GET_COBJ );
    MomfbdJob::maybeOverride( tree.get<bool>( "GET_DIVERSITY", flags & MFBD_GET_DIVERSITY ), flags, MFBD_GET_DIVERSITY );
    MomfbdJob::maybeOverride( tree.get<bool>( "GET_MODES", flags & MFBD_GET_MODES ), flags, MFBD_GET_MODES );
    MomfbdJob::maybeOverride( tree.get<bool>( "GET_PSF", flags & MFBD_GET_PSF ), flags, MFBD_GET_PSF );
    MomfbdJob::maybeOverride( tree.get<bool>( "GET_PSF_AVG", flags & MFBD_GET_PSF_AVG ), flags, MFBD_GET_PSF_AVG );
    MomfbdJob::maybeOverride( tree.get<bool>( "GET_RESIDUAL", flags & MFBD_GET_RESIDUAL ), flags, MFBD_GET_RESIDUAL );
    MomfbdJob::maybeOverride( tree.get<bool>( "NO_FILTER", flags & MFBD_NO_FILTER ), flags, MFBD_NO_FILTER );
    MomfbdJob::maybeOverride( tree.get<bool>( "SAVE_FFDATA", flags & MFBD_SAVE_FFDATA ), flags, MFBD_SAVE_FFDATA );

    if( ( flags & MFBD_GET_PSF ) && ( flags & MFBD_GET_PSF_AVG ) ) {
        LOG_WARN << "both GET_PSF and GET_PSF_AVG mode requested";
    }

    stokesWeights = tree.get<vector<double>>( "VECTOR", myJob.stokesWeights );

    for( auto & it : tree ) {
        if( iequals( it.first, "CHANNEL" ) ) {
            ChannelCfg* tmpCh = new ChannelCfg( *this, myJob );
            tmpCh->parseProperties( it.second );
            channels.push_back( shared_ptr<ChannelCfg>( tmpCh ) );
        }
    }


    if( fn.length() > 0 ) { // command line arg override
        outputFileName = fn;
    }
    else {
        outputFileName = tree.get<string>( "OUTPUT_FILE", "" );
    }

    if( outputFileName.empty() ) {
        if( channels.size() > 0 ) {
            string tpl = channels[0]->imageTemplate;
            size_t q, p = tpl.find_first_of( '%' );
            q = tpl.find_first_not_of( '%', p );
            if( ( q = tpl.find_first_not_of( '%', p ) ) != string::npos ) {
                tmpString = boost::str( boost::format( tpl ) % 1 );
                auto it = tmpString.begin();
                auto it2 = tpl.begin();
                p = 0;
                size_t r = std::min( tmpString.length(), tpl.length() );
                while( p < r && tmpString[p] == tpl[p] ) p++;
                it = tmpString.end();
                it2 = tpl.end();
                q = tmpString.length() - 1;
                r = tpl.length() - 1;
                while( q && r && tmpString[q] == tpl[r] ) {
                    q--;
                    r--;
                }
                tmpString.replace( p, q - p + 1, "%d..%d" );
                if( count( tmpString.begin(), tmpString.end(), '%' ) != 2 ) {
                    LOG_CRITICAL << boost::format( "failed to generate output file from \"%s\"  (->\"%s\")." ) % tpl % tmpString;
                }
                else outputFileName = boost::str( boost::format( tmpString ) % channels[0]->imageNumbers[0] % *channels[0]->imageNumbers.rbegin() );
            }
            else LOG_CRITICAL << boost::format( "first filename template \"%s\" does not contain valid format specifier." ) % tpl;
        }
    }

    bfs::path tmpOF( outputFileName );
    for( int i = 1; i & MFBD_FT_MASK; i <<= 1 ) {
        if( i & flags ) { // this filetype is specified.
            tmpOF.replace_extension( FileTypeExtensions.at( ( FileType )i ) );
            if( bfs::exists( tmpOF ) ) LOG_CRITICAL << boost::format( "output file %s already exists! Use -f to overwrite." ) % tmpOF;
            else LOG << boost::format( "output file is %s" ) % tmpOF;
        }
    }

    LOG_DEBUG << "Object::parseProperties() done.";

}


bpt::ptree ObjectCfg::getPropertyTree( bpt::ptree* root ) {

    bpt::ptree tree;

    if( reg_gamma != myJob.reg_gamma ) tree.put( "REG_GAMMA", reg_gamma );
    if( weight != DEF_WEIGHT ) tree.put( "WEIGHT", weight );
    if( angle != DEF_ANGLE ) tree.put( "ANGLE", angle );
    tree.put( "WAVELENGTH", wavelength );
    if( imageNumbers != myJob.imageNumbers ) tree.put( "IMAGE_NUM", imageNumbers );
    if( sequenceNumber != myJob.sequenceNumber ) tree.put( "SEQUENCE_NUM", sequenceNumber );
    if( darkNumbers != myJob.darkNumbers ) tree.put( "DARK_NUM", darkNumbers );
    if( !wf_num.empty() ) tree.put( "WFINDEX", wf_num );
    if( imageDataDir != myJob.imageDataDir ) tree.put( "IMAGE_DATA_DIR", imageDataDir );
    if( fillpix_method != myJob.fillpix_method ) tree.put( "FPMETHOD", fillpix_method );
    // TODO FILE_TYPE/DATA_TYPE
    uint32_t dflags = flags ^ myJob.flags;
    if( dflags & MFBD_FIT_PLANE ) tree.put( "FIT_PLANE", ( bool )( flags & MFBD_FIT_PLANE ) );
    if( dflags & MFBD_GET_ALPHA ) tree.put( "GET_ALPHA", ( bool )( flags & MFBD_GET_ALPHA ) );
    if( dflags & MFBD_GET_COBJ ) tree.put( "GET_COBJ", ( bool )( flags & MFBD_GET_COBJ ) );
    if( dflags & MFBD_GET_DIVERSITY ) tree.put( "GET_DIVERSITY", ( bool )( flags & MFBD_GET_DIVERSITY ) );
    if( dflags & MFBD_GET_MODES ) tree.put( "GET_MODES", ( bool )( flags & MFBD_GET_MODES ) );
    if( dflags & MFBD_GET_PSF ) tree.put( "GET_PSF", ( bool )( flags & MFBD_GET_PSF ) );
    if( dflags & MFBD_GET_PSF_AVG ) tree.put( "GET_PSF_AVG", ( bool )( flags & MFBD_GET_PSF_AVG ) );
    if( dflags & MFBD_GET_RESIDUAL ) tree.put( "GET_RESIDUAL", ( bool )( flags & MFBD_GET_RESIDUAL ) );
    if( dflags & MFBD_NO_FILTER ) tree.put( "NO_FILTER", ( bool )( flags & MFBD_NO_FILTER ) );
    if( dflags & MFBD_SAVE_FFDATA ) tree.put( "SAVE_FFDATA", ( bool )( flags & MFBD_SAVE_FFDATA ) );
    if( stokesWeights != myJob.stokesWeights ) tree.put( "VECTOR", stokesWeights );
    if( !outputFileName.empty() ) tree.put( "OUTPUT_FILE", outputFileName );

    for( auto & it : channels ) {
        it->getPropertyTree( &tree );
    }

    if( root ) {
        root->push_back( bpt::ptree::value_type( "object", tree ) );
    }
    
    return tree;
    
}


size_t ObjectCfg::size(void) const {
    size_t sz = 2;                          // fillpix_method, output_data_type
    sz += 4*sizeof(uint32_t);               // objectSize, sequenceNumber, objectPupilSize, flags
    sz += 9*sizeof(double);                 // reg_gamma, weight, angle, wavelength, lim_freq, r_c, cf2pix, pix2cf, cf2def
    sz += imageNumbers.size()*sizeof(uint32_t)+sizeof(uint64_t);
    sz += sequenceNumbers.size()*sizeof(uint32_t)+sizeof(uint64_t);
    sz += darkNumbers.size()*sizeof(uint32_t)+sizeof(uint64_t);
    sz += wf_num.size()*sizeof(uint32_t)+sizeof(uint64_t);
    sz += stokesWeights.size()*sizeof(double)+sizeof(uint64_t);
    sz += imageDataDir.length() + outputFileName.length() + pupilFile.length() + 3;
    sz += pupil.size();
    sz += 2*sizeof(uint32_t);                   // channels.size() & modes.size()
    for(auto& it: channels) {
        sz += it->size();
    }
    for(auto& it: modes) {
        sz += sizeof(uint32_t);
        sz += it.second->size();
    }
    return sz;
}


uint64_t ObjectCfg::pack(char* ptr) const {
    using redux::util::pack;
    uint64_t count = pack(ptr, fillpix_method);
    count += pack(ptr+count, output_data_type);
    count += pack(ptr+count, objectSize);
    count += pack(ptr+count, sequenceNumber);
    count += pack(ptr+count, objectPupilSize);
    count += pack(ptr+count, flags);
    count += pack(ptr+count, reg_gamma);
    count += pack(ptr+count, weight);
    count += pack(ptr+count, angle);
    count += pack(ptr+count, wavelength);
    count += pack(ptr+count, lim_freq);
    count += pack(ptr+count, r_c);
    count += pack(ptr+count, cf2pix);
    count += pack(ptr+count, pix2cf);
    count += pack(ptr+count, cf2def);
    count += pack(ptr+count, imageNumbers);
    count += pack(ptr+count, sequenceNumbers);
    count += pack(ptr+count, darkNumbers);
    count += pack(ptr+count, wf_num);
    count += pack(ptr+count, stokesWeights);
    count += pack(ptr+count, imageDataDir);
    count += pack(ptr+count, outputFileName);
    count += pack(ptr+count, pupilFile);
    count += pupil.pack(ptr+count);
    count += pack(ptr+count, (uint32_t)modes.size());
    for(auto& it: modes) {
        count += pack(ptr+count, it.first);
        count += it.second->pack(ptr+count);
    }
    count += pack(ptr+count, (uint32_t)channels.size());
    for(auto& it: channels) {
        count += it->pack(ptr+count);
    }
    if(count != size()) cout << "Obj " << hexString(this) << " has a size mismatch: " << count << "  sz = " << size() << "  diff = " << (size()-count) <<endl;
    return count;
}


uint64_t ObjectCfg::unpack(const char* ptr, bool swap_endian) {
    using redux::util::unpack;

    uint64_t count = unpack(ptr, fillpix_method);
    count += unpack(ptr+count, output_data_type);
    count += unpack(ptr+count, objectSize, swap_endian);
    count += unpack(ptr+count, sequenceNumber, swap_endian);
    count += unpack(ptr+count, objectPupilSize, swap_endian);
    count += unpack(ptr+count, flags, swap_endian);
    count += unpack(ptr+count, reg_gamma, swap_endian);
    count += unpack(ptr+count, weight, swap_endian);
    count += unpack(ptr+count, angle, swap_endian);
    count += unpack(ptr+count, wavelength, swap_endian);
    count += unpack(ptr+count, lim_freq, swap_endian);
    count += unpack(ptr+count, r_c, swap_endian);
    count += unpack(ptr+count, cf2pix, swap_endian);
    count += unpack(ptr+count, pix2cf, swap_endian);
    count += unpack(ptr+count, cf2def, swap_endian);
    count += unpack(ptr+count, imageNumbers, swap_endian);
    count += unpack(ptr+count, sequenceNumbers, swap_endian);
    count += unpack(ptr+count, darkNumbers, swap_endian);
    count += unpack(ptr+count, wf_num, swap_endian);
    count += unpack(ptr+count, stokesWeights, swap_endian);
    count += unpack(ptr+count, imageDataDir);
    count += unpack(ptr+count, outputFileName);
    count += unpack(ptr+count, pupilFile);
    count += pupil.unpack(ptr+count, swap_endian);
    uint32_t tmp;
    count += unpack(ptr+count, tmp, swap_endian);
    if(tmp) {
        ModeCache& cache = ModeCache::getCache();
        while(tmp--) {
            uint32_t modeNumber;
            PupilMode::Ptr md( new PupilMode );
            count += unpack(ptr+count, modeNumber, swap_endian);
            count += md->unpack(ptr+count, swap_endian);
            modes.emplace( modeNumber, cache.mode( modeNumber, objectPupilSize, r_c, wavelength, angle ) );
        }
    }
    count += unpack(ptr+count, tmp, swap_endian);
    channels.resize(tmp);
    for(auto& it: channels) {
        it.reset(new ChannelCfg(*this,myJob));
        count += it->unpack(ptr+count, swap_endian);
    }
   return count;
}


size_t ObjectCfg::nImages(size_t offset) {
    size_t n(0);
    for( auto & it : channels ) n += it->nImages(offset+n);
    return n;
}


void ObjectCfg::collectImages(redux::util::Array<float>& stack) const {
    for( auto & it : channels ) it->collectImages(stack);
}


void ObjectCfg::initWorkSpace( WorkSpace& ws ) {
    for( auto & it : channels ) it->initWorkSpace(ws);
}


bool ObjectCfg::checkCfg(void) {
    
    if( channels.empty() ) return false;     // need at least 1 channel
    
    for( auto & it : channels ) {
        if( !it->checkCfg() ) return false;
    }
    
    return true;
}


bool ObjectCfg::checkData(void) {
    
 //   fn = bfs::path( outputFileName );
/*    if( !bfs::exists( fn ) ) {
        LOG_CRITICAL << "Not found !!! \"" << fn.string() << "\"";
        imageNumbers.erase( imageNumbers.begin() + i );
        continue;
    }
  */  

    for( auto & it : channels ) {
        if( !it->checkData() ) return false;
    }
    
    return true;
}


void ObjectCfg::init( void ) {

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


void ObjectCfg::cleanup( void ) {

}


void ObjectCfg::loadData( boost::asio::io_service& service ) {
    
    LOG << "Object::loadData()";
    if( pupil.nDimensions() < 2 && pupilFile != "" ) {
        service.post( std::bind( util::loadPupil, pupilFile, std::ref(pupil), objectPupilSize ) );
    }

    for( auto & it : channels ) {
        it->loadData( service );
    }
    
}


void ObjectCfg::preprocessData(boost::asio::io_service& service ) {
    
    LOG_TRACE << "Object::preprocessData()";
    ModeCache& cache = ModeCache::getCache();
    if( pupil.nDimensions() < 2 ) {
        LOG_DETAIL << "Object::preprocessData()  Generating pupil.";
        const std::pair<Array<double>, double>& ppair = cache.pupil(objectPupilSize,r_c);
        pupil = ppair.first;
    }
    
    if ( modes.size() < myJob.modes.size() ) {
        LOG_DETAIL << "Object::preprocessData()  Generating modes.";
        for( auto & modeNumber: myJob.modes ) {
            if( myJob.basis == CFG_KARHUNEN_LOEVE ) {
                modes.emplace( modeNumber, cache.mode( myJob.klMinMode, myJob.klMaxMode, modeNumber,objectPupilSize, r_c, wavelength,angle ) );
            } else {
                modes.emplace( modeNumber, cache.mode( modeNumber, objectPupilSize, r_c, wavelength, angle ) );
            }
        }
    }
    
    cf2pix = util::cf2pix(myJob.arcSecsPerPixel, myJob.telescopeDiameter);
    cf2def = util::cf2def(1.0,myJob.telescopeDiameter/myJob.telescopeFocalLength);
    auto it = modes.find(2);
    if(it != modes.end()) {
        LOG_DETAIL << "Object::preprocessData()  adjusting pixel <-> coefficient conversion for this object.";
        double tmp = it->second->at(objectPupilSize/2,objectPupilSize/2+1) - it->second->at(objectPupilSize/2,objectPupilSize/2);
        tmp *= 0.5 * wavelength * lim_freq;
        cf2pix *= tmp;
    }
    pix2cf = 1.0/cf2pix;
    
    for(auto& it: channels) {
        it->preprocessData(service);
    }
}


void ObjectCfg::normalize(boost::asio::io_service& service ) {
    
    LOG_TRACE << "Object::normalize()";
    double maxMean = std::numeric_limits<double>::lowest();
    for(auto it: channels) {
        double mM = it->getMaxMean();
        if( mM > maxMean ) maxMean = mM;
    }
    for(auto it: channels) {
        it->normalizeData(service, maxMean);
    }
}


void ObjectCfg::prepareStorage(void) {

    bfs::path fn = bfs::path( outputFileName );
    fn.replace_extension( "momfbd" );
    
    LOG_DEBUG << "Preparing file " << fn << " for temporary, and possibly final, storage.";
    
    FileMomfbd* infoPtr = new FileMomfbd();
    std::shared_ptr<FileMomfbd> info ( infoPtr );

    // Extract date/time from the git commit.
    int day, month, year, hour;
    char buffer [15];
    sscanf(reduxCommitTime, "%4d-%2d-%2d %2d", &year, &month, &day, &hour);
    sprintf (buffer, "%4d%02d%02d.%02d", year, month, day, hour);
    infoPtr->versionString = buffer;
    infoPtr->version = atof ( infoPtr->versionString.c_str() );

    infoPtr->dateString = "FIXME";
    infoPtr->timeString = "FIXME";
//     if(false) {
//         for( auto& it: channels ) {
//             infoPtr->fileNames.push_back ( "FIXME" );
//         }
//         infoPtr->dataMask |= MOMFBD_NAMES;
//     }
    infoPtr->nFileNames = infoPtr->fileNames.size();
    
    int32_t n_img = nImages();
    int32_t nChannels = infoPtr->nChannels = channels.size();
    infoPtr->clipStartX = sharedArray<int16_t>(nChannels);
    infoPtr->clipEndX = sharedArray<int16_t>(nChannels);
    infoPtr->clipStartY = sharedArray<int16_t>(nChannels);
    infoPtr->clipEndY = sharedArray<int16_t>(nChannels);
    for(int i=0; i<nChannels; ++i ) {
        infoPtr->clipStartX.get() [ i ] = channels[i]->alignClip[0];
        infoPtr->clipEndX.get() [ i ] = channels[i]->alignClip[1];
        infoPtr->clipStartY.get() [ i ] = channels[i]->alignClip[2];
        infoPtr->clipEndY.get() [ i ] = channels[i]->alignClip[3];
    }
    
    infoPtr->nPH = objectPupilSize;
    
    Array<float> tmp;
    if(flags&MFBD_GET_MODES && (infoPtr->nPH>0)) {
        tmp.resize(modes.size()+1,infoPtr->nPH,infoPtr->nPH);
        tmp.zero();
        Array<float> tmp_slice(tmp, 0, 0, 0, infoPtr->nPH-1, 0, infoPtr->nPH-1);
        tmp_slice = pupil;
        infoPtr->phOffset = 0;
        if(modes.size()) {
            infoPtr->nModes = modes.size();
            infoPtr->modesOffset = objectPupilSize*objectPupilSize*sizeof(float);
            for( auto& it: modes ) {
                tmp_slice.shift(0,1);
                tmp_slice = *it.second;
            }
        }
    }
    infoPtr->pix2cf = pix2cf;
    infoPtr->cf2pix = cf2pix;
    
    infoPtr->nPatchesX = myJob.nPatchesX;
    infoPtr->nPatchesY = myJob.nPatchesY;
    infoPtr->patches.resize ( infoPtr->nPatchesY, infoPtr->nPatchesX );
    
    auto dummy = sharedArray<int32_t>(nChannels);
    for ( int x = 0; x < infoPtr->nPatchesX; ++x ) {
        for ( int y = 0; y < infoPtr->nPatchesY; ++y ) {
            infoPtr->patches(x,y).region[0] = infoPtr->patches(x,y).region[2] = 1;
            infoPtr->patches(x,y).region[1] = infoPtr->patches(x,y).region[3] = objectSize;
            infoPtr->patches(x,y).nChannels = nChannels;
            infoPtr->patches(x,y).nim = sharedArray<int32_t>(nChannels); //dummy;
            infoPtr->patches(x,y).dx = sharedArray<int32_t>(nChannels); //dummy;
            infoPtr->patches(x,y).dy = sharedArray<int32_t>(nChannels);; //dummy;
            for(int i=0; i<nChannels; ++i) {
                infoPtr->patches(x,y).nim.get()[i] = 1000+x*100+y*10+i;
                infoPtr->patches(x,y).dx.get()[i] = 2000+x*100+y*10+i;
                infoPtr->patches(x,y).dy.get()[i] = 3000+x*100+y*10+i;
            }
            infoPtr->patches(x,y).npsf = n_img;
            infoPtr->patches(x,y).nobj = n_img;
            infoPtr->patches(x,y).nres = n_img;
            infoPtr->patches(x,y).nalpha = n_img;
            infoPtr->patches(x,y).ndiv = n_img;
            infoPtr->patches(x,y).nm = infoPtr->nModes;
            infoPtr->patches(x,y).nphx = infoPtr->nPH;
            infoPtr->patches(x,y).nphy = infoPtr->nPH;

        }   // y-loop
    }   // x-loop

    
    
    
    uint8_t writeMask = MOMFBD_IMG;                                                 // always output image
    if(flags&MFBD_GET_PSF || flags&MFBD_GET_PSF_AVG)    writeMask |= MOMFBD_PSF;
    if(flags&MFBD_GET_COBJ)    writeMask |= MOMFBD_OBJ;
    if(flags&MFBD_GET_RESIDUAL)    writeMask |= MOMFBD_RES;
    if(flags&MFBD_GET_ALPHA)    writeMask |= MOMFBD_ALPHA;
    if(flags&MFBD_GET_DIVERSITY)    writeMask |= MOMFBD_DIV;
    if(flags&MFBD_GET_MODES)    writeMask |= MOMFBD_MODES;
    
    //cout << "prepareStorage: " << bitString(writeMask) << endl;
    infoPtr->write ( fn.string(), reinterpret_cast<char*>(tmp.ptr()), writeMask );
    //cout << "prepareStorage done."  << endl;
 
}


size_t ObjectCfg::sizeOfPatch(uint32_t npixels) const {
    size_t sz(0);
    for( auto & it : channels ) {
        sz += it->sizeOfPatch(npixels);
    }
    return sz;
}


void  ObjectCfg::applyLocalOffsets(PatchData::Ptr patch) const {
    //LOG_TRACE << "ObjectCfg::applyLocalOffsets #" << patch->id;
    for( auto & it : channels ) {
        it->applyLocalOffsets(patch);
    }
}
 

Point ObjectCfg::clipImages(void) {
    Point sizes;
    for(auto& it: channels) {
        Point tmp = it->clipImages();
        if(sizes.x == 0) {
            sizes = tmp;
        } else if( tmp != sizes ) {
            throw std::logic_error("The clipped images have different sizes for the different channels, please verify the ALIGN_CLIP values.");
        }
    }
    return sizes;
}

