#include "redux/momfbd/object.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/defines.hpp"
#include "redux/momfbd/channel.hpp"

#include "redux/translators.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/constants.hpp"
#include "redux/logger.hpp"

#include <limits>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/info_parser.hpp>

namespace bfs = boost::filesystem;
using namespace redux::momfbd;
using namespace redux::util;
using namespace redux;
using namespace std;
using boost::algorithm::iequals;

#define lg Logger::mlg
namespace {
    const string thisChannel = "momfbdobj";

    void pupcfg( double &lim_freq, double &r_c, uint32_t &nph, double lambda, uint32_t np, double telescope_d, double telescope_f, double arcsecperpix ) {
        double rad2deg = 360.0 / ( 2.0 * redux::PI );
        double scale_angle2CCD = arcsecperpix / ( rad2deg * 3600.0 );
        double q_number = lambda / ( scale_angle2CCD * telescope_d );
        lim_freq = ( double )np / q_number;
        nph = np>>2;
        r_c = lim_freq / 2.0;       // telescope radius in pupil pixels...
        if( nph < r_c ) {           // this should only be needed for oversampled images
            uint32_t goodsizes[] = { 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45, 48, 50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128, 135, 144 };
            for( int i = 0; ( nph = max( goodsizes[i], nph ) ) < r_c; ++i ); // find right size
        }
        nph <<= 1;
    }

}


Object::Object( const MomfbdJob& j ) : imageNumbers(j.imageNumbers), darkNumbers(j.darkNumbers), reg_gamma(j.reg_gamma),
                                 weight(1), angle(0), lambda(0), nPoints(j.patchSize), sequenceNumber(j.sequenceNumber),
                                 nph(0), stokesWeights(j.stokesWeights), flags(j.flags), imageDataDir(j.imageDataDir),
                                 fillpix_method(j.fillpix_method), output_data_type(j.output_data_type),
                                 lim_freq(), r_c(), myJob( j ), pupil(j.pupil) {
    
    
}

Object::~Object() {
}

void Object::parseProperties( bpt::ptree& tree, const string& fn ) {

    LOG_DEBUG << "Object::parseProperties()";


    reg_gamma = tree.get<double>( "REG_GAMMA", myJob.reg_gamma );
    weight = tree.get<double>( "WEIGHT", DEF_WEIGHT );
    angle = tree.get<double>( "ANGLE", DEF_ANGLE );
    lambda = tree.get<double>( "WAVELENGTH" );

    imageNumbers = tree.get<vector<uint32_t>>( "IMAGE_NUM", myJob.imageNumbers );
    sequenceNumber = tree.get<uint32_t>( "SEQUENCE_NUM", myJob.sequenceNumber );
    darkNumbers = tree.get<vector<uint32_t>>( "DARK_NUM", myJob.darkNumbers );
    wf_num = tree.get<vector<uint32_t>>( "WFINDEX", vector<uint32_t>() );

    imageDataDir = cleanPath( tree.get<string>( "IMAGE_DATA_DIR", myJob.imageDataDir ) );

    string tmpString;
    fillpix_method = MomfbdJob::getFromMap( tmpString = tree.get<string>( "FPMETHOD", "" ), fillpixMap );
    if( fillpix_method == 0 ) {
        if( tmpString.length() > 0 ) {
            LOG_ERR << "unknown fillpix method \"" << tmpString << "\"\n  Valid entries currently are: "
                    "\"median\", \"invdistweight\" or \"horint\"";
        }
        fillpix_method = myJob.fillpix_method;
    }

    pupcfg( lim_freq, r_c, nph, lambda, nPoints, myJob.telescopeDiameter, myJob.telescopeFocalLength, myJob.arcSecsPerPixel );

    if( myJob.pupil.get() ) { // pupil?
        if( myJob.pupilSize > nph ) {
            LOG_WARN << "MomfbdObject: the pupil image needs to be rebinned...";
            // pupil=rebin(myJob.pupil,1,myJob.pupilSize,1,myJob.pupilSize,r_c,1,nph,1,nph);
        }
        else if( myJob.pupilSize < nph ) {
            LOG_WARN << "MomfbdObject: the pupil image is " << myJob.pupilSize << " but nph is " << nph;
            LOG_WARN << "MomfbdObject: the pupil image needs to be interpolated, this is probably not very accurate.";
            // pupil=resample(myJob.pupil,1,myJob.pupilSize,1,myJob.pupilSize,r_c,1,nph,1,nph);
        }
        else {
            LOG_DETAIL << "MomfbdObject: matching pupil found.";
            pupil = myJob.pupil;    // NB: shared data.
        }
    }
    else pupil.resize();

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
            Channel* tmpCh = new Channel( *this, myJob );
            tmpCh->parseProperties( it.second );
            channels.push_back( shared_ptr<Channel>( tmpCh ) );
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


bpt::ptree Object::getPropertyTree( bpt::ptree* root ) {

    bpt::ptree tree;

    if( reg_gamma != myJob.reg_gamma ) tree.put( "REG_GAMMA", reg_gamma );
    if( weight != DEF_WEIGHT ) tree.put( "WEIGHT", weight );
    if( angle != DEF_ANGLE ) tree.put( "ANGLE", angle );
    tree.put( "WAVELENGTH", lambda );
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


size_t Object::size(void) const {
    size_t sz = 2;
    sz += 4*sizeof(uint32_t);
    sz += 6*sizeof(double);
    sz += imageNumbers.size()*sizeof(uint32_t)+sizeof(size_t);
    sz += sequenceNumbers.size()*sizeof(uint32_t)+sizeof(size_t);
    sz += darkNumbers.size()*sizeof(uint32_t)+sizeof(size_t);
    sz += wf_num.size()*sizeof(uint32_t)+sizeof(size_t);
    sz += stokesWeights.size()*sizeof(double)+sizeof(size_t);
    sz += imageDataDir.length() + outputFileName.length() + 2;
    sz += sizeof(size_t);                   // channels.size()
    for(auto& it: channels) {
        sz += it->size();
    }
    sz += pupil.size();
    return sz;
}


uint64_t Object::pack(char* ptr) const {
    using redux::util::pack;

    uint64_t count = pack(ptr, fillpix_method);
    count += pack(ptr+count, output_data_type);
    count += pack(ptr+count, nPoints);
    count += pack(ptr+count, sequenceNumber);
    count += pack(ptr+count, nph);
    count += pack(ptr+count, flags);
    count += pack(ptr+count, reg_gamma);
    count += pack(ptr+count, weight);
    count += pack(ptr+count, angle);
    count += pack(ptr+count, lambda);
    count += pack(ptr+count, lim_freq);
    count += pack(ptr+count, r_c);
    count += pack(ptr+count, imageNumbers);
    count += pack(ptr+count, sequenceNumbers);
    count += pack(ptr+count, darkNumbers);
    count += pack(ptr+count, wf_num);
    count += pack(ptr+count, stokesWeights);
    count += pack(ptr+count, imageDataDir);
    count += pack(ptr+count, outputFileName);
    count += pack(ptr+count, channels.size());
    for(auto& it: channels) {
        count += it->pack(ptr+count);
    }
    count += pupil.pack(ptr+count);
    return count;
}


uint64_t Object::unpack(const char* ptr, bool swap_endian) {
    using redux::util::unpack;

    uint64_t count = unpack(ptr, fillpix_method);
    count += unpack(ptr+count, output_data_type);
    count += unpack(ptr+count, nPoints, swap_endian);
    count += unpack(ptr+count, sequenceNumber, swap_endian);
    count += unpack(ptr+count, nph, swap_endian);
    count += unpack(ptr+count, flags, swap_endian);
    count += unpack(ptr+count, reg_gamma, swap_endian);
    count += unpack(ptr+count, weight, swap_endian);
    count += unpack(ptr+count, angle, swap_endian);
    count += unpack(ptr+count, lambda, swap_endian);
    count += unpack(ptr+count, lim_freq, swap_endian);
    count += unpack(ptr+count, r_c, swap_endian);
    count += unpack(ptr+count, imageNumbers, swap_endian);
    count += unpack(ptr+count, sequenceNumbers, swap_endian);
    count += unpack(ptr+count, darkNumbers, swap_endian);
    count += unpack(ptr+count, wf_num, swap_endian);
    count += unpack(ptr+count, stokesWeights, swap_endian);
    count += unpack(ptr+count, imageDataDir);
    count += unpack(ptr+count, outputFileName);
    size_t tmp;
    count += unpack(ptr+count, tmp, swap_endian);
    channels.resize(tmp);
    for(auto& it: channels) {
        it.reset(new Channel(*this,myJob));
        count += it->unpack(ptr+count, swap_endian);
    }
    count += pupil.unpack(ptr+count, swap_endian);
    
   return count;
}


bool Object::checkCfg(void) {
    
    if( channels.empty() ) return false;     // need at least 1 channel
    
    for( auto & it : channels ) {
        if( !it->checkCfg() ) return false;
    }
    
    return true;
}


bool Object::checkData(void) {

    for( auto & it : channels ) {
        if( !it->checkData() ) return false;
    }
    
    return true;
}


void Object::init( void ) {

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


void Object::cleanup( void ) {

}


void Object::loadData( boost::asio::io_service& service, boost::thread_group& pool ) {
    LOG << "Object::loadData()";
    for( auto & it : channels ) {
        it->loadData( service, pool );
    }
    
}


void Object::preprocessData(boost::asio::io_service& service, boost::thread_group& pool) {
    LOG << "Object::preprocessData()";
    for(auto& it: channels) {
        it->preprocessData(service, pool);
    }
}


void Object::normalize(boost::asio::io_service& service, boost::thread_group& pool) {
    LOG << "Object::normalize()";
    double maxMean = std::numeric_limits<double>::min();
    for(auto it: channels) {
        double mM = it->getMaxMean();
        if( mM > maxMean ) maxMean = mM;
    }
    for(auto it: channels) {
        it->normalizeData(service, pool, maxMean);
    }
}


size_t Object::sizeOfPatch(uint32_t npixels) const {
    size_t sz(0);
    for( auto & it : channels ) {
        sz += it->sizeOfPatch(npixels);
    }
    return sz;
}


uint64_t Object::packPatch( Patch::Ptr patch, char* ptr ) const {
    
    uint64_t count(0);
    for( auto & it : channels ) {
        count += it->packPatch(patch,ptr);
    }
    return count;
}
 

Point Object::clipImages(void) {
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

