#include "redux/momfbd/momfbdjob.hpp"

#include "redux/momfbd/defines.hpp"
#include "redux/momfbd/object.hpp"

#include "redux/translators.hpp"
#include "redux/file/fileio.hpp"
#include "redux/constants.hpp"
#include "redux/logger.hpp"
#include "redux/util/bitoperations.hpp"

#include <boost/algorithm/string.hpp>

using namespace redux::momfbd;
using namespace redux::file;
using namespace redux::util;
using namespace redux;
using namespace std;
using boost::algorithm::iequals;


const map<string, int> redux::momfbd::getstepMap( { {"getstep_steepest_descent", CFG_GETSTEP_SDSC},
    {"getstep_conjugate_gradient", CFG_GETSTEP_CNJG},
    {"getstep_BFGS_inv", CFG_GETSTEP_BFGS_inv}
} );

const map<string, int> redux::momfbd::gradientMap = { {"gradient_diff", CFG_GRADIENT_DIFF},
    {"gradient_Vogel", CFG_GRADIENT_VOGEL}
};

const map<string, int> redux::momfbd::fillpixMap = { {"median", CFG_FPM_MEDIAN},
    {"invdistweight", CFG_FPM_INVDISTWEIGHT},
    {"horint", CFG_FPM_HORINT}
};


#define lg Logger::mlg
namespace {

    const string thisChannel = "momfbdjob";
    static Job* createMomfbdJob( void ) {
        return new MomfbdJob();
    }

    void checkFAP( double& F, double& A, double& P ) {
        double rad2asec = 180.0 * 3600.0 / redux::PI;
        size_t count = F > 0 ? 1 : 0;
        count += A > 0 ? 1 : 0;
        count += P > 0 ? 1 : 0;
        if( count == 0 ) {
            F = DEF_TELESCOPE_F;
            A = DEF_ARCSECPERPIX;
            P = F * A / rad2asec;
        }
        else if( count == 3 ) {
            LOG_WARN << "Too many parameters given: replacing telescope focal length (" << F << ") with computed value = " << ( P / A * rad2asec );
            F = P * rad2asec / A;
        }
        else if( count == 2 ) {
            if( !( F > 0 ) ) F = P * rad2asec / A;
            else if( A > 0 ) P = F * A / rad2asec;
            else A = P / F * rad2asec;
        }
        else {
            if( !( F > 0 ) ) {
                F = DEF_TELESCOPE_F;
                if( A > 0 ) P = F * A / rad2asec;
                else A = P / F * rad2asec;
            }
            else {
                A = DEF_ARCSECPERPIX;
                P = F * A / rad2asec;
            }

        }
    }

}
size_t MomfbdJob::jobType = Job::registerJob( "momfbd", createMomfbdJob );

int MomfbdJob::getFromMap( string str, const map<string, int>& m ) {
    auto res = m.find( str );
    if( res != m.end() ) {
        return res->second;
    }
    return 0;
}

void MomfbdJob::maybeOverride( bool value, uint32_t& flagset, uint32_t flag ) {
    if( value != bool( flagset & flag ) ) {
        flagset ^= flag;
    }
}

MomfbdJob::MomfbdJob( void ) : basis( 0 ), fillpix_method( 0 ), output_data_type( 0 ), flags( 0 ), patchSize( 0 ), sequenceNumber( 0 ),
    klMinMode( 0 ), klMaxMode( 0 ), borderClip( 0 ), minIterations( 0 ), maxIterations( 0 ), nDoneMask( 0 ),
    gradient_method( 0 ), getstep_method( 0 ), max_local_shift( 0 ), mstart( 0 ), mstep( 0 ), pupilSize( 0 ),
    nPatchesX( 0 ), nPatchesY( 0 ), ncal( 0 ), telescopeFocalLength( 0 ), telescopeDiameter( 0 ), arcSecsPerPixel( 0 ),
    pixelSize( 0 ), reg_gamma( 0 ), FTOL( 0 ), EPS( 0 ), svd_reg( 0 ) {

    //LOG_DEBUG << "MomfbdJob::MomfbdJob()   (jobType = " << jobType << ")";
    info.typeString = "momfbd";
}

MomfbdJob::~MomfbdJob( void ) {
    //LOG_DEBUG << "MomfbdJob::~MomfbdJob()";
}

uint64_t MomfbdJob::unpackParts( const char* ptr, std::vector<Part::Ptr>& parts, bool swap_endian ) {

    using redux::util::unpack;
    size_t nParts;
    uint64_t count = unpack( ptr, nParts, swap_endian );
    parts.resize( nParts );
    for( auto & it : parts ) {
        it.reset( new Patch );
        count += it->unpack( ptr+count, swap_endian );
    }
    return count;
}

void MomfbdJob::parseProperties( po::variables_map& vm, bpt::ptree& tree ) {

    Job::parseProperties( vm, tree );
    LOG_DEBUG << "MomfbdJob::parseProperties()";

    bpt::ptree cmdTree;      // just to be able to use the VectorTranslator

    if( vm.count( "simx" ) ) cmdTree.put( "simx", vm["simx"].as<string>() );
    if( vm.count( "simy" ) ) cmdTree.put( "simy", vm["simy"].as<string>() );
    if( vm.count( "imgn" ) ) cmdTree.put( "imgn", vm["imgn"].as<string>() );

    if( vm.count( "output-file" ) ) cmdTree.put( "output-file", vm["output-file"].as<string>() );
    if( vm.count( "sequence" ) ) sequenceNumber = vm["sequence"].as<int>();

    string tmpString = tree.get<string>( "BASIS", "Zernike" );
    basis = iequals( tmpString, "Karhunen-Loeve" ) ? CFG_KARHUNEN_LOEVE : CFG_ZERNIKE;

    modes = tree.get<vector<uint32_t>>( "MODES", vector<uint32_t>( { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                                     20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35
                                                                   } ) );
    imageNumbers = cmdTree.get<vector<uint32_t>>( "imgn", {} );
    if( imageNumbers.size() == 0 ) imageNumbers = tree.get<vector<uint32_t>>( "IMAGE_NUM", {} );
    darkNumbers = tree.get<vector<uint32_t>>( "DARK_NUM", {} );
    sequenceNumber = tree.get<uint32_t>( "SEQUENCE_NUM", sequenceNumber );

    klMinMode  = tree.get<uint32_t>( "KL_MIN_MODE", DEF_KL_MIN_MODE );
    klMaxMode  = tree.get<uint32_t>( "KL_MAX_MODE", DEF_KL_MAX_MODE );
    patchSize    = tree.get<uint32_t>( "NUM_POINTS", DEF_NUM_POINTS );
    borderClip = tree.get<uint32_t>( "BORDER_CLIP", DEF_BORDER_CLIP );

    telescopeDiameter = tree.get<double>( "TELESCOPE_D", DEF_TELESCOPE_D );
    telescopeFocalLength = tree.get<double>( "TELESCOPE_F", 0 );
    arcSecsPerPixel = tree.get<double>( "ARCSECPERPIX", 0 );
    pixelSize = tree.get<double>( "PIXELSIZE", 0 );

    checkFAP( telescopeFocalLength, arcSecsPerPixel, pixelSize );

    LOG << "Pixel size = " << pixelSize;
    LOG << "Arcseconds per pixel = " << arcSecsPerPixel;
    LOG << "Telescope focal length = " << telescopeFocalLength;

    minIterations = tree.get<uint32_t>( "MIN_ITER", DEF_MIN_ITER );
    maxIterations = tree.get<uint32_t>( "MAX_ITER", DEF_MAX_ITER );
    nDoneMask = tree.get<uint32_t>( "N_DONE_ITER", DEF_N_DONE_ITER );
    nDoneMask = ( 1 << nDoneMask ) - 1;

    FTOL = tree.get<double>( "FTOL", DEF_FTOL );
    EPS = tree.get<double>( "EPS", DEF_EPS );
    svd_reg = tree.get<double>( "SVD_REG", DEF_SVD_REG );
    reg_gamma = tree.get<double>( "REG_GAMMA", DEF_REG_GAMMA );

    programDataDir = cleanPath( tree.get<string>( "PROG_DATA_DIR", DEF_PROG_DATA_DIR ) );
    imageDataDir = cleanPath( tree.get<string>( "IMAGE_DATA_DIR", "" ) );

    stokesWeights = tree.get<vector<double>>( "VECTOR", vector<double>() );

    getstep_method = getFromMap( tree.get<string>( "GETSTEP", DEF_GETSTEP ), getstepMap );
    if( getstep_method == 0 ) {
        LOG_ERR << "unknown getstep method \"" << tree.get<string>( "GETSTEP", DEF_GETSTEP ) << "\"\n  Valid entries currently are: "
                "\"getstep_steepest_descent\", \"getstep_conjugate_gradient\" or \"getstep_BFGS_inv\"";
    }
    gradient_method = getFromMap( tree.get<string>( "GRADIENT", DEF_GRADIENT ), gradientMap );
    if( gradient_method == 0 ) {
        LOG_ERR << "unknown gradient method \"" << tree.get<string>( "GRADIENT", DEF_GRADIENT ) << "\"\n  Valid entries currently are: "
                "\"gradient_diff\" or \"gradient_Vogel\"";
    }
    fillpix_method = getFromMap( tree.get<string>( "FPMETHOD", DEF_FPMETHOD ), fillpixMap );
    if( fillpix_method == 0 ) {
        LOG_ERR << "unknown fillpix method \"" << tree.get<string>( "FPMETHOD", DEF_FPMETHOD ) << "\"\n  Valid entries currently are: "
                "\"median\", \"invdistweight\" or \"horint\"";
    }

    time_obs = tree.get<string>( "TIME_OBS", "" );
    date_obs = tree.get<string>( "DATE_OBS", DEF_DATE_OBS );

    max_local_shift = tree.get<uint32_t>( "MAX_LOCAL_SHIFT", DEF_MAX_LOCAL_SHIFT );
    mstart = tree.get<uint32_t>( "MODE_START", DEF_MODE_START );
    mstep = tree.get<uint32_t>( "MODE_STEP", DEF_MODE_STEP );

    if( tree.get<bool>( "FLATFIELD", false ) && tree.get<bool>( "CALIBRATE", false ) ) {
        LOG_ERR << "both FLATFIELD and CALIBRATE mode requested";
    }

    tmpString = cleanPath( tree.get<string>( "PUPIL", "" ), imageDataDir );
    if( tmpString != "" ) {
        readFile( tmpString, pupil );
        if( pupil.nDimensions() != 2 ) {
            LOG_ERR << "pupil file \"" << tmpString << "\" not 2 dimensional.";
        }
        pupilSize = pupil.dimSize(0);
        uint32_t tmpInt = pupil.dimSize(1);
        if( pupilSize != tmpInt ) {
            LOG_ERR << "pupil file \"" << tmpString << "\" not square (" << pupilSize << "x" << tmpInt << ")";
        }
        pupil = Array<double>( pupilSize, pupilSize );
    }

    flags = 0;
    if( tree.get<bool>( "CALIBRATE", false ) )             flags |= MFBD_CALIBRATE;
    if( tree.get<bool>( "DONT_MATCH_IMAGE_NUMS", false ) ) flags |= MFBD_DONT_MATCH_IMAGE_NUMS;
    if( tree.get<bool>( "FAST_QR", false ) )               flags |= MFBD_FAST_QR;
    if( tree.get<bool>( "FIT_PLANE", false ) )             flags |= MFBD_FIT_PLANE;
    if( tree.get<bool>( "FLATFIELD", false ) )             flags |= MFBD_FLATFIELD;
    if( tree.get<bool>( "GET_ALPHA", false ) || ( flags & MFBD_CALIBRATE ) ) flags |= MFBD_GET_ALPHA;
    if( tree.get<bool>( "GET_COBJ", false ) )              flags |= MFBD_GET_COBJ;
    if( tree.get<bool>( "GET_DIVERSITY", false ) )         flags |= MFBD_GET_DIVERSITY;
    if( tree.get<bool>( "GET_METRIC", false ) )            flags |= MFBD_GET_METRIC;
    if( tree.get<bool>( "GET_MODES", false ) )             flags |= MFBD_GET_MODES;
    if( tree.get<bool>( "GET_PSF", false ) )               flags |= MFBD_GET_PSF;
    if( tree.get<bool>( "GET_PSF_AVG", false ) )           flags |= MFBD_GET_PSF_AVG;
    if( tree.get<bool>( "GET_RESIDUAL", false ) )          flags |= MFBD_GET_RESIDUAL;
    if( tree.get<bool>( "GLOBAL_NOISE", false ) )          flags |= MFBD_GLOBAL_NOISE;
    if( tree.get<bool>( "NEW_CONSTRAINTS", false ) )       flags |= MFBD_NEW_CONSTRAINTS;
    if( tree.get<bool>( "NO_CLIP", false ) )               flags |= MFBD_NO_CLIP;
    if( tree.get<bool>( "NO_CONSTRAINTS", false ) )        flags |= MFBD_NO_CONSTRAINTS;
    if( tree.get<bool>( "NO_FILTER", false ) )             flags |= MFBD_NO_FILTER;
    if( tree.get<bool>( "SAVE_FFDATA", false ) )           flags |= MFBD_SAVE_FFDATA;

    if( vm.count( "force" ) ) flags |= MFBD_FORCE_WRITE;
    if( vm.count( "swap" ) )  flags |= MFBD_SWAP;

    if( ( flags & MFBD_CALIBRATE ) && ( flags & MFBD_FLATFIELD ) ) {
        LOG_WARN << "both FLATFIELD and CALIBRATE mode requested, forcing CALIBRATE";
        flags &= ~MFBD_FLATFIELD;
    }
    if( ( flags & MFBD_GET_PSF ) && ( flags & MFBD_GET_PSF_AVG ) ) {
        LOG_WARN << "both GET_PSF and GET_PSF_AVG mode requested";
    }
    if( ( flags & MFBD_CALIBRATE ) && ( flags & MFBD_NEW_CONSTRAINTS ) ) {
        LOG_WARN << "calibration mode uses old style constraints, ignoring NEW_CONSTRAINTS";
        flags &= ~MFBD_NEW_CONSTRAINTS;
    }

    subImagePosX = cmdTree.get<vector<uint32_t>>( "simx", {} );
    if( subImagePosX.size() == 0 ) subImagePosX = tree.get<vector<uint32_t>>( "SIM_X", {} );
    subImagePosY = cmdTree.get<vector<uint32_t>>( "simy", {} );
    if( subImagePosY.size() == 0 ) subImagePosY = tree.get<vector<uint32_t>>( "SIM_Y", {} );

    if( tree.get<bool>( "CAL_X", false ) ) {
        if( tree.get<bool>( "CAL_Y", false ) ) {
            if( subImagePosX.size() || subImagePosY.size() ) LOG << "Note: SIM_X/SIM_Y replaced by CAL_X/CAL_Y";
            subImagePosX = tree.get<vector<uint32_t>>( "CAL_X", {} );
            subImagePosY = tree.get<vector<uint32_t>>( "CAL_Y", {} );
            if( subImagePosX.empty() || ( subImagePosX.size() != subImagePosY.size() ) ) {
                LOG_ERR << "CAL_X and CAL_Y must have the same number of elements!";
            }
            else ncal = subImagePosX.size();
        }
        else LOG_ERR << "CAL_Y must be provided if CAL_X is!";
    }
    LOG << "MomfbdJob::parseProps()  " << printArray(subImagePosX,"SubX");
    LOG << "MomfbdJob::parseProps()  " << printArray(subImagePosY,"SubY");

    auto filetypes = tree.get<vector<FileType>>( "FILE_TYPE", vector<FileType>( 1, flags & MFBD_CALIBRATE ? FT_ANA : FT_FITS ) );
    for( auto & it : filetypes ) flags |= it;
    if( ( flags & MFBD_FT_MASK ) == 0 ) LOG_ERR << "\"FILE_TYPE\" has to be one of ANA/FITS/MOMFBD.";

    tmpString = tree.get<string>( "DATA_TYPE", DEF_DATA_TYPE );
    if( iequals( tmpString, "FLOAT" ) ) output_data_type = MFBD_F32T;
    else if( iequals( tmpString, "SHORT" ) ) output_data_type = MFBD_I16T;
    else {
        LOG_WARN << "\"DATA_TYPE\" unrecognized data type \"" << tmpString << "\", using INT16";
        output_data_type = MFBD_I16T;
    }

    tmpString = tree.get<string>( "OUTPUT_FILES", "" );
    boost::split( outputFiles, tmpString, boost::is_any_of( "," ) );

    size_t nObj( 0 );
    for( auto & it : tree ) {
        if( iequals( it.first, "OBJECT" ) ) {
            tmpString = ( nObj >= outputFiles.size() ) ? "" : outputFiles[nObj++];
            Object* tmpObj = new Object( *this );
            tmpObj->parseProperties( it.second, tmpString );
            objects.push_back( shared_ptr<Object>( tmpObj ) );
        }
    }
    if( outputFiles.size() > objects.size() ) {
        LOG_WARN << outputFiles.size() << " output file names specified but only " << objects.size() << " objects found.";
    }
    LOG_DEBUG << "MomfbdJob::parseProperties() done.";

}

bpt::ptree MomfbdJob::getPropertyTree( bpt::ptree* root ) {

    bpt::ptree tree = Job::getPropertyTree();         // get common properties

    for( auto & it : objects ) {
        it->getPropertyTree( &tree );
    }

    // TODO BASIS
    tree.put( "MODES", modes );
    if( !imageNumbers.empty() ) tree.put( "IMAGE_NUM", imageNumbers );
    if( !darkNumbers.empty() ) tree.put( "DARK_NUM", darkNumbers );
    tree.put( "SEQUENCE_NUM", sequenceNumber );
    if( klMinMode != DEF_KL_MIN_MODE ) tree.put( "KL_MIN_MODE", klMinMode );
    if( klMaxMode != DEF_KL_MAX_MODE ) tree.put( "KL_MAX_MODE", klMaxMode );
    if( patchSize != DEF_NUM_POINTS ) tree.put( "NUM_POINTS", patchSize );
    if( borderClip != DEF_BORDER_CLIP ) tree.put( "BORDER_CLIP", borderClip );
    if( telescopeDiameter != DEF_TELESCOPE_D ) tree.put( "TELESCOPE_D", telescopeDiameter );
    if( telescopeFocalLength != 0 ) tree.put( "TELESCOPE_F", telescopeFocalLength );
    if( arcSecsPerPixel != 0 ) tree.put( "ARCSECPERPIX", arcSecsPerPixel );
    if( pixelSize != 0 ) tree.put( "PIXELSIZE", pixelSize );
    if( minIterations != DEF_MIN_ITER ) tree.put( "MIN_ITER", minIterations );
    if( maxIterations != DEF_MAX_ITER ) tree.put( "MAX_ITER", maxIterations );
    if( nDoneMask != DEF_N_DONE_ITER ) tree.put( "N_DONE_ITER", nDoneMask );
    if( FTOL != DEF_FTOL ) tree.put( "FTOL", FTOL );
    if( EPS != DEF_EPS ) tree.put( "EPS", EPS );
    if( svd_reg != DEF_SVD_REG ) tree.put( "SVD_REG", svd_reg );
    if( reg_gamma != DEF_REG_GAMMA ) tree.put( "REG_GAMMA", reg_gamma );
    if( programDataDir != DEF_PROG_DATA_DIR ) tree.put( "PROG_DATA_DIR", programDataDir );
    if( !imageDataDir.empty() ) tree.put( "IMAGE_DATA_DIR", imageDataDir );
    if( !stokesWeights.empty() ) tree.put( "VECTOR", stokesWeights );
    //if(getstep_method!=DEF_GETSTEP) tree.put("GETSTEP",getstep_method);
    //if(gradient_method!=DEF_GRADIENT) tree.put("GRADIENT",gradient_method);
    //if(fillpix_method!=DEF_FPMETHOD) tree.put("FPMETHOD",fillpix_method);
    if( !time_obs.empty() ) tree.put( "TIME_OBS", time_obs );
    if( date_obs != DEF_DATE_OBS ) tree.put( "DATE_OBS", date_obs );
    if( max_local_shift != DEF_MAX_LOCAL_SHIFT ) tree.put( "MAX_LOCAL_SHIFT", max_local_shift );
    if( mstart != DEF_MODE_START ) tree.put( "MODE_START", mstart );
    if( mstep != DEF_MODE_STEP ) tree.put( "MODE_STEP", mstep );
    if( !outputFiles.empty() ) tree.put( "OUTPUT_FILES", outputFiles );
    // TODO PUPIL/SIM_X/SIM_Y/CAL_X/CAL_Y/FILE_TYPE/DATA_TYPE

    if( flags & MFBD_CALIBRATE ) tree.put( "CALIBRATE", ( bool )( flags & MFBD_CALIBRATE ) );
    if( flags & MFBD_DONT_MATCH_IMAGE_NUMS ) tree.put( "DONT_MATCH_IMAGE_NUMS", ( bool )( flags & MFBD_DONT_MATCH_IMAGE_NUMS ) );
    if( flags & MFBD_FAST_QR ) tree.put( "FAST_QR", ( bool )( flags & MFBD_FAST_QR ) );
    if( flags & MFBD_FIT_PLANE ) tree.put( "FIT_PLANE", ( bool )( flags & MFBD_FIT_PLANE ) );
    if( flags & MFBD_FLATFIELD ) tree.put( "FLATFIELD", ( bool )( flags & MFBD_FLATFIELD ) );
    if( flags & MFBD_GET_ALPHA ) tree.put( "GET_ALPHA", ( bool )( flags & MFBD_GET_ALPHA ) );
    if( flags & MFBD_GET_COBJ ) tree.put( "GET_COBJ", ( bool )( flags & MFBD_GET_COBJ ) );
    if( flags & MFBD_GET_DIVERSITY ) tree.put( "GET_DIVERSITY", ( bool )( flags & MFBD_GET_DIVERSITY ) );
    if( flags & MFBD_GET_METRIC ) tree.put( "GET_METRIC", ( bool )( flags & MFBD_GET_METRIC ) );
    if( flags & MFBD_GET_MODES ) tree.put( "GET_MODES", ( bool )( flags & MFBD_GET_MODES ) );
    if( flags & MFBD_GET_PSF ) tree.put( "GET_PSF", ( bool )( flags & MFBD_GET_PSF ) );
    if( flags & MFBD_GET_PSF_AVG ) tree.put( "GET_PSF_AVG", ( bool )( flags & MFBD_GET_PSF_AVG ) );
    if( flags & MFBD_GET_RESIDUAL ) tree.put( "GET_RESIDUAL", ( bool )( flags & MFBD_GET_RESIDUAL ) );
    if( flags & MFBD_GLOBAL_NOISE ) tree.put( "GLOBAL_NOISE", ( bool )( flags & MFBD_GLOBAL_NOISE ) );
    if( flags & MFBD_NEW_CONSTRAINTS ) tree.put( "NEW_CONSTRAINTS", ( bool )( flags & MFBD_NEW_CONSTRAINTS ) );
    if( flags & MFBD_NO_CLIP ) tree.put( "NO_CLIP", ( bool )( flags & MFBD_NO_CLIP ) );
    if( flags & MFBD_NO_CONSTRAINTS ) tree.put( "NO_CONSTRAINTS", ( bool )( flags & MFBD_NO_CONSTRAINTS ) );
    if( flags & MFBD_NO_FILTER ) tree.put( "NO_FILTER", ( bool )( flags & MFBD_NO_FILTER ) );
    if( flags & MFBD_SAVE_FFDATA ) tree.put( "SAVE_FFDATA", ( bool )( flags & MFBD_SAVE_FFDATA ) );


    if( root ) {
        root->push_back( bpt::ptree::value_type( "momfbd", tree ) );
    }

    return tree;
}

size_t MomfbdJob::size( void ) const {
    size_t sz = Job::size();
    sz += 3;                                    // basis, fillpix_method, output_data_type;
    sz += 18 * sizeof( uint32_t );
    sz += 8 * sizeof( double );
    sz += modes.size() * sizeof( uint32_t ) + sizeof( size_t );
    sz += imageNumbers.size() * sizeof( uint32_t ) + sizeof( size_t );
    sz += darkNumbers.size() * sizeof( uint32_t ) + sizeof( size_t );
    sz += subImagePosX.size() * sizeof( uint32_t ) + sizeof( size_t );
    sz += subImagePosY.size() * sizeof( uint32_t ) + sizeof( size_t );
    sz += stokesWeights.size() * sizeof( double ) + sizeof( size_t );
    sz += imageDataDir.length() + programDataDir.length() + time_obs.length() + date_obs.length() + 4;        // strings + \0
    sz += 2 * sizeof( size_t );           // outputFiles.size() + objects.size()
    for( auto & it : outputFiles ) {
        sz += it.length() + 1;
    }
    for( auto & it : objects ) {
        sz += it->size();
    }
    sz += pupil.size();
    return sz;
}

uint64_t MomfbdJob::pack( char* ptr ) const {
    
    using redux::util::pack;
    
    uint64_t count = Job::pack( ptr );
    count += pack( ptr+count, basis );
    count += pack( ptr+count, fillpix_method );
    count += pack( ptr+count, output_data_type );
    count += pack( ptr+count, flags );
    count += pack( ptr+count, patchSize );
    count += pack( ptr+count, sequenceNumber );
    count += pack( ptr+count, klMinMode );
    count += pack( ptr+count, klMaxMode );
    count += pack( ptr+count, borderClip );
    count += pack( ptr+count, minIterations );
    count += pack( ptr+count, maxIterations );
    count += pack( ptr+count, nDoneMask );
    count += pack( ptr+count, gradient_method );
    count += pack( ptr+count, getstep_method );
    count += pack( ptr+count, max_local_shift );
    count += pack( ptr+count, mstart );
    count += pack( ptr+count, mstep );
    count += pack( ptr+count, pupilSize );
    count += pack( ptr+count, nPatchesX );
    count += pack( ptr+count, nPatchesY );
    count += pack( ptr+count, ncal );
    count += pack( ptr+count, telescopeFocalLength );
    count += pack( ptr+count, telescopeDiameter );
    count += pack( ptr+count, arcSecsPerPixel );
    count += pack( ptr+count, pixelSize );
    count += pack( ptr+count, reg_gamma );
    count += pack( ptr+count, FTOL );
    count += pack( ptr+count, EPS );
    count += pack( ptr+count, svd_reg );
    count += pack( ptr+count, modes );
    count += pack( ptr+count, imageNumbers );
    count += pack( ptr+count, darkNumbers );
    count += pack( ptr+count, subImagePosX );
    count += pack( ptr+count, subImagePosY );
    count += pack( ptr+count, stokesWeights );
    count += pack( ptr+count, imageDataDir );
    count += pack( ptr+count, programDataDir );
    count += pack( ptr+count, time_obs );
    count += pack( ptr+count, date_obs );
    count += pack( ptr+count, outputFiles.size() );
    for( auto & it : outputFiles ) {
        count += pack( ptr+count, it );
    }
    count += pack( ptr+count, objects.size() );
    for( auto & it : objects ) {
        count += it->pack( ptr+count );
    }
    count += pupil.pack( ptr+count);
    
    return count;
    
}

uint64_t MomfbdJob::unpack( const char* ptr, bool swap_endian ) {
    
    using redux::util::unpack;
    
    uint64_t count = Job::unpack( ptr, swap_endian );
    count += unpack( ptr+count, basis );
    count += unpack( ptr+count, fillpix_method );
    count += unpack( ptr+count, output_data_type );
    count += unpack( ptr+count, flags, swap_endian );
    count += unpack( ptr+count, patchSize, swap_endian );
    count += unpack( ptr+count, sequenceNumber, swap_endian );
    count += unpack( ptr+count, klMinMode, swap_endian );
    count += unpack( ptr+count, klMaxMode, swap_endian );
    count += unpack( ptr+count, borderClip, swap_endian );
    count += unpack( ptr+count, minIterations, swap_endian );
    count += unpack( ptr+count, maxIterations, swap_endian );
    count += unpack( ptr+count, nDoneMask, swap_endian );
    count += unpack( ptr+count, gradient_method, swap_endian );
    count += unpack( ptr+count, getstep_method, swap_endian );
    count += unpack( ptr+count, max_local_shift, swap_endian );
    count += unpack( ptr+count, mstart, swap_endian );
    count += unpack( ptr+count, mstep, swap_endian );
    count += unpack( ptr+count, pupilSize, swap_endian );
    count += unpack( ptr+count, nPatchesX, swap_endian );
    count += unpack( ptr+count, nPatchesY, swap_endian );
    count += unpack( ptr+count, ncal, swap_endian );
    count += unpack( ptr+count, telescopeFocalLength, swap_endian );
    count += unpack( ptr+count, telescopeDiameter, swap_endian );
    count += unpack( ptr+count, arcSecsPerPixel, swap_endian );
    count += unpack( ptr+count, pixelSize, swap_endian );
    count += unpack( ptr+count, reg_gamma, swap_endian );
    count += unpack( ptr+count, FTOL, swap_endian );
    count += unpack( ptr+count, EPS, swap_endian );
    count += unpack( ptr+count, svd_reg, swap_endian );
    count += unpack( ptr+count, modes, swap_endian );
    count += unpack( ptr+count, imageNumbers, swap_endian );
    count += unpack( ptr+count, darkNumbers, swap_endian );
    count += unpack( ptr+count, subImagePosX, swap_endian );
    count += unpack( ptr+count, subImagePosY, swap_endian );
    count += unpack( ptr+count, stokesWeights, swap_endian );
    count += unpack( ptr+count, imageDataDir );
    count += unpack( ptr+count, programDataDir );
    count += unpack( ptr+count, time_obs );
    count += unpack( ptr+count, date_obs );
    size_t tmp;
    count += unpack( ptr+count, tmp, swap_endian );
    outputFiles.resize( tmp );
    for( auto & it : outputFiles ) {
        count += unpack( ptr+count, it, swap_endian );
    }
    count += unpack( ptr+count, tmp, swap_endian );
    objects.resize( tmp );
    for( auto & it : objects ) {
        it.reset(new Object(*this));
        count += it->unpack( ptr+count, swap_endian );
    }
    count += pupil.unpack( ptr+count, swap_endian );
    
    return count;
    
}



void MomfbdJob::checkParts( void ) {

    uint8_t mask = 0;
    for( auto & it : patches ) {
        /*if( it.second->step & JSTEP_ERR && (it.second->nRetries<info.maxPartRetries)) {    // TODO: handle failed parts.
            it.second->nRetries++;
            it.second->step &= ~JSTEP_ERR;
        }*/
        mask |= it.second->step;
    }

    if( mask & JSTEP_ERR ) {    // TODO: handle failed parts.

    }

    if( countBits( mask ) == 1 ) {  // if all parts have the same "step", set the whole job to that step.
        info.step.store( mask );
    }

}


size_t MomfbdJob::getParts( WorkInProgress& wip, uint8_t nThreads ) {

    uint8_t step = info.step.load();
    wip.parts.clear();
    if( step == JSTEP_QUEUED || step == JSTEP_RUNNING ) {
        unique_lock<mutex> lock( jobMutex );
//         size_t nParts = wip.peer->status.nThreads;
//         if( info.nThreads ) nParts = std::min( wip.peer->status.nThreads, info.nThreads );
        for( auto & it : patches ) {
            if( it.second->step == JSTEP_QUEUED ) {
                it.second->step = JSTEP_RUNNING;
                wip.parts.push_back( it.second );
                info.step.store( JSTEP_RUNNING );
                info.state.store( JSTATE_ACTIVE );
                if( wip.parts.size() > 0 ) break;
            }
        }
        checkParts();
    }
    return wip.parts.size();
}


void MomfbdJob::ungetParts( WorkInProgress& wip ) {
    unique_lock<mutex> lock( jobMutex );
    for( auto & it : wip.parts ) {
        it->step = JSTEP_QUEUED;
    }
    wip.parts.clear();
}


void MomfbdJob::returnParts( WorkInProgress& wip ) {
    unique_lock<mutex> lock( jobMutex );
    checkParts();
    for( auto & it : wip.parts ) {
        auto patch = static_pointer_cast<Patch>( it );
        patches[it->id]->step = patch->step;
        //patches[it->id]->result = patch->result;
    }
    wip.parts.clear();
    checkParts();
}


void MomfbdJob::init(void) {
    
    //redux::image::KL_cfg* coeff = nullptr;
    if( klMinMode || klMaxMode ) {
        LOG << "calculating Karhunen-Loeve coefficients";
        //coeff = legacy::klConfig(klMinMode,klMaxMode);
    }
    
    for( auto& it : objects ) {
        it->init(); //coeff);
    }
}


void MomfbdJob::cleanup(void) {
    for( auto& it : objects ) {
        it->cleanup();
    }
}


bool MomfbdJob::run( WorkInProgress& wip, boost::asio::io_service& service, boost::thread_group& pool, uint8_t maxThreads ) {
    
    uint8_t step = info.step.load();
    if( step < JSTEP_SUBMIT ) {
        info.step.store( JSTEP_SUBMIT );        // do nothing before submitting
        return true;                            // run again
    }
    else if( step == JSTEP_RECEIVED ) {
        preProcess(service, pool);                           // preprocess on master: load, flatfield, split in patches
    }
    else if( step == JSTEP_RUNNING || step == JSTEP_QUEUED ) {          // main processing
        size_t nThreads = std::min( maxThreads, info.maxThreads);
        service.reset();

        for( auto & it : wip.parts ) {
            service.post( boost::bind( &MomfbdJob::runMain, this, boost::ref( it ) ) );
        }
        for( size_t t = 0; t < nThreads; ++t ) {
            pool.create_thread( boost::bind( &boost::asio::io_service::run, &service ) );
        }

        pool.join_all();

    }
    else if( step == JSTEP_POSTPROCESS ) {
        postProcess(service, pool);                          // postprocess on master, collect results, save...
    }
    else {
        LOG << "MomfbdJob::run()  unrecognized step = " << ( int )info.step.load();
        info.step.store( JSTEP_ERR );
    }
    return false;
    
}



void MomfbdJob::preProcess( boost::asio::io_service& service, boost::thread_group& pool ) {

    // TODO: start logging (to file)

    LOG_TRACE << "MomfbdJob::preProcess()";
    
    if ( !checkData() ) {
        LOG_ERR << "MomfbdJob::preProcess(): sanity check failed.";
        info.step.store( JSTEP_ERR );
        info.state.store( JSTATE_IDLE );
        return;
    }
    
    // load shared files synchronously (dark,gain,psf,offset...)
    service.reset();
    for( auto & it : objects ) {
        it->loadData(service, pool);
    }

    info.maxThreads = 12;
    // load remaining files asynchronously  (images)
    for( size_t t = 0; t < info.maxThreads; ++t ) {
        LOG_TRACE << "Adding load-thread #" << (t+1);
        pool.create_thread( boost::bind( &boost::asio::io_service::run, &service ) );
    }
    pool.join_all();

    // Done loading files -> start the preprocessing (flatfielding etc.)
    
    service.reset();
    Point imageSizes;
    for( auto & it : objects ) {
        Point tmp = it->clipImages();
        if(imageSizes.x == 0) {
            imageSizes = tmp;
        } else if( tmp != imageSizes ) {
            throw std::logic_error("The clipped images have different sizes for the different objects, please verify the ALIGN_CLIP values.");
        }
        it->preprocessData(service, pool);
    }

    for( size_t t = 0; t < info.maxThreads; ++t ) {
        pool.create_thread( boost::bind( &boost::asio::io_service::run, &service ) );
    }
    pool.join_all();

    // Done pre-processing -> normalize within each object
    
    service.reset();
    for( auto & it : objects ) {
        it->normalize(service, pool);
    }

    for( size_t t = 0; t < info.maxThreads; ++t ) {
        pool.create_thread( boost::bind( &boost::asio::io_service::run, &service ) );
    }
    pool.join_all();

    // Done normalizing -> split in patches
    
    int minimumOverlap = 16;                // desired width of blending zone in pixels
    //int patchSeparation = 3 * patchSize / 4 - minimumOverlap; // target separation
    if( subImagePosX.empty() ) {
        // TODO: do we need to handle single patches ?? (this only supports nPatches >= 2)
        // TODO: verify michiel's splitting method
        int firstPos = max_local_shift + patchSize/2;
        int lastPos = imageSizes.x - firstPos - 1;
        nPatchesX = 2;
        double separation = (lastPos-firstPos)/static_cast<double>(nPatchesX-1);
        double overlap = std::max(patchSize-separation,0.0);
        while(overlap < minimumOverlap) {
            ++nPatchesX;
            separation = (lastPos-firstPos)/static_cast<double>(nPatchesX-1);
            overlap = std::max(patchSize-separation,0.0);
        }
        for( size_t i = 0; i < nPatchesX; ++i ) {
            subImagePosX.push_back(static_cast<uint32_t>( i*separation + firstPos ) );
        }
        LOG << "MomfbdJob::preProcess(): Generated patch positions  " << printArray(subImagePosX,"x-pos");
    }

    if( subImagePosY.empty() ) {
        int firstPos = max_local_shift + patchSize/2;
        int lastPos = imageSizes.y - firstPos - 1;
        nPatchesY = 2;
        double separation = (lastPos-firstPos)/static_cast<double>(nPatchesY-1);
        double overlap = std::max(patchSize-separation,0.0);
        while(overlap < minimumOverlap) {
            ++nPatchesY;
            separation = (lastPos-firstPos)/static_cast<double>(nPatchesY-1);
            overlap = std::max(patchSize-separation,0.0);
        }
        for( size_t i = 0; i < nPatchesY; ++i ) {
            subImagePosY.push_back(static_cast<uint32_t>( i*separation + firstPos ) );
        }
        LOG << "MomfbdJob::preProcess(): Generated patch positions  " << printArray(subImagePosY,"y-pos");
    }
 
    if( subImagePosX.empty() || subImagePosY.empty() ) {
        LOG_ERR << "MomfbdJob::preProcess(): No patches specified or generated, can't continue.";
        info.step.store( JSTEP_ERR );
        info.state.store( JSTATE_IDLE );
        return;
    }
    

   
    size_t count = 0;
    uint32_t yid=0,xid=0;
    service.reset();
    uint32_t span = patchSize/2 + max_local_shift;
    for( auto posY : subImagePosY ) {
        uint32_t trimmedPosY = std::min(std::max(span,posY),imageSizes.y-span);
        if( trimmedPosY != posY ) LOG_WARN << "MomfbdJob::preProcess() y-position of patch was outside the image area and was trimmed: " << posY << " -> " << trimmedPosY;
        for( auto posX : subImagePosX ) {
            uint32_t trimmedPosX = std::min(std::max(span,posX),imageSizes.x);
            if( trimmedPosX != posX ) LOG_WARN << "MomfbdJob::preProcess() x-position of patch was outside the image area and was trimmed: " << posX << " -> " << trimmedPosX;
            Patch::Ptr patch( new Patch() );
            patch->first.x = trimmedPosX-span;
            patch->first.y = trimmedPosY-span;
            patch->last.x = patch->first.x+2*span-1;
            patch->last.y = patch->first.y+2*span-1;
            patch->id = ++count;
            //cout << "PATCH: " << count << "  first.x = " << patch->first.x << "  first.y = " << patch->first.y << "  last.x = " << patch->last.x << "  last.y. = " << patch->last.y;
            //cout << "  szX = " << (patch->last.x-patch->first.x+1) << "  szY = " << (patch->last.y-patch->first.y+1) << "  npix = " << patch->nPixels()  << endl;
            //cout << "  posX = " << posX << "  posY = " << posY << "  trimmedPosX = " << trimmedPosX;
            //cout << "  trimmedPosY = " << trimmedPosY << "  span = " << span << "  patchSize = " << patchSize << endl;
            patch->setIndex( yid, xid++ );
            service.post( boost::bind( &MomfbdJob::packPatch, this, patch ) );
            patches.insert( make_pair( patch->id, patch ) );
        }
        yid++;
    }

    LOG_DETAIL << "MomfbdJob::preProcess()  nPatches = " << patches.size();

    for( size_t t = 0; t < 1/*info.nThreads*/; ++t ) {
        pool.create_thread( boost::bind( &boost::asio::io_service::run, &service ) );
    }
    pool.join_all();

    info.step.store( JSTEP_QUEUED );

}

void MomfbdJob::packPatch( Patch::Ptr patch ) {
    
    LOG_TRACE << "MomfbdJob::packPatch() #" << patch->id;
    size_t totalPatchSize(0);
    for( auto & it : objects ) {
        totalPatchSize += it->sizeOfPatch(patch->nPixels());
    }

    patch->dataSize = totalPatchSize;
    patch->data = sharedArray<char>(totalPatchSize);
    char* ptr = patch->data.get();
    uint64_t count(0);
    for( auto & it : objects ) {
        count += it->packPatch(patch,ptr+count);
    }
    
    if(count != totalPatchSize) {
        LOG_WARN << "Estimation of patch data-size was wrong:  est = " << totalPatchSize << "  real = " << ptrdiff_t(ptr-patch->data.get());
    }
    // TODO: compress and store in swapfile
}
 

void MomfbdJob::runMain( Part::Ptr& part ) {

    auto pptr = static_pointer_cast<Patch>( part );


//    LOG << "MomfbdJob::runMain()";
//     // temporaries, to avoid cache collisions.
//     uint32_t sizeX = pptr->xPixelH - pptr->xPixelL + 1;
//     uint32_t sizeY = pptr->yPixelH - pptr->yPixelL + 1;
//     double stepX = ( pptr->endX - pptr->beginX ) / ( sizeX - 1 );
//     double stepY = ( pptr->endY - pptr->beginY ) / ( sizeY - 1 );
//     double beginX = pptr->beginX;
//     double beginY = pptr->beginY;
//
//     size_t id = pptr->id;
//     size_t sid = pptr->sortedID;
//     uint32_t max_iters = maxIterations;
//
//     int32_t pid = getpid();
//
//     auto tmp = sharedArray<int64_t>( sizeY, sizeX );
//     auto ptr = tmp.get();
//     double x, y;
//     for( uint32_t ix = 0; ix < sizeX; ++ix ) {
//         x = beginX + ix * stepX;
//         for( uint32_t iy = 0; iy < sizeY; ++iy ) {
//             y = beginY + iy * stepY;
//
//             ptr[iy][ix] = mandelbrot( complex<double>( x, y ), max_iters );
//
//             if( ptr[iy][ix] < 0 ) continue;
//
//             if( ix < iy ) {                                 // top-left triangle showing the real part-ID (should increase upwards and to the right)
//                 ptr[iy][ix] = sid;
//             }
//             else if( ix > ( sizeY - iy ) ) {                // right triangle: the unsorted part-ID (=processing order)
//                 ptr[iy][ix] = id;
//             }
//             else  {                                         // bottom left triangle: pid, to distinguish parts processed on different machines or instances.
//                 ptr[iy][ix] = pid;
//             }
//         }
//     }
//
//     pptr->result.reset( sizeY, sizeX );
//     memcpy( pptr->result.ptr(), tmp.get()[0], sizeY * sizeX * sizeof( int64_t ) );

    //sleep(1);
    part->step = JSTEP_POSTPROCESS;

}



void MomfbdJob::postProcess( boost::asio::io_service& service, boost::thread_group& pool ) {

    LOG << "MomfbdJob::postProcess()";
    
//     auto image = sharedArray<int16_t>( ySize, xSize );
//     int16_t** img = image.get();
//
//     int64_t minPID, maxPID, minID, maxID, minSID, maxSID;
//     minPID = minID = minSID = UINT32_MAX;
//     maxPID = maxID = maxSID = 0;
//     for( auto & it : jobParts ) {
//
//         auto ptr = static_pointer_cast<DebugPart>( it.second );
//
//         uint32_t sizeX = ptr->xPixelH - ptr->xPixelL + 1;
//         uint32_t sizeY = ptr->yPixelH - ptr->yPixelL + 1;
//
//         auto blaha = reshapeArray( ptr->result.ptr( 0 ), sizeY, sizeX );
//         auto res = blaha.get();
//
//         for( uint32_t ix = 0; ix < sizeX; ++ix ) {
//             for( uint32_t iy = 0; iy < sizeY; ++iy ) {
//                 int64_t tmp = res[iy][ix];
//                 if( tmp < 0 ) {
//                     continue;      // to skip the contour for the normalization
//                 }
//                 if( ix < iy ) {
//                     if( tmp > maxSID ) maxSID = tmp;
//                     if( tmp < minSID ) minSID = tmp;
//                 }
//                 else if( ix > ( sizeY - iy ) ) {
//                     if( tmp > maxID ) maxID = tmp;
//                     if( tmp < minID ) minID = tmp;
//                 }
//                 else {
//                     if( tmp > maxPID ) maxPID = tmp;
//                     if( tmp < minPID ) minPID = tmp;
//                 }
//             }
//         }
//     }
//
//     for( auto & it : jobParts ) {
//
//         auto ptr = static_pointer_cast<DebugPart>( it.second );
//
//         uint32_t sizeX = ptr->xPixelH - ptr->xPixelL + 1;
//         uint32_t sizeY = ptr->yPixelH - ptr->yPixelL + 1;
//
//         auto blaha = reshapeArray( ptr->result.ptr( 0 ), sizeY, sizeX );
//         auto res = blaha.get();
//
//         for( uint32_t ix = 0; ix < sizeX; ++ix ) {
//             for( uint32_t iy = 0; iy < sizeY; ++iy ) {
//                 size_t tmp = res[iy][ix];
//
//                 if( tmp < 0 ) {
//                     img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
//                     continue;
//                 }
//
//                 if( ix < iy ) {
//                     if( maxSID == minSID ) img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
//                     else img[ptr->yPixelL + iy][ptr->xPixelL + ix] = ( tmp - minSID + 1 ) * 1.0 / ( maxSID - minSID + 1 ) * INT16_MAX;
//                 }
//                 else if( ix > ( sizeY - iy ) ) {
//                     if( maxID == minID ) img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
//                     else img[ptr->yPixelL + iy][ptr->xPixelL + ix] = ( tmp - minID + 1 ) * 1.0 / ( maxID - minID + 1 ) * INT16_MAX;
//                 }
//                 else {
//                     if( maxPID == minPID ) img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
//                     else img[ptr->yPixelL + iy][ptr->xPixelL + ix] = ( tmp - minPID + 1 ) * 1.0 / ( maxPID - minPID + 1 ) * INT16_MAX;
//                 }
//
//             }
//         }
//
//     }
//
//
//     Ana::Ptr hdr( new Ana() );
//
//     hdr->m_ExtendedHeader = "DebugJob";
//     hdr->m_Header.datyp = Ana::ANA_WORD;
//
//     hdr->m_Header.ndim = 2;
//     hdr->m_Header.dim[0] = xSize;
//     hdr->m_Header.dim[1] = ySize;
//
//     std::ofstream file( "debugjob_output.f0" );
//
//     Ana::write( file, reinterpret_cast<char*>( *img ), hdr );

    info.step.store( JSTEP_COMPLETED );
    info.state.store( JSTATE_IDLE );

}


bool MomfbdJob::checkCfg(void) {
    
    if( objects.empty() ) return false;     // nothing to do
    
    for( auto & it : objects ) {
        if( !it->checkCfg() ) return false;
    }
    
    return true;
}


bool MomfbdJob::checkData(void) {

    for( auto & it : objects ) {
        if( !it->checkData() ) return false;
    }
    
    return true;
}
        


