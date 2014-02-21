#include "redux/momfbd/momfbdjob.hpp"

#include "redux/momfbd/defines.hpp"
#include "redux/momfbd/object.hpp"

#include "redux/translators.hpp"
#include "redux/file/fileana.hpp"
#include "redux/constants.hpp"
#include "redux/logger.hpp"

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

void MomfbdJob::maybeOverride( bool value, uint32_t& set, uint32_t flag ) {
    if( value != ( set & flag ) ) {
        set ^= flag;
    }
}

MomfbdJob::MomfbdJob( void ) : sequenceNumber( 0 ) {
    //LOG_DEBUG << "MomfbdJob::MomfbdJob()   (jobType = " << jobType << ")";
    info.typeString = "momfbd";
}

MomfbdJob::~MomfbdJob( void ) {
    //LOG_DEBUG << "MomfbdJob::~MomfbdJob()";
}

void MomfbdJob::parseProperties( po::variables_map& vm, bpt::ptree& tree ) {

    LOG_DEBUG << "MomfbdJob::parseProperties()";

    bpt::ptree cmdTree;      // just to be able to use the VectorTranslator

    if( vm.count( "simx" ) ) cmdTree.put( "simx", vm["simx"].as<string>() );
    if( vm.count( "simy" ) ) cmdTree.put( "simy", vm["simy"].as<string>() );
    if( vm.count( "imgn" ) ) cmdTree.put( "imgn", vm["imgn"].as<string>() );

    if( vm.count( "output-file" ) ) cmdTree.put( "output-file", vm["output-file"].as<string>() );
    if( vm.count( "sequence" ) ) sequenceNumber = vm["sequence"].as<int>();

    string tmpString = tree.get<string>( "BASIS", "Zernike" );
    basis = iequals( tmpString, "Karhunen-Loeve" ) ? CFG_KARHUNEN_LOEVE : CFG_ZERNIKE;

    modes = tree.get<vector<uint32_t>>( "MODES", vector<uint32_t>( {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                                    20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35
                                                                   } ) );
    imageNumbers = tree.get<vector<uint32_t>>( "imgn", vector<uint32_t>() );
    if( imageNumbers.size() == 0 ) imageNumbers = tree.get<vector<uint32_t>>( "IMAGE_NUM", vector<uint32_t>() );
    darkNumbers = tree.get<vector<uint32_t>>( "DARK_NUM", vector<uint32_t>() );
    sequenceNumber = tree.get<uint32_t>( "SEQUENCE_NUM", sequenceNumber );

    klMinMode  = tree.get<uint32_t>( "KL_MIN_MODE", DEF_KL_MIN_MODE );
    klMaxMode  = tree.get<uint32_t>( "KL_MAX_MODE", DEF_KL_MAX_MODE );
    nPoints    = tree.get<uint32_t>( "NUM_POINTS", DEF_NUM_POINTS );
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
                "\"median\", \"invdistweight\" or \"horizontal interpolation\"";
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
        ifstream file( tmpString.c_str(), ifstream::binary );
        std::shared_ptr<Ana> header( new Ana() );
        header->read( file );

        if( header->m_Header.ndim != 2 ) {
            LOG_ERR << "pupil file \"" << tmpString << "\" not 2 dimensional.";
        }
        pupilSize = header->m_Header.dim[0];
        uint32_t tmpInt = header->m_Header.dim[1];
        if( pupilSize != tmpInt ) {
            LOG_ERR << "pupil file \"" << tmpString << "\" not square (" << pupilSize << "x" << tmpInt << ")";
        }
        pupil = Array<double>( pupilSize, tmpInt );
        Ana::read( file, pupil, header );

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

    xl = tree.get<vector<uint32_t>>( "simx", vector<uint32_t>() );
    if( xl.size() == 0 ) xl = tree.get<vector<uint32_t>>( "SIM_X", vector<uint32_t>() );
    yl = tree.get<vector<uint32_t>>( "simy", vector<uint32_t>() );
    if( yl.size() == 0 ) yl = tree.get<vector<uint32_t>>( "SIM_Y", vector<uint32_t>() );
    xh.resize( xl.size() );
    yh.resize( yl.size() );
    for( size_t i = 0; i < xl.size(); ++i ) {
        xl[i] -= nPoints / 2;
        xh[i] = xl[i] + nPoints - 1;
    }
    for( size_t i = 0; i < yl.size(); ++i ) {
        yl[i] -= nPoints / 2;
        yh[i] = yl[i] + nPoints - 1;
    }

    if( tree.get<bool>( "CAL_X", false ) ) {
        if( tree.get<bool>( "CAL_Y", false ) ) {
            vector<uint32_t> tmpx = tree.get<vector<uint32_t>>( "CAL_X", vector<uint32_t>() );
            vector<uint32_t> tmpy = tree.get<vector<uint32_t>>( "CAL_Y", vector<uint32_t>() );
            if( tmpx.size() && ( tmpx.size() == tmpy.size() ) ) {
                ncal = tmpx.size();
                if( xl.size() ) LOG << "Note: SIM_X/SIM_Y replaced by CAL_X/CAL_Y";
                for( size_t i = 0; i < tmpx.size(); ++i ) {
                    xl[i] = tmpx[i] - nPoints / 2;
                    xh[i] = xl[i] + nPoints - 1;
                    yl[i] = tmpy[i] - nPoints / 2;
                    yh[i] = yl[i] + nPoints - 1;
                }
            }
            else LOG_ERR << "CAL_X and CAL_Y must have the same number of elements!";
        }
        else LOG_ERR << "CAL_Y must be provided if CAL_X is!";
    }

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
    boost::split(outputFiles, tmpString, boost::is_any_of(",") );

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
    if( nPoints != DEF_NUM_POINTS ) tree.put( "NUM_POINTS", nPoints );
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

size_t MomfbdJob::size(void) const {
    size_t sz = Job::size();
    
    return sz;
}

char* MomfbdJob::pack(char* ptr) const {
    ptr = Job::pack(ptr);
    return ptr;
}

const char* MomfbdJob::unpack(const char* ptr, bool swap_endian) {
    ptr = Job::unpack(ptr, swap_endian);
    return ptr;
}
