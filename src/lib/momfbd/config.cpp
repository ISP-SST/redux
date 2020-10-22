#include "redux/momfbd/config.hpp"

#include "redux/file/fileio.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/constants.hpp"
#include "redux/image/zernike.hpp"
#include "redux/logging/logger.hpp"
#include "redux/translators.hpp"

#include <boost/algorithm/string.hpp>
//#include <boost/range/algorithm.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
namespace bpo = boost::program_options;
namespace bpt = boost::property_tree;
namespace bfs = boost::filesystem;

using namespace redux::file;
using namespace redux::momfbd;
using namespace redux::image;
using namespace redux::logging;
using namespace redux::util;
using namespace std;
using boost::algorithm::iequals;


namespace {

    const char* basisTags[] = { "", "Zernike", "Karhunen-Loeve" };
    const char* fpmTags[] = { "", "median", "invdistweight", "horint" };
    const char* gmTags[] = { "", "gradient_diff", "gradient_Vogel" };
    const char* gsmTags[] = { "", "getstep_steepest_descent", "getstep_conjugate_gradient", "getstep_BFGS", "getstep_BFGS_inv" };
    const char* ftTags[] = { "", "ANA", "FITS", "ANA,FITS", "MOMFBD", "ANA,MOMFBD", "FITS,MOMFBD", "ANA,FITS,MOMFBD" };
    const char* ftExt[] = { "", "f0", "fits", "", "momfbd" };
    const char* dtTags[] = { "byte", "short", "int", "int64", "float", "double" };

    template <typename T>
    int getFromMap( const string& str, const map<string, int, T>& m ) {
        auto res = m.find( str );

        if( res != m.end() ) {
            return res->second;
        }

        return 0;
    }
    
    template<typename T>
    void setFlag( T& flagset, T flag, bool value ) {
        if( value != bool( flagset&flag ) ) {
            flagset ^= flag;
        }
    }

    
    template<typename T>
    T getValue( const bpt::ptree& tree, string name, const T& defaultValue ) {
        T ret;
        try {
            ret = tree.get<T>( name, defaultValue );
        } catch( exception& e ) {
            string msg = "Failed to convert entry \"" + name + "\" to number(s).\n"
            + "Value: \"" +  tree.get<string>( name )
            + "\"\nReason: " + e.what();
            throw logic_error( msg );
        }
        return ret;
    }

}

const map<FileType, string> redux::momfbd::FileTypeNames = {
    { FT_ANA, ftTags[FT_ANA] },
    { FT_FITS, ftTags[FT_FITS] },
    { FT_MOMFBD, ftTags[FT_MOMFBD] }
};

const map<FileType, string> redux::momfbd::FileTypeExtensions = {
    { FT_ANA, ftExt[FT_ANA] },
    { FT_FITS, ftExt[FT_FITS] },
    { FT_MOMFBD, ftExt[FT_MOMFBD] }
};

const map<string, int, cicomp> redux::momfbd::fillpixMap = {
    { fpmTags[FPM_MEDIAN], FPM_MEDIAN },
    { fpmTags[FPM_INVDISTWEIGHT], FPM_INVDISTWEIGHT },
    { fpmTags[FPM_HORINT], FPM_HORINT }
};

const map<string, int, cicomp> redux::momfbd::gradientMap = {
    { gmTags[GM_DIFF], GM_DIFF },
    { gmTags[GM_VOGEL], GM_VOGEL },
    { "diff", GM_DIFF },            // alternative spellings
    { "vogel", GM_VOGEL}
};

const map<string, int, cicomp> redux::momfbd::getstepMap = {
    { gsmTags[GSM_SDSC], GSM_SDSC },
    { gsmTags[GSM_CNJG], GSM_CNJG },
    { gsmTags[GSM_BFGS], GSM_BFGS },
    { gsmTags[GSM_BFGS_inv], GSM_BFGS_inv },
    { "steepest", GSM_SDSC },       // alternative spellings
    { "steepest_descent", GSM_SDSC },
    { "conjugate", GSM_CNJG },
    { "conjugate_gradient", GSM_CNJG },
    { "bfgs", GSM_BFGS },
    { "bfgsinv", GSM_BFGS_inv },
    { "bfgs_inv", GSM_BFGS_inv }
};







/********************  Channel  ********************/

ChannelCfg::ChannelCfg() : rotationAngle(0), noiseFudge(1), weight(1), diversityBasis(ZERNIKE), noRestore(false),
        borderClip(100), incomplete(0), discard(2,0),
        mmRow(0), mmWidth(0), imageNumberOffset(0) {

}


ChannelCfg::~ChannelCfg() {

}


ChannelCfg::operator std::string() const {
    bpt::ptree dump;
    std::stringstream ss;
    this->getProperties( dump );
    bpt::write_info( ss, dump );
    return ss.str();
}


void ChannelCfg::parseProperties( bpt::ptree& tree, Logger& logger, const ChannelCfg& defaults ) {
    
    rotationAngle = getValue( tree, "ANGLE", defaults.rotationAngle );

    noiseFudge = tree.get<double>("NF", defaults.noiseFudge);
    weight = tree.get<double>("WEIGHT", defaults.weight);
    noRestore = tree.get<bool>("NO_RESTORE", defaults.noRestore);
    
    diversityBasis = tree.get<ModeBase>( "DIV_BASIS", defaults.diversityBasis );
    diversityValues = tree.get<vector<DiversityValue>>("DIVERSITY", defaults.diversityValues);
    diversityModes = tree.get<ModeList>("DIV_ORDERS", defaults.diversityModes);
    diversityModes.setDefaultModeType( diversityBasis );
    
    // if there is only 1 coefficient, we assume it is a PD defocus term, and don't require DIV_ORDERS
    if( diversityValues.size() == 1 && diversityModes.empty() ) { 
        diversityModes.push_back( ModeID(4,ZERNIKE) );
    }
    for( auto& d: diversityModes ) {
        if( d.mode == 2 || d.mode == 3 ) d.type = ZERNIKE;      // Force Zernike for all tilt modes.
    }


    alignMap = getValue( tree, "ALIGN_MAP", defaults.alignMap );
    alignClip = getValue( tree, "ALIGN_CLIP", defaults.alignClip );
    discard = getValue( tree, "DISCARD", defaults.discard );
    borderClip = getValue( tree, "BORDER_CLIP", defaults.borderClip );
    incomplete = getValue<bool>( tree, "INCOMPLETE", defaults.incomplete );
    
    subImagePosXY = getValue( tree, "SIM_XY", defaults.subImagePosXY );
    subImagePosX = getValue( tree, "SIM_X", defaults.subImagePosX );
    subImagePosY = getValue( tree, "SIM_Y", defaults.subImagePosY );
    if( getValue<bool>( tree, "CAL_X", false ) ) {
        if( getValue<bool>( tree, "CAL_Y", false ) ) {
            if( subImagePosX.size() || subImagePosY.size() ) LOG_WARN << "Note: SIM_X/SIM_Y replaced by CAL_X/CAL_Y" << ende;
            subImagePosX = getValue( tree, "CAL_X", defaults.subImagePosX );
            subImagePosY = getValue( tree, "CAL_Y", defaults.subImagePosY );
            if( subImagePosX.empty() || ( subImagePosX.size() != subImagePosY.size() ) ) {
                LOG_ERR << "CAL_X and CAL_Y must have the same number of elements!" << ende;
            }
        } else LOG_ERR << "CAL_Y must be provided if CAL_X is!" << ende;
    }

    
    imageDataDir = getValue<string>( tree, "IMAGE_DATA_DIR", defaults.imageDataDir );
    imageTemplate = getValue<string>( tree, "FILENAME_TEMPLATE", defaults.imageTemplate );
    darkTemplate = getValue<string>( tree, "DARK_TEMPLATE", defaults.darkTemplate );
    gainFile = getValue<string>( tree, "GAIN_FILE", defaults.gainFile );
    responseFile = getValue<string>( tree, "CCD_RESPONSE", defaults.responseFile );
    backgainFile = getValue<string>( tree, "BACK_GAIN", defaults.backgainFile );
    psfFile = getValue<string>( tree, "PSF", defaults.psfFile );
    mmFile = getValue<string>( tree, "MODMAT", defaults.mmFile );
    mmRow = getValue( tree, "MMROW", defaults.mmRow );
    mmWidth = getValue( tree, "MMWIDTH",defaults.mmWidth );
    xOffsetFile = getValue<string>( tree, "XOFFSET", "" );
    yOffsetFile = getValue<string>( tree, "YOFFSET", "" );
    imageNumberOffset = getValue( tree, "DT", defaults.imageNumberOffset );
    fileNumbers = getValue( tree, "IMAGE_NUM", defaults.fileNumbers );
    if( fileNumbers.empty() ) {
        fileNumbers = getValue( tree, "IMAGE_NUMS", defaults.fileNumbers );
    }
    waveFrontList = getValue( tree, "WFINDEX", defaults.waveFrontList );
    darkNumbers = getValue( tree, "DARK_NUM", defaults.darkNumbers );
    stokesWeights = getValue( tree, "VECTOR", defaults.stokesWeights );
    if( mmFile.length() > 0 ) {
        if( !mmRow ) {
            LOG_ERR << "a modulation matrix was provided but no row specified (MMROW)." << ende;
        }
        if( !mmWidth ) {
            LOG_ERR << "modulation matrix dimensions cannot be autodetected (yet): you must provide the matrix width (MMWIDTH)!" << ende;
        }
        if( stokesWeights.size() == 0 ) {
            LOG_ERR << "modulation matrix specified but no VECTOR input given!" << ende;
        } else if( stokesWeights.size() != mmWidth ) {
            LOG_ERR << "VECTOR input has " << stokesWeights.size() << " elements, but MMWIDTH=" << (int)mmWidth << ende;
        }
    }
    else {  // TODO: don't modify cfg values!! ...make the main code use weight 1 as default instead.
        mmRow = mmWidth = 1;
        stokesWeights.resize( 1, 1.0 );
    }

    // Resolve relative paths etc.
    if( imageDataDir.empty() ) imageDataDir = "./";
    imageDataDir = cleanPath( imageDataDir );
    if( !gainFile.empty() ) gainFile = cleanPath( gainFile );
    if( !responseFile.empty() ) responseFile = cleanPath( responseFile );
    if( !backgainFile.empty() ) backgainFile = cleanPath( backgainFile );
    if( !psfFile.empty() ) psfFile = cleanPath( psfFile );
    if( !mmFile.empty() ) mmFile = cleanPath( mmFile );
    if( !xOffsetFile.empty() ) xOffsetFile = cleanPath( xOffsetFile );
    if( !yOffsetFile.empty() ) yOffsetFile = cleanPath( yOffsetFile );

}


void ChannelCfg::getProperties( bpt::ptree& tree, const ChannelCfg& defaults, bool showAll ) const {

    if( showAll || rotationAngle != defaults.rotationAngle ) tree.put( "ANGLE", rotationAngle );

    if( showAll || noiseFudge != defaults.noiseFudge ) tree.put("NF", noiseFudge);
    if( showAll || weight != defaults.weight ) tree.put("WEIGHT", weight);
    if( showAll || noRestore != defaults.noRestore ) tree.put("NO_RESTORE", noRestore);
    
    if( showAll || diversityBasis != defaults.diversityBasis ) tree.put( "DIV_BASIS", diversityBasis );
    if( showAll || diversityModes != defaults.diversityModes ) tree.put( "DIV_ORDERS", diversityModes );
    if( showAll || diversityValues != defaults.diversityValues ) tree.put( "DIVERSITY", diversityValues );
    
    if( showAll || alignMap != defaults.alignMap ) tree.put( "ALIGN_MAP", alignMap );
    if( showAll || alignClip != defaults.alignClip ) tree.put( "ALIGN_CLIP", alignClip );
    if( showAll || discard != defaults.discard ) tree.put( "DISCARD", discard );
    if( showAll || borderClip != defaults.borderClip ) tree.put( "BORDER_CLIP", borderClip );
    if( showAll || incomplete != defaults.incomplete ) tree.put( "INCOMPLETE", ( bool)incomplete );

    if( showAll || subImagePosXY != defaults.subImagePosXY ) tree.put( "SIM_XY", subImagePosXY );
    if( showAll || subImagePosXY.empty() ) {
        if( subImagePosX != defaults.subImagePosX ) tree.put( "SIM_X", subImagePosX );
        if( subImagePosY != defaults.subImagePosY ) tree.put( "SIM_Y", subImagePosY );
    }
    
    if( showAll || imageDataDir != defaults.imageDataDir ) tree.put( "IMAGE_DATA_DIR", imageDataDir );
    if( showAll || imageTemplate != defaults.imageTemplate ) tree.put( "FILENAME_TEMPLATE", imageTemplate );
    if( showAll || darkTemplate != defaults.darkTemplate ) tree.put( "DARK_TEMPLATE", darkTemplate );
    if( showAll || gainFile != defaults.gainFile ) tree.put( "GAIN_FILE", gainFile );
    if( showAll || responseFile != defaults.responseFile ) tree.put( "CCD_RESPONSE", responseFile );
    if( showAll || backgainFile != defaults.backgainFile ) tree.put( "BACK_GAIN", backgainFile );
    if( showAll || psfFile != defaults.psfFile ) tree.put( "PSF", psfFile );
    if( showAll || mmFile != defaults.mmFile ) tree.put( "MODMAT", mmFile );
    if( showAll || mmRow != defaults.mmRow ) tree.put( "MMROW", mmRow );
    if( showAll || mmWidth != defaults.mmWidth ) tree.put( "MMWIDTH", mmWidth );
    
    if( showAll || xOffsetFile != defaults.xOffsetFile ) tree.put( "XOFFSET", xOffsetFile );
    if( showAll || yOffsetFile != defaults.yOffsetFile ) tree.put( "YOFFSET", yOffsetFile );

    if( showAll || imageNumberOffset != defaults.imageNumberOffset ) tree.put( "DT", imageNumberOffset );
    if( showAll || fileNumbers != defaults.fileNumbers ) tree.put( "IMAGE_NUM", fileNumbers );
    if( showAll || waveFrontList != defaults.waveFrontList ) tree.put( "WFINDEX", waveFrontList );
    if( showAll || darkNumbers != defaults.darkNumbers ) tree.put( "DARK_NUM", darkNumbers );
    
    if( showAll || stokesWeights != defaults.stokesWeights ) tree.put( "VECTOR", stokesWeights );

}


uint64_t ChannelCfg::size( void ) const {
    // static sizes (PoD types)
    static uint64_t ssz = sizeof( borderClip ) + sizeof( diversityBasis ) + sizeof( imageNumberOffset )
        + sizeof( incomplete ) + sizeof( mmRow ) + sizeof( mmWidth ) + sizeof( noiseFudge ) + sizeof( noRestore )
        + sizeof( rotationAngle ) + sizeof( weight );
    uint64_t sz = ssz;
    // strings
    sz += backgainFile.length() + darkTemplate.length() + gainFile.length() + 3;
    sz += imageDataDir.length() + imageTemplate.length() + mmFile.length() + psfFile.length() + 4;
    sz += responseFile.length() + xOffsetFile.length() + yOffsetFile.length() + 3;
    // vectors etc.
    sz += alignClip.size()*sizeof( int16_t ) + sizeof( uint64_t );
    sz += alignMap.size()*sizeof( float ) + sizeof( uint64_t );
    sz += darkNumbers.size()*sizeof( uint32_t ) + sizeof( uint64_t );
    sz += discard.size()*sizeof( uint16_t ) + sizeof( uint64_t );
    sz += redux::util::size( diversityModes );
    sz += redux::util::size( diversityValues );
    sz += fileNumbers.size()*sizeof( uint32_t ) + sizeof( uint64_t );
    sz += stokesWeights.size()*sizeof( float ) + sizeof( uint64_t );
    sz += subImagePosX.size()*sizeof( uint16_t ) + sizeof( uint64_t );
    sz += subImagePosXY.size()*sizeof( uint16_t ) + sizeof( uint64_t );
    sz += subImagePosY.size()*sizeof( uint16_t ) + sizeof( uint64_t );
    sz += waveFrontList.size()*sizeof( uint32_t ) + sizeof( uint64_t );
    return sz;
}


uint64_t ChannelCfg::pack( char* ptr ) const {
    using redux::util::pack;
    // scalar values
    uint64_t count = pack( ptr, borderClip );
    count += pack( ptr+count, diversityBasis );
    count += pack( ptr+count, imageNumberOffset );
    count += pack( ptr+count, incomplete );
    count += pack( ptr+count, mmRow );
    count += pack( ptr+count, mmWidth );
    count += pack( ptr+count, noiseFudge );
    count += pack( ptr+count, noRestore );    //b
    count += pack( ptr+count, rotationAngle );
    count += pack( ptr+count, weight );
    // strings
    count += pack( ptr+count, backgainFile );
    count += pack( ptr+count, darkTemplate );
    count += pack( ptr+count, gainFile );
    count += pack( ptr+count, imageDataDir );
    count += pack( ptr+count, imageTemplate );
    count += pack( ptr+count, mmFile );
    count += pack( ptr+count, psfFile );
    count += pack( ptr+count, responseFile );
    count += pack( ptr+count, xOffsetFile );
    count += pack( ptr+count, yOffsetFile );
    // vectors etc.
    count += pack( ptr+count, alignClip );
    count += pack( ptr+count, alignMap );
    count += pack( ptr+count, darkNumbers );
    count += pack( ptr+count, discard );
    count += pack( ptr+count, diversityModes );
    count += pack( ptr+count, diversityValues );
    count += pack( ptr+count, fileNumbers );
    count += pack( ptr+count, stokesWeights );
    count += pack( ptr+count, subImagePosX );
    count += pack( ptr+count, subImagePosXY );
    count += pack( ptr+count, subImagePosY );
    count += pack( ptr+count, waveFrontList );
    return count;
}


uint64_t ChannelCfg::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    // scalar values
    uint64_t count = unpack( ptr, borderClip, swap_endian );
    count += unpack( ptr+count, diversityBasis, swap_endian );
    count += unpack( ptr+count, imageNumberOffset, swap_endian );
    count += unpack( ptr+count, incomplete );
    count += unpack( ptr+count, mmRow );
    count += unpack( ptr+count, mmWidth );
    count += unpack( ptr+count, noiseFudge, swap_endian );
    count += unpack( ptr+count, noRestore );
    count += unpack( ptr+count, rotationAngle, swap_endian );
    count += unpack( ptr+count, weight, swap_endian );
    // strings
    count += unpack( ptr+count, backgainFile );
    count += unpack( ptr+count, darkTemplate );
    count += unpack( ptr+count, gainFile );
    count += unpack( ptr+count, imageDataDir );
    count += unpack( ptr+count, imageTemplate );
    count += unpack( ptr+count, mmFile );
    count += unpack( ptr+count, psfFile );
    count += unpack( ptr+count, responseFile );
    count += unpack( ptr+count, xOffsetFile );
    count += unpack( ptr+count, yOffsetFile );
    // vectoruns etc.
    count += unpack( ptr+count, alignClip, swap_endian );
    count += unpack( ptr+count, alignMap, swap_endian );
    count += unpack( ptr+count, darkNumbers, swap_endian );
    count += unpack( ptr+count, discard, swap_endian );
    count += unpack( ptr+count, diversityModes, swap_endian );
    count += unpack( ptr+count, diversityValues, swap_endian );
    count += unpack( ptr+count, fileNumbers, swap_endian );
    count += unpack( ptr+count, stokesWeights, swap_endian );
    count += unpack( ptr+count, subImagePosX, swap_endian );
    count += unpack( ptr+count, subImagePosXY, swap_endian );
    count += unpack( ptr+count, subImagePosY, swap_endian );
    count += unpack( ptr+count, waveFrontList, swap_endian );
    return count;
}


bool ChannelCfg::operator==( const ChannelCfg& rhs ) const {
    return ( rotationAngle == rhs.rotationAngle ) &&
           ( noiseFudge == rhs.noiseFudge ) &&
           ( weight == rhs.weight ) &&
           ( borderClip == rhs.borderClip ) &&
           ( incomplete == rhs.incomplete ) &&
           ( subImagePosX == rhs.subImagePosX ) &&
           ( subImagePosY == rhs.subImagePosY ) &&
           ( imageDataDir == rhs.imageDataDir ) &&
           ( darkTemplate == rhs.darkTemplate ) &&
           ( imageTemplate == rhs.imageTemplate ) &&
           ( fileNumbers == rhs.fileNumbers ) &&
           ( waveFrontList == rhs.waveFrontList ) &&
           ( darkNumbers == rhs.darkNumbers );
}



/********************   Object  ********************/

ObjectCfg::ObjectCfg() : telescopeF(0), arcSecsPerPixel(0), pixelSize(1E-5),
                         alphaToPixels(0), pixelsToAlpha(0),
                         alphaToDefocus(0), defocusToAlpha(0),
                         maxLocalShift(5), minimumOverlap(16), 
                         patchSize(128), pupilPixels(64), saveMask(0), wavelength(0), traceObject(false) {

}


ObjectCfg::~ObjectCfg() {

}


void ObjectCfg::parseProperties( bpt::ptree& tree, Logger& logger, const ChannelCfg& def ) {
    
    const ObjectCfg& defaults = reinterpret_cast<const ObjectCfg&>( def );

    telescopeF = getValue( tree, "TELESCOPE_F", defaults.telescopeF );
    arcSecsPerPixel = getValue( tree, "ARCSECPERPIX", defaults.arcSecsPerPixel );
    pixelSize = getValue( tree, "PIXELSIZE", defaults.pixelSize );

    maxLocalShift = getValue( tree, "MAX_LOCAL_SHIFT", defaults.maxLocalShift );
    minimumOverlap = getValue( tree, "MINIMUM_OVERLAP", defaults.minimumOverlap );
    patchSize  = getValue( tree, "NUM_POINTS", defaults.patchSize );
    pupilPixels  = getValue( tree, "PUPIL_POINTS", defaults.pupilPixels );
    
    saveMask = 0;
    if( getValue<bool>( tree, "GET_ALPHA", defaults.saveMask&SF_SAVE_ALPHA ) ) saveMask |= SF_SAVE_ALPHA;
    if( getValue<bool>( tree, "GET_COBJ", defaults.saveMask&SF_SAVE_COBJ ) ) saveMask |= SF_SAVE_COBJ;
    if( getValue<bool>( tree, "GET_DIVERSITY", defaults.saveMask&SF_SAVE_DIVERSITY ) ) saveMask |= SF_SAVE_DIVERSITY;
    if( getValue<bool>( tree, "GET_METRIC", defaults.saveMask&SF_SAVE_METRIC ) ) saveMask |= SF_SAVE_METRIC;
    if( getValue<bool>( tree, "GET_MODES", defaults.saveMask&SF_SAVE_MODES ) ) saveMask |= SF_SAVE_MODES;
    if( getValue<bool>( tree, "GET_PSF", defaults.saveMask&SF_SAVE_PSF ) ) saveMask |= SF_SAVE_PSF;
    if( getValue<bool>( tree, "GET_PSF_AVG", defaults.saveMask&SF_SAVE_PSF_AVG ) ) saveMask |= SF_SAVE_PSF_AVG;
    if( getValue<bool>( tree, "GET_RESIDUAL", defaults.saveMask&SF_SAVE_RESIDUAL ) ) saveMask |= SF_SAVE_RESIDUAL;
    if( getValue<bool>( tree, "SAVE_FFDATA", defaults.saveMask&SF_SAVE_FFDATA ) ) saveMask |= SF_SAVE_FFDATA;
    traceObject = tree.get<bool>( "TRACE_REF", defaults.traceObject );
    outputFileName = getValue<string>( tree, "OUTPUT_FILE", defaults.outputFileName );
    initFile = getValue<string>( tree, "INIT_FILE", defaults.initFile );
    modeFile = getValue<string>( tree, "MODE_FILE", defaults.modeFile );
    pupilFile = getValue<string>( tree, "PUPIL", defaults.pupilFile );
    wavelength = getValue( tree, "WAVELENGTH", defaults.wavelength );

    if( ( saveMask & SF_SAVE_PSF ) && ( saveMask & SF_SAVE_PSF_AVG ) ) {
        saveMask &= ~SF_SAVE_PSF;
        //LOG_WARN << "both GET_PSF and GET_PSF_AVG mode requested" << ende;
    }

    // Resolve relative paths etc.
    if( !initFile.empty() ) initFile = cleanPath( initFile );
    if( !modeFile.empty() ) modeFile = cleanPath( modeFile );
    if( !pupilFile.empty() ) pupilFile = cleanPath( pupilFile );
    
    ChannelCfg::parseProperties( tree, logger, defaults );

}


void ObjectCfg::getProperties( bpt::ptree& tree, const ChannelCfg& def, bool showAll ) const {

    const ObjectCfg& defaults = reinterpret_cast<const ObjectCfg&>( def );
    
    if( showAll || telescopeF != defaults.telescopeF ) tree.put( "TELESCOPE_F", telescopeF );
    if( showAll || arcSecsPerPixel != defaults.arcSecsPerPixel ) tree.put( "ARCSECPERPIX", arcSecsPerPixel );
    if( showAll || pixelSize != defaults.pixelSize ) tree.put( "PIXELSIZE", pixelSize );
    if( showAll || maxLocalShift != defaults.maxLocalShift ) tree.put( "MAX_LOCAL_SHIFT", maxLocalShift );
    if( showAll || minimumOverlap != defaults.minimumOverlap ) tree.put( "MINIMUM_OVERLAP", minimumOverlap );
    if( showAll || patchSize != defaults.patchSize ) tree.put( "NUM_POINTS", patchSize );
    if( showAll || pupilPixels != defaults.pupilPixels ) tree.put( "PUPIL_POINTS", pupilPixels );

    uint16_t diff = saveMask ^ defaults.saveMask;
    if( diff & SF_SAVE_ALPHA ) tree.put( "GET_ALPHA", bool( saveMask & SF_SAVE_ALPHA ) );
    if( diff & SF_SAVE_COBJ ) tree.put( "GET_COBJ", bool( saveMask & SF_SAVE_COBJ ) );
    if( diff & SF_SAVE_DIVERSITY ) tree.put( "GET_DIVERSITY", bool( saveMask & SF_SAVE_DIVERSITY ) );
    if( diff & SF_SAVE_METRIC ) tree.put( "GET_METRIC", bool( saveMask & SF_SAVE_METRIC ) );
    if( diff & SF_SAVE_MODES ) tree.put( "GET_MODES", bool( saveMask & SF_SAVE_MODES ) );
    if( diff & SF_SAVE_PSF ) tree.put( "GET_PSF", bool( saveMask & SF_SAVE_PSF ) );
    if( diff & SF_SAVE_PSF_AVG ) tree.put( "GET_PSF_AVG", bool( saveMask & SF_SAVE_PSF_AVG ) );
    if( diff & SF_SAVE_RESIDUAL ) tree.put( "GET_RESIDUAL", bool( saveMask & SF_SAVE_RESIDUAL ) );
    if( diff & SF_SAVE_FFDATA ) tree.put( "SAVE_FFDATA", bool( saveMask & SF_SAVE_FFDATA ) );
    if( showAll || traceObject != defaults.traceObject ) tree.put( "TRACE_REF", traceObject );
    if( showAll || outputFileName != defaults.outputFileName ) tree.put( "OUTPUT_FILE", outputFileName );
    if( showAll || initFile != defaults.initFile ) tree.put( "INIT_FILE", initFile );
    if( showAll || modeFile != defaults.modeFile ) tree.put( "MODE_FILE", modeFile );
    if( showAll || pupilFile != defaults.pupilFile ) tree.put( "PUPIL", pupilFile );
    if( showAll || wavelength != defaults.wavelength ) tree.put( "WAVELENGTH", wavelength );

    ChannelCfg::getProperties( tree, defaults );

}


uint64_t ObjectCfg::size( void ) const {
    // static sizes (PoD types)
    uint64_t ssz = sizeof( alphaToPixels ) + sizeof( alphaToDefocus ) + sizeof( arcSecsPerPixel ) + sizeof( defocusToAlpha )
                 + sizeof( maxLocalShift ) + sizeof( minimumOverlap ) + sizeof( patchSize ) + sizeof( pixelSize )
                 + sizeof( pixelsToAlpha ) + sizeof( pupilPixels ) + sizeof( saveMask ) + sizeof( telescopeF )
                 + sizeof( traceObject ) + sizeof( wavelength );
    uint64_t sz = ssz + ChannelCfg::size();
    sz += initFile.length() + modeFile.length() + 2;
    sz += outputFileName.length() + pupilFile.length() + 2;
    return sz;
}


uint64_t ObjectCfg::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = ChannelCfg::pack( ptr );
    // scalar values
    count += pack( ptr+count, alphaToPixels );
    count += pack( ptr+count, alphaToDefocus );
    count += pack( ptr+count, arcSecsPerPixel );
    count += pack( ptr+count, defocusToAlpha );
    count += pack( ptr+count, maxLocalShift );
    count += pack( ptr+count, minimumOverlap );
    count += pack( ptr+count, patchSize );
    count += pack( ptr+count, pixelSize );
    count += pack( ptr+count, pixelsToAlpha );
    count += pack( ptr+count, pupilPixels );
    count += pack( ptr+count, saveMask );
    count += pack( ptr+count, telescopeF );
    count += pack( ptr+count, traceObject );
    count += pack( ptr+count, wavelength );
    // strings
    count += pack( ptr+count, initFile );
    count += pack( ptr+count, modeFile );
    count += pack( ptr+count, outputFileName );
    count += pack( ptr+count, pupilFile );
    return count;
}


uint64_t ObjectCfg::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = ChannelCfg::unpack( ptr, swap_endian );
    // scalar values
    count += unpack( ptr+count, alphaToPixels, swap_endian );
    count += unpack( ptr+count, alphaToDefocus, swap_endian );
    count += unpack( ptr+count, arcSecsPerPixel, swap_endian );
    count += unpack( ptr+count, defocusToAlpha, swap_endian );
    count += unpack( ptr+count, maxLocalShift, swap_endian );
    count += unpack( ptr+count, minimumOverlap, swap_endian );
    count += unpack( ptr+count, patchSize, swap_endian );
    count += unpack( ptr+count, pixelSize, swap_endian );
    count += unpack( ptr+count, pixelsToAlpha, swap_endian );
    count += unpack( ptr+count, pupilPixels, swap_endian );
    count += unpack( ptr+count, saveMask, swap_endian );
    count += unpack( ptr+count, telescopeF, swap_endian );
    count += unpack( ptr+count, traceObject );
    count += unpack( ptr+count, wavelength, swap_endian );
    // strings
    count += unpack( ptr+count, initFile );
    count += unpack( ptr+count, modeFile );
    count += unpack( ptr+count, outputFileName );
    count += unpack( ptr+count, pupilFile );
    return count;
}


const ObjectCfg& ObjectCfg::operator=( const ChannelCfg& rhs ) {
    ChannelCfg::operator=( rhs );
    return *this;
}


bool ObjectCfg::operator==( const ObjectCfg& rhs ) const {
    return ( saveMask == rhs.saveMask ) &&
           ( patchSize == rhs.patchSize ) &&
           ( pupilPixels == rhs.pupilPixels ) &&
           ( maxLocalShift == rhs.maxLocalShift ) &&
           ( wavelength == rhs.wavelength ) &&
           ( outputFileName == rhs.outputFileName ) &&
           ChannelCfg::operator==( rhs );
}



/********************   Global   ********************/

GlobalCfg::GlobalCfg() : runFlags( 0), modeBasis(ZERNIKE), klMinMode( 2), klMaxMode( 2000), klCutoff( 1E-3),
    nInitialModes( 5), nModeIncrement(5), nModes(0),
    telescopeD(0), telescopeCO(0), minIterations(5), maxIterations(500), targetIterations(3),
    fillpixMethod(FPM_INVDISTWEIGHT), gradientMethod(GM_DIFF), getstepMethod(GSM_BFGS_inv),
    badPixelThreshold(1E-5), FTOL(1E-3), EPS(1E-10), reg_alpha(0), graddiff_step(1E-2), trace(false),
    outputFileType(FT_NONE), outputDataType(DT_I16T), sequenceNumber(0),
    observationTime(""), observationDate("N/A"), tmpDataDir("./data") {


}


GlobalCfg::~GlobalCfg() {

}


void GlobalCfg::parseProperties( bpt::ptree& tree, Logger& logger, const ChannelCfg& def ) {
    
    const GlobalCfg& defaults = reinterpret_cast<const GlobalCfg&>(def);

    if( getValue<bool>( tree, "CALIBRATE", false ) )            runFlags |= RF_CALIBRATE;
    if( getValue<bool>( tree, "DONT_MATCH_IMAGE_NUMS", false )) runFlags |= RF_DONT_MATCH_IMAGE_NUMS;
    if( getValue<bool>( tree, "FAST_QR", false ) )              runFlags |= RF_FAST_QR;
    if( getValue<bool>( tree, "FIT_PLANE", false ) )            runFlags |= RF_FIT_PLANE;
    if( getValue<bool>( tree, "FLATFIELD", false ) )            runFlags |= RF_FLATFIELD;
    if( getValue<bool>( tree, "GLOBAL_NOISE", false ) )         runFlags |= RF_GLOBAL_NOISE;
    if( getValue<bool>( tree, "NEW_CONSTRAINTS", false ) )      runFlags |= RF_NEW_CONSTRAINTS;
    if( getValue<bool>( tree, "NOSWAP", false ) )               runFlags |= RF_NOSWAP;
    if( getValue<bool>( tree, "NO_CLIP", false ) )              runFlags |= RF_NO_CLIP;
    if( getValue<bool>( tree, "NO_CONSTRAINTS", false ) )       runFlags |= RF_NO_CONSTRAINTS;
    if( getValue<bool>( tree, "NO_FILTER", false ) )            runFlags |= RF_NO_FILTER;
    if( getValue<bool>( tree, "OLD_NS", false ) )               runFlags |= RF_OLD_NS;
    if( getValue<bool>( tree, "OVERWRITE", false ) )            runFlags |= RF_FORCE_WRITE;
    if( getValue<bool>( tree, "SORT_MODES", false ) )           runFlags |= RF_SORT_MODES;
    
    trace = tree.get<bool>( "TRACE", defaults.trace );
    
/*    if( ( runFlags & RF_CALIBRATE ) && ( runFlags & RF_FLATFIELD ) ) {
        LOG_WARN << "both FLATFIELD and CALIBRATE mode requested, forcing CALIBRATE";
        runFlags &= ~RF_FLATFIELD;
    }
    if( ( runFlags & RF_CALIBRATE ) && ( runFlags & RF_NEW_CONSTRAINTS ) ) {
        LOG_WARN << "calibration mode uses old style constraints, ignoring NEW_CONSTRAINTS";
        runFlags &= ~RF_NEW_CONSTRAINTS;
    }
*/
    modeBasis = tree.get<ModeBase>( "BASIS", defaults.modeBasis );


    klMinMode  = getValue( tree, "KL_MIN_MODE", defaults.klMinMode );
    klMaxMode  = getValue( tree, "KL_MAX_MODE", defaults.klMaxMode );
    klCutoff = getValue( tree, "SVD_REG", defaults.klCutoff );
    nInitialModes  = getValue( tree, "MODE_START", defaults.nInitialModes );
    nModeIncrement  = getValue( tree, "MODE_STEP", defaults.nModeIncrement );
    
    modeList = getValue( tree, "MODES", defaults.modeList );
    modeList.setDefaultModeType( modeBasis );
    if( (runFlags&RF_SORT_MODES) && (modeBasis == KARHUNEN_LOEVE) ) {
        const std::map<uint16_t, Zernike::KLPtr>& kle = Zernike::karhunenLoeveExpansion( klMinMode, klMaxMode );
        ModeList tmpL;
        vector<Zernike::KLPtr> tmp;
        for( const auto& kl: kle ) tmp.push_back( kl.second );
        std::sort( tmp.begin(), tmp.end(), [](const Zernike::KLPtr& a, const Zernike::KLPtr& b){
            if( a->covariance == b->covariance ) return a->id < b->id;
            return a->covariance > b->covariance;
            
        });
        for( auto& id: modeList ) {
            if( id.mode<klMinMode || id.mode >= tmp.size() ) continue;
            ModeBase tp = KARHUNEN_LOEVE;
            uint16_t klID = tmp[id.mode-klMinMode]->id;
            if( (klID == 2) || (klID == 3) ) tp = ZERNIKE;
            tmpL.push_back( ModeID( klID, tp ) );
        }
        std::swap( tmpL, modeList );
    } else {
        for( auto& d: diversityModes ) {
            if( d.mode == 2 || d.mode == 3 ) d.type = ZERNIKE;      // Force Zernike for all tilt modes.
        }
    }
    nModes = modeList.size();

    telescopeD = getValue( tree, "TELESCOPE_D", defaults.telescopeD );
    telescopeCO = getValue( tree, "TELESCOPE_CO", defaults.telescopeCO );

    minIterations = getValue( tree, "MIN_ITER", defaults.minIterations );
    maxIterations = getValue( tree, "MAX_ITER", defaults.maxIterations );
    targetIterations = getValue( tree, "N_DONE_ITER", defaults.targetIterations );
    fillpixMethod = defaults.fillpixMethod;
    string tmpString = getValue<string>( tree, "FPMETHOD", "" );
    int tmpInt;
    if( tmpString.length() ) {
        tmpInt = getFromMap( tmpString, fillpixMap );
        if( tmpInt ) {
            fillpixMethod = tmpInt;
        } else {
            string msg = "Unrecognized FPMETHOD value \"" + tmpString + "\"\n  Valid entries are: ";
            for( const auto& entry: fillpixMap ) msg += "\"" + entry.first + "\" ";
            LOG_ERR << msg << ende;
        }
    }
    gradientMethod = defaults.gradientMethod;
    tmpString = getValue<string>( tree, "GRADIENT", "" );
    if( tmpString.length() ) {
        tmpInt = getFromMap( tmpString, gradientMap );
        if( tmpInt ) {
            gradientMethod = tmpInt;
        } else {
            string msg = "Unrecognized GRADIENT value \"" + tmpString + "\"\n  Valid entries are: ";
            for( const auto& entry: gradientMap ) msg += "\"" + entry.first + "\" ";
            LOG_ERR << msg << ende;
        }
    }
    getstepMethod = defaults.getstepMethod;
    tmpString = getValue<string>( tree, "GETSTEP", "" );
    if( tmpString.length() ) {
        tmpInt = getFromMap( tmpString, getstepMap );
        if( tmpInt ) {
            getstepMethod = tmpInt;
        } else {
            string msg = "Unrecognized GETSTEP value \"" + tmpString + "\"\n  Valid entries are: ";
            for( const auto& entry: getstepMap ) msg += "\"" + entry.first + "\" ";
            LOG_ERR << msg << ende;
        }
    }
    badPixelThreshold = getValue( tree, "BADPIXEL", defaults.badPixelThreshold );
    FTOL = getValue( tree, "FTOL", defaults.FTOL );
    EPS = getValue( tree, "EPS", defaults.EPS );
    reg_alpha = getValue( tree, "REG_ALPHA", defaults.reg_alpha );
    graddiff_step = getValue( tree, "GRADDIFF_STEP", defaults.graddiff_step );

    vector<FileType> filetypes = getValue( tree, "FILE_TYPE", vector<FileType>( 1, ( runFlags & RF_CALIBRATE ) ? FT_ANA : FT_FITS ) );
    for( const FileType& it : filetypes ) outputFileType |= it;
    if( ( outputFileType & FT_MASK ) == 0 ) {
        string msg = "\"FILE_TYPE\" has to be specified as ANA, FITS or MOMFBD. Multiple values supported, but not fully implemented yet.";
        throw logic_error( msg );
    }

    outputDataType = defaults.outputDataType;
    tmpString = getValue<string>( tree, "DATA_TYPE", dtTags[defaults.outputDataType] );
    if( iequals( tmpString, "FLOAT") ) outputDataType = DT_F32T;
    else if( iequals( tmpString, "SHORT") ) outputDataType = DT_I16T;
    else {
        LOG_WARN << "\"DATA_TYPE\" unrecognized data type \"" << tmpString
                 << "\", using default ( "  << dtTags[outputDataType] << ")" << ende;
    }

    sequenceNumber = getValue( tree, "SEQUENCE_NUM", defaults.sequenceNumber );
    observationTime = getValue<string>( tree, "TIME_OBS", defaults.observationTime );
    observationDate = getValue<string>( tree, "DATE_OBS", defaults.observationDate );
    //tmpDataDir = cleanPath( getValue<string>( tree, "PROG_DATA_DIR", defaults.tmpDataDir ) );
    tmpDataDir = getValue<string>( tree, "PROG_DATA_DIR", defaults.tmpDataDir );
    tmpString = getValue<string>( tree, "OUTPUT_FILES", "" );
    outputFiles = defaults.outputFiles;
    if( tmpString != "" ) {
        boost::split( outputFiles, tmpString, boost::is_any_of( "," ) );
    }
    
    if( tree.count( "INIT_FILES" ) ) {
        tmpString = getValue<string>( tree, "INIT_FILES", "" );
        initFiles = defaults.initFiles;
        if( tmpString != "" ) {
            boost::split( initFiles, tmpString, boost::is_any_of( "," ) );
        } else initFiles.resize( 1,"OUTPUT" );
    }

    if( runFlags & RF_CALIBRATE ) {
        //saveMask |= SF_SAVE_ALPHA; // necessary for calibration runs.
        outputFileType |= FT_ANA;
    }
    
    // Resolve relative paths etc.
    if( !tmpDataDir.empty() ) tmpDataDir = cleanPath( tmpDataDir );
    for( auto& f: outputFiles ) if( !f.empty() ) f = cleanPath( f );
    for( auto& f: initFiles ) if( !f.empty() ) f = cleanPath( f );
    
    ObjectCfg::parseProperties( tree, logger, defaults );

}


void GlobalCfg::getProperties( bpt::ptree& tree, const ChannelCfg& def, bool showAll ) const {
    
    const GlobalCfg& defaults = reinterpret_cast<const GlobalCfg&>(def);

    uint16_t diff = runFlags ^ defaults.runFlags;
    if( diff & RF_CALIBRATE ) tree.put( "CALIBRATE", bool( runFlags & RF_CALIBRATE ) );
    if( diff & RF_DONT_MATCH_IMAGE_NUMS ) tree.put( "DONT_MATCH_IMAGE_NUMS", bool( runFlags & RF_DONT_MATCH_IMAGE_NUMS ) );
    if( diff & RF_FAST_QR ) tree.put( "FAST_QR", bool( runFlags & RF_FAST_QR ) );
    if( diff & RF_FIT_PLANE ) tree.put( "FIT_PLANE", bool( runFlags & RF_FIT_PLANE ) );
    if( diff & RF_FLATFIELD ) tree.put( "FLATFIELD", bool( runFlags & RF_FLATFIELD ) );
    if( diff & RF_GLOBAL_NOISE ) tree.put( "GLOBAL_NOISE", bool( runFlags & RF_GLOBAL_NOISE ) );
    if( diff & RF_NEW_CONSTRAINTS ) tree.put( "NEW_CONSTRAINTS", bool( runFlags & RF_NEW_CONSTRAINTS ) );
    if( diff & RF_NO_CLIP ) tree.put( "NO_CLIP", bool( runFlags & RF_NO_CLIP ) );
    if( diff & RF_NO_CONSTRAINTS ) tree.put( "NO_CONSTRAINTS", bool( runFlags & RF_NO_CONSTRAINTS ) );
    if( diff & RF_NO_FILTER ) tree.put( "NO_FILTER", bool( runFlags & RF_NO_FILTER ) );
    if( diff & RF_FORCE_WRITE ) tree.put( "OVERWRITE", bool( runFlags & RF_FORCE_WRITE ) );
    if( diff & RF_NOSWAP ) tree.put( "NOSWAP", bool( runFlags & RF_NOSWAP ) );
    if( diff & RF_OLD_NS ) tree.put( "OLD_NS", bool( runFlags & RF_OLD_NS ) );
    if( diff & RF_SORT_MODES ) tree.put( "SORT_MODES", bool( runFlags & RF_SORT_MODES ) );

    if( showAll || trace != defaults.trace ) tree.put( "TRACE", trace );
    
    if( showAll || (modeBasis && (modeBasis != defaults.modeBasis)) ) tree.put( "BASIS", basisTags[modeBasis%3] );
    if( showAll || klMinMode != defaults.klMinMode ) tree.put( "KL_MIN_MODE", klMinMode );
    if( showAll || klMaxMode != defaults.klMaxMode ) tree.put( "KL_MAX_MODE", klMaxMode );
    if( showAll || klCutoff != defaults.klCutoff ) tree.put( "SVD_REG", klCutoff );
    if( showAll || nInitialModes != defaults.nInitialModes ) tree.put( "MODE_START", nInitialModes );
    if( showAll || nModeIncrement != defaults.nModeIncrement ) tree.put( "MODE_STEP", nModeIncrement );
    if( showAll || modeList != defaults.modeList ) tree.put( "MODES", modeList );
    
    if( showAll || telescopeD != defaults.telescopeD ) tree.put( "TELESCOPE_D", telescopeD );
    if( showAll || telescopeCO != defaults.telescopeCO ) tree.put( "TELESCOPE_CO", telescopeCO );
    
    if( showAll || minIterations != defaults.minIterations ) tree.put( "MIN_ITER", minIterations );
    if( showAll || maxIterations != defaults.maxIterations ) tree.put( "MAX_ITER", maxIterations );
    if( showAll || targetIterations != defaults.targetIterations ) tree.put( "N_DONE_ITER", targetIterations );
    if( showAll || fillpixMethod != defaults.fillpixMethod ) tree.put( "FPMETHOD", fpmTags[fillpixMethod%4] );
    if( showAll || gradientMethod != defaults.gradientMethod ) tree.put( "GRADIENT", gmTags[gradientMethod%3] );
    if( showAll || getstepMethod != defaults.getstepMethod ) tree.put( "GETSTEP", gsmTags[getstepMethod%5] );
    if( showAll || badPixelThreshold != defaults.badPixelThreshold ) tree.put( "BADPIXEL", badPixelThreshold );
    if( showAll || FTOL != defaults.FTOL ) tree.put( "FTOL", FTOL );
    if( showAll || EPS != defaults.EPS ) tree.put( "EPS", EPS );
    if( showAll || reg_alpha != defaults.reg_alpha ) tree.put( "REG_ALPHA", reg_alpha );
    if( showAll || graddiff_step != defaults.graddiff_step ) tree.put( "GRADDIFF_STEP", graddiff_step );
 
    if( showAll || outputFileType != ( ( runFlags & RF_CALIBRATE ) ? FT_ANA : FT_FITS) ) tree.put( "FILE_TYPE", ftTags[outputFileType%8] );
    if( showAll || outputDataType != defaults.outputDataType ) tree.put( "DATA_TYPE", dtTags[outputDataType%5] );
    if( showAll || sequenceNumber != defaults.sequenceNumber ) tree.put( "SEQUENCE_NUM", sequenceNumber );
    if( showAll || observationTime != defaults.observationTime ) tree.put( "TIME_OBS", observationTime );
    if( showAll || observationDate != defaults.observationDate ) tree.put( "DATE_OBS", observationDate );
    if( showAll || tmpDataDir != defaults.tmpDataDir ) tree.put( "PROG_DATA_DIR", tmpDataDir );
    if( showAll || outputFiles != defaults.outputFiles ) tree.put( "OUTPUT_FILES", outputFiles );
    if( showAll || initFiles != defaults.initFiles ) tree.put( "INIT_FILES", initFiles );

    ObjectCfg::getProperties( tree, defaults );

}


uint64_t GlobalCfg::size( void ) const {
    // static sizes (PoD types)
    uint64_t ssz = sizeof( badPixelThreshold ) + sizeof( EPS ) + sizeof( fillpixMethod ) + sizeof( FTOL )
                 + sizeof( getstepMethod ) + sizeof( graddiff_step ) + sizeof( gradientMethod ) + sizeof( klCutoff )
                 + sizeof( klMaxMode ) + sizeof( klMinMode ) + sizeof( maxIterations ) + sizeof( minIterations )
                 + sizeof( modeBasis ) + sizeof( nInitialModes ) + sizeof( nModeIncrement ) + sizeof( nModes )
                 + sizeof( outputFileType ) + sizeof( outputDataType ) + sizeof( reg_alpha ) + sizeof( runFlags )
                 + sizeof( sequenceNumber ) + sizeof( targetIterations ) + sizeof( telescopeCO )
                 + sizeof( telescopeD ) + sizeof( trace );
    uint64_t sz = ssz + ObjectCfg::size();
    // strings
    sz += observationTime.length() + observationDate.length() + tmpDataDir.length() + 3;
    // vectors
    sz += redux::util::size( modeList );
    sz += redux::util::size( outputFiles );
    return sz;
}


uint64_t GlobalCfg::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = ObjectCfg::pack( ptr );
    // scalar values
    count += pack( ptr+count, badPixelThreshold );
    count += pack( ptr+count, EPS );
    count += pack( ptr+count, fillpixMethod );
    count += pack( ptr+count, FTOL );
    count += pack( ptr+count, getstepMethod );
    count += pack( ptr+count, graddiff_step );
    count += pack( ptr+count, gradientMethod );
    count += pack( ptr+count, klCutoff );
    count += pack( ptr+count, klMaxMode );
    count += pack( ptr+count, klMinMode );
    count += pack( ptr+count, maxIterations );
    count += pack( ptr+count, minIterations );
    count += pack( ptr+count, modeBasis );
    count += pack( ptr+count, nInitialModes );
    count += pack( ptr+count, nModeIncrement );
    count += pack( ptr+count, nModes );
    count += pack( ptr+count, outputFileType );
    count += pack( ptr+count, outputDataType );
    count += pack( ptr+count, reg_alpha );
    count += pack( ptr+count, runFlags );
    count += pack( ptr+count, sequenceNumber );
    count += pack( ptr+count, targetIterations );
    count += pack( ptr+count, telescopeCO );
    count += pack( ptr+count, telescopeD );
    count += pack( ptr+count, trace );
    // strings
    count += pack( ptr+count, observationDate );
    count += pack( ptr+count, observationTime );
    count += pack( ptr+count, tmpDataDir );
    // vectors
    count += pack( ptr+count, modeList );
    count += pack( ptr+count, outputFiles );

    return count;
}


uint64_t GlobalCfg::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;

    uint64_t count = ObjectCfg::unpack( ptr, swap_endian );
    // scalar values
    count += unpack( ptr+count, badPixelThreshold, swap_endian );
    count += unpack( ptr+count, EPS, swap_endian );
    count += unpack( ptr+count, fillpixMethod );//b
    count += unpack( ptr+count, FTOL, swap_endian );
    count += unpack( ptr+count, getstepMethod );//b
    count += unpack( ptr+count, graddiff_step, swap_endian );
    count += unpack( ptr+count, gradientMethod );//b
    count += unpack( ptr+count, klCutoff, swap_endian );
    count += unpack( ptr+count, klMaxMode, swap_endian );
    count += unpack( ptr+count, klMinMode, swap_endian );
    count += unpack( ptr+count, maxIterations, swap_endian );
    count += unpack( ptr+count, minIterations, swap_endian );
    count += unpack( ptr+count, modeBasis, swap_endian );
    count += unpack( ptr+count, nInitialModes, swap_endian );
    count += unpack( ptr+count, nModeIncrement, swap_endian );
    count += unpack( ptr+count, nModes, swap_endian );
    count += unpack( ptr+count, outputFileType );//b
    count += unpack( ptr+count, outputDataType );//b
    count += unpack( ptr+count, reg_alpha, swap_endian );
    count += unpack( ptr+count, runFlags, swap_endian );
    count += unpack( ptr+count, sequenceNumber, swap_endian );
    count += unpack( ptr+count, targetIterations, swap_endian );
    count += unpack( ptr+count, telescopeCO, swap_endian );
    count += unpack( ptr+count, telescopeD, swap_endian );
    count += unpack( ptr+count, trace );//b
    // strings
    count += unpack( ptr+count, observationDate );
    count += unpack( ptr+count, observationTime );
    count += unpack( ptr+count, tmpDataDir );
    // vectors
    count += unpack( ptr+count, modeList, swap_endian );
    count += unpack( ptr+count, outputFiles, swap_endian );
    return count;
}


const GlobalCfg& GlobalCfg::operator=( const ObjectCfg& rhs ) {
    ObjectCfg::operator=( rhs );
    return *this;
}


const GlobalCfg& GlobalCfg::operator=( const ChannelCfg& rhs ) {
    ChannelCfg::operator=( rhs );
    return *this;
}


bool GlobalCfg::operator==( const GlobalCfg& rhs ) const {
    return ( runFlags == rhs.runFlags ) &&
           ( modeBasis == rhs.modeBasis ) &&
           ( klMinMode == rhs.klMinMode ) &&
           ( klMaxMode == rhs.klMaxMode ) &&
           ( klCutoff == rhs.klCutoff ) &&
           ( nInitialModes == rhs.nInitialModes ) &&
           ( nModeIncrement == rhs.nModeIncrement ) &&
           ( telescopeD == rhs.telescopeD ) &&
           ( telescopeCO == rhs.telescopeCO ) &&
           ( minIterations == rhs.minIterations ) &&
           ( maxIterations == rhs.maxIterations ) &&
           ( fillpixMethod == rhs.fillpixMethod ) &&
           ( gradientMethod == rhs.gradientMethod ) &&
           ( getstepMethod == rhs.getstepMethod ) &&
           ( badPixelThreshold == rhs.badPixelThreshold ) &&
           ( FTOL == rhs.FTOL ) &&
           ( EPS == rhs.EPS ) &&
           ( reg_alpha == rhs.reg_alpha ) &&
           ( outputFileType == rhs.outputFileType ) &&
           ( outputDataType == rhs.outputDataType ) &&
           ( sequenceNumber == rhs.sequenceNumber ) &&
           ( observationTime == rhs.observationTime ) &&
           ( observationDate == rhs.observationDate ) &&
           ( tmpDataDir == rhs.tmpDataDir ) &&
           ( modeList == rhs.modeList ) &&
           ObjectCfg::operator==( rhs ) &&
           ( outputFiles == rhs.outputFiles );
}
