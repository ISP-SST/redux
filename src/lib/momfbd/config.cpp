#include "redux/momfbd/config.hpp"

#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/constants.hpp"
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

using namespace redux::momfbd;
using namespace redux::util;
using namespace std;
using boost::algorithm::iequals;


namespace {

    boost::system::error_code ec;
    /*const*/ GlobalCfg globalDefaults;
    const string basisTags[] = {"","Zernike","Karhunen-Loeve"};
    const string fpmTags[] = {"","median","invdistweight","horint"};
    const string gmTags[] = {"","gradient_diff","gradient_Vogel"};
    const string gsmTags[] = {"","getstep_steepest_descent","getstep_conjugate_gradient","getstep_BFGS","getstep_BFGS_inv"};
    const string ftTags[] = {"","ANA","FITS","ANA,FITS","MOMFBD", "ANA,MOMFBD", "FITS,MOMFBD", "ANA,FITS,MOMFBD"};
    const string ftExt[] = {"","f0","fits","","momfbd"};
    const string dtTags[] = {"byte","short","int","int64","float","double"};

    int getFromMap(const string& str, const map<string, int>& m) {
        auto res = m.find(str);

        if(res != m.end()) {
            return res->second;
        }

        return 0;
    }
    
    template<typename T>
    void setFlag(T& flagset, T flag, bool value ) {
        if( value != bool(flagset&flag) ) {
            flagset ^= flag;
        }
    }


    void parseSegment( vector<uint16_t>& divs, vector<uint16_t>& types, string elem ) {
        size_t n = std::count( elem.begin(), elem.end(), '-' );
        uint16_t tp = 0;
        if( elem.find_first_of( "Zz" ) != string::npos ) tp |= ZERNIKE;
        if( elem.find_first_of( "Kk" ) != string::npos ) tp |= KARHUNEN_LOEVE;
        if( tp == 3 ) {
            //LLOG_ERR("channel") << "Different mode-types in specified mode range \"" << elem << "\"";
        } else if( tp == 0 ) tp = ZERNIKE;
        
        elem.erase( remove_if( elem.begin(), elem.end(),
                               [](const char&a ){ return (a == 'Z' || a == 'z' || a == 'K' || a == 'k'); } ), elem.end() );
        if( n == 0 ) {
            divs.push_back( boost::lexical_cast<uint16_t>( elem ) );
            types.push_back( tp );
            return;
        }
        else if( n == 1 ) {
            n = elem.find_first_of( '-' );
            uint16_t first = boost::lexical_cast<uint16_t>( elem.substr( 0, n ) );
            uint16_t last = boost::lexical_cast<uint16_t>( elem.substr( n + 1 ) );
            while( first <= last ) {
                divs.push_back( first++ );
                types.push_back( tp );
            }
        }
    }

    /*double def2cf( double pd_defocus, double telescope_r ) { // defocus distance in meters
        static const double tmp = ( 8.0 * sqrt( 3.0 ) );
        return -pd_defocus * redux::PI * telescope_r * telescope_r * tmp;
    }*/

}

const map<FileType, string> redux::momfbd::FileTypeNames = {
    {FT_ANA, ftTags[FT_ANA]},
    {FT_FITS, ftTags[FT_FITS]},
    {FT_MOMFBD, ftTags[FT_MOMFBD]}
};

const map<FileType, string> redux::momfbd::FileTypeExtensions = {
    {FT_ANA, ftExt[FT_ANA]},
    {FT_FITS, ftExt[FT_FITS]},
    {FT_MOMFBD, ftExt[FT_MOMFBD]}
};

const map<string, int> redux::momfbd::fillpixMap = {
    {fpmTags[FPM_MEDIAN], FPM_MEDIAN},
    {fpmTags[FPM_INVDISTWEIGHT], FPM_INVDISTWEIGHT},
    {fpmTags[FPM_HORINT], FPM_HORINT}
};

const map<string, int> redux::momfbd::gradientMap = {
    {gmTags[GM_DIFF], GM_DIFF},
    {gmTags[GM_VOGEL], GM_VOGEL}
};

const map<string, int> redux::momfbd::getstepMap = {
    {gsmTags[GSM_SDSC], GSM_SDSC},
    {gsmTags[GSM_CNJG], GSM_CNJG},
    {gsmTags[GSM_BFGS], GSM_BFGS},
    {gsmTags[GSM_BFGS_inv], GSM_BFGS_inv}
};







/********************  Channel  ********************/

ChannelCfg::ChannelCfg() : rotationAngle(0), noiseFudge(1), weight(1), borderClip(100), incomplete(0),
        mmRow(0), mmWidth(0), imageNumberOffset(0), logChannel("config") {

}


ChannelCfg::~ChannelCfg() {

}


ChannelCfg::operator std::string() const {
    bpt::ptree dump;
    std::stringstream ss;
    this->getProperties(dump);
    bpt::write_info( ss, dump );
    return ss.str();
}


void ChannelCfg::parseProperties(bpt::ptree& tree, const ChannelCfg& defaults) {
    
    rotationAngle = tree.get<double>("ANGLE", defaults.rotationAngle);

    noiseFudge = tree.get<double>("NF", defaults.noiseFudge);
    weight = tree.get<double>("WEIGHT", defaults.weight);
    
    // TODO: collect diversity settings in a struct and write a translator
    string tmpString = tree.get<string>( "DIVERSITY", "" );
    diversity.clear();
    diversityModes.clear();
    diversityTypes.clear();
    if( ! tmpString.empty() ) {
        double tmpD = 1.0;
        if( tmpString.find( "mm" ) != string::npos ) tmpD = 1.00E-03;
        else if( tmpString.find( "cm" ) != string::npos ) tmpD = 1.00E-02;
        //tmpString.erase( boost::remove_if( tmpString, boost::is_any_of( "cm\" " ) ), tmpString.end() );
        tmpString.erase( remove_if( tmpString.begin(), tmpString.end(),
                               [](const char&a ){ return (a == 'c' || a == 'm' || a == ' ' || a == '"'); } ), tmpString.end() );
        bpt::ptree tmpTree;                         // just to be able to use the VectorTranslator
        tmpTree.put( "tmp", tmpString );
        diversity = tmpTree.get<vector<double>>( "tmp", vector<double>() );
        tmpString = tree.get<string>( "DIV_ORDERS", "" );
        if( tmpString.empty() ) {
            if( diversity.size() > 1 ) {
                //LOG_ERR << "Multiple diversity coefficients specified, but no orders provided!";
            } else if( diversity.size() == 1 ) {  // A single diversity value is interpreted as (de)focus.
                diversityModes.resize( 1, 4 );
                diversityTypes.resize( 1,  ZERNIKE );
            }
        } else {
            vector<string> tokens;
            boost::split( tokens, tmpString, boost::is_any_of( "," ) );
            for( auto & token : tokens ) {
                parseSegment( diversityModes, diversityTypes, token );
            }
        }

        if( diversity.size() == diversityModes.size() ) {
            for(unsigned int i=0; i<diversity.size(); ++i) {
                if(diversityModes[i] == 4) {   // focus term, convert from physical length (including mm/cm) to coefficient
                    diversity[i] = tmpD*diversity[i]; // TODO: verify conversion: def2cf( tmpD*diversity[i], globalDefaults.telescopeD / telescopeF );
                }
            }
        } else {
            //LOG_ERR << "Number of diversity orders does not match number of diversity coefficients!";
        }
        
        if( diversity.size() ) {
//             LOG << "Diversity: " << printArray(diversityModes, "modes")
//                 << "    " << printArray(diversityTypes, "types")
//                 << "    " << printArray(diversity, "values");
        }
        
    }

    alignMap = tree.get<vector<float>>( "ALIGN_MAP", defaults.alignMap );
    alignClip = tree.get<vector<int16_t>>( "ALIGN_CLIP", defaults.alignClip );
    borderClip = tree.get<uint16_t>("BORDER_CLIP", defaults.borderClip);
    incomplete = tree.get<bool>( "INCOMPLETE", defaults.incomplete );
    
    subImagePosX = tree.get<vector<uint16_t>>( "SIM_X", defaults.subImagePosX );
    subImagePosY = tree.get<vector<uint16_t>>( "SIM_Y", defaults.subImagePosY );
    if( tree.get<bool>( "CAL_X", false ) ) {
        if( tree.get<bool>( "CAL_Y", false ) ) {
            //if( subImagePosX.size() || subImagePosY.size() ) LOG << "Note: SIM_X/SIM_Y replaced by CAL_X/CAL_Y";
            subImagePosX = tree.get<vector<uint16_t>>( "CAL_X", defaults.subImagePosX );
            subImagePosY = tree.get<vector<uint16_t>>( "CAL_Y", defaults.subImagePosY );
            if( subImagePosX.empty() || ( subImagePosX.size() != subImagePosY.size() ) ) {
                //LOG_ERR << "CAL_X and CAL_Y must have the same number of elements!";
            }
        }
        //else LOG_ERR << "CAL_Y must be provided if CAL_X is!";
    }

    
    imageDataDir = tree.get<string>("IMAGE_DATA_DIR", defaults.imageDataDir);
    imageTemplate = tree.get<string>( "FILENAME_TEMPLATE", defaults.imageTemplate );
    darkTemplate = tree.get<string>( "DARK_TEMPLATE", defaults.darkTemplate );
    gainFile = tree.get<string>( "GAIN_FILE", defaults.gainFile );
    responseFile = tree.get<string>( "CCD_RESPONSE", defaults.responseFile );
    backgainFile = tree.get<string>( "BACK_GAIN", defaults.backgainFile );
    psfFile = tree.get<string>( "PSF", defaults.psfFile );
    mmFile = tree.get<string>( "MODMAT", defaults.mmFile );
    mmRow = tree.get<uint8_t>( "MMROW", defaults.mmRow );
    mmWidth = tree.get<uint8_t>( "MMWIDTH",defaults.mmWidth);
    xOffsetFile = tree.get<string>( "XOFFSET", "" );
    yOffsetFile = tree.get<string>( "YOFFSET", "" );
    imageNumberOffset = tree.get<uint32_t>( "DT", defaults.imageNumberOffset );
    fileNumbers = tree.get<vector<uint32_t>>("IMAGE_NUM", defaults.fileNumbers);
    if( fileNumbers.empty() ) {
        fileNumbers = tree.get<vector<uint32_t>>("IMAGE_NUMS", defaults.fileNumbers);
    }
    waveFronts = tree.get<vector<uint32_t>>( "WFINDEX", defaults.waveFronts );
    darkNumbers = tree.get<vector<uint32_t>>("DARK_NUM", defaults.darkNumbers);
    stokesWeights = tree.get<vector<float>>( "VECTOR", defaults.stokesWeights );
    if( mmFile.length() > 0 ) {
        if( !mmRow ) {
            //LOG_FATAL << "a modulation matrix was provided but no row specified (MMROW).";
        }
        if( !mmWidth ) {
            //LOG_FATAL << "modulation matrix dimensions cannot be autodetected (yet): you must provide the matrix width (MMWIDTH)!";
        }
        if( stokesWeights.size() == 0 ) {
            //LOG_ERR << "modulation matrix specified but no VECTOR input given!";
        } else if( stokesWeights.size() != mmWidth ) {
            //LOG_ERR << "VECTOR input has " << stokesWeights.size() << " elements, but MMWIDTH=" << (int)mmWidth;
        }
    }
    else {  // TODO: don't modify cfg values!! ...make the main code use weight 1 as default instead.
        mmRow = mmWidth = 1;
        stokesWeights.resize( 1, 1.0 );
    }


}


void ChannelCfg::getProperties(bpt::ptree& tree, const ChannelCfg& defaults) const {

    if(rotationAngle != defaults.rotationAngle) tree.put("ANGLE", rotationAngle);

    if(noiseFudge != defaults.noiseFudge) tree.put("NF", noiseFudge);
    if(weight != defaults.weight) tree.put("WEIGHT", weight);
    
    if(alignMap != defaults.alignMap) tree.put("ALIGN_MAP", alignMap);
    if(alignClip != defaults.alignClip) tree.put("ALIGN_CLIP", alignClip);
    if(borderClip != defaults.borderClip) tree.put("BORDER_CLIP", borderClip);
    if(incomplete != defaults.incomplete) tree.put("INCOMPLETE", (bool)incomplete);

    if(subImagePosX != defaults.subImagePosX) tree.put("SIM_X", subImagePosX);
    if(subImagePosY != defaults.subImagePosY) tree.put("SIM_Y", subImagePosY);

    if(imageDataDir != defaults.imageDataDir) tree.put("IMAGE_DATA_DIR", imageDataDir);
    if(imageTemplate != defaults.imageTemplate) tree.put("FILENAME_TEMPLATE", imageTemplate);
    if(darkTemplate != defaults.darkTemplate) tree.put("DARK_TEMPLATE", darkTemplate);
    if(gainFile != defaults.gainFile) tree.put("GAIN_FILE", gainFile);
    if(responseFile != defaults.responseFile) tree.put("CCD_RESPONSE", responseFile);
    if(backgainFile != defaults.backgainFile) tree.put("BACK_GAIN", backgainFile);
    if(psfFile != defaults.psfFile) tree.put("PSF", psfFile);
    if(mmFile != defaults.mmFile) tree.put("MODMAT", mmFile);
    if(mmRow != defaults.mmRow) tree.put("MMROW", mmRow);
    if(mmWidth != defaults.mmWidth) tree.put("MMWIDTH", mmWidth);
    
    if(xOffsetFile != defaults.xOffsetFile) tree.put("XOFFSET", xOffsetFile);
    if(yOffsetFile != defaults.yOffsetFile) tree.put("YOFFSET", yOffsetFile);

    if(imageNumberOffset != defaults.imageNumberOffset) tree.put("DT", imageNumberOffset);
    if(fileNumbers != defaults.fileNumbers) tree.put("IMAGE_NUM", fileNumbers);
    if(waveFronts != defaults.waveFronts) tree.put("WFINDEX", waveFronts);
    if(darkNumbers != defaults.darkNumbers) tree.put("DARK_NUM", darkNumbers);
    
    if(stokesWeights != defaults.stokesWeights) tree.put("VECTOR", stokesWeights);

}


uint64_t ChannelCfg::size(void) const {
    uint64_t sz = 3*sizeof(double);         // rotationAngle, weight, noiseFudge
    sz += sizeof(uint16_t)+3;             // borderClip, incomplete, mmRow, mmWidth
    sz += sizeof(uint32_t);                 // imageNumberOffset
    sz += subImagePosX.size()*sizeof(uint16_t) + sizeof(uint64_t);
    sz += subImagePosY.size()*sizeof(uint16_t) + sizeof(uint64_t);
    sz += diversity.size() * sizeof( double ) + sizeof( uint64_t );
    sz += diversityModes.size() * sizeof( uint16_t ) + sizeof( uint64_t );
    sz += diversityTypes.size() * sizeof( uint16_t ) + sizeof( uint64_t );
    sz += alignMap.size()*sizeof(float) + sizeof(uint64_t);
    sz += alignClip.size()*sizeof(int16_t) + sizeof(uint64_t);
    sz += imageDataDir.length() + 1;
    sz += imageTemplate.length() + darkTemplate.length() + gainFile.length() + 3;
    sz += responseFile.length() + backgainFile.length() + psfFile.length() + mmFile.length() + 4;
    sz += xOffsetFile.length() + yOffsetFile.length() + 2;
    sz += fileNumbers.size()*sizeof(uint32_t) + sizeof(uint64_t);
    sz += waveFronts.size()*sizeof(uint32_t) + sizeof(uint64_t);
    sz += darkNumbers.size()*sizeof(uint32_t) + sizeof(uint64_t);
    sz += stokesWeights.size()*sizeof(float) + sizeof(uint64_t);
    return sz;
}


uint64_t ChannelCfg::pack(char* ptr) const {
    using redux::util::pack;
    uint64_t count = pack(ptr, rotationAngle);
    count += pack(ptr+count, noiseFudge);
    count += pack(ptr+count, weight);
    count += pack(ptr+count, diversity);
    count += pack(ptr+count, diversityModes);
    count += pack(ptr+count, diversityTypes);
    count += pack(ptr+count, alignMap);
    count += pack(ptr+count, alignClip);
    count += pack(ptr+count, borderClip);
    count += pack(ptr+count, incomplete);
    count += pack(ptr+count, subImagePosX);
    count += pack(ptr+count, subImagePosY);
    count += pack(ptr+count, imageDataDir);
    count += pack(ptr+count, imageTemplate);
    count += pack(ptr+count, darkTemplate);
    count += pack(ptr+count, gainFile);
    count += pack(ptr+count, responseFile);
    count += pack(ptr+count, backgainFile);
    count += pack(ptr+count, psfFile);
    count += pack(ptr+count, mmFile);
    count += pack(ptr+count, mmRow);
    count += pack(ptr+count, mmWidth);
    count += pack(ptr+count, xOffsetFile);
    count += pack(ptr+count, yOffsetFile);
    count += pack(ptr+count, imageNumberOffset);
    count += pack(ptr+count, fileNumbers);
    count += pack(ptr+count, waveFronts);
    count += pack(ptr+count, darkNumbers);
    count += pack(ptr+count, stokesWeights);
    return count;
}


uint64_t ChannelCfg::unpack(const char* ptr, bool swap_endian) {
    using redux::util::unpack;
    uint64_t count = unpack(ptr, rotationAngle, swap_endian);
    count += unpack(ptr+count, noiseFudge, swap_endian);
    count += unpack(ptr+count, weight, swap_endian);
    count += unpack(ptr+count, diversity, swap_endian);
    count += unpack(ptr+count, diversityModes, swap_endian);
    count += unpack(ptr+count, diversityTypes, swap_endian);
    count += unpack(ptr+count, alignMap, swap_endian);
    count += unpack(ptr+count, alignClip, swap_endian);
    count += unpack(ptr+count, borderClip, swap_endian);
    count += unpack(ptr+count, incomplete);
    count += unpack(ptr+count, subImagePosX, swap_endian);
    count += unpack(ptr+count, subImagePosY, swap_endian);
    count += unpack(ptr+count, imageDataDir);
    count += unpack(ptr+count, imageTemplate);
    count += unpack(ptr+count, darkTemplate);
    count += unpack(ptr+count, gainFile);
    count += unpack(ptr+count, responseFile);
    count += unpack(ptr+count, backgainFile);
    count += unpack(ptr+count, psfFile);
    count += unpack(ptr+count, mmFile);
    count += unpack(ptr+count, mmRow);
    count += unpack(ptr+count, mmWidth);
    count += unpack(ptr+count, xOffsetFile);
    count += unpack(ptr+count, yOffsetFile);
    count += unpack(ptr+count, imageNumberOffset, swap_endian);
    count += unpack(ptr+count, fileNumbers, swap_endian);
    count += unpack(ptr+count, waveFronts, swap_endian);
    count += unpack(ptr+count, darkNumbers, swap_endian);
    count += unpack(ptr+count, stokesWeights, swap_endian);
    return count;
}


bool ChannelCfg::operator==(const ChannelCfg& rhs) const {
    return (rotationAngle == rhs.rotationAngle) &&
           (noiseFudge == rhs.noiseFudge) &&
           (weight == rhs.weight) &&
           (borderClip == rhs.borderClip) &&
           (incomplete == rhs.incomplete) &&
           (subImagePosX == rhs.subImagePosX) &&
           (subImagePosY == rhs.subImagePosY) &&
           (imageDataDir == rhs.imageDataDir) &&
           (darkTemplate == rhs.darkTemplate) &&
           (imageTemplate == rhs.imageTemplate) &&
           (fileNumbers == rhs.fileNumbers) &&
           (waveFronts == rhs.waveFronts) &&
           (darkNumbers == rhs.darkNumbers);
}



/********************   Object  ********************/

ObjectCfg::ObjectCfg() : telescopeF(0), arcSecsPerPixel(0), pixelSize(1E-5),
                         alphaToPixels(0), pixelsToAlpha(0),
                         alphaToDefocus(0), defocusToAlpha(0),
                         maxLocalShift(5), minimumOverlap(16), 
                         patchSize(128), pupilPixels(64), saveMask(0), wavelength(0) {

}


ObjectCfg::~ObjectCfg() {

}


void ObjectCfg::parseProperties(bpt::ptree& tree, const ChannelCfg& def) {
    
    const ObjectCfg& defaults = reinterpret_cast<const ObjectCfg&>(def);

    telescopeF = tree.get<double>("TELESCOPE_F", defaults.telescopeF);
    arcSecsPerPixel = tree.get<double>("ARCSECPERPIX", defaults.arcSecsPerPixel);
    pixelSize = tree.get<double>("PIXELSIZE", defaults.pixelSize);

    maxLocalShift = tree.get<uint16_t>("MAX_LOCAL_SHIFT", defaults.maxLocalShift);
    minimumOverlap = tree.get<uint16_t>("MINIMUM_OVERLAP", defaults.minimumOverlap);
    patchSize  = tree.get<uint16_t>("NUM_POINTS", defaults.patchSize);
    pupilPixels  = tree.get<uint16_t>("PUPIL_POINTS", defaults.pupilPixels);
    
    saveMask = 0;
    if( tree.get<bool>( "GET_ALPHA", defaults.saveMask&SF_SAVE_ALPHA ) ) saveMask |= SF_SAVE_ALPHA;
    if( tree.get<bool>( "GET_COBJ", defaults.saveMask&SF_SAVE_COBJ ) ) saveMask |= SF_SAVE_COBJ;
    if( tree.get<bool>( "GET_DIVERSITY", defaults.saveMask&SF_SAVE_DIVERSITY ) ) saveMask |= SF_SAVE_DIVERSITY;
    if( tree.get<bool>( "GET_METRIC", defaults.saveMask&SF_SAVE_METRIC ) ) saveMask |= SF_SAVE_METRIC;
    if( tree.get<bool>( "GET_MODES", defaults.saveMask&SF_SAVE_MODES ) ) saveMask |= SF_SAVE_MODES;
    if( tree.get<bool>( "GET_PSF", defaults.saveMask&SF_SAVE_PSF ) ) saveMask |= SF_SAVE_PSF;
    if( tree.get<bool>( "GET_PSF_AVG", defaults.saveMask&SF_SAVE_PSF_AVG ) ) saveMask |= SF_SAVE_PSF_AVG;
    if( tree.get<bool>( "GET_RESIDUAL", defaults.saveMask&SF_SAVE_RESIDUAL ) ) saveMask |= SF_SAVE_RESIDUAL;
    if( tree.get<bool>( "SAVE_FFDATA", defaults.saveMask&SF_SAVE_FFDATA ) ) saveMask |= SF_SAVE_FFDATA;
    outputFileName = tree.get<string>("OUTPUT_FILE", defaults.outputFileName);
    modeFile = tree.get<string>("MODE_FILE", defaults.modeFile);
    pupilFile = tree.get<string>("PUPIL", defaults.pupilFile);
    wavelength = tree.get<double>("WAVELENGTH", defaults.wavelength);

    if( ( saveMask & SF_SAVE_PSF ) && ( saveMask & SF_SAVE_PSF_AVG ) ) {
        //LOG_WARN << "both GET_PSF and GET_PSF_AVG mode requested";
    }

    ChannelCfg::parseProperties(tree, defaults);

}


void ObjectCfg::getProperties(bpt::ptree& tree, const ChannelCfg& def) const {

    const ObjectCfg& defaults = reinterpret_cast<const ObjectCfg&>(def);
    
    if(telescopeF != defaults.telescopeF) tree.put("TELESCOPE_F", telescopeF);
    if(arcSecsPerPixel != defaults.arcSecsPerPixel) tree.put("ARCSECPERPIX", arcSecsPerPixel);
    if(pixelSize != defaults.pixelSize) tree.put("PIXELSIZE", pixelSize);
    if(maxLocalShift != defaults.maxLocalShift) tree.put("MAX_LOCAL_SHIFT", maxLocalShift);
    if(minimumOverlap != defaults.minimumOverlap) tree.put("MINIMUM_OVERLAP", minimumOverlap);
    if(patchSize != defaults.patchSize) tree.put("NUM_POINTS", patchSize);
    if(pupilPixels != defaults.pupilPixels) tree.put("PUPIL_POINTS", pupilPixels);

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
    if(outputFileName != defaults.outputFileName) tree.put("OUTPUT_FILE", outputFileName);
    if(modeFile != defaults.modeFile) tree.put("MODE_FILE", modeFile);
    if(pupilFile != defaults.pupilFile) tree.put("PUPIL", pupilFile);
    if(wavelength != defaults.wavelength) tree.put("WAVELENGTH", wavelength);

    ChannelCfg::getProperties(tree, defaults);

}


uint64_t ObjectCfg::size(void) const {
    uint64_t sz = ChannelCfg::size();
    sz += 5*sizeof(uint16_t);
    sz += outputFileName.length() + 1;
    sz += modeFile.length() + pupilFile.length() + 2;
    sz += 8*sizeof(double);
    return sz;
}


uint64_t ObjectCfg::pack(char* ptr) const {
    using redux::util::pack;
    uint64_t count = ChannelCfg::pack(ptr);
    count += pack(ptr+count, telescopeF);
    count += pack(ptr+count, arcSecsPerPixel);
    count += pack(ptr+count, pixelSize);
    count += pack(ptr+count, alphaToPixels);
    count += pack(ptr+count, pixelsToAlpha);
    count += pack(ptr+count, alphaToDefocus);
    count += pack(ptr+count, defocusToAlpha);
    count += pack(ptr+count, maxLocalShift);
    count += pack(ptr+count, minimumOverlap);
    count += pack(ptr+count, patchSize);
    count += pack(ptr+count, pupilPixels);
    count += pack(ptr+count, saveMask);
    count += pack(ptr+count, outputFileName);
    count += pack(ptr+count, modeFile);
    count += pack(ptr+count, pupilFile);
    count += pack(ptr+count, wavelength);
    return count;
}


uint64_t ObjectCfg::unpack(const char* ptr, bool swap_endian) {
    using redux::util::unpack;
    uint64_t count = ChannelCfg::unpack(ptr, swap_endian);
    count += unpack(ptr+count, telescopeF, swap_endian);
    count += unpack(ptr+count, arcSecsPerPixel, swap_endian);
    count += unpack(ptr+count, pixelSize, swap_endian);
    count += unpack(ptr+count, alphaToPixels, swap_endian);
    count += unpack(ptr+count, pixelsToAlpha, swap_endian);
    count += unpack(ptr+count, alphaToDefocus, swap_endian);
    count += unpack(ptr+count, defocusToAlpha, swap_endian);
    count += unpack(ptr+count, maxLocalShift, swap_endian);
    count += unpack(ptr+count, minimumOverlap, swap_endian);
    count += unpack(ptr+count, patchSize, swap_endian);
    count += unpack(ptr+count, pupilPixels, swap_endian);
    count += unpack(ptr+count, saveMask, swap_endian);
    count += unpack(ptr+count, outputFileName);
    count += unpack(ptr+count, modeFile);
    count += unpack(ptr+count, pupilFile);
    count += unpack(ptr+count, wavelength, swap_endian);
    return count;
}


const ObjectCfg& ObjectCfg::operator=(const ChannelCfg& rhs) {
    ChannelCfg::operator=(rhs);
    return *this;
}


bool ObjectCfg::operator==(const ObjectCfg& rhs) const {
    return (saveMask == rhs.saveMask) &&
           (patchSize == rhs.patchSize) &&
           (pupilPixels == rhs.pupilPixels) &&
           (maxLocalShift == rhs.maxLocalShift) &&
           (wavelength == rhs.wavelength) &&
           (outputFileName == rhs.outputFileName) &&
           ChannelCfg::operator==(rhs);
}



/********************   Global   ********************/

GlobalCfg::GlobalCfg() : runFlags(0), modeBasis(ZERNIKE), klMinMode(2), klMaxMode(2000), klCutoff(1E-3),
    nInitialModes(5), nModeIncrement(5),
    modeNumbers({ 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35 }),
    telescopeD(0), minIterations(5), maxIterations(500), targetIterations(3),
    fillpixMethod(FPM_INVDISTWEIGHT), gradientMethod(GM_DIFF), getstepMethod(GSM_BFGS_inv),
    badPixelThreshold(1E-5), FTOL(1E-03), EPS(1E-10), reg_alpha(0), graddiff_step(1E-2),
    outputFileType(FT_NONE), outputDataType(DT_F32T), sequenceNumber(0),
    observationTime(""), observationDate("N/A"), tmpDataDir("./data") {


}


GlobalCfg::~GlobalCfg() {

}


void GlobalCfg::parseProperties(bpt::ptree& tree, const ChannelCfg& def) {
    
    const GlobalCfg& defaults = reinterpret_cast<const GlobalCfg&>(def);

    if( tree.get<bool>( "CALIBRATE", false ) )             runFlags |= RF_CALIBRATE;
    if( tree.get<bool>( "DONT_MATCH_IMAGE_NUMS", false ) ) runFlags |= RF_DONT_MATCH_IMAGE_NUMS;
    if( tree.get<bool>( "FAST_QR", false ) )               runFlags |= RF_FAST_QR;
    if( tree.get<bool>( "FIT_PLANE", false ) )             runFlags |= RF_FIT_PLANE;
    if( tree.get<bool>( "FLATFIELD", false ) )             runFlags |= RF_FLATFIELD;
    if( tree.get<bool>( "GLOBAL_NOISE", false ) )          runFlags |= RF_GLOBAL_NOISE;
    if( tree.get<bool>( "NEW_CONSTRAINTS", false ) )       runFlags |= RF_NEW_CONSTRAINTS;
    if( tree.get<bool>( "NO_CLIP", false ) )               runFlags |= RF_NO_CLIP;
    if( tree.get<bool>( "NO_CONSTRAINTS", false ) )        runFlags |= RF_NO_CONSTRAINTS;
    if( tree.get<bool>( "NO_FILTER", false ) )             runFlags |= RF_NO_FILTER;
    if( tree.get<bool>( "OVERWRITE", false ) )             runFlags |= RF_FORCE_WRITE;
    if( tree.get<bool>( "SWAP", false ) )                  runFlags |= RF_SWAP;
    
/*    if( ( runFlags & RF_CALIBRATE ) && ( runFlags & RF_FLATFIELD ) ) {
        LOG_WARN << "both FLATFIELD and CALIBRATE mode requested, forcing CALIBRATE";
        runFlags &= ~RF_FLATFIELD;
    }
    if( ( runFlags & RF_CALIBRATE ) && ( runFlags & RF_NEW_CONSTRAINTS ) ) {
        LOG_WARN << "calibration mode uses old style constraints, ignoring NEW_CONSTRAINTS";
        runFlags &= ~RF_NEW_CONSTRAINTS;
    }
*/
    string tmpString = tree.get<string>("BASIS", "");
    modeBasis = defaults.modeBasis;
    if(tmpString.length()) {
        if(iequals(tmpString, "Karhunen-Loeve")) {
            modeBasis = KARHUNEN_LOEVE;
        } else
            if(iequals(tmpString, "Zernike")) {
                modeBasis = ZERNIKE;
            } else {
                //LOG_ERR << "Unrecognized BASIS value \"" << tmpString << "\", using default \"" << basisTags[globalDefaults.modeBasis] << "\"";
                //modeBasis = defaults.modeBasis;
            }
    }

    klMinMode  = tree.get<uint16_t>("KL_MIN_MODE", defaults.klMinMode);
    klMaxMode  = tree.get<uint16_t>("KL_MAX_MODE", defaults.klMaxMode);
    klCutoff = tree.get<float>("SVD_REG", defaults.klCutoff);
    nInitialModes  = tree.get<uint16_t>("MODE_START", defaults.nInitialModes);
    nModeIncrement  = tree.get<uint16_t>("MODE_STEP", defaults.nModeIncrement);
    modeNumbers = tree.get<vector<uint16_t>>("MODES", defaults.modeNumbers);

    telescopeD = tree.get<double>("TELESCOPE_D", defaults.telescopeD);

    minIterations = tree.get<uint16_t>("MIN_ITER", defaults.minIterations);
    maxIterations = tree.get<uint16_t>("MAX_ITER", defaults.maxIterations);
    targetIterations = tree.get<uint16_t>("N_DONE_ITER", defaults.targetIterations);
    fillpixMethod = defaults.fillpixMethod;
    tmpString = tree.get<string>("FPMETHOD", "");
    int tmpInt;
    if(tmpString.length()) {
        tmpInt = getFromMap(tmpString, fillpixMap);
        if(tmpInt) {
            fillpixMethod = tmpInt;
        } else {
            string msg = "Unrecognized FPMETHOD value \"" + tmpString + "\"\n  Valid entries are: ";
            for(const auto& entry: fillpixMap) msg += "\"" + entry.first + "\" ";
            //LOG_ERR << msg;
        }
    }
    gradientMethod = defaults.gradientMethod;
    tmpString = tree.get<string>("GRADIENT", "");
    if(tmpString.length()) {
        tmpInt = getFromMap(tmpString, gradientMap);
        if(tmpInt) {
            gradientMethod = tmpInt;
        } else {
            string msg = "Unrecognized GRADIENT value \"" + tmpString + "\"\n  Valid entries are: ";
            for(const auto& entry: gradientMap) msg += "\"" + entry.first + "\" ";
            //LOG_ERR << msg;
        }
    }
    getstepMethod = defaults.getstepMethod;
    tmpString = tree.get<string>("GETSTEP", "");
    if(tmpString.length()) {
        tmpInt = getFromMap(tmpString, getstepMap);
        if(tmpInt) {
            getstepMethod = tmpInt;
        } else {
            string msg = "Unrecognized GETSTEP value \"" + tmpString + "\"\n  Valid entries are: ";
            for(const auto& entry: getstepMap) msg += "\"" + entry.first + "\" ";
            //LOG_ERR << msg;
        }
    }
    badPixelThreshold = tree.get<float>("BADPIXEL", defaults.badPixelThreshold);
    FTOL = tree.get<float>("FTOL", defaults.FTOL);
    EPS = tree.get<float>("EPS", defaults.EPS);
    reg_alpha = tree.get<float>("REG_ALPHA", defaults.reg_alpha);
    graddiff_step = tree.get<float>("GRADDIFF_STEP", defaults.graddiff_step);

    vector<FileType> filetypes = tree.get<vector<FileType>>("FILE_TYPE", vector<FileType>(1, (runFlags & RF_CALIBRATE) ? FT_ANA : FT_FITS));
    for(const FileType& it : filetypes) outputFileType |= it;
    if((outputFileType & FT_MASK) == 0) cout << "\"FILE_TYPE\" has to be one of ANA/FITS/MOMFBD." << endl;

    tmpString = tree.get<string>("DATA_TYPE", dtTags[defaults.outputDataType]);
    if(iequals(tmpString, "FLOAT")) outputDataType = DT_F32T;
    else if(iequals(tmpString, "SHORT")) outputDataType = DT_I16T;
    else {
        //LOG_WARN << "\"DATA_TYPE\" unrecognized data type \"" << tmpString << "\", using default (" +dtTags[defaults.outputDataType]+ ")";
    }

    sequenceNumber = tree.get<uint32_t>("SEQUENCE_NUM", defaults.sequenceNumber);
    observationTime = tree.get<string>( "TIME_OBS", defaults.observationTime );
    observationDate = tree.get<string>( "DATE_OBS", defaults.observationDate );
    //tmpDataDir = cleanPath( tree.get<string>( "PROG_DATA_DIR", defaults.tmpDataDir ) );
    tmpDataDir = tree.get<string>( "PROG_DATA_DIR", defaults.tmpDataDir );
    tmpString = tree.get<string>( "OUTPUT_FILES", "" );
    outputFiles = defaults.outputFiles;
    if( tmpString != "" ) {
        boost::split( outputFiles, tmpString, boost::is_any_of( "," ) );
    }

    if( runFlags & RF_CALIBRATE ) {
        //saveMask |= SF_SAVE_ALPHA; // necessary for calibration runs.
        outputFileType |= FT_ANA;
    }
    
    ObjectCfg::parseProperties(tree, defaults);

}


void GlobalCfg::getProperties(bpt::ptree& tree, const ChannelCfg& def) const {
    
    const GlobalCfg& defaults = reinterpret_cast<const GlobalCfg&>(def);
    
/*
    defaults.modeBasis = 123;
    defaults.klMinMode = 123;
    defaults.klMaxMode = 123;
    defaults.klCutoff = 123;
    defaults.nInitialModes = 123;
    defaults.nModeIncrement = 123;
    defaults.modeNumbers = {17,66};
    defaults.telescopeD = 123;
    defaults.telescopeF = 123;
    defaults.minIterations = 123;
    defaults.maxIterations = 123;
    defaults.fillpixMethod = 123;
    defaults.gradientMethod = 123;
    defaults.getstepMethod = 123;
    defaults.badPixelThreshold = 123;
    defaults.FTOL = 123;
    defaults.EPS = 123;
    defaults.reg_gamma = 123;
    defaults.outputFileType = 123;
    defaults.outputDataType = 123;
    defaults.sequenceNumber = 123;
    defaults.observationTime = "time";
    defaults.observationDate = "date";
    defaults.tmpDataDir = "ddir";
    defaults.outputFiles = {"file1", "file2"};

    defaults.wavelength = 123;

    defaults.arcSecsPerPixel = 123;
    defaults.pixelSize = 123;
    defaults.rotationAngle = 123;
    defaults.weight = 123;
*/
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
    if( diff & RF_SWAP ) tree.put( "SWAP", bool( runFlags & RF_SWAP ) );

    if(modeBasis && modeBasis != defaults.modeBasis) tree.put("BASIS", basisTags[modeBasis%3]);
    if(klMinMode != defaults.klMinMode) tree.put("KL_MIN_MODE", klMinMode);
    if(klMaxMode != defaults.klMaxMode) tree.put("KL_MAX_MODE", klMaxMode);
    if(klCutoff != defaults.klCutoff) tree.put("SVD_REG", klCutoff);
    if(nInitialModes != defaults.nInitialModes) tree.put("MODE_START", nInitialModes);
    if(nModeIncrement != defaults.nModeIncrement) tree.put("MODE_STEP", nModeIncrement);
    if(modeNumbers != defaults.modeNumbers) tree.put("MODES", modeNumbers);
    
    if(telescopeD != defaults.telescopeD) tree.put("TELESCOPE_D", telescopeD);
    
    if(minIterations != defaults.minIterations) tree.put("MIN_ITER", minIterations);
    if(maxIterations != defaults.maxIterations) tree.put("MAX_ITER", maxIterations);
    if(targetIterations != defaults.targetIterations) tree.put("N_DONE_ITER", targetIterations);
    if(fillpixMethod != defaults.fillpixMethod) tree.put("FPMETHOD", fpmTags[fillpixMethod%4]);
    if(gradientMethod != defaults.gradientMethod) tree.put("GRADIENT", gmTags[gradientMethod%3]);
    if(getstepMethod != defaults.getstepMethod) tree.put("GETSTEP", gsmTags[getstepMethod%5]);
    if(badPixelThreshold != defaults.badPixelThreshold) tree.put("BADPIXEL", badPixelThreshold);
    if(FTOL != defaults.FTOL) tree.put("FTOL", FTOL);
    if(EPS != defaults.EPS) tree.put("EPS", EPS);
    if(reg_alpha != defaults.reg_alpha) tree.put("REG_ALPHA", reg_alpha);
    if(graddiff_step != defaults.graddiff_step) tree.put("GRADDIFF_STEP", graddiff_step);
 
    if(outputFileType != ((runFlags & RF_CALIBRATE) ? FT_ANA : FT_FITS)) tree.put("FILE_TYPE", ftTags[outputFileType%8]);
    if(outputDataType != defaults.outputDataType) tree.put("DATA_TYPE", dtTags[outputDataType%5]);
    if(sequenceNumber != defaults.sequenceNumber) tree.put("SEQUENCE_NUM", sequenceNumber);
    if(observationTime != defaults.observationTime) tree.put("TIME_OBS", observationTime);
    if(observationDate != defaults.observationDate) tree.put("DATE_OBS", observationDate);
    if(tmpDataDir != defaults.tmpDataDir) tree.put("PROG_DATA_DIR", tmpDataDir);
    if(outputFiles != defaults.outputFiles) tree.put("OUTPUT_FILES", outputFiles);

    ObjectCfg::getProperties(tree, defaults);

}


uint64_t GlobalCfg::size(void) const {
    uint64_t sz = ObjectCfg::size();
    sz += 6*sizeof(uint8_t);                 // modeBasis, fillpixMethod, gradientMethod, getstepMethod, outputFileType, outputDataType
    sz += 9*sizeof(uint16_t);                // runFlags, klMinMode, klMaxMode, nInitialModes, nModeIncrement, minIterations, maxIterations, targetIterations, outputFiles.size()
    sz += modeNumbers.size()*sizeof(uint16_t) + sizeof(uint64_t);
    sz += 6*sizeof(float);                   // klCutoff, badPixelThreshold, FTOL, EPS, reg_gamma, graddiff_step
    sz += sizeof(double);                    // telescopeD
    sz += sizeof(uint32_t);                  // sequenceNumber
    sz += observationTime.length() + 1;
    sz += observationDate.length() + 1;
    sz += tmpDataDir.length() + 1;
    for( const auto& filename : outputFiles ) {
        sz += filename.length() + 1;
    }
    return sz;
}


uint64_t GlobalCfg::pack(char* ptr) const {
    using redux::util::pack;
    uint64_t count = ObjectCfg::pack(ptr);
    count += pack(ptr+count, runFlags);
    count += pack(ptr+count, modeBasis);
    count += pack(ptr+count, klMinMode);
    count += pack(ptr+count, klMaxMode);
    count += pack(ptr+count, klCutoff);
    count += pack(ptr+count, nInitialModes);
    count += pack(ptr+count, nModeIncrement);
    count += pack(ptr+count, modeNumbers);
    count += pack(ptr+count, telescopeD);
    count += pack(ptr+count, minIterations);
    count += pack(ptr+count, maxIterations);
    count += pack(ptr+count, targetIterations);
    count += pack(ptr+count, fillpixMethod);
    count += pack(ptr+count, gradientMethod);
    count += pack(ptr+count, getstepMethod);
    count += pack(ptr+count, badPixelThreshold);
    count += pack(ptr+count, FTOL);
    count += pack(ptr+count, EPS);
    count += pack(ptr+count, reg_alpha);
    count += pack(ptr+count, graddiff_step);
    count += pack(ptr+count, outputFileType);
    count += pack(ptr+count, outputDataType);
    count += pack(ptr+count, sequenceNumber);
    count += pack(ptr+count, observationTime);
    count += pack(ptr+count, observationDate);
    count += pack(ptr+count, tmpDataDir);
    count += pack( ptr+count, (uint16_t)outputFiles.size() );
    for( auto & filename : outputFiles ) {
        count += pack( ptr+count, filename );
    }

    return count;
}


uint64_t GlobalCfg::unpack(const char* ptr, bool swap_endian) {
    using redux::util::unpack;

    uint64_t count = ObjectCfg::unpack(ptr, swap_endian);
    count += unpack(ptr+count, runFlags, swap_endian);
    count += unpack(ptr+count, modeBasis);
    count += unpack(ptr+count, klMinMode, swap_endian);
    count += unpack(ptr+count, klMaxMode, swap_endian);
    count += unpack(ptr+count, klCutoff, swap_endian);
    count += unpack(ptr+count, nInitialModes, swap_endian);
    count += unpack(ptr+count, nModeIncrement, swap_endian);
    count += unpack(ptr+count, modeNumbers, swap_endian);
    count += unpack(ptr+count, telescopeD, swap_endian);
    count += unpack(ptr+count, minIterations, swap_endian);
    count += unpack(ptr+count, maxIterations, swap_endian);
    count += unpack(ptr+count, targetIterations, swap_endian);
    count += unpack(ptr+count, fillpixMethod);
    count += unpack(ptr+count, gradientMethod);
    count += unpack(ptr+count, getstepMethod);
    count += unpack(ptr+count, badPixelThreshold, swap_endian);
    count += unpack(ptr+count, FTOL, swap_endian);
    count += unpack(ptr+count, EPS, swap_endian);
    count += unpack(ptr+count, reg_alpha, swap_endian);
    count += unpack(ptr+count, graddiff_step, swap_endian);
    count += unpack(ptr+count, outputFileType);
    count += unpack(ptr+count, outputDataType);
    count += unpack(ptr+count, sequenceNumber, swap_endian);
    count += unpack(ptr+count, observationTime, swap_endian);
    count += unpack(ptr+count, observationDate, swap_endian);
    count += unpack(ptr+count, tmpDataDir, swap_endian);
    uint16_t tmp;
    count += unpack( ptr+count, tmp, swap_endian );
    outputFiles.resize( tmp );
    for( auto& filename : outputFiles ) {
        count += unpack( ptr+count, filename, swap_endian );
    }

    return count;
}


const GlobalCfg& GlobalCfg::operator=(const ObjectCfg& rhs) {
    ObjectCfg::operator=(rhs);
    return *this;
}


const GlobalCfg& GlobalCfg::operator=(const ChannelCfg& rhs) {
    ChannelCfg::operator=(rhs);
    return *this;
}


bool GlobalCfg::operator==(const GlobalCfg& rhs) const {
    return (runFlags == rhs.runFlags) &&
           (modeBasis == rhs.modeBasis) &&
           (klMinMode == rhs.klMinMode) &&
           (klMaxMode == rhs.klMaxMode) &&
           (klCutoff == rhs.klCutoff) &&
           (nInitialModes == rhs.nInitialModes) &&
           (nModeIncrement == rhs.nModeIncrement) &&
           (telescopeD == rhs.telescopeD) &&
           (minIterations == rhs.minIterations) &&
           (maxIterations == rhs.maxIterations) &&
           (fillpixMethod == rhs.fillpixMethod) &&
           (gradientMethod == rhs.gradientMethod) &&
           (getstepMethod == rhs.getstepMethod) &&
           (badPixelThreshold == rhs.badPixelThreshold) &&
           (FTOL == rhs.FTOL) &&
           (EPS == rhs.EPS) &&
           (reg_alpha == rhs.reg_alpha) &&
           (outputFileType == rhs.outputFileType) &&
           (outputDataType == rhs.outputDataType) &&
           (sequenceNumber == rhs.sequenceNumber) &&
           (observationTime == rhs.observationTime) &&
           (observationDate == rhs.observationDate) &&
           (tmpDataDir == rhs.tmpDataDir) &&
           (modeNumbers == rhs.modeNumbers) &&
           ObjectCfg::operator==(rhs) &&
           (outputFiles == rhs.outputFiles);
}
