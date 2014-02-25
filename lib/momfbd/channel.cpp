#include "redux/momfbd/channel.hpp"

#include "redux/momfbd/defines.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/object.hpp"

#include "redux/constants.hpp"
#include "redux/logger.hpp"
#include "redux/translators.hpp"
#include "redux/util/stringutil.hpp"

#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;
using namespace redux::momfbd;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

    const string thisChannel = "momfbdch";

    void parseSegment( vector<uint32_t>& divs, vector<uint32_t>& types, const string& elem ) {
        size_t n = std::count( elem.begin(), elem.end(), '-' );
        uint32_t tp = CFG_DEFAULT;
        if( elem.find_first_of( "Zz" ) != string::npos ) tp = CFG_ZERNIKE;
        if( elem.find_first_of( "Kk" ) != string::npos ) {
            if( tp == CFG_DEFAULT ) {
                tp = CFG_KARHUNEN_LOEVE;
            }
            else {
                LOG_CRITICAL << "error: different types in specified mode range \"" << elem << "\"";
            }
        }
        string tmp = elem;
        tmp.erase(boost::remove_if(tmp, boost::is_any_of("ZzKk")), tmp.end());
        if( n == 0 ) {
            divs.push_back( boost::lexical_cast<uint32_t>( tmp ) );
            types.push_back( tp );
            return;
        }
        else if( n == 1 ) {
            n = tmp.find_first_of( '-' );
            uint32_t first = boost::lexical_cast<uint32_t>( tmp.substr( 0, n ) );
            uint32_t last = boost::lexical_cast<uint32_t>( tmp.substr( n + 1 ) );
            while( first <= last ) {
                divs.push_back( first++ );
                types.push_back( tp );
            }
        }
    }

    double def2cf( double pd_defocus, double telescope_r ) { // defocus distance in meters
        static const double tmp = ( 8.0 * sqrt( 3.0 ) );
        return -pd_defocus * redux::PI * telescope_r * telescope_r * tmp;
    }

}


Channel::Channel( Object& o, MomfbdJob& j ) : flags(o.flags), mmRow(0), mmWidth(0), fillpix_method(o.fillpix_method),
                                              image_num_offs(0), sequenceNumber(o.sequenceNumber), nf(0), myObject( o ), myJob( j ) {
    
}

Channel::~Channel() {

}

void Channel::parseProperties( bpt::ptree& tree ) {

    LOG_DEBUG << "Channel::parseProperties()";

    imageNumbers = tree.get<vector<uint32_t>>( "IMAGE_NUM", myObject.imageNumbers );

    if( imageNumbers.size() == 0 ) LOG_CRITICAL << "Still no image sequence numbers at channel level.";
    sequenceNumber = tree.get<uint32_t>( "SEQUENCE_NUM", myObject.sequenceNumber );
    darkNumbers = tree.get<vector<uint32_t>>( "DARK_NUM", myObject.darkNumbers );

    if( tree.get<bool>( "ALIGN_CLIP", false ) ) {
        alignClip = tree.get<vector<int16_t>>( "ALIGN_CLIP", vector<int16_t>() );
        if( alignClip.size() != 4 ) LOG_ERR << "argument to ALIGN_CLIP could not be translated to 4 numbers.";
    }

    wf_num = tree.get<vector<uint32_t>>( "WFINDEX", myObject.wf_num );

    imageDataDir = cleanPath( tree.get<string>( "IMAGE_DATA_DIR", myJob.imageDataDir ) );
    if( imageDataDir.length() == 0 ) {
        imageDataDir = cleanPath( "./" );      // Nothing specified in cfg file, use current directory.
    }

    filenameTemplate = cleanPath( tree.get<string>( "FILENAME_TEMPLATE", "" ) );
    if( filenameTemplate.length() == 0 ) {
        LOG_ERR << "no filename template specified.";
    }

    darkTemplate = tree.get<string>( "DARK_TEMPLATE", "" );
    if( darkTemplate.length() > 0 ) {
        if( darkNumbers.size() == 0 ) {
            LOG_ERR << "darkfield template specified but no dark numbers.";
        }
    }
    else if( darkNumbers.size() > 0 ) {
        LOG_ERR << "darkfield dark numbers specified but no darkfield template.";
    }

    gainFile = tree.get<string>( "GAIN_FILE", "" );
    if( gainFile.length() > 0 ) {
        if( darkTemplate.length() == 0 ) {
            LOG_ERR << "a gain file name but no dark field was specified.";
        }
    }
    else if( darkTemplate.length() > 0 ) {
        LOG_ERR << "a dark field name but no gain file was specified.";
    }

    responseFile = tree.get<string>( "CCD_RESPONSE", "" );
    if( ( responseFile.length() > 0 ) && ( gainFile.length() == 0 ) ) {
        LOG_ERR << "detector response correction only possible when flatfielding.";
    }

    backgainFile = cleanPath( tree.get<string>( "BACK_GAIN", "" ), imageDataDir );
    psfFile = cleanPath( tree.get<string>( "PSF", "" ), imageDataDir );
    mmFile = cleanPath( tree.get<string>( "MODMAT", "" ), imageDataDir );

    if( mmFile.length() > 0 ) {
        mmRow = tree.get<uint8_t>( "MMROW", 0 );
        if( tree.get<bool>( "MMROW", false ) ) {
            LOG_CRITICAL << "a modulation matrix was provided but no row specified (MMROW).";
        }
        mmWidth = tree.get<uint8_t>( "MMWIDTH", 0 );
        if( tree.get<bool>( "MMWIDTH", false ) ) {
            LOG_CRITICAL << "modulation matrix dimensions cannot be autodetected (yet): you must provide the matrix width (MMWIDTH)!";
        }
        stokesWeights = tree.get<vector<double>>( "VECTOR", myObject.stokesWeights );
        if( stokesWeights.size() > 0 ) {
            if( stokesWeights.size() != mmWidth ) {
                LOG_ERR << "VECTOR input has " << stokesWeights.size() << " elements, but MMWIDTH=" << mmWidth;
            }
        }
        else {
            stokesWeights = myObject.stokesWeights;
            if( stokesWeights.size() == 0 ) {
                LOG_ERR << "modulation matrix specified but no VECTOR input given!";
            }
        }
    }
    else {
        mmRow = mmWidth = 1;
        stokesWeights.resize( 1, 1.0 );
    }

    offxFile = tree.get<string>( "XOFFSET", "" );
    offyFile = tree.get<string>( "YOFFSET", "" );

    string tmpString = tree.get<string>( "DIVERSITY", "" );
    if( tmpString.size() == 0 ) {
        LOG_WARN << "no diversity specified (assuming zero).";
        diversity.resize( 1, 0.0 );
        diversityOrders.resize( 1, 4 );
        diversityTypes.resize( 1, CFG_ZERNIKE );
    }
    else {
        double tmpD;
        if( tmpString.find( "mm" ) != string::npos ) tmpD = 1.00E-03;
        else if( tmpString.find( "cm" ) != string::npos ) tmpD = 1.00E-02;
        else tmpD = 1.0;
        tmpString.erase(boost::remove_if(tmpString, boost::is_any_of("cm\" ")), tmpString.end());
        bpt::ptree tmpTree;                         // just to be able to use the VectorTranslator
        tmpTree.put( "tmp", tmpString );
        diversity = tmpTree.get<vector<double>>( "tmp", vector<double>() );
        tmpString = tree.get<string>( "DIV_ORDERS", "" );
        if( tmpString.empty() ) {
            if( diversity.size() > 1 ) {
                LOG_ERR << "multiple coefficients found but no diversity orders specified!";
            }
            else {
                diversityOrders.resize( 1, 4 );
                diversityTypes.resize( 1, CFG_ZERNIKE );
                diversity[1] = def2cf( tmpD * diversity[1], myJob.telescopeDiameter / myJob.telescopeFocalLength );
            }
        }
        else {
            vector<string> tmp;
            boost::split( tmp, tmpString, boost::is_any_of(",") );
            for( auto & it : tmp ) {
                parseSegment( diversityOrders, diversityTypes, it );
            }
            if( diversity.size() != diversityOrders.size() ) {
                LOG_ERR << "number of diversity orders does not match number of diversity coefficients!";
            }
        }
    }

    fillpix_method = MomfbdJob::getFromMap( tmpString = tree.get<string>( "FPMETHOD", "" ), fillpixMap );
    if( fillpix_method == 0 ) {
        if( tmpString.length() > 0 ) {
            LOG_ERR << "unknown fillpix method \"" << tmpString << "\"\n  Valid entries currently are: "
                    "\"median\", \"invdistweight\" or \"horizontal interpolation\"";
        }
        fillpix_method = myObject.fillpix_method;
    }

    nf = tree.get<double>( "NF", DEF_NF );
    image_num_offs = tree.get<int>( "DT", 0 );

    flags = myObject.flags;
    MomfbdJob::maybeOverride( tree.get<bool>( "NO_RESTORE", flags & MFBD_NO_RESTORE ), flags, MFBD_NO_RESTORE );
    MomfbdJob::maybeOverride( tree.get<bool>( "SAVE_FFDATA", flags & MFBD_SAVE_FFDATA ), flags, MFBD_SAVE_FFDATA );

    size_t p;
    if( ( p = filenameTemplate.find_first_of( '%' ) ) != string::npos ) {
        if( sequenceNumber > 0 ) {
            size_t q;
            if( ( q = filenameTemplate.find_first_of( '%', p + 1 ) ) != string::npos ) {
                tmpString = boost::str( boost::format( filenameTemplate.substr( 0, q ) ) % sequenceNumber );
                filenameTemplate = tmpString + filenameTemplate.substr( q );
            }
            else  LOG_WARN << boost::format( "file name template \"%s\" does not contain a 2nd format specifier (needs 2)" ) % filenameTemplate;
        }
    }
    else {
        LOG_WARN << boost::format( "file name template \"%s\" does not contain a format specifier (needs %d)" ) % filenameTemplate % ( 1 + ( sequenceNumber >= 0 ) );
    }

    if( tree.get<bool>( "INCOMPLETE", false ) ) {
        for( size_t i( 0 ); i < imageNumbers.size(); ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( filenameTemplate ) % ( image_num_offs + imageNumbers[i] ) ) );
            if( !bfs::exists( fn ) ) {
                imageNumbers.erase( imageNumbers.begin() + i );
                continue;
            }
            ++i;
        }
        if( imageNumbers.empty() ) {
            LOG_CRITICAL << boost::format( "no files found for incomplete object with filename template \"%s\" in directory \"%s\"" ) % filenameTemplate % imageDataDir;
        }
    }
    else {
        for( size_t i( 0 ); i < imageNumbers.size(); ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( filenameTemplate ) % ( image_num_offs + imageNumbers[i] ) ) );
            if( !bfs::exists( fn ) ) {
                LOG_CRITICAL << boost::format( "file %s not found!" ) % fn;
            }
            ++i;
        }

    }


    LOG_DEBUG << "Channel::parseProperties() done.";
    
}



bpt::ptree Channel::getPropertyTree( bpt::ptree* root ) {

    bpt::ptree tree;

    if( imageNumbers != myObject.imageNumbers ) tree.put( "IMAGE_NUM", imageNumbers );
    if( sequenceNumber != myObject.sequenceNumber ) tree.put( "SEQUENCE_NUM", sequenceNumber );
    if( darkNumbers != myObject.darkNumbers ) tree.put( "DARK_NUM", darkNumbers );
    if( alignClip.size() == 4 ) tree.put( "ALIGN_CLIP", alignClip );
    if( wf_num != myObject.wf_num ) tree.put( "WFINDEX", wf_num );
    if( imageDataDir != myObject.imageDataDir ) tree.put( "IMAGE_DATA_DIR", imageDataDir );
    if( !filenameTemplate.empty() ) tree.put( "FILENAME_TEMPLATE", filenameTemplate );
    if( !darkTemplate.empty() ) tree.put( "DARK_TEMPLATE", darkTemplate );
    if( !gainFile.empty() ) tree.put( "GAIN_FILE", gainFile );
    if( !responseFile.empty() ) tree.put( "CCD_RESPONSE", responseFile );
    if( !backgainFile.empty() ) tree.put( "BACK_GAIN", backgainFile );
    if( !psfFile.empty() ) tree.put( "PSF", psfFile );
    if( !mmFile.empty() ) tree.put( "MODMAT", mmFile );
    if( mmRow ) tree.put( "MMROW", mmRow );
    if( mmWidth ) tree.put( "MMWIDTH", mmWidth );
    if( stokesWeights != myObject.stokesWeights ) tree.put( "VECTOR", stokesWeights );
    if( !offxFile.empty() ) tree.put( "XOFFSET", offxFile );
    if( !offyFile.empty() ) tree.put( "YOFFSET", offyFile );
    if( !diversity.empty() ) tree.put( "DIVERSITY", diversity );
    if( !diversityOrders.empty() ) tree.put( "DIV_ORDERS", diversityOrders ); // TODO types missing
    if( fillpix_method != myObject.fillpix_method ) tree.put( "FPMETHOD", fillpix_method );
    if( nf != DEF_NF ) tree.put( "NF", nf );
    if( image_num_offs != 0 ) tree.put( "DT", image_num_offs );
    uint32_t dflags = flags ^ myObject.flags;
    if( dflags & MFBD_NO_RESTORE ) tree.put( "NO_RESTORE", ( bool )( flags & MFBD_NO_RESTORE ) );
    if( dflags & MFBD_SAVE_FFDATA ) tree.put( "SAVE_FFDATA", ( bool )( flags & MFBD_SAVE_FFDATA ) );
    // TODO "INCOMPLETE"

    if( root ) {
        root->push_back( bpt::ptree::value_type( "channel", tree ) );
    }
    
    return tree;

}


size_t Channel::size(void) const {
    
    size_t sz = 3;
    sz += 3*sizeof(uint32_t);
    sz += sizeof(double);
    sz += imageNumbers.size()*sizeof(uint32_t)+sizeof(size_t);
    sz += darkNumbers.size()*sizeof(uint32_t)+sizeof(size_t);
    sz += alignClip.size()*sizeof(int16_t)+sizeof(size_t);
    sz += wf_num.size()*sizeof(uint32_t)+sizeof(size_t);
    sz += stokesWeights.size()*sizeof(double)+sizeof(size_t);
    sz += diversity.size()*sizeof(double)+sizeof(size_t);
    sz += diversityOrders.size()*sizeof(uint32_t)+sizeof(size_t);
    sz += diversityTypes.size()*sizeof(uint32_t)+sizeof(size_t);
    sz += imageDataDir.length() + filenameTemplate.length() + darkTemplate.length() + gainFile.length() + 4;
    sz += responseFile.length() + backgainFile.length() + psfFile.length() + mmFile.length() + 4;
    sz += offxFile.length() + offyFile.length() + 2;

    return sz;
}


char* Channel::pack(char* ptr) const {
    using redux::util::pack;

    ptr = pack(ptr,fillpix_method);
    ptr = pack(ptr,mmRow);
    ptr = pack(ptr,mmWidth);
    ptr = pack(ptr,sequenceNumber);
    ptr = pack(ptr,image_num_offs);
    ptr = pack(ptr,flags);
    ptr = pack(ptr,nf);
    ptr = pack(ptr,imageNumbers);
    ptr = pack(ptr,darkNumbers);
    ptr = pack(ptr,alignClip);
    ptr = pack(ptr,wf_num);
    ptr = pack(ptr,stokesWeights);
    ptr = pack(ptr,diversity);
    ptr = pack(ptr,diversityOrders);
    ptr = pack(ptr,diversityTypes);
    ptr = pack(ptr,imageDataDir);
    ptr = pack(ptr,filenameTemplate);
    ptr = pack(ptr,darkTemplate);
    ptr = pack(ptr,gainFile);
    ptr = pack(ptr,responseFile);
    ptr = pack(ptr,backgainFile);
    ptr = pack(ptr,psfFile);
    ptr = pack(ptr,mmFile);
    ptr = pack(ptr,offxFile);
    ptr = pack(ptr,offyFile);

    return ptr;
}


const char* Channel::unpack(const char* ptr, bool swap_endian) {
    using redux::util::unpack;

    ptr = unpack(ptr,fillpix_method);
    ptr = unpack(ptr,mmRow);
    ptr = unpack(ptr,mmWidth);
    ptr = unpack(ptr,sequenceNumber, swap_endian);
    ptr = unpack(ptr,image_num_offs, swap_endian);
    ptr = unpack(ptr,flags, swap_endian);
    ptr = unpack(ptr,nf, swap_endian);
    ptr = unpack(ptr,imageNumbers, swap_endian);
    ptr = unpack(ptr,darkNumbers, swap_endian);
    ptr = unpack(ptr,alignClip, swap_endian);
    ptr = unpack(ptr,wf_num, swap_endian);
    ptr = unpack(ptr,stokesWeights, swap_endian);
    ptr = unpack(ptr,diversity, swap_endian);
    ptr = unpack(ptr,diversityOrders, swap_endian);
    ptr = unpack(ptr,diversityTypes, swap_endian);
    ptr = unpack(ptr,imageDataDir);
    ptr = unpack(ptr,filenameTemplate);
    ptr = unpack(ptr,darkTemplate);
    ptr = unpack(ptr,gainFile);
    ptr = unpack(ptr,responseFile);
    ptr = unpack(ptr,backgainFile);
    ptr = unpack(ptr,psfFile);
    ptr = unpack(ptr,mmFile);
    ptr = unpack(ptr,offxFile);
    ptr = unpack(ptr,offyFile);
    
    return ptr;
}
