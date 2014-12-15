#include "redux/momfbd/channelcfg.hpp"

#include "redux/momfbd/defines.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/objectcfg.hpp"

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
using namespace redux::image;
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
        tmp.erase( boost::remove_if( tmp, boost::is_any_of( "ZzKk" ) ), tmp.end() );
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


ChannelCfg::ChannelCfg( const ObjectCfg& o, const MomfbdJob& j ) : fillpix_method( o.fillpix_method ), mmRow( 0 ), mmWidth( 0 ), incomplete(0), flags( o.flags ),
    image_num_offs( 0 ), sequenceNumber( o.sequenceNumber ), noiseFudge( 0 ), myObject( o ), myJob( j ) {

}

ChannelCfg::~ChannelCfg() {

}

void ChannelCfg::parseProperties( bpt::ptree& tree ) {

    LOG_ERR << "myObject." << printArray(myObject.imageNumbers,"imageNumbers" );
    imageNumbers = tree.get<vector<uint32_t>>( "IMAGE_NUM", myObject.imageNumbers );
    LOG_ERR << "channel." << printArray(myObject.imageNumbers,"imageNumbers" );

    if( imageNumbers.size() == 0 ) LOG_CRITICAL << "Still no image sequence numbers at channel level.";
    sequenceNumber = tree.get<uint32_t>( "SEQUENCE_NUM", myObject.sequenceNumber );
    darkNumbers = tree.get<vector<uint32_t>>( "DARK_NUM", myObject.darkNumbers );

    alignClip = tree.get<vector<int16_t>>( "ALIGN_CLIP", vector<int16_t>() );
    if( !alignClip.empty() && alignClip.size() != 4 ) {
        LOG_ERR << "argument to ALIGN_CLIP could not be translated to 4 integers. Whole image area will be used !!";
        alignClip.clear();
    }
    wf_num = tree.get<vector<uint32_t>>( "WFINDEX", myObject.wf_num );

    imageDataDir = cleanPath( tree.get<string>( "IMAGE_DATA_DIR", myJob.imageDataDir ) );
    if( imageDataDir.length() == 0 ) {
        imageDataDir = cleanPath( "./" );      // Nothing specified in cfg file, use current directory.
    }

    imageTemplate = cleanPath( tree.get<string>( "FILENAME_TEMPLATE", "" ) );
    if( imageTemplate.length() == 0 ) {
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
        tmpString.erase( boost::remove_if( tmpString, boost::is_any_of( "cm\" " ) ), tmpString.end() );
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
            boost::split( tmp, tmpString, boost::is_any_of( "," ) );
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
                    "\"median\", \"invdistweight\" or \"horint\"";
        }
        fillpix_method = myObject.fillpix_method;
    }

    noiseFudge = tree.get<double>( "NF", DEF_NF );
    image_num_offs = tree.get<int>( "DT", 0 );

    flags = myObject.flags;
    MomfbdJob::maybeOverride( tree.get<bool>( "NO_RESTORE", flags & MFBD_NO_RESTORE ), flags, MFBD_NO_RESTORE );
    MomfbdJob::maybeOverride( tree.get<bool>( "SAVE_FFDATA", flags & MFBD_SAVE_FFDATA ), flags, MFBD_SAVE_FFDATA );

    size_t p;
    if( ( p = imageTemplate.find_first_of( '%' ) ) != string::npos ) {
        if( sequenceNumber > 0 ) {
            size_t q;
            if( ( q = imageTemplate.find_first_of( '%', p + 1 ) ) != string::npos ) {
                tmpString = boost::str( boost::format( imageTemplate.substr( 0, q ) ) % sequenceNumber );
                imageTemplate = tmpString + imageTemplate.substr( q );
            }
            else  LOG_WARN << boost::format( "file name template %s does not contain a 2nd format specifier (needs 2)" ) % imageTemplate;
        }
    }
    else {
        LOG_WARN << boost::format( "file name template %s does not contain a format specifier (needs %d)" ) % imageTemplate % ( 1 + ( sequenceNumber >= 0 ) );
    }

    incomplete = tree.get<bool>( "INCOMPLETE", false );

    LOG_DEBUG << "Channel::parseProperties() done.";

}



bpt::ptree ChannelCfg::getPropertyTree( bpt::ptree* root ) {

    bpt::ptree tree;

    if( imageNumbers != myObject.imageNumbers ) tree.put( "IMAGE_NUM", imageNumbers );
    if( sequenceNumber != myObject.sequenceNumber ) tree.put( "SEQUENCE_NUM", sequenceNumber );
    if( darkNumbers != myObject.darkNumbers ) tree.put( "DARK_NUM", darkNumbers );
    if( alignClip.size() == 4 ) tree.put( "ALIGN_CLIP", alignClip );
    if( wf_num != myObject.wf_num ) tree.put( "WFINDEX", wf_num );
    if( imageDataDir != myObject.imageDataDir ) tree.put( "IMAGE_DATA_DIR", imageDataDir );
    if( !imageTemplate.empty() ) tree.put( "FILENAME_TEMPLATE", imageTemplate );
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
    if( noiseFudge != DEF_NF ) tree.put( "NF", noiseFudge );
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


size_t ChannelCfg::size( void ) const {

    size_t sz = 4;                      // fillpix_method, mmRow, mmWidth, incomplete;
    sz += 4 * sizeof( uint32_t );       // flags, image_num_offs, sequenceNumber, imageOffset;
    sz += sizeof( double );             // noiseFudge
    sz += imageNumbers.size() * sizeof( uint32_t ) + sizeof( uint64_t );
    sz += wf_num.size() * sizeof( uint32_t ) + sizeof( uint64_t );
    sz += stokesWeights.size() * sizeof( double ) + sizeof( uint64_t );
    sz += diversity.size() * sizeof( double ) + sizeof( uint64_t );
    sz += diversityOrders.size() * sizeof( uint32_t ) + sizeof( uint64_t );
    sz += diversityTypes.size() * sizeof( uint32_t ) + sizeof( uint64_t );
    sz += darkNumbers.size() * sizeof( uint32_t ) + sizeof( uint64_t );
    sz += alignClip.size() * sizeof( int16_t ) + sizeof( uint64_t );
    sz += imageDataDir.length() + imageTemplate.length() + darkTemplate.length() + gainFile.length() + 4;
    sz += responseFile.length() + backgainFile.length() + psfFile.length() + mmFile.length() + 4;
    sz += offxFile.length() + offyFile.length() + 2;
    sz += dark.size();
    sz += imageStats.size() * Statistics::size() + sizeof(uint32_t);

    return sz;
}


uint64_t ChannelCfg::pack( char* ptr ) const {
    using redux::util::pack;

    uint64_t count = pack( ptr, fillpix_method );
    count += pack( ptr+count, mmRow );
    count += pack( ptr+count, mmWidth );
    count += pack( ptr+count, incomplete );
    count += pack( ptr+count, flags );
    count += pack( ptr+count, image_num_offs );
    count += pack( ptr+count, sequenceNumber );
    count += pack( ptr+count, imageOffset );
    count += pack( ptr+count, noiseFudge );
    count += pack( ptr+count, imageNumbers );
    count += pack( ptr+count, wf_num );
    count += pack( ptr+count, stokesWeights );
    count += pack( ptr+count, diversity );
    count += pack( ptr+count, diversityOrders );
    count += pack( ptr+count, diversityTypes );
    count += pack( ptr+count, darkNumbers );
    count += pack( ptr+count, alignClip );
    count += pack( ptr+count, imageDataDir );
    count += pack( ptr+count, imageTemplate );
    count += pack( ptr+count, darkTemplate );
    count += pack( ptr+count, gainFile );
    count += pack( ptr+count, responseFile );
    count += pack( ptr+count, backgainFile );
    count += pack( ptr+count, psfFile );
    count += pack( ptr+count, mmFile );
    count += pack( ptr+count, offxFile );
    count += pack( ptr+count, offyFile );
    
    count += dark.pack( ptr+count );
    
    uint32_t statSize = imageStats.size();
    count += pack( ptr+count, statSize );
    for( auto &it : imageStats ) count += it->pack(ptr+count);
    

    if(count != size()) cout << "Ch " << hexString(this) << " has a size mismatch: " << count << "  sz = " << size() << "  diff = " << (size()-count) <<endl;
    return count;
}


uint64_t ChannelCfg::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;

    uint64_t count = unpack( ptr, fillpix_method );
    count += unpack( ptr+count, mmRow );
    count += unpack( ptr+count, mmWidth );
    count += unpack( ptr+count, incomplete );
    count += unpack( ptr+count, flags, swap_endian );
    count += unpack( ptr+count, image_num_offs, swap_endian );
    count += unpack( ptr+count, sequenceNumber, swap_endian );
    count += unpack( ptr+count, imageOffset, swap_endian );
    count += unpack( ptr+count, noiseFudge, swap_endian );
    count += unpack( ptr+count, imageNumbers, swap_endian );
    count += unpack( ptr+count, wf_num, swap_endian );
    count += unpack( ptr+count, stokesWeights, swap_endian );
    count += unpack( ptr+count, diversity, swap_endian );
    count += unpack( ptr+count, diversityOrders, swap_endian );
    count += unpack( ptr+count, diversityTypes, swap_endian );
    count += unpack( ptr+count, darkNumbers, swap_endian );
    count += unpack( ptr+count, alignClip, swap_endian );
    count += unpack( ptr+count, imageDataDir );
    count += unpack( ptr+count, imageTemplate );
    count += unpack( ptr+count, darkTemplate );
    count += unpack( ptr+count, gainFile );
    count += unpack( ptr+count, responseFile );
    count += unpack( ptr+count, backgainFile );
    count += unpack( ptr+count, psfFile );
    count += unpack( ptr+count, mmFile );
    count += unpack( ptr+count, offxFile );
    count += unpack( ptr+count, offyFile );

    count += dark.unpack( ptr+count, swap_endian );

    uint32_t statSize;
    count += unpack( ptr+count, statSize, swap_endian);
    imageStats.resize(statSize);
    for( auto &it : imageStats ) {
        it.reset( new Statistics() );
        count += it->unpack(ptr+count, swap_endian);
        //cout << "STATS: " << it->noise << endl; 
    }

    return count;
}


bool ChannelCfg::checkCfg(void) {
    
    LOG_TRACE << "Channel::checkCfg()";

    // Do we have a correct filename template ?
    if( imageTemplate.empty() ) {
        LOG_ERR << "No filename template specified.";
        return false;
    }
    size_t nWild = std::count(imageTemplate.begin(), imageTemplate.end(), '%');
    if( nWild > 2 ) {
        LOG_ERR << "Filename template contains too many wildcards: \"" << imageTemplate << "\"";
        return false;
    } else if( nWild == 1 && imageNumbers.empty() ) {
        LOG_ERR << "Filename template contains wildcard and no image-numbers given (with IMAGE_NUM)";
        return false;
    } else if( nWild == 2 && sequenceNumber == 0 ) {
        LOG_ERR << "Filename template contains 2 wildcards and no sequence-number given (with SEQUENCE_NUM)";
        return false;
    }
    
    // Do we have a correct dark template ?
    if( darkTemplate.empty() ) {
        LOG_ERR << "No filename template specified.";
        return false;
    }
    nWild = std::count(darkTemplate.begin(), darkTemplate.end(), '%');
    if( nWild > 1 ) {
        LOG_ERR << "Dark template contains too many wildcards: \"" << darkTemplate << "\"";
        return false;
    } else if( nWild == 1 && darkNumbers.empty() ) {
        LOG_ERR << "Dark template contains wildcard and no dark-numbers given (with DARK_NUM)";
        return false;
    } else if( nWild == 0 && darkNumbers.size() ) {
        LOG_WARN << "Dark template contains no wildcard AND dark-numbers specified. Numbers will be ignored and the dark-template used as a single filename.";
        darkNumbers.clear();    // TODO: fix this properly, numbers might reappear after transfer (because of inheritance)
    }

    return true;
    
}


bool ChannelCfg::checkData(void) {
    
    LOG_TRACE << "Channel::checkData()";
    
    // Images
    if( incomplete ) {  // check if files are present
        for( size_t i( 0 ); i < imageNumbers.size(); ) {
            bfs::path fn = bfs::path( boost::str( boost::format( imageTemplate ) % (image_num_offs + imageNumbers[i]) ) );
            if( !bfs::exists( fn ) ) {
                fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( imageTemplate ) % (image_num_offs + imageNumbers[i]) ) );
                if( !bfs::exists( fn ) ) {
                    LOG_CRITICAL << "Not found !!! \"" << fn.string() << "\"";
                    imageNumbers.erase( imageNumbers.begin() + i );
                    continue;
                }
            }
            ++i;
        }
        if( imageNumbers.empty() ) {
            LOG_CRITICAL << boost::format( "No files found for incomplete object with filename template \"%s\" in directory \"%s\"" ) % imageTemplate % imageDataDir;
            return false;
        }
    }
    if( imageNumbers.empty() ) {        // single file
        bfs::path fn = bfs::path( imageDataDir ) / bfs::path( imageTemplate );
        if( ! bfs::exists( fn ) ) {
            LOG_ERR << boost::format( "Image-file %s not found!" ) % fn;
            return false;
        }
    } else {                            // template + numbers 
        for( auto & it : imageNumbers ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( imageTemplate ) % (image_num_offs + it) ) );
            if( !bfs::exists( fn ) ) {
                LOG_ERR << boost::format( "Image-file %s not found!" ) % boost::str( boost::format( imageTemplate ) % (image_num_offs + it) );
                return false;
            }
        }
    }
    

    // Dark(s)
    size_t nWild = std::count(darkTemplate.begin(), darkTemplate.end(), '%');
    if( nWild == 0 || darkNumbers.empty() ) {         // single file, DARK_NUM will be ignored if no wildcard in the template
        if( ! bfs::exists( bfs::path( darkTemplate ) ) ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( darkTemplate );
            if( ! bfs::exists( fn ) ) {
               LOG_ERR << boost::format( "Dark-file %s not found!" ) % darkTemplate;
                return false;
            } else darkTemplate = fn.c_str();
        }
    } else {                            // template
        for( auto & it : darkNumbers ) {
            bfs::path fn = bfs::path( boost::str( boost::format( darkTemplate ) % it ) );
            if( !bfs::exists( fn ) ) {
                fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( darkTemplate ) % it ) );
                if( !bfs::exists( fn ) ) {
                    LOG_ERR << boost::format( "Dark-file %s not found!" ) % boost::str( boost::format( darkTemplate ) % it );
                    return false;
                } else darkTemplate = fn.c_str();
            }
        }
    }
    

    // Gain
    if( !gainFile.empty() ) {
        if( ! bfs::exists( bfs::path( gainFile ) ) ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( gainFile );
            if( ! bfs::exists( fn ) ) {
                LOG_ERR << boost::format( "Gain-file %s not found!" ) % gainFile;
                return false;
            } else gainFile = fn.c_str();
        }
    }

    if( !responseFile.empty() ) {
        if( ! bfs::exists( bfs::path( responseFile ) ) ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( responseFile );
            if( ! bfs::exists( fn ) ) {
                LOG_ERR << boost::format( "Response-file %s not found!" ) % responseFile;
                return false;
            } else responseFile = fn.c_str();
        }
    }

    if( !backgainFile.empty() ) {
        if( ! bfs::exists( bfs::path( backgainFile ) ) ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( backgainFile );
            if( ! bfs::exists( fn ) ) {
                LOG_ERR << boost::format( "Backgain-file %s not found!" ) % backgainFile;
                return false;
            } else backgainFile = fn.c_str();
        }
    }

    if( !psfFile.empty() ) {
        if( ! bfs::exists( bfs::path( psfFile ) ) ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( psfFile );
            if( ! bfs::exists( fn ) ) {
                LOG_ERR << boost::format( "PSF-file %s not found!" ) % psfFile;
                return false;
            } else psfFile = fn.c_str();
        }
    }

    if( !mmFile.empty() ) {
        if( ! bfs::exists( bfs::path( mmFile ) ) ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( mmFile );
            if( ! bfs::exists( fn ) ) {
                LOG_ERR << boost::format( "Modulation-matrix file %s not found!" ) % mmFile;
                return false;
            } else mmFile = fn.c_str();
        }
    }

    if( !offxFile.empty() ) {
        if( ! bfs::exists( bfs::path( offxFile ) ) ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( offxFile );
            if( ! bfs::exists( fn ) ) {
                LOG_ERR << boost::format( "Offset-file %s not found!" ) % offxFile;
                return false;
            } else offxFile = fn.c_str();
        }
    }

    if( !offyFile.empty() ) {
        if( ! bfs::exists( bfs::path( offyFile ) ) ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( offyFile );
            if( ! bfs::exists( fn ) ) {
                LOG_ERR << boost::format( "Offset-file %s not found!" ) % offyFile;
                return false;
            } else offyFile = fn.c_str();
        }
    }
    
    return true;
}


namespace {
    template <typename T>
    void loadWrapper(const string& fn, T& img) {
        redux::file::readFile( fn, img );
        LOG_DETAIL << boost::format( "Loaded file \"%s\"" ) % fn;
    }
}


void ChannelCfg::loadData( boost::asio::io_service& service ) {

    LOG_TRACE << "Channel::loadData()";
    // TODO: absolute/relative paths
    // TODO: cache files and just fetch shared_ptr

    if( !darkTemplate.empty() ) {       // needs to be read synchronously because of adding/normalization
        size_t nWild = std::count(darkTemplate.begin(), darkTemplate.end(), '%');
        if( nWild == 0 || darkNumbers.empty() ) {
            LOG_DETAIL << boost::format( "Loading file %s" ) % darkTemplate;
            redux::file::readFile( darkTemplate, dark );
            checkIfMultiFrames(dark);
        }
        else {
            Image<float> tmp;
            for( size_t di = 0; di < darkNumbers.size(); ++di ) {
                bfs::path fn = bfs::path( boost::str( boost::format( darkTemplate ) % darkNumbers[di] ) );
                LOG_DETAIL << boost::format( "Loading file %s" ) % fn;
                if( !di ) {
                    redux::file::readFile( fn.string(), dark );
                    checkIfMultiFrames(dark);
                }
                else {
                    redux::file::readFile( fn.string(), tmp );
                    checkIfMultiFrames(tmp);
                    dark += tmp;
                }
            }
        }
        dark.normalize();
    }


    if( !gainFile.empty() ) {
        service.post( std::bind( loadWrapper< Image<float> >, gainFile, std::ref(gain) ) );
    }


    if( !responseFile.empty() ) {
        service.post( std::bind( loadWrapper< Image<float> >, responseFile, std::ref(ccdResponse) ) );
    }

    if( !backgainFile.empty() ) {
        service.post( std::bind( loadWrapper< Image<float> >, backgainFile, std::ref(ccdScattering) ) );
    }

    if( !psfFile.empty() ) {
        service.post( std::bind( loadWrapper< Image<float> >, psfFile, std::ref(psf) ) );
    }

    if( !mmFile.empty() ) {
        service.post( std::bind( loadWrapper< Image<float> >, mmFile, std::ref(modulationMatrix) ) );
    }

    if( !offxFile.empty() ) {
        service.post( std::bind( loadWrapper< Image<int16_t> >, offxFile, std::ref(xOffset) ) );
    }

    if( !offyFile.empty() ) {
        service.post( std::bind( loadWrapper< Image<int16_t> >, offyFile, std::ref(yOffset) ) );
    }

    size_t nImages = imageNumbers.size();
    if( nImages ) {
        imageStats.resize( nImages );
        Image<float> tmp;
        bfs::path fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( imageTemplate ) % imageNumbers[0] ) );
        redux::file::readFile( fn.string(), tmp );
        images.resize( nImages, tmp.dimSize( 0 ), tmp.dimSize( 1 ) );
        for( size_t i = 0; i < nImages; ++i ) {
            imageStats[i].reset(new Statistics());
            service.post( std::bind( &ChannelCfg::loadImage, this, i ) );
        }
    } else  {
        bfs::path fn = bfs::path( imageDataDir ) / bfs::path( imageTemplate );
        redux::file::readFile( fn.string(), images );
        loadWrapper(fn.string(), images );                              // Single image read synchronously so stats can be done below
        imageStats.resize( 1 );
        imageStats[0].reset(new Statistics());
        imageStats[0]->getStats(myJob.borderClip, images, ST_VALUES);   // only get min/max/mean
    }


}


void ChannelCfg::preprocessData( boost::asio::io_service& service ) {

    size_t nImages = imageNumbers.size();
    double avgMean = 0.0;
    for( size_t i = 0; i < nImages; ++i ) {
        avgMean += imageStats[i]->mean;
    }
    avgMean /= static_cast<double>( nImages );

    for( size_t i = 0; i < nImages; ++i ) {
        service.post( std::bind( &ChannelCfg::preprocessImage, this, i, avgMean ) );
    }

}


double ChannelCfg::getMaxMean(void) const {
    double maxMean = std::numeric_limits<double>::lowest();
    for(auto it: imageStats ) {
        if( it->mean > maxMean ) maxMean = it->mean;
    }
    return maxMean;
}


void ChannelCfg::collectImages(redux::util::Array<float>& stack) const {
    size_t n = images.dimSize(0);
    if ( n ) {
        Array<float> block( stack, imageOffset, imageOffset+n-1, 0, images.dimSize(1)-1, 0, images.dimSize(2)-1 );       // sub-array of stack, starting at offset.
        images.copy(block);
    }
}

            
void ChannelCfg::initWorkSpace( WorkSpace& ws ) {
    
}


void ChannelCfg::normalizeData(boost::asio::io_service& service, double value) {
    size_t nImages = imageNumbers.size();
    for( size_t i = 0; i < nImages; ++i ) {
        service.post( std::bind( &ChannelCfg::normalizeImage, this, i, value ) );
    }
}


void ChannelCfg::loadImage( size_t index ) {
    Image<float> subimg( images, index, index, 0, images.dimSize( 1 ) - 1, 0, images.dimSize( 2 ) - 1 );
    bfs::path fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( imageTemplate ) % imageNumbers[index] ) );
    redux::file::readFile( fn.string(), subimg );
    LOG_DETAIL << boost::format( "Loaded file %s" ) % fn;
    imageStats[index]->getStats(myJob.borderClip,subimg,ST_VALUES);                          // only get min/max/mean
}


void ChannelCfg::preprocessImage( size_t index, double avgMean ) {

    double& mean = imageStats[index]->mean;
    bool modified = false;

    Image<float> subimg( images, index, index, 0, images.dimSize( 1 ) - 1, 0, images.dimSize( 2 ) - 1 );
    bfs::path fn = bfs::path( boost::str( boost::format( imageTemplate ) % imageNumbers[index] ) );
    LOG_DETAIL << boost::format( "Pre-processing image %s" ) % fn;
    
    
    bfs::path fn2 = bfs::path( fn.leaf().string() + ".orig" );
    LOG_DETAIL << boost::format( "Saving RAW file %s" ) % fn2.string();
    redux::file::Ana::write( fn2.string(), subimg );

    
    // Michiel's method for detecting bitshifted Sarnoff images.
    if( mean > 4 * avgMean ) {
        LOG_WARN << boost::format( "Image bit shift detected for image %d (mean > 4*avgMean). adjust factor=0.625 (keep your fingers crossed)!" ) % index;
        subimg *= 0.625;
        modified = true;
    }
    else if( mean < 0.25 * avgMean ) {
        LOG_WARN << boost::format( "Image bit shift detected for image %d (mean < 0.25*avgMean). adjust factor=16 (keep your fingers crossed)!" ) % index;
        subimg *= 16;
        modified = true;
    }
    
    fn2 = bfs::path( fn.leaf().string() + ".1.orig" );
    redux::file::Ana::write( fn2.string(), subimg );

    if( dark.valid() && gain.valid() ) {
        if( ! subimg.sameSize( dark ) ) {
            LOG_ERR << boost::format( "Dimensions of dark (%s) does not match this image (%s), skipping flatfielding !!" ) % printArray( dark.dimensions(), "" ) % printArray( subimg.dimensions(), "" );
            return;
        }
        if( ! subimg.sameSize( gain ) ) {
            LOG_ERR << boost::format( "Dimensions of gain (%s) does not match this image (%s), skipping flatfielding !!" ) % printArray( gain.dimensions(), "" ) % printArray( subimg.dimensions(), "" );
            return;
        }
        if( ccdResponse.valid() && !subimg.sameSize( ccdResponse ) ) {
            LOG_WARN << boost::format( "Dimensions of ccd-response (%s) does not match this image (%s), will not be used !!" ) % printArray( ccdResponse.dimensions(), "" ) % printArray( subimg.dimensions(), "" );
            ccdResponse.resize();
        }

        subimg -= dark;
        modified = true;
    fn2 = bfs::path( fn.leaf().string() + ".2.orig" );
    redux::file::Ana::write( fn2.string(), subimg );

        if( ccdResponse.valid() ) {  // correct for the detector response (this should not contain the gain correction and must be done before descattering)
            subimg *= ccdResponse;
        }
    fn2 = bfs::path( fn.leaf().string() + ".3.orig" );
    redux::file::Ana::write( fn2.string(), subimg );

        if( ccdScattering.valid() && psf.valid() ) {          // apply backscatter correction
            if( subimg.sameSize( ccdScattering ) && subimg.sameSize( psf ) ) {
                LOG_TRACE << "Applying correction for CCD transparency.";
                redux::image::descatter( subimg, ccdScattering, psf);
            }
            else {
                LOG_ERR << boost::format( "Dimensions of ccdScattering (%s) or psf (%s) does not match this image (%s), skipping flatfielding !!" )
                        % printArray( ccdScattering.dimensions(), "" ) % printArray( psf.dimensions(), "" ) % printArray( subimg.dimensions(), "" );
            }
        }

    fn2 = bfs::path( fn.leaf().string() + ".4.orig" );
    redux::file::Ana::write( fn2.string(), subimg );
        subimg *= gain;
        namespace sp = std::placeholders;
        size_t ni = images.dimSize(0);
        size_t sy = images.dimSize(1);
        size_t sx = images.dimSize(2);
        shared_ptr<float**> array = images.get(ni,sy,sx);
        float*** arrayPtr = array.get();
        float badPixelThreshold = 1E-5;     // TODO: configurable threshold.
        switch(fillpix_method) {
            case CFG_FPM_HORINT: {
                LOG_TRACE << "Filling bad pixels using horizontal interpolation.";
                function<float(size_t,size_t)> func = bind(horizontalInterpolation<float>, arrayPtr[index], sy, sx, sp::_1, sp::_2);
                fillPixels(arrayPtr[index], sy, sx, func, std::bind2nd(std::less_equal<float>(),badPixelThreshold) );
                break;
            }
            case CFG_FPM_MEDIAN: {
                // TODO: median method
                break;
            }
            case CFG_FPM_INVDISTWEIGHT:       // inverse distance weighting is the default method, so fall through
            default: {
                function<double(size_t,size_t)> func = bind(inverseDistanceWeight<float>, arrayPtr[index], sy, sx, sp::_1, sp::_2);
                //function<double(size_t,size_t)> func = bind(inv_dist_wght, arrayPtr[index], sy, sx, sp::_1, sp::_2);
                fillPixels(arrayPtr[index], sy, sx, func, std::bind2nd(std::less_equal<float>(),badPixelThreshold) );
            }
        }

    fn2 = bfs::path( fn.leaf().string() + ".5.orig" );
    redux::file::Ana::write( fn2.string(), subimg );
    }

    imageStats[index]->getStats(myJob.borderClip, subimg);      // get all stats for the cleaned up data

    if( modified && flags & MFBD_SAVE_FFDATA ) {
        //string fnT = fn.leaf().string();
        fn = bfs::path( fn.leaf().string() + ".cor" );

        LOG_DETAIL << boost::format( "Saving flat/dark corrected file %s" ) % fn.string();
        redux::file::Ana::write( fn.string(), subimg );
        
      /*  
        //FourierTransform bla(subimg);
        Image<float> apod = subimg.copy();
        apod.trim();
        cout << printArray(apod.dimensions(),"blaDims") << endl;
        apodize(apod,45);
        //subimg *= apod;
        Image<double> psf(subimg.dimensions(true));
        size_t sy = psf.dimSize(0);
        size_t sx = psf.dimSize(1);
        cout << printArray(psf.dimensions(),"psfDims") << endl;
        redux::math::gauss(psf.ptr(),sy,sx,5,5,sy/2,sx/2);
        //FourierTransform::reorder(psf);
        //redux::file::Ana::write( fnT+".rpsf", psf );
        FourierTransform psft(psf, FT_REORDER|FT_NORMALIZE);
        //blat.reorder();
        auto conv = psft.convolve(subimg);
        cout << printArray(conv.dimensions(),"convDims") << endl;
        Image<float> power(psft.dimensions());
        auto rit = psft.begin();
        for(auto& it: power) {
           //it = rit++->real();
           it = std::abs(*rit++);
        }
        //bla.copy(bla2);
        redux::file::Ana::write( fnT+".psf", psf );
        redux::file::Ana::write( fnT+".pwr", power );
        redux::file::Ana::write( fnT+".apod", apod );
        redux::file::Ana::write( fnT+".conv", conv );
        */
        
    }



}


void ChannelCfg::normalizeImage( size_t index, double value ) {
    LOG_TRACE << boost::format( "Normalizing image %d" ) % (imageOffset+index);
    Image<float> subimg( images, index, index, 0, images.dimSize( 1 ) - 1, 0, images.dimSize( 2 ) - 1 );
    subimg *= (value/imageStats[index]->mean);
    double noise1 = imageStats[index]->noise;
    imageStats[index]->getStats(subimg);
    LOG_TRACE << "  image #" << (imageOffset+index) << "  noise1 = " << noise1 << "  noise2 = " << imageStats[index]->noise;
}


size_t ChannelCfg::sizeOfPatch(uint32_t npixels) const {
    size_t sz = sizeof(size_t) + imageStats.size()*sizeof(float);
    sz += npixels*images.dimSize(0)*sizeof(float);
    return sz;
}


void ChannelCfg::applyLocalOffsets(PatchData::Ptr patch) const {
   
    //LOG_TRACE << "ChannelCfg::applyLocalOffsets #" << patch->id;
    int patchOffsetX, patchOffsetY;
    float residualOffsetX, residualOffsetY;
    residualOffsetX = residualOffsetY = patchOffsetX = patchOffsetY = 0;
    
    if(xOffset.valid()) {
        Image<int16_t> tmpOff(xOffset, patch->first.y, patch->last.y, patch->first.x, patch->last.x);
        Statistics stats;
        stats.getStats(tmpOff,ST_VALUES);
        double tmp;
        residualOffsetX = modf(stats.mean/100 , &tmp);
        patchOffsetX = lrint(tmp);
    }
    
    if(yOffset.valid()) {
        Image<int16_t> tmpOff(xOffset, patch->first.y, patch->last.y, patch->first.x, patch->last.x);
        Statistics stats;
        stats.getStats(Image<int16_t>(yOffset, patch->first.y, patch->last.y, patch->first.x, patch->last.x),ST_VALUES);
        double tmp;
        residualOffsetY = modf(stats.mean/100 , &tmp);
        patchOffsetY = lrint(tmp);
    }
   
    if( patchOffsetX ) {
        int shift = patch->images.shift(2,patchOffsetX);
        if( shift != patchOffsetX) {
            residualOffsetX += (patchOffsetX-shift);
        }
    }

    if( patchOffsetY ) {
        int shift = patch->images.shift(1,patchOffsetY);
        if( shift != patchOffsetY) {
            residualOffsetY += (patchOffsetY-shift);
        }
    }
    
    patch->residualOffset.y = residualOffsetY;
    patch->residualOffset.x = residualOffsetX;
   // LOG_TRACE << "ChannelCfg::applyLocalOffsets #" << patch->id << "  ...end";

}
 

Point ChannelCfg::clipImages(void) {

    if( alignClip.size() != 4 ) {
        LOG_WARN << "No (or faulty) align-clip supplied, using full images: " << printArray(alignClip,"clip"); 
        return Point(images.dimSize(1),images.dimSize(2));
    }

    LOG_DETAIL << "Clipping images using " << printArray(alignClip,"alignClip");
    bool flipX=false,flipY=false;
    if( alignClip[0] > alignClip[1]) {     // we have the y (row/slow) dimension first, momfbd cfg-files (and thus alignClip) has x first. TBD: how should this be ?
        std::swap(alignClip[0], alignClip[1]);
        flipX = true;
        LOG_DETAIL << "Reversing x-coordinate for this channel.";
    }
    if( alignClip[2] > alignClip[3]) {
        std::swap(alignClip[2], alignClip[3]);
        flipY = true;
        LOG_DETAIL << "Reversing y-coordinate for this channel.";
    }
    for( auto& it: alignClip ) --it;        // NOTE: momfbd uses 1-based indexes, we start with 0...
    
    string tmp =  printArray(images.dimensions(),"original");
    images.setLimits( 0, imageNumbers.size()-1, alignClip[2], alignClip[3], alignClip[0], alignClip[1] );
    images.trim(false);
    LOG_DEBUG << "          image stack: " << tmp << printArray(images.dimensions(),"  clipped");
    
    if( dark.valid() ) {
        LOG_DEBUG << "                 dark: " << printArray(dark.dimensions(),"original");
        dark.setLimits( alignClip[2], alignClip[3], alignClip[0], alignClip[1] );
        dark.trim();
    }

    if( gain.valid() ) {
        LOG_DEBUG << "                 gain: " << printArray(gain.dimensions(),"original");
        gain.setLimits( alignClip[2], alignClip[3], alignClip[0], alignClip[1] );
        gain.trim();
    }

    if( ccdResponse.valid() ) {
        LOG_DEBUG << "          ccdResponse: " << printArray(ccdResponse.dimensions(),"original");
        ccdResponse.setLimits( alignClip[2], alignClip[3], alignClip[0], alignClip[1] );
        dark.trim();
    }

    if( ccdScattering.valid() ) {
        LOG_DEBUG << "        ccdScattering: " << printArray(ccdScattering.dimensions(),"original");
        ccdScattering.setLimits( alignClip[2], alignClip[3], alignClip[0], alignClip[1] );
        ccdScattering.trim();
    }

    if( psf.valid() ) {     // The PSF has to be symmetrically clipped, otherwise the convolution will be skewed.
        const std::vector<size_t>& dims = psf.dimensions();
        auto tmp = alignClip;
        size_t sy = alignClip[3] - alignClip[2] + 1;
        size_t sx = alignClip[1] - alignClip[0] + 1;
        if (dims.size() != 2 || dims[0] < sy || dims[1] < sx) throw std::logic_error("PSF has wrong dimensions: " + printArray(dims,"dims"));
        int skewY = (dims[0] - sy)/2  - alignClip[2];
        int skewX = (dims[1] - sx)/2  - alignClip[0];
        tmp[0] += skewX;
        tmp[1] += skewX;
        tmp[2] += skewY;
        tmp[3] += skewY;
        LOG_DEBUG << "                  psf: " << printArray(dims,"original") << printArray(tmp,"  symmetric clip");
        psf.setLimits( tmp[2], tmp[3], tmp[0], tmp[1] );
        psf.trim();
    }

    if( flipX || flipY ) {
        size_t sy = alignClip[3] - alignClip[2] + 1;
        size_t sx = alignClip[1] - alignClip[0] + 1;
        size_t ni = imageNumbers.size();
        shared_ptr<float**> imgShared = images.get(ni,sy,sx);
        float*** imgPtr = imgShared.get();
        for(size_t i=0; i<ni; ++i) {
            if (flipX) reverseX(imgPtr[i],sy,sx);
            if (flipY) reverseY(imgPtr[i],sy,sx);
        }
        shared_ptr<float*> tmpShared;
        if( dark.valid() ) {
            tmpShared = dark.get(sy,sx);
            if (flipX) reverseX(tmpShared.get(),sy,sx);
            if (flipY) reverseY(tmpShared.get(),sy,sx);
        }

        if( gain.valid() ) {
            tmpShared = gain.get(sy,sx);
            if (flipX) reverseX(tmpShared.get(),sy,sx);
            if (flipY) reverseY(tmpShared.get(),sy,sx);
        }

        if( ccdResponse.valid() ) {
            tmpShared = ccdResponse.get(sy,sx);
            if (flipX) reverseX(tmpShared.get(),sy,sx);
            if (flipY) reverseY(tmpShared.get(),sy,sx);
        }

        if( ccdScattering.valid() ) {
            tmpShared = ccdScattering.get(sy,sx);
            if (flipX) reverseX(tmpShared.get(),sy,sx);
            if (flipY) reverseY(tmpShared.get(),sy,sx);
        }

        if( psf.valid() ) {
            tmpShared = psf.get(sy,sx);
            if (flipX) reverseX(tmpShared.get(),sy,sx);
            if (flipY) reverseY(tmpShared.get(),sy,sx);
        }
    }
    
    return Point(images.dimSize(1),images.dimSize(2));
    
}


