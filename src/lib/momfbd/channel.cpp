#include "redux/momfbd/channel.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/object.hpp"
#include "redux/momfbd/subimage.hpp"
#include "redux/momfbd/wavefront.hpp"

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

    const string thisChannel = "channel";
        
    bool checkImageScale( float& F, float& A, float& P ) {
        
        double rad2asec = 180.0 * 3600.0 / redux::PI;
        size_t count = F > 0 ? 1 : 0;
        count += A > 0 ? 1 : 0;
        count += P > 0 ? 1 : 0;
        if( count > 2 ) {
            LOG_WARN << "Too many parameters specified: replacing telescope focal length (" << F
                     << ") with computed value (" << (P * rad2asec / A) << ")";
            F = P * rad2asec / A;
            return true;
        } else if ( count < 2 ) {
            LOG_ERR << "At least two of the parameters \"TELESCOPE_F\", \"ARCSECPERPIX\" and \"PIXELSIZE\" has to be provided.";
        } else {    // count == 2
            if( F <= 0 ) {
                F = P * rad2asec / A;
            } else if ( A <= 0 ) {
                A = P * rad2asec / F;
            } else if ( P <= 0 ) {
                P = F * A / rad2asec;
            }
            return true;
        }
        return false;
    }

    void calculatePupilSize( double &frequencyCutoff, double &pupilRadiusInPixels, uint16_t &nPupilPixels, double wavelength, uint32_t nPixels, double telescopeDiameter, double arcSecsPerPixel ) {
        static double radians_per_arcsec = redux::PI/(180.0*3600.0);             // (2.0*redux::PI)/(360.0*3600.0)
        double radians_per_pixel = arcSecsPerPixel * radians_per_arcsec;
        double q_number = wavelength / ( radians_per_pixel * telescopeDiameter );
        frequencyCutoff = ( double )nPixels / q_number;
        nPupilPixels = nPixels>>2;
        pupilRadiusInPixels = frequencyCutoff / 2.0;                   // telescope radius in pupil pixels...
        if( nPupilPixels < pupilRadiusInPixels ) {           // this should only be needed for oversampled images
            uint16_t goodsizes[] = { 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45, 48, 50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128, 135, 144 };
            for( int i = 0; ( nPupilPixels = max( goodsizes[i], nPupilPixels ) ) < pupilRadiusInPixels; ++i ); // find right size
        }
        nPupilPixels <<= 1;
    }


}


Channel::Channel( Object& o, MomfbdJob& j, uint16_t id ) : ID( id ), myObject( o ), myJob( j ) {

}

Channel::~Channel() {

}

void Channel::parsePropertyTree( bpt::ptree& tree ) {
    
    ChannelCfg::parseProperties(tree, myObject);

    /*if( !alignClip.empty() && alignClip.size() != 4 ) {
        LOG_ERR << "argument to ALIGN_CLIP could not be translated to 4 integers. Whole image area will be used !!";
        alignClip.clear();
    }*/

    /*if( imageDataDir.length() == 0 ) {
        imageDataDir = cleanPath( "./" );      // Nothing specified in cfg file, use current directory.
    }*/

    if( imageTemplate.length() == 0 ) {
        LOG_ERR << "no filename template specified.";
    }

/*    if( darkTemplate.length() > 0 ) {
        if( darkNumbers.size() == 0 ) {
            LOG_ERR << "darkfield template specified but no dark numbers.";
        }
    }
    else if( darkNumbers.size() > 0 ) {
        LOG_ERR << "darkfield dark numbers specified but no darkfield template.";
    }
*/
    if( gainFile.length() > 0 ) {
        if( darkTemplate.length() == 0 ) {
            LOG_ERR << "a gain file name but no dark field was specified.";
        }
    }
    else if( darkTemplate.length() > 0 ) {
        LOG_ERR << "a dark field name but no gain file was specified.";
    }

    if( ( responseFile.length() > 0 ) && ( gainFile.length() == 0 ) ) {
        LOG_ERR << "detector response correction only possible when flatfielding.";
    }



    //MomfbdJob::maybeOverride( tree.get<bool>( "NO_RESTORE", flags & MFBD_NO_RESTORE ), flags, MFBD_NO_RESTORE );

    size_t p;
    if( ( p = imageTemplate.find_first_of( '%' ) ) != string::npos ) {
   /*     if( sequenceNumber > 0 ) {
            size_t q;
            if( ( q = imageTemplate.find_first_of( '%', p + 1 ) ) != string::npos ) {
                string tmpString = boost::str( boost::format( imageTemplate.substr( 0, q ) ) % sequenceNumber );
                imageTemplate = tmpString + imageTemplate.substr( q );
            }
            else  LOG_WARN << boost::format( "file name template %s does not contain a 2nd format specifier (needs 2)" ) % imageTemplate;
        }*/
    }
    else {
        //LOG_WARN << boost::format( "file name template %s does not contain a format specifier (needs %d)" ) % imageTemplate % ( 1 + ( sequenceNumber >= 0 ) );
    }


    //LOG_DEBUG << "Channel::parseProperties() done.";

}



bpt::ptree Channel::getPropertyTree( bpt::ptree& tree ) {

    bpt::ptree node;

/*
    if( sequenceNumber != myObject.sequenceNumber ) node.put( "SEQUENCE_NUM", sequenceNumber );
    if( mmRow ) node.put( "MMROW", mmRow );
    if( mmWidth ) node.put( "MMWIDTH", mmWidth );
    if( stokesWeights != myObject.stokesWeights ) node.put( "VECTOR", stokesWeights );
    if( !diversity.empty() ) node.put( "DIVERSITY", diversity );
    if( !diversityOrders.empty() ) node.put( "DIV_ORDERS", diversityOrders ); // TODO types missing
    if( image_num_offs != 0 ) node.put( "DT", image_num_offs );
    uint32_t dflags = flags ^ myObject.flags;
    if( dflags & MFBD_NO_RESTORE ) node.put( "NO_RESTORE", ( bool )( flags & MFBD_NO_RESTORE ) );
    if( dflags & MFBD_SAVE_FFDATA ) node.put( "SAVE_FFDATA", ( bool )( flags & MFBD_SAVE_FFDATA ) );
    // TODO "INCOMPLETE"
*/
    ChannelCfg::getProperties(node,myObject);

    tree.push_back( bpt::ptree::value_type( "channel", node ) );

    return node;

}


size_t Channel::size( void ) const {

    size_t sz = ChannelCfg::size();
    sz += sizeof( uint16_t );       // ID;
    sz += sizeof( uint32_t );       // dataOffset;
    sz += dark.size();
    sz += imageStats.size() * Statistics::size() + sizeof(uint16_t);
    return sz;
}


uint64_t Channel::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = ChannelCfg::pack( ptr );
    count += pack( ptr + count, ID );
    count += pack( ptr + count, dataOffset );
    count += dark.pack( ptr + count );
    uint16_t statSize = imageStats.size();
    count += pack( ptr+count, statSize );
    for( auto &it : imageStats ) count += it->pack(ptr+count);
    if(count != size()) {
        LOG_ERR << "(" << hexString(this) << "): Packing failed, there is a size mismatch:  count = " << count << "  sz = " << size();
    }
    return count;
}


uint64_t Channel::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;

    uint64_t count = ChannelCfg::unpack( ptr, swap_endian );
    count += unpack( ptr + count, ID, swap_endian );
    count += unpack( ptr + count, dataOffset, swap_endian );
    count += dark.unpack( ptr + count, swap_endian );
    uint16_t statSize;
    count += unpack(ptr+count, statSize, swap_endian);
    imageStats.resize(statSize);
    for( auto &it : imageStats ) {
        it.reset( new Statistics() );
        count += it->unpack(ptr+count, swap_endian);
    }
    return count;
}


bool Channel::checkCfg(void) {
    
    if( !checkImageScale(telescopeF, arcSecsPerPixel, pixelSize) ) {
        return false;
    }

    //calculatePupilSize( lim_freq, r_c, pupilSize, wavelength, patchSize, myJob.telescopeD, arcSecsPerPixel );

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
    } /*else if( nWild == 2 && sequenceNumber == 0 ) {
        LOG_ERR << "Filename template contains 2 wildcards and no sequence-number given (with SEQUENCE_NUM)";
        return false;
    }*/
    
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


bool Channel::checkData(void) {
    
    // Images
    if( incomplete ) {  // check if files are present
        for( size_t i( 0 ); i < imageNumbers.size(); ) {
            bfs::path fn = bfs::path( boost::str( boost::format( imageTemplate ) % (imageNumberOffset + imageNumbers[i]) ) );
            if( !bfs::exists( fn ) ) {
                fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( imageTemplate ) % (imageNumberOffset + imageNumbers[i]) ) );
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
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( imageTemplate ) % (imageNumberOffset + it) ) );
            if( !bfs::exists( fn ) ) {
                LOG_ERR << boost::format( "Image-file %s not found!" ) % boost::str( boost::format( imageTemplate ) % (imageNumberOffset + it) );
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

    if( !xOffsetFile.empty() ) {
        if( ! bfs::exists( bfs::path( xOffsetFile ) ) ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( xOffsetFile );
            if( ! bfs::exists( fn ) ) {
                LOG_ERR << boost::format( "Offset-file %s not found!" ) % xOffsetFile;
                return false;
            } else xOffsetFile = fn.c_str();
        }
    }

    if( !yOffsetFile.empty() ) {
        if( ! bfs::exists( bfs::path( yOffsetFile ) ) ) {
            bfs::path fn = bfs::path( imageDataDir ) / bfs::path( yOffsetFile );
            if( ! bfs::exists( fn ) ) {
                LOG_ERR << boost::format( "Offset-file %s not found!" ) % yOffsetFile;
                return false;
            } else yOffsetFile = fn.c_str();
        }
    }
    
    return true;
}

void Channel::init( void ) {

}

#define INDEX_THRESHOLD  0
//#define INDEX_THRESHOLD  1E-12
// TBD: this should be a parameter somewhere...
void Channel::initCache( void ) {

    LOG_DETAIL << "wavelength = " << myObject.wavelength << "   patchSize = " << patchSize << "  telescopeD = " << myJob.telescopeD << "  arcSecsPerPixel = " << arcSecsPerPixel;
    calculatePupilSize( frequencyCutoff, pupilRadiusInPixels, pupilPixels, myObject.wavelength, patchSize, myJob.telescopeD, arcSecsPerPixel );
    myJob.patchSize = myObject.patchSize = patchSize;     // TODO: fulhack until per-channel sizes is implemented
    myJob.pupilPixels = myObject.pupilPixels = pupilPixels;
    size_t otfPixels = 2 * pupilPixels;
    LOG_DETAIL << "frequencyCutoff = " << frequencyCutoff << "  pupilSize = " << pupilPixels << "  pupilRadiusInPixels = " << pupilRadiusInPixels;
    pupil = myJob.globalData->fetch(pupilPixels,pupilRadiusInPixels);
    // Create a temporary OTF and store the indices where the OTF/pupil are non-zero. This will be used in loops to skip irrelevant evaluations.
    Array<double> tmpImg( 2 * pupilPixels, 2 * pupilPixels );
    tmpImg.zero();
    Array<double> tmpSubImg( tmpImg, 0, pupilPixels - 1, 0, pupilPixels - 1 );
    pupil.first.copy( tmpSubImg );
    FourierTransform::autocorrelate( tmpImg );
    double* tmpPtr = pupil.first.get();

    for( size_t index = 0; index < pupil.first.nElements(); ++index ) {     // map indices where the pupil-mask is non-zero.
        if( tmpPtr[index] > INDEX_THRESHOLD ) {
            pupilIndices.insert( index );
            size_t otfOffset = ( index / pupilPixels ) * otfPixels + ( index % pupilPixels ) + ( pupilPixels / 2 ) + ( pupilPixels / 2 ) * otfPixels;
            pupilInOTF.insert( make_pair( index, otfOffset ) );
        }
    }

    LOG_DETAIL << "Generated pupilIndices with " << pupilIndices.size() << " elements. (full size = " << pupil.first.nElements() << ", area = " << pupil.second << " )";
    tmpPtr = tmpImg.get();

    for( size_t index = 0; index < tmpImg.nElements(); ++index ) {          // map indices where the OTF-mask (auto-correlated pupil-mask) is non-zero.
        if( tmpPtr[index] > INDEX_THRESHOLD ) {
            otfIndices.insert( index );
        }
    }

    LOG_DETAIL << "Generated otfIndices with " << otfIndices.size() << " elements. (full size = " << tmpImg.nElements() << ")";

    Cache::ModeID id( myJob.klMinMode, myJob.klMaxMode, 0, pupilPixels, pupilRadiusInPixels, rotationAngle );

    for( uint i = 0; i < diversityModes.size(); ++i ) {
        uint16_t modeNumber = diversityModes[i];
        Cache::ModeID id2 = id;
        if(diversityTypes[i] == ZERNIKE) {
            id2.firstMode = id2.lastMode = 0;
        }
        id2.modeNumber = modeNumber;
        myJob.globalData->fetch(id2);
    }

    for( uint16_t& it: myJob.modeNumbers ) {
        if( myJob.modeBasis == ZERNIKE || it == 2 || it == 3 ) {    // force use of Zernike modes for all tilts
            id.firstMode = id.lastMode = 0;
        }
        id.modeNumber = it;
        const PupilMode::Ptr mode = myJob.globalData->fetch(id);
        modes.emplace( it, myJob.globalData->fetch(id) );
    }

}


void Channel::cleanup( void ) {

}



namespace {
    template <typename T>
    void loadWrapper(const string& fn, T& img) {
        redux::file::readFile( fn, img );
        LOG_DETAIL << boost::format( "Loaded file \"%s\"" ) % fn;
    }
}


void Channel::loadData( boost::asio::io_service& service ) {

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

    if( !pupilFile.empty() ) {
        //service.post( std::bind( util::loadPupil, pupilFile, std::ref(pupil), pupilSize ) );
    }
    
    if( !mmFile.empty() ) {
        service.post( std::bind( loadWrapper< Image<float> >, mmFile, std::ref(modulationMatrix) ) );
    }

    if( !xOffsetFile.empty() ) {
        service.post( std::bind( loadWrapper< Image<int16_t> >, xOffsetFile, std::ref(xOffset) ) );
    }

    if( !yOffsetFile.empty() ) {
        service.post( std::bind( loadWrapper< Image<int16_t> >, yOffsetFile, std::ref(yOffset) ) );
    }

    size_t nImages = imageNumbers.size();
    if( nImages ) {
        imageStats.resize( nImages );
        Image<float> tmp;
        bfs::path fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( imageTemplate ) % imageNumbers[0] ) );
        redux::file::readFile( fn.string(), tmp );                      // read first file to get dimensions
        images.resize( nImages, tmp.dimSize( 0 ), tmp.dimSize( 1 ) );
        for( size_t i = 0; i < nImages; ++i ) {
            imageStats[i].reset(new Statistics());
            service.post( std::bind( &Channel::loadImage, this, i ) );
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


void Channel::preprocessData( boost::asio::io_service& service ) {

    size_t nImages = imageNumbers.size();
    double avgMean = 0.0;
    for( size_t i = 0; i < nImages; ++i ) {
        avgMean += imageStats[i]->mean;
    }
    avgMean /= static_cast<double>( nImages );

    for( size_t i = 0; i < nImages; ++i ) {
        service.post( std::bind( &Channel::preprocessImage, this, i, avgMean ) );
    }
    /*
    Cache& cache = Cache::getCache();
    if( pupil.nDimensions() < 2 ) {
        LOG_DETAIL << "Object::preprocessData()  Generating pupil.";
        const std::pair<Array<double>, double>& ppair = cache.pupil(pupilSize,r_c);
        pupil = ppair.first;
    }
    
    if ( modes.size() < myJob.modeNumbers.size() ) {
        LOG_DETAIL << "Object::preprocessData()  Generating modes.";
        for( auto & modeNumber: myJob.modeNumbers ) {
            if( myJob.modeBasis == KARHUNEN_LOEVE ) {
                modes.emplace( modeNumber, cache.mode( myJob.klMinMode, myJob.klMaxMode, modeNumber, pupilSize, r_c, wavelength, rotationAngle ) );
            } else {
                modes.emplace( modeNumber, cache.mode( modeNumber, pupilSize, r_c, wavelength, rotationAngle ) );
            }
        }
    }
    
    cf2pix = util::cf2pix(myJob.arcSecsPerPixel, myJob.telescopeD);
    cf2def = util::cf2def(1.0,myJob.telescopeD/myJob.telescopeF);
    auto it = modes.find(2);
    if(it != modes.end()) {
        LOG_DETAIL << "Object::preprocessData()  adjusting pixel <-> coefficient conversion for this object.";
        double tmp = it->second->at(pupilSize/2,pupilSize/2+1) - it->second->at(pupilSize/2,pupilSize/2);
        tmp *= 0.5 * wavelength * lim_freq;
        cf2pix *= tmp;
    }
    pix2cf = 1.0/cf2pix;
    */

}


double Channel::getMaxMean(void) const {
    double maxMean = std::numeric_limits<double>::lowest();
    for(Statistics::Ptr imStat: imageStats ) {
        if( imStat->mean > maxMean ) maxMean = imStat->mean;
    }
    return maxMean;
}


void Channel::collectImages(redux::util::Array<float>& stack) const {
    size_t n = images.dimSize(0);
    if ( n ) {
        Array<float> block( stack, dataOffset, dataOffset+n-1, 0, images.dimSize(1)-1, 0, images.dimSize(2)-1 );       // sub-array of stack, starting at offset.
        images.copy(block);
    }
}

            
void Channel::initProcessing( WorkSpace::Ptr ws ) {
    workspace = ws;
    initCache();        // this will initialize modes & pupil for this channel
    initPhiFixed();
    for( auto& m : modes ) {           // check if the tilts are present
        if( m.first == 2 || m.first == 3 ) {
            shared_ptr<Tilts> tilts( new Tilts(*this,m.first));
            auto ret = workspace->tilts.emplace(m.first,tilts);
            if ( ret.second ) {     // new tiltmode, i.e. this will be the reference channel
                tilts->nFreeAlpha = imageNumbers.size();
            } else {                // existing tiltmode, this channel is added as a "relative tilt"
                ret.first->second->addRelativeTilt(tilts);
            }
            tilts->init();
        }
    }
    
    for ( uint16_t i=0; i<imageNumbers.size(); ++i ) {
        uint32_t imageNumber = imageNumbers[i];
        std::shared_ptr<WaveFront>& wf = workspace->wavefronts[imageNumber];
        if(!wf) wf.reset(new WaveFront());
        for(auto& m: modes) {
            if( m.first > 3 ) {     // tilts are treated separately
                wf->addWeight( m.first, 1.0 / ( m.second->atm_rms * myObject.wavelength*myObject.wavelength) );
            }
        }
    }

}


void Channel::initPatch( ChannelData& cd ) {
    
    if(imageNumbers.size() != cd.images.dimSize()) {
        LOG_ERR << "Number of images in stack does not match the imageNumbers.";
    }
    cout << "Channel::initPatch()" << endl;
    subImages.clear();
    for ( uint16_t i=0; i<imageNumbers.size(); ++i ) {
        uint32_t imageNumber = imageNumbers[i];
        std::shared_ptr<SubImage> simg( new SubImage( myObject, *this, workspace->window, cd.images, i, maxLocalShift, maxLocalShift,
                                        patchSize, pupilPixels ) ); // TODO: fix offsets
        subImages.push_back( simg );
        simg->init();
        std::shared_ptr<WaveFront>& wf = workspace->wavefronts[imageNumber];
        if(!wf) cout << "Channel::initPatch(): wf = NULL for imageNumber = " << imageNumber << endl;
        wf->addImage(simg);
    }
}


void Channel::initPhiFixed(void) {
    phi_fixed.resize(pupilPixels,pupilPixels);
    phi_fixed.zero();
    Cache::ModeID id( myJob.klMinMode, myJob.klMaxMode, 0, pupilPixels, pupilRadiusInPixels, rotationAngle );
    uint16_t modeNumber;
    for( uint i = 0; i<diversityModes.size(); ++i ) {
        Cache::ModeID id2 = id;
        modeNumber = diversityModes[i];
        cout << "Channel::initPhiFixed()  i=" << i << endl;
        if( modeNumber==2 || modeNumber==3 || diversityTypes[i] == ZERNIKE ) {
            id2.firstMode = id2.lastMode = 0;
        }
        id2.modeNumber = modeNumber;
        const PupilMode::Ptr mode = myJob.globalData->fetch(id2);
        redux::file::Ana::write( "mode_" + to_string( modeNumber ) + "_" + to_string( i ) + ".f0", *mode );
        phi_fixed.add(*mode, diversity[i]);
        redux::file::Ana::write( "phi-mode_" + to_string( modeNumber ) + "_" + to_string( i ) + ".f0", phi_fixed );
    }
    computePhi();   // no tilts for now, just initialize once
}


void Channel::computePhi(void) {
    //cout << "Channel::computePhi()" << endl;
    phi_channel = phi_fixed;
    static int bla(0);
    if(diversityModes.size()) {
        redux::file::Ana::write( "phi_" + to_string( bla++ ) + ".f0", phi_channel );
    }
    // TODO: add tilt corrections
}


void Channel::addMode(redux::util::Array<double>& phi, uint16_t modenumber, double weight) const {
    const PupilMode::Ptr mode = modes.at(modenumber);
   // cout << "Channel::addMode()  mode = " << modenumber << "  weight = " << weight << endl;
 redux::file::Ana::write("mode_"+to_string(modenumber)+".f0", *mode);
 redux::file::Ana::write("pupil.f0", pupil.first);
      if( mode ) {
        phi.add(*mode, weight);
    }
}


void Channel::getPhi(redux::util::Array<double>& phi, const WaveFront& wf) const {
    phi = phi_channel;
    //return;
    //cout << "Channel::getPhi()" << endl;
    for( auto& it: wf.modes ) {
        //cout << "Channel::getPhi()  it.first = " << it.first << endl;
        const PupilMode::Ptr mode = modes.at(it.first);
        if( mode ) { //&& it.second.second ) { // TODO: possibility to enable/disable modes
        //cout << "Channel::getPhi()  mode = " << hexString(mode.get()) << endl;
            phi.add(*mode, *it.second.value);
        }
    }
}


void Channel::addAllFT(redux::util::Array<double>& ftsum) {
    for ( shared_ptr<SubImage>& it: subImages ) {
        it->addFT(ftsum);
    }
}

double Channel::metric(void) {

  double sum = 0.0;
//   for(shared_ptr<SubImage> &im: subImages) {
//       for( auto& a: im->wf->alpha) {
//           double coeff = a.second.first;
//           sum += coeff*coeff * modes.at(a.first)->inv_atm_rms;
//       }
//   }
  return sum;

}


void Channel::normalizeData(boost::asio::io_service& service, double value) {
    size_t nImages = imageNumbers.size();
    for( size_t i = 0; i < nImages; ++i ) {
        service.post( std::bind( &Channel::normalizeImage, this, i, value ) );
    }
}


void Channel::loadImage( size_t index ) {
    Image<float> subimg( images, index, index, 0, images.dimSize( 1 ) - 1, 0, images.dimSize( 2 ) - 1 );
    bfs::path fn = bfs::path( imageDataDir ) / bfs::path( boost::str( boost::format( imageTemplate ) % imageNumbers[index] ) );
    redux::file::readFile( fn.string(), subimg );
    LOG_DETAIL << boost::format( "Loaded file %s" ) % fn;
    imageStats[index]->getStats(myJob.borderClip,subimg,ST_VALUES);                          // only get min/max/mean
}


void Channel::preprocessImage( size_t index, double avgMean ) {

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
        switch(myJob.fillpixMethod) {
            case FPM_HORINT: {
                LOG_TRACE << "Filling bad pixels using horizontal interpolation.";
                function<float(size_t,size_t)> func = bind(horizontalInterpolation<float>, arrayPtr[index], sy, sx, sp::_1, sp::_2);
                fillPixels(arrayPtr[index], sy, sx, func, std::bind2nd(std::less_equal<float>(),myJob.badPixelThreshold) );
                break;
            }
            case FPM_MEDIAN: {
                // TODO: median method
                break;
            }
            case FPM_INVDISTWEIGHT:       // inverse distance weighting is the default method, so fall through
            default: {
                function<double(size_t,size_t)> func = bind(inverseDistanceWeight<float>, arrayPtr[index], sy, sx, sp::_1, sp::_2);
                //function<double(size_t,size_t)> func = bind(inv_dist_wght, arrayPtr[index], sy, sx, sp::_1, sp::_2);
                fillPixels(arrayPtr[index], sy, sx, func, std::bind2nd(std::less_equal<float>(),myJob.badPixelThreshold) );
            }
        }

    fn2 = bfs::path( fn.leaf().string() + ".5.orig" );
    redux::file::Ana::write( fn2.string(), subimg );
    }

    imageStats[index]->getStats(myJob.borderClip, subimg);      // get all stats for the cleaned up data

    if( modified /*&& saveMask & SF_SAVE_FFDATA*/ ) {
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


void Channel::normalizeImage( size_t index, double value ) {
    LOG_TRACE << boost::format( "Normalizing image %d" ) % (dataOffset+index);
    Image<float> subimg( images, index, index, 0, images.dimSize( 1 ) - 1, 0, images.dimSize( 2 ) - 1 );
    subimg *= (value/imageStats[index]->mean);
    double noise1 = imageStats[index]->noise;
    imageStats[index]->getStats(subimg);
    LOG_TRACE << "  image #" << (dataOffset+index) << "  noise1 = " << noise1 << "  noise2 = " << imageStats[index]->noise;
}


size_t Channel::sizeOfPatch(uint32_t npixels) const {
    size_t sz = sizeof(size_t) + imageStats.size()*sizeof(float);
    sz += npixels*images.dimSize(0)*sizeof(float);
    return sz;
}


void Channel::getPatchData(ChannelData& chData, Point16 patchID) const {
   
    if( imageNumbers.empty() ) return;
    
    chData.offset.x = chData.offset.y = chData.residualOffset.x = chData.residualOffset.y = 0;
    
    uint16_t blockSize = patchSize + 2*maxLocalShift;
    uint16_t halfBlockSize = blockSize/2;
    
    Point16 first(subImagePosY[patchID.y]-halfBlockSize, subImagePosX[patchID.x]-halfBlockSize);
    Point16 last(first.y+blockSize-1, first.x+blockSize-1);
    
    double tmpD;
    Image<float> tmpImages(images, 0, imageNumbers.size()-1, first.y, last.y, first.x, last.x);
    
    if(xOffset.valid()) {
        Statistics stats;
        stats.getStats(Image<int16_t>(xOffset, first.y, last.y, first.x, last.x), ST_VALUES);
        chData.residualOffset.x = modf(stats.mean/100 , &tmpD);
        chData.offset.x = lrint(tmpD);
        tmpD = stats.mean/100;
    }
    
    if(yOffset.valid()) {
        Statistics stats;
        stats.getStats(Image<int16_t>(yOffset, first.y, last.y, first.x, last.x), ST_VALUES);
        chData.residualOffset.y = modf(stats.mean/100 , &tmpD);
        chData.offset.y = lrint(tmpD);
        tmpD = stats.mean/100;
    }

    if( chData.offset.x ) {
        int shift = tmpImages.shift(2,chData.offset.x);
        if( shift != chData.offset.x) {
            chData.residualOffset.x += (chData.offset.x-shift);
            chData.offset.x = shift;
        }
    }

    if( chData.offset.y ) {
        int shift = tmpImages.shift(1,chData.offset.y);
        if( shift != chData.offset.y) {
            chData.residualOffset.y += (chData.offset.y-shift);
            chData.offset.y = shift;
        }
    }
    
    chData.images = tmpImages;

}


void Channel::calcPatchPositions(const std::vector<uint16_t>& y, const std::vector<uint16_t>& x) {
    // For now just copy the anchor positions.
    // TODO: use calibration to map the pixelcoordinates for each channel to allow for non-identical hardware & image scales in different objects/channels.
    subImagePosY = y;
    subImagePosX = x;
}


Point16 Channel::clipImages(void) {

/*    if( alignClip.size() != 4 ) {
        LOG_WARN << "No (or faulty) align-clip supplied, using full images: " << printArray(alignClip,"clip"); 
        return Point16(images.dimSize(1),images.dimSize(2));
    }
*/
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
    
    return Point16(images.dimSize(1),images.dimSize(2));
    
}


