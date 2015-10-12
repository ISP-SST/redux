#include "redux/momfbd/modes.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/util.hpp"

#include "redux/constants.hpp"
#include "redux/logger.hpp"
#include "redux/file/fileio.hpp"
#include "redux/image/utils.hpp"
#include "redux/image/zernike.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/cache.hpp"
#include "redux/util/stringutil.hpp"

#include <cmath>
#include <iostream>

#include <boost/filesystem.hpp>

using namespace redux::file;
using namespace redux::image;
using namespace redux::momfbd;
using namespace redux::util;
using namespace std;

namespace bfs = boost::filesystem;


#define lg Logger::mlg
namespace {

    const string thisChannel = "mode";
    
}


ModeInfo::ModeInfo( string filename, uint16_t nPixels )
    : firstMode(0), lastMode(0), modeNumber(0), nPupilPixels(nPixels),
      pupilRadius(0), angle(0), cutoff(0), filename(filename) {

}


ModeInfo::ModeInfo(uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double angle)
    : firstMode(0), lastMode(0), modeNumber(modeNumber), nPupilPixels(nPoints),
      pupilRadius(pupilRadius), angle(angle), cutoff(0), filename("") {

}


ModeInfo::ModeInfo(uint16_t firstMode, uint16_t lastMode, uint16_t modeNumber, uint16_t nPoints, double pupilRadius, double angle, double cutoff)
    : firstMode(firstMode), lastMode(lastMode), modeNumber(modeNumber), nPupilPixels(nPoints),
      pupilRadius(pupilRadius), angle(angle), cutoff(cutoff), filename("") {
          
}

      
uint64_t ModeInfo::size( void ) const {
    static uint64_t sz = 4*sizeof(uint16_t) + 3*sizeof(double) + 1;
    sz += filename.length();
    return sz;
}


uint64_t ModeInfo::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack(ptr,firstMode);
    count += pack(ptr+count,lastMode);
    count += pack(ptr+count,modeNumber);
    count += pack(ptr+count,nPupilPixels);
    count += pack(ptr+count,pupilRadius);
    count += pack(ptr+count,angle);
    count += pack(ptr+count,cutoff);
    count += pack(ptr+count,filename);
    return count;
}


uint64_t ModeInfo::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = unpack(ptr,firstMode,swap_endian);
    count += unpack(ptr+count,lastMode,swap_endian);
    count += unpack(ptr+count,modeNumber,swap_endian);
    count += unpack(ptr+count,nPupilPixels,swap_endian);
    count += unpack(ptr+count,pupilRadius,swap_endian);
    count += unpack(ptr+count,angle,swap_endian);
    count += unpack(ptr+count,cutoff,swap_endian);
    count += unpack(ptr+count,filename,swap_endian);
    return count;
}


bool ModeInfo::operator<(const ModeInfo& rhs) const {
    if(filename != rhs.filename)
        return filename < rhs.filename;
    if(firstMode != rhs.firstMode)
        return firstMode < rhs.firstMode;
    if(lastMode != rhs.lastMode)
        return lastMode < rhs.lastMode;
    if(modeNumber != rhs.modeNumber)
        return modeNumber < rhs.modeNumber;
    if(nPupilPixels != rhs.nPupilPixels)
        return nPupilPixels < rhs.nPupilPixels;
    if(pupilRadius != rhs.pupilRadius)
        return pupilRadius < rhs.pupilRadius;
    if(angle != rhs.angle)
        return angle < rhs.angle;
    return cutoff < rhs.cutoff;
}


PupilMode::PupilMode(uint16_t modeNumber, uint16_t nPoints, double r_c, double angle) :
    Array<double> (nPoints, nPoints), atm_rms(0)  {      // Zernike

    if(modeNumber == 1) {
        Array<double>::operator=(1.0);
    } else {
        double** modePtr = makePointers(get(), nPoints, nPoints);
        //makeZernike_thi(modePtr,modeNumber,nPoints,r_c,angle);   //FIXME: using MvN's Zernike-generator for comparisons
        makeZernike_mvn(modePtr,modeNumber,nPoints,r_c,angle);
        delPointers(modePtr);
    }

    atm_rms = sqrt(Zernike::covariance(modeNumber,modeNumber));


}


PupilMode::PupilMode(uint16_t firstMode, uint16_t lastMode, uint16_t klModeNumber, uint16_t nPoints, double r_c, double angle, double cutoff) :
     Array<double> (nPoints, nPoints), atm_rms(0) {

    if(firstMode > lastMode) swap(firstMode, lastMode);

    if(klModeNumber < firstMode || klModeNumber > lastMode) {
        throw invalid_argument("klModeNumber (" + to_string(klModeNumber) +
                               ") is not in the range [ firstMode (" + to_string(firstMode) +
                               "), lastMode (" + to_string(lastMode) + ")]");
    }

    zero();
        
    const Zernike::KLPtr& kle = Zernike::karhunenLoeveExpansion(firstMode, lastMode).at(klModeNumber);
    double c;
    
    ModeInfo z_info(0, 0, 0, nPoints, r_c, angle, cutoff);
    for(auto & weight : kle->zernikeWeights) {
        if(fabs(c = weight.second) >= cutoff) {
            z_info.modeNumber = weight.first;
            auto& mode = redux::util::Cache::get< ModeInfo, PupilMode::Ptr >( z_info, PupilMode::Ptr() );
            if( !mode ) {
                mode.reset( new PupilMode( weight.first, nPoints, r_c, angle ) );    // generate Zernike
            }
            this->add(*mode, c);
        }
    }

    atm_rms = sqrt(kle->covariance);

}


ModeSet::ModeSet() : info(""), xTiltIndex(-1), yTiltIndex(-1) {
    
}


ModeSet::ModeSet(ModeSet&& rhs) : redux::util::Array<double>(std::move(reinterpret_cast<redux::util::Array<double>&>(rhs))),
    info(std::move(rhs.info)), xTiltIndex(std::move(rhs.xTiltIndex)), yTiltIndex(std::move(rhs.yTiltIndex)),
    modeNumbers(std::move(rhs.modeNumbers)), modePointers(std::move(rhs.modePointers)) {

}


ModeSet::ModeSet(const ModeSet& rhs) : redux::util::Array<double>(reinterpret_cast<const redux::util::Array<double>&>(rhs)),
    info(rhs.info), xTiltIndex(rhs.xTiltIndex), yTiltIndex(rhs.yTiltIndex),
    modeNumbers(rhs.modeNumbers), modePointers(rhs.modePointers) {
    
}


uint64_t ModeSet::size( void ) const {
    uint64_t sz = Array<double>::size();
    sz += info.size() + 2*sizeof(int32_t);
    sz += modeNumbers.size()*sizeof(uint16_t) + sizeof(uint64_t);
    return sz;
}


uint64_t ModeSet::pack( char* data ) const {
    using redux::util::pack;
    uint64_t count = Array<double>::pack(data);
    count += info.pack(data+count);
    count += pack(data+count,xTiltIndex);
    count += pack(data+count,yTiltIndex);
    count += pack(data+count,modeNumbers);
    return count;
}


uint64_t ModeSet::unpack( const char* data, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Array<double>::unpack(data,swap_endian);
    count += info.unpack(data+count,swap_endian);
    count += unpack(data+count,xTiltIndex,swap_endian);
    count += unpack(data+count,yTiltIndex,swap_endian);
    count += unpack(data+count,modeNumbers,swap_endian);
    modePointers.clear();
    for(unsigned int i=0; i<dimSize(0); ++i) {
        modePointers.push_back( ptr(i,0,0) );
    }
    return count;
}


ModeSet& ModeSet::operator=( const ModeSet& rhs ) {
    redux::util::Array<double>::operator=( reinterpret_cast<const redux::util::Array<double>&>(rhs) );
    modeNumbers = rhs.modeNumbers;
    modePointers = rhs.modePointers;
    xTiltIndex = rhs.xTiltIndex;
    yTiltIndex = rhs.yTiltIndex;
    info = rhs.info;
    return *this;
}


bool ModeSet::load( const string& filename, uint16_t pixels ) {
    
    if ( bfs::exists( bfs::path( filename ) ) ) {
        redux::file::readFile( filename, *this );
        if( nDimensions() != 3 || dimSize(1) != pixels || dimSize(2) != pixels ) {    // mismatch
            LOG_WARN << "File: " << filename << " does not match this " << pixels << "x" << pixels << " ModeSet."
            << printArray(dimensions(),"  filedims");
            clear();
        } else {
            // TODO rescale file to right size
            info.nPupilPixels = pixels;
            info.pupilRadius = info.angle = 0;
            modeNumbers.resize( dimSize(0) );
            std::iota(modeNumbers.begin(), modeNumbers.end(), 0);
            modePointers.clear();
            for(unsigned int i=0; i<dimSize(0); ++i) modePointers.push_back( ptr(i,0,0) );
            xTiltIndex = 1;
            yTiltIndex = 0;     // FIXME  properly detect tilts, this is hardcoded for the mode-file I am testing stuff with !!!
            return true;
        }
    }
    return false;
    
}



void ModeSet::init( const MomfbdJob& job, const Object& obj ) {
    
    // From file
    bfs::path fn = bfs::path(obj.modeFile);
    if ( bfs::exists(fn) ) {
        ModeInfo info( obj.modeFile, obj.pupilPixels );
        ModeSet& ret = job.globalData->get(info);
        unique_lock<mutex> lock(ret.mtx);
        if( ret.empty() ) {    // this set was inserted, so it is not loaded yet.
            LOG << "Loading mode-file: " << obj.modeFile;
            redux::file::readFile( fn.string(), ret );
            if( ret.nDimensions() != 3 || ret.dimSize(1) != obj.pupilPixels || ret.dimSize(2) != obj.pupilPixels ) {    // mismatch
                LOG_WARN << "File: " << obj.modeFile << " does not match this " << obj.pupilPixels << "x" << obj.pupilPixels << " ModeSet."
                << printArray(ret.dimensions(),"  retdims");
            } else {
                // TODO rescale file to right size and detect tilt-modes
                ret.info.nPupilPixels = ret.dimSize(1);
                ret.modeNumbers.resize(ret.dimSize(0));
                ret.info.pupilRadius = info.angle = 0;
                std::iota(ret.modeNumbers.begin(), ret.modeNumbers.end(), 0);
                ret.modePointers.clear();
                for(unsigned int i=0; i<ret.dimSize(0); ++i) ret.modePointers.push_back( ret.ptr(i,0,0) );
                *this = ret;
                return;
            }
        } else {
            if( ret.info.nPupilPixels && ret.info.nPupilPixels == obj.pupilPixels ) {    // matching modes
                LOG_TRACE << "Cloning ModeSet";
                *this = ret;
                return;
            } else {
                LOG_ERR << "The Cache returned a non-matching ModeSet. This might happen if a loaded ModeSet was rescaled (which is not implemented yet).";
            }
        }
    }
    
    // From configuration
    ModeInfo info(job.klMinMode, job.klMaxMode, 0, obj.pupilPixels, obj.pupilRadiusInPixels, obj.rotationAngle, job.klCutoff);
    if (job.modeBasis == ZERNIKE) {     // force use of Zernike modes for all tilts
        info.firstMode = info.lastMode = 0;
    }
    ModeSet& ret = job.globalData->get(info);
    unique_lock<mutex> lock(ret.mtx);
    if( ret.empty() ) {    // this set was inserted, so it is not generated yet.
        if(job.modeBasis == ZERNIKE) {
            ret.generate( obj.pupilPixels, obj.pupilRadiusInPixels, obj.rotationAngle, job.modeNumbers );
        } else {
            ret.generate( obj.pupilPixels, obj.pupilRadiusInPixels, obj.rotationAngle, job.klMinMode, job.klMaxMode, job.modeNumbers, job.klCutoff );
        }

        if( ret.nDimensions() != 3 || ret.dimSize(1) != obj.pupilPixels || ret.dimSize(2) != obj.pupilPixels ) {    // mismatch
            LOG_ERR << "Generated ModeSet does not match. This should NOT happen!!";
        } else {
           *this = ret;
        }
    } else {
        if( ret.info.nPupilPixels && ret.info.nPupilPixels == obj.pupilPixels ) {    // matching modes
            LOG_TRACE << "Cloning ModeSet";
            *this = ret;
        } else {
            LOG_ERR << "The Cache returned a non-matching ModeSet. This should NOT happen!!";
        }
    }


}


void ModeSet::generate( uint16_t pixels, double radius, double angle, const vector<uint16_t>& modes ) {
    
    if ( modes.empty() ) return;
    
    info.nPupilPixels = pixels;
    info.pupilRadius = radius;
    info.angle = angle;
    
    LOG_TRACE << "Generating Zernike ModeSet.";

    ModeInfo base_info(0, 0, 0, pixels, radius, angle, 0);
 
    resize( modes.size(), pixels, pixels );
    
    Array<double> view( reinterpret_cast<const redux::util::Array<double>&>(*this), 0, 0, 0, pixels-1, 0, pixels-1 );
    
    for( auto& it : modes ) {
        ModeInfo info = base_info;
        if ( it == 2 || it == 3 ) {     // force use of Zernike modes for all tilts
            info.firstMode = info.lastMode = 0;
        }
        info.modeNumber = it;
        auto& mode = redux::util::Cache::get< ModeInfo, PupilMode::Ptr >( info );
        if( !mode ) {
            mode.reset( new PupilMode( it, pixels, radius, angle ) );    // Zernike
        }

        view.assign( reinterpret_cast<const redux::util::Array<double>&>(*mode) );
        modePointers.push_back(view.ptr());
        view.shift(0,1);
        modeNumbers.push_back(it);
    }
    
    vector<uint16_t>::const_iterator it = std::find(modeNumbers.begin(), modeNumbers.end(), 2);
    if( it != modeNumbers.end() ) xTiltIndex = (it-modeNumbers.begin());
    it = std::find(modeNumbers.begin(), modeNumbers.end(), 3);
    if( it != modeNumbers.end() ) yTiltIndex = (it-modeNumbers.begin());

}


void ModeSet::generate( uint16_t pixels, double radius, double angle, uint16_t firstZernike, uint16_t lastZernike, const vector<uint16_t>& modes, double cutoff ) {
    
    if ( modes.empty() ) return;
    
    info.nPupilPixels = pixels;
    info.pupilRadius = radius;
    info.angle = angle;
    
    LOG_TRACE << "Generating Karhunen-Loeve ModeSet.";

    ModeInfo base_info(firstZernike, lastZernike, 0, pixels, radius, angle, cutoff);
 
    resize( modes.size(), pixels, pixels );
    
    Array<double> view( reinterpret_cast<const redux::util::Array<double>&>(*this), 0, 0, 0, pixels-1, 0, pixels-1 );
    
    for( auto& it : modes ) {
        ModeInfo info = base_info;
        if ( it == 2 || it == 3 ) {     // force use of Zernike modes for all tilts
            info.firstMode = info.lastMode = 0;
        }
        info.modeNumber = it;
        auto& mode = redux::util::Cache::get< ModeInfo, PupilMode::Ptr >( info );
        if( !mode ) {
            if( it == 2 || it == 3 ) {     // force use of Zernike modes for all tilts
                mode.reset( new PupilMode( it, pixels, radius, angle ) );    // Zernike
            } else {
                mode.reset( new PupilMode( firstZernike, lastZernike, it, pixels, radius, angle, cutoff ) );    // K-L
            }
        }

        view.assign( reinterpret_cast<const redux::util::Array<double>&>(*mode) );
        modePointers.push_back(view.ptr());
        view.shift(0,1);
        modeNumbers.push_back(it);
    }
    
    vector<uint16_t>::const_iterator it = std::find(modeNumbers.begin(), modeNumbers.end(), 2);
    if( it != modeNumbers.end() ) xTiltIndex = (it-modeNumbers.begin());
    it = std::find(modeNumbers.begin(), modeNumbers.end(), 3);
    if( it != modeNumbers.end() ) yTiltIndex = (it-modeNumbers.begin());
   
}


void ModeSet::normalize( const redux::image::Pupil& pup ) {

    Array<double> view( reinterpret_cast<const redux::util::Array<double>&>(*this), 0, 0, 0, info.nPupilPixels-1, 0, info.nPupilPixels-1 );
    ArrayStats stats;

    for( unsigned int i=0; i<dimSize(0); ++i ) {
        stats.getMinMaxMean( view * pup );
        view *= (pup.area/stats.norm);
        view.shift(0,1);
    }
   
   
    
}


void ModeSet::setPupilSize( uint16_t nPixels, double radiusInPixels , double rot ) {
    
    info.nPupilPixels = nPixels;
    info.pupilRadius = radiusInPixels;
    info.angle = rot;
        
}

