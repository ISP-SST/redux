#include "redux/momfbd/modes.hpp"


#include "redux/momfbd/util.hpp"

#include "redux/file/fileio.hpp"
#include "redux/image/utils.hpp"
#include "redux/image/zernike.hpp"
#include "redux/util/cache.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

#include <boost/filesystem.hpp>

using namespace redux::file;
using namespace redux::image;
using namespace redux::momfbd;
using namespace redux::util;
using namespace std;

namespace bfs = boost::filesystem;


ModeInfo::ModeInfo( const string& filename, uint16_t nPixels )
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


ModeInfo::ModeInfo(uint16_t firstMode, uint16_t lastMode, const ModeList& modeList, uint16_t nPoints, double pupilRadius, double angle, double cutoff)
    : firstMode(firstMode), lastMode(lastMode), modeNumber(0), nPupilPixels(nPoints), modeList(modeList),
      pupilRadius(pupilRadius), angle(angle), cutoff(cutoff), filename("") {
          
}

      
uint64_t ModeInfo::size( void ) const {
    static uint64_t sz = 4*sizeof(uint16_t) + 3*sizeof(double) + 1;
    sz += redux::util::size(modeList);
    sz += filename.length();
    return sz;
}


uint64_t ModeInfo::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack(ptr,firstMode);
    count += pack(ptr+count,lastMode);
    count += pack(ptr+count,modeNumber);
    count += pack(ptr+count,modeList);
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
    count += unpack(ptr+count,modeList,swap_endian);
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
    if(modeList != rhs.modeList)
        return modeList < rhs.modeList;
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


ModeInfo::operator string() const {
    string ret = to_string(nPupilPixels)+":"+to_string(pupilRadius)+":";
    if( firstMode == lastMode ) {
        ret += "Z";
    } else {
        ret += "KL";
    }
    if( !modeList.empty() ) {
        //ret += "["+uIntsToString( modeNumbers )+"]";  TODO FIXME
        ret += "[]";
    } else  {
        ret += "("+ to_string(modeNumber) + ")";
    }
    if( firstMode != lastMode ) {
        ret += ":(Z-basis=[" + to_string(firstMode) + "," +to_string(lastMode)+ "])";
    }
    if( !filename.empty() ) {
        ret += ":\"" + filename + "\"";
    }
    return ret;
}


PupilMode::PupilMode(uint16_t modeNumber, uint16_t nPoints, double r_c, double angle, int flags ) :
    Array<double> (nPoints, nPoints), atm_rms(0)  {      // Zernike

    if( modeNumber < 2 ) {
        Array<double>::operator=(modeNumber);
    } else {
        Zernike::getZernike( get(), nPoints, r_c, angle, modeNumber, flags );
    }

    atm_rms = sqrt(Zernike::getCovariance (modeNumber,modeNumber));


}


PupilMode::PupilMode(uint16_t firstZernike, uint16_t lastZernike, uint16_t klModeNumber,
                     uint16_t nPoints, double r_c, double angle, double cutoff, int flags ) :
     Array<double> (nPoints, nPoints), atm_rms(0) {

    if( firstZernike > lastZernike ) swap( firstZernike, lastZernike );

    if(klModeNumber < firstZernike || klModeNumber > lastZernike ) {
        throw invalid_argument("klModeNumber (" + to_string(klModeNumber) +
                               ") is not in the range [ firstMode (" + to_string( firstZernike ) +
                               "), lastMode (" + to_string( lastZernike ) + ")]");
    }

    zero();
    
        
    const Zernike::KLPtr& kle = Zernike::karhunenLoeveExpansion( firstZernike, lastZernike ).at(klModeNumber);

    ModeInfo z_info(0, 0, 0, nPoints, r_c, angle, cutoff);
    for(auto & weight : kle->zernikeWeights) {
        double c = weight.second;
        if(fabs(c) >= cutoff) {
            z_info.modeNumber = weight.first;
            auto& mode = redux::util::Cache::get< ModeInfo, PupilMode::Ptr >( z_info, PupilMode::Ptr() );
            if( !mode ) {
                mode.reset( new PupilMode( weight.first, nPoints, r_c, angle, flags ) );    // generate Zernike
            }
            this->add(*mode, c);
        }
    }
    
    if( flags&Zernike::NORMALIZE ) {     // numerical normalization
        size_t blockSize = nPoints*nPoints;
        long double normalization(0.0);
        redux::image::Pupil pupil = Pupil::fetch( nPoints, r_c );
        double* pupPtr = pupil.ptr();
        double* mPtr = ptr();
        for( size_t i(0); i<blockSize; ++i ) {
            if( pupPtr[i] > 0 ) {
                normalization += mPtr[i]*mPtr[i]*pupPtr[i];
            }
        }
        normalization = sqrtl( pupil.area/normalization );
        std::transform( mPtr, mPtr+blockSize, mPtr, std::bind(std::multiplies<long double>(), std::placeholders::_1, normalization));

    }
    atm_rms = sqrt(kle->covariance);

}


ModeSet::ModeSet() : info(""), tiltMode(-1,-1) {
    
}


ModeSet::ModeSet(ModeSet&& rhs) : redux::util::Array<double>(std::move(reinterpret_cast<redux::util::Array<double>&>(rhs))),
    info(std::move(rhs.info)), tiltMode(std::move(rhs.tiltMode)), modeList(std::move(rhs.modeList)),
    modePointers(std::move(rhs.modePointers)), atm_rms(std::move(rhs.atm_rms)) {

}


ModeSet::ModeSet(const ModeSet& rhs) : redux::util::Array<double>(reinterpret_cast<const redux::util::Array<double>&>(rhs)),
    info(rhs.info), tiltMode(rhs.tiltMode), modeList(rhs.modeList), modePointers(rhs.modePointers),
    atm_rms(rhs.atm_rms) {
    
}


uint64_t ModeSet::size( void ) const {
    uint64_t sz = Array<double>::size();
    sz += info.size() + tiltMode.size() + shiftToAlpha.size();
    sz += redux::util::size(modeList);
    sz += atm_rms.size()*sizeof(double) + sizeof(uint64_t);
    sz += norms.size()*sizeof(double) + sizeof(uint64_t);
    return sz;
}


uint64_t ModeSet::pack( char* data ) const {
    using redux::util::pack;
    uint64_t count = Array<double>::pack(data);
    count += info.pack(data+count);
    count += tiltMode.pack(data+count);
    count += shiftToAlpha.pack(data+count);
    count += pack(data+count,modeList);
    count += pack(data+count,atm_rms);
    count += pack(data+count,norms);
    return count;
}


uint64_t ModeSet::unpack( const char* data, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Array<double>::unpack(data,swap_endian);
    count += info.unpack(data+count,swap_endian);
    count += tiltMode.unpack(data+count,swap_endian);
    count += shiftToAlpha.unpack(data+count,swap_endian);
    count += unpack(data+count,modeList,swap_endian);
    count += unpack(data+count,atm_rms,swap_endian);
    count += unpack(data+count,norms,swap_endian);
    modePointers.clear();
    for(unsigned int i=0; i<dimSize(0); ++i) {
        modePointers.push_back( ptr(i,0,0) );
    }
    return count;
}


ModeSet& ModeSet::operator=( const ModeSet& rhs ) {
    redux::util::Array<double>::operator=( reinterpret_cast<const redux::util::Array<double>&>(rhs) );
    modeList = rhs.modeList;
    atm_rms = rhs.atm_rms;
    norms = rhs.norms;
    modePointers = rhs.modePointers;
    tiltMode = rhs.tiltMode;
    shiftToAlpha = rhs.shiftToAlpha;
    info = rhs.info;
    return *this;
}


ModeSet ModeSet::clone( void ) const {

    ModeSet tmp(*this);
    copy( reinterpret_cast<redux::util::Array<double>&>(tmp) );

    tmp.modePointers.clear();
    if( tmp.nDimensions() == 2 ) {
        tmp.modePointers.push_back( tmp.ptr(0,0) );
    } else {
        for(unsigned int i=0; i<tmp.dimSize(0); ++i) {
            tmp.modePointers.push_back( tmp.ptr(i,0,0) );
        }
    }

    return tmp;
    
}


bool ModeSet::load( const string& filename, uint16_t pixels ) {
    
    if ( bfs::is_regular_file(filename) ) {
        redux::file::readFile( filename, *this );
        if( nDimensions() != 3 || dimSize(1) != pixels || dimSize(2) != pixels ) {    // mismatch
            clear();
        } else {
            // TODO rescale file to right size
            info.nPupilPixels = pixels;
            info.pupilRadius = info.angle = 0;
            modeList.resize( dimSize(0) );
            uint16_t cnt(0);
            for( auto& m: modeList ) {
                m.type = MB_NONE;
                m.mode = cnt++;
            }
            modePointers.clear();
            for( unsigned int i=0; i<dimSize(0); ++i) {
                modePointers.push_back( ptr(i,0,0) );
            }
            // FIXME  properly detect tilts, this is hardcoded for the mode-files with tilts as the old MOMFBD code !!!
            tiltMode.x = 1;
            tiltMode.y = 0;
            return true;
        }
    }
    return false;
    
}


void ModeSet::generate( uint16_t pixels, double radius, double angle, const ModeList& modes, int flags ) {
    
    resize();   // clear
    modeList = modes;
    modePointers.clear();

    if ( modes.empty() ) return;
    
    info.nPupilPixels = pixels;
    info.pupilRadius = radius;
    info.angle = angle;

    ModeInfo base_info(0, 0, 0, pixels, radius, angle, 0);
 
    resize( modeList.size(), pixels, pixels );
    
    Array<double> view( reinterpret_cast<const redux::util::Array<double>&>(*this), 0, 0, 0, pixels-1, 0, pixels-1 );
    set<Point16> forced;
    uint16_t n;
    int16_t m;
    for( auto& it : modeList ) {
        int tmp_flags = flags;
        ModeInfo minfo = base_info;
        minfo.modeNumber = it.mode;
        if( minfo.modeNumber == 2 || minfo.modeNumber == 3 || (it.type == ZERNIKE) ) {     // force use of Zernike modes for all tilts
            minfo.firstMode = minfo.lastMode = 0;
            if ( minfo.modeNumber == 2 ) tiltMode.x = modePointers.size();
            else if ( minfo.modeNumber == 3 ) tiltMode.y = modePointers.size();
            it.type = ZERNIKE;
        }
        Zernike::NollToNM( minfo.modeNumber, n, m );
        Point16 nm(n,abs(m));
        if( forced.count(nm) ) {  // already forced.
            tmp_flags &= ~Zernike::FORCE;
        } else {
            forced.insert(nm);
        }
        auto& mode = redux::util::Cache::get< ModeInfo, PupilMode::Ptr >( minfo );
        if( !mode || (tmp_flags&Zernike::FORCE) ) {
            mode.reset( new PupilMode( minfo.modeNumber, pixels, radius, angle, tmp_flags ) );    // Zernike
        }

        view.assign( reinterpret_cast<const redux::util::Array<double>&>(*mode) );
        modePointers.push_back(view.ptr(0,0,0));
        view.shift(0,1);
        atm_rms.push_back(mode->atm_rms);
    }

}


void ModeSet::generate( uint16_t pixels, double radius, double angle, uint16_t firstZernike, uint16_t lastZernike, const ModeList& modes, double cutoff, int flags ) {
    
    resize();       // clear
    modeList = modes;
    modePointers.clear();
    
    if ( modes.empty() ) return;
    
    info.nPupilPixels = pixels;
    info.pupilRadius = radius;
    info.angle = angle;
    
    ModeInfo base_info(firstZernike, lastZernike, 0, pixels, radius, angle, cutoff);
 
    resize( modeList.size(), pixels, pixels );
    
    Array<double> view( reinterpret_cast<const redux::util::Array<double>&>(*this), 0, 0, 0, pixels-1, 0, pixels-1 );
    
    for( auto& it : modeList ) {
        ModeInfo info = base_info;
        info.modeNumber = it.mode;
        if ( info.modeNumber == 2 || info.modeNumber == 3 || (it.type == ZERNIKE) ) {     // force use of Zernike modes for all tilts
            info.firstMode = info.lastMode = 0;
            if ( info.modeNumber == 2 ) tiltMode.x = modePointers.size();
            else if ( info.modeNumber == 3 ) tiltMode.y = modePointers.size();
            it.type = ZERNIKE;
        }
        auto& mode = redux::util::Cache::get< ModeInfo, PupilMode::Ptr >( info );
        if( !mode || (flags&Zernike::FORCE) ) {
            if( it.type == ZERNIKE ) {     // force use of Zernike modes for all tilts
                mode.reset( new PupilMode( info.modeNumber, pixels, radius, angle, flags ) );    // Zernike
            } else {
                mode.reset( new PupilMode( firstZernike, lastZernike, info.modeNumber, pixels, radius, angle, cutoff, flags ) );    // K-L
            }
        }

        view.assign( reinterpret_cast<const redux::util::Array<double>&>(*mode) );
        modePointers.push_back(view.ptr(0,0,0));
        view.shift(0,1);
        atm_rms.push_back(mode->atm_rms);
    }
   
}


// Find the normalizations such that the integral of each mode over the pupil equals the pupil-area.
void ModeSet::getNorms( const redux::image::Pupil& pup ) {

    const double* pupPtr = pup.get();
    size_t nPixels = info.nPupilPixels*info.nPupilPixels;

    norms.resize(modePointers.size());
    for ( uint16_t i=0; i<modePointers.size(); ++i ) {
        double* ptr = modePointers[i];
        if( !ptr ) continue;
        double norm(0);
        double mx = std::numeric_limits<double>::min();
        double mn = std::numeric_limits<double>::max();
        bool isTilt = (i == tiltMode.x) || (i == tiltMode.y);
        for( size_t ind=0; ind<nPixels; ++ind) {
            if( pupPtr[ind] > 0 ) {
                double tmp = ptr[ind]*pupPtr[ind];
                norm += tmp*ptr[ind];               // i.e. = mode² * pupil
                if( isTilt ) {
                    mx = std::max( mx, tmp );
                    mn = std::min( mn, tmp );
                }
            }
        }
        norms[i] = sqrtl( norm/pup.area );

        if( i == tiltMode.x ) {
            shiftToAlpha.x = 2 * M_PI / (mx-mn); // * norms[i];     // A shift of 1 pixel corresponds to an introduced phase-shift across the pupil of 1 period.
        } else if ( i == tiltMode.y ) {
            shiftToAlpha.y = 2 * M_PI / (mx-mn); // * norms[i];
        }

    }

}


void ModeSet::setPupilSize( uint16_t nPixels, double radiusInPixels , double rot ) {
    
    info.nPupilPixels = nPixels;
    info.pupilRadius = radiusInPixels;
    info.angle = rot;
        
}


void ModeSet::normalize( double scale ) {
    
    size_t nPixels = info.nPupilPixels*info.nPupilPixels;

    for ( uint16_t i=0; i<modePointers.size(); ++i ) {
        double* ptr = modePointers[i];
        double mode_scale = scale/norms[i];
        std::transform( ptr, ptr+nPixels, ptr, std::bind1st(std::multiplies<double>(), mode_scale) );
        if( i == tiltMode.x ) {
            shiftToAlpha.x /= mode_scale;
        } else if (i == tiltMode.y) {
            shiftToAlpha.y /= mode_scale;
        }
        norms[i] = mode_scale;
    }
    
}
