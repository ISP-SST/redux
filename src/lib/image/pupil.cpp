#include "redux/image/pupil.hpp"

#include "redux/file/fileana.hpp"
#include "redux/file/fileio.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/grid.hpp"
#include "redux/image/utils.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/util/arraystats.hpp"

#include <boost/filesystem.hpp>

using namespace redux;
using namespace redux::file;
using namespace redux::image;
using namespace redux::util;

using namespace std;

namespace bfs = boost::filesystem;


PupilInfo::PupilInfo( string filename, double pupilRadius, uint16_t pixels )
    : nPixels(pixels), pupilRadius(pupilRadius), coRadius(0.0), filename(filename) {

}


PupilInfo::PupilInfo( uint16_t pupilPixels, double pupilRadius, double coRadius )
    : nPixels(pupilPixels), pupilRadius(pupilRadius), coRadius(coRadius), filename("") {

}

      
uint64_t PupilInfo::size( void ) const {
    static uint64_t sz = sizeof(uint16_t) + 2*sizeof(double) + 1;
    sz += filename.length();
    return sz;
}


uint64_t PupilInfo::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack(ptr,nPixels);
    count += pack(ptr+count,pupilRadius);
    count += pack(ptr+count,coRadius);
    count += pack(ptr+count,filename);
    return count;
}


uint64_t PupilInfo::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = unpack(ptr,nPixels,swap_endian);
    count += unpack(ptr+count,pupilRadius,swap_endian);
    count += unpack(ptr+count,coRadius,swap_endian);
    count += unpack(ptr+count,filename,swap_endian);
    return count;
}


bool PupilInfo::operator<(const PupilInfo& rhs) const {
    if(filename != rhs.filename) return (filename < rhs.filename);
    if(nPixels != rhs.nPixels) return (nPixels < rhs.nPixels);
    if(pupilRadius != rhs.pupilRadius) return (pupilRadius < rhs.pupilRadius);
    return (coRadius < rhs.coRadius);
}


PupilInfo::operator string() const {
    string ret = to_string(nPixels)+":"+to_string(pupilRadius);
    if( coRadius > 0 ) ret +=  ":" + bitString(coRadius);
    if( !filename.empty() ) ret += ":\"" + filename + "\"";
    return ret;
}


void Pupil::calculatePupilSize (double &frequencyCutoff, double &pupilRadiusInPixels, uint16_t &nPupilPixels, double wavelength, uint32_t nPixels, double telescopeDiameter, double arcSecsPerPixel) {
    static double radians_per_arcsec = M_PI / (180.0 * 3600.0);     // (2.0*PI)/(360.0*3600.0)
    static const vector<uint16_t> goodsizes = { 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45, 48, 50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128, 135, 144 };
    if( nPixels < 4 ) throw logic_error("calculatePupilSize: nPixels (="+to_string(nPixels)+") < 4, which doesn't make much sense.");
    if( wavelength <= 0.0 ) throw logic_error("calculatePupilSize: wavelength (="+to_string(wavelength)+") <= 0.0, which doesn't make much sense.");
    if( telescopeDiameter <= 0.0 ) throw logic_error("calculatePupilSize: telescopeDiameter (="+to_string(telescopeDiameter)+") <= 0.0, which doesn't make much sense.");
    if( arcSecsPerPixel <= 0.0 ) throw logic_error("calculatePupilSize: arcSecsPerPixel (="+to_string(arcSecsPerPixel)+") <= 0.0, which doesn't make much sense.");
    double radians_per_pixel = arcSecsPerPixel * radians_per_arcsec;
    double q_number = wavelength / (radians_per_pixel * telescopeDiameter);
    frequencyCutoff = (double) nPixels / q_number;                  // Diffraction limit (radius) in Fourier space.
    nPupilPixels = nPixels >> 2;                                    // Divide nPixels by 4 (the loop below actually works with half pupil-sizes).
    pupilRadiusInPixels = frequencyCutoff / 2.0;                    // Telescope radius in pupil pixels.
    if( nPupilPixels < pupilRadiusInPixels ) {                      // Only increase pupil-pixels for oversampled data
        for( auto& gs: goodsizes ) {                                // Find a usable size in the list specified above
            nPupilPixels = max( gs, nPupilPixels );
            if( nPupilPixels >= pupilRadiusInPixels ) {             // Now the pupil fits -> break
                break;
            }
        }
    }
    nPupilPixels <<= 1;

}


Pupil::Pupil( uint16_t pixels, double pupilRadius )
    : info( pixels, pupilRadius, 0 ), nPixels(pixels), radius(pupilRadius), co_radius(0), area(0) {

    generate(pixels,pupilRadius);
    
}

Pupil::Pupil(Pupil&& rhs) : redux::util::Array<double>(std::move(reinterpret_cast<redux::util::Array<double>&>(rhs))),
    info( std::move(rhs.info) ), nPixels(std::move(rhs.nPixels)), radius(std::move(rhs.radius)),
    co_radius(std::move(rhs.co_radius)), area(std::move(rhs.area)), pupilSupport(std::move(rhs.pupilSupport)),
    otfSupport(std::move(rhs.otfSupport)), pupilInOTF(std::move(rhs.pupilInOTF)) {

}


Pupil::Pupil(const Pupil& rhs) : redux::util::Array<double>(reinterpret_cast<const redux::util::Array<double>&>(rhs)),
    info(rhs.info), nPixels(rhs.nPixels), radius(rhs.radius), co_radius(rhs.co_radius), area(rhs.area),
    pupilSupport(rhs.pupilSupport), otfSupport(rhs.otfSupport), pupilInOTF(rhs.pupilInOTF)  {
    
}


uint64_t Pupil::size( void ) const {
    uint64_t sz = Array<double>::size();
    sz += info.size();
    sz += sizeof(nPixels) + sizeof(radius) + sizeof(area);
    sz += pupilSupport.size()*sizeof(size_t) + sizeof(uint64_t);
    sz += otfSupport.size()*sizeof(size_t) + sizeof(uint64_t);
    sz += pupilInOTF.size()*2*sizeof(size_t) + sizeof(uint64_t);
    return sz;
}


uint64_t Pupil::pack( char* data ) const {
    using redux::util::pack;
    uint64_t count = Array<double>::pack(data);
    count += info.pack(data+count);
    count += pack(data+count,nPixels);
    count += pack(data+count,radius);
    count += pack(data+count,area);
    count += pack(data+count,pupilSupport);
    count += pack(data+count,otfSupport);
    count += pack(data+count,(uint64_t)pupilInOTF.size());
    for( const auto& index: pupilInOTF ) {
        count += pack(data+count,index.first);
        count += pack(data+count,index.second);
    }
    return count;
}


uint64_t Pupil::unpack( const char* data, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Array<double>::unpack(data,swap_endian);
    count += info.unpack(data+count,swap_endian);
    count += unpack(data+count,nPixels,swap_endian);
    count += unpack(data+count,radius,swap_endian);
    count += unpack(data+count,area,swap_endian);
    count += unpack(data+count,pupilSupport,swap_endian);
    count += unpack(data+count,otfSupport,swap_endian);
    uint64_t tmp;
    count += unpack(data+count,tmp,swap_endian);
    pupilInOTF.resize(tmp);
    for( auto& index: pupilInOTF ) {
        count += unpack(data+count,index.first,swap_endian);
        count += unpack(data+count,index.second,swap_endian);
    }
    return count;
}


bool Pupil::load( const string& filename, uint16_t pixels, double pupilRadius ) {
    
    if ( bfs::is_regular_file(filename) ) {
        redux::file::readFile( filename, *this );
        if( nDimensions() != 2 ) {    // not a 2D image
            clear();
        } else {
            size_t inSz = dimSize(0);
            if( inSz != dimSize(1)) {    // not a square image
                clear();
                return false;
            } 
            if( inSz != pixels ) {    // size mismatch
                Array<double> tmp = Array<double>::copy(true);
                resize( pixels, pixels );
                cv::Mat src( inSz, inSz, cv::cvType<double>(), tmp.get() );
                cv::Mat dst( pixels, pixels, cv::cvType<double>(), get() );
                cv::Mat map = cv::Mat::eye( 2, 3, CV_32F );         // unit affine transformation
                double scale = 2.0*pupilRadius/inSz;                // input image is assumed to have a pupil-diameter = inSz.
                map *= scale;
                map.at<float>(0,2) = (pixels-1-inSz*scale)/2.0;     // center pupil image on the output image
                map.at<float>(1,2) = (pixels-1-inSz*scale)/2.0;
                if( !(pixels&1) ) {                                 // center pupil on the (N/2,N/2) pixel for even pixels
                    map.at<float>(0,2) += 0.5;
                    map.at<float>(1,2) += 0.5;
                }
                int flags = cv::INTER_AREA;                  // INTER_LINEAR, INTER_CUBIC, INTER_AREA, INTER_LANCZOS4
                int borderMode = cv::BORDER_CONSTANT;
                const cv::Scalar borderValue = cv::Scalar();
                cv::Size dsize = dst.size();
                cv::warpAffine( src, dst, map, dsize, flags, borderMode, borderValue );
            }
            info.filename = filename;
            info.nPixels = pixels;
            info.pupilRadius = pupilRadius;
            nPixels = pixels;
            radius = pupilRadius;
            normalize();
            generateSupport(1E-9);                         // TODO: tweak or make into a config parameter?
            return true;
        }
    }
    return false;
    
}


void Pupil::generate( uint16_t pixels, double pupilRadius, double coRadius ) {
    
    info = PupilInfo( pixels, pupilRadius, coRadius );
    nPixels = pixels;
    radius = pupilRadius;
    co_radius = coRadius;
    resize( nPixels, nPixels );
    auto ptr = reshape( nPixels, nPixels );        // returns a 2D shared_ptr
    area = makePupil( ptr.get(), nPixels, radius, co_radius );

    normalize();

    generateSupport(1E-9);                         // TODO: tweak or make into a config parameter?
    
}



void Pupil::generateSupport(double threshold){
    
    if(nDimensions() != 2 || dimSize(0) != nPixels || dimSize(1) != nPixels) return;            // init needs to be called first.

    size_t otfPixels = 2*nPixels;
    Array<double> OTF(otfPixels, otfPixels);
    OTF.zero();
    
    Array<double> subOTF(OTF, 0, nPixels-1, 0, nPixels-1);              // bottom-left quadrant of the OTF images
    copy(subOTF);                                                       // copy the pupil
    
    size_t cnt(0);
    area = 0;
    for (auto & value: subOTF) {                                           // find the indices where the pupil is > threshold
        if( value > threshold ) {
            pupilSupport.push_back(cnt);
            size_t otfOffset = (cnt / nPixels) * otfPixels + (cnt % nPixels) + (nPixels / 2) + (nPixels / 2) * otfPixels;
            pupilInOTF.push_back(make_pair (cnt, otfOffset));
            area += value;
        }
        cnt++;
    }
    
    FourierTransform::autocorrelate(OTF,true);                           // auto-correlate the pupil to generate the support of the OTF.

    double* tmpPtr = OTF.get();
    for (size_t index = 0; index < OTF.nElements(); ++index) {           // map indices where the OTF-mask (auto-correlated pupil-mask) is non-zero.
        if (fabs(tmpPtr[index]) > threshold) {
            otfSupport.push_back(index);
        }
    }
    
}


void Pupil::normalize( void ) {
    
    ArrayStats stats;
    stats.getMinMaxMean(*this);

    if( stats.min != 0.0 || stats.max != 1.0 ) {
        //cerr << "The pupil will be naively re-scaled to the interval [0,1].\n";
    }

    // FIXME:  decide better normalization scheme to allow max values < 1 and min values > 0 (i.e. realistic transmission of an aperture)
    if( stats.min == stats.max ) {
        //cerr << "The pupil is a constant value.\n";
    } else {
        *this -= stats.min;
        *this *= 1.0/(stats.max-stats.min);
    }
    
    

}


void Pupil::dump( string tag ) const {

    if( nElements() ) {
        Ana::write( tag + ".f0", *this );
        vector<size_t> dims = dimensions();
        Array<uint8_t> support(dims);
        support *= 0;
        uint8_t* ptr = support.get();
        for( const size_t& ind: pupilSupport ) {
            ptr[ind] = 1;
        }
        Ana::write( tag + "_support.f0", support );
        for( size_t& d: dims ) d *= 2;
        support.resize(dims);
        support.zero();
        ptr = support.get();
        for( const size_t& ind: otfSupport ) {
            ptr[ind] = 1;
        }
        Ana::write( tag + "_otfsupport.f0", support );
        support.zero();
        ptr = support.get();
        for( const auto& ind: pupilInOTF ) {
            ptr[ind.second] = 1;
        }
        Ana::write( tag + "_pupilinotf.f0", support );
    }


}


Pupil& Pupil::fetch( uint16_t pupilPixels, double pupilRadius, double coRadius ) {
    
    PupilInfo pi( pupilPixels, pupilRadius, coRadius );
    Pupil& pupil = redux::util::Cache::get<PupilInfo,Pupil>(pi);
    if( pupil.empty() ) {    // pupil does not exist in the cache, so it has to be generated
        pupil.generate( pupilPixels, pupilRadius, coRadius );
    }
    
    return pupil;
    
}


Pupil& Pupil::operator=( const Pupil& rhs ) {
    redux::util::Array<double>::operator=( reinterpret_cast<const redux::util::Array<double>&>(rhs) );
    nPixels = rhs.nPixels;
    radius = rhs.radius;
    area = rhs.area;
    pupilSupport = rhs.pupilSupport;
    otfSupport = rhs.otfSupport;
    pupilInOTF = rhs.pupilInOTF;
    return *this;
}


bool Pupil::operator<(const Pupil& rhs) const {
    if( nPixels == rhs.nPixels ) {
        return (radius < rhs.radius);
    }
    return (nPixels<rhs.nPixels);
}

