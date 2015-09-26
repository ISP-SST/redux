#include "redux/image/pupil.hpp"

#include "redux/logger.hpp"
#include "redux/file/fileio.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/image/utils.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/channel.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/util/arraystats.hpp"

#include <boost/filesystem.hpp>

using namespace redux;
using namespace redux::file;
using namespace redux::image;
using namespace redux::momfbd;
using namespace redux::util;

using namespace std;

namespace bfs = boost::filesystem;


#define lg Logger::mlg

namespace {
        
    const string thisChannel = "pupil";
        
    double makePupil(double** pupil, uint16_t nPoints, double radius) {

        double area = 0.0, origin = 0.5;
        memset (*pupil, 0, nPoints * nPoints * sizeof (double));
        uint16_t mid = nPoints / 2;                         // N.B: Don't use odd number of points for pupils !!
                                                            // Pupil should be centered ON pixel (mid,mid), to match the location of the origin for Fourier-transforms.

        if (nPoints % 2) {
            // TODO: warn or throw if nPoints is odd.
        }

        Grid grid (mid + 1, origin, origin);
        float** distPtr = grid.distance.get();              // distance(i,j) is the distance from the centre of the pupil, to the inner boundary of pixel (i,j)
                                                            // i.e. dist(0,0) = dist(0,1) = dist(1,0) = dist(1,1) = sqrt(2)/2  (it is centered on that pixel)
        double val;
        for (uint x = 0; x < mid; ++x) {
            for (uint y = 0; y <= x; ++y) {                 // We only generate the first octant, then copy.
                val = 0;
                if (distPtr[y + 1][x + 1] < radius) {
                    val = 1;
                } else if (distPtr[y][x] < radius) {        // partial pixel
                    if (x == 0 && y == 0) {                 // central pixel = 1 for all practical cases
                        if (radius < 0.5) val = M_PI*radius*radius;    // a pupil of size < sqrt(2) pixel is a bit absurd...
                        else val = M_PI*radius*radius + (radius - 0.5) / (sqrt (0.5) - 0.5) * (1 - M_PI*radius*radius);
                    } else {
                        // TBD: better approximation of pixel fill-factor ??
                        val = (radius - distPtr[y][x]) / (distPtr[y + 1][x + 1] - distPtr[y][x]); // linear fill-factor from radial ratio
                    }
                }
                if (val > 0) {
                    pupil[mid + y][mid + x] = val;
                    area += val;
                    if (x != y) {
                        pupil[mid + x][mid + y] = val; // Use symmetry to fill the second octant
                        area += val;
                    }
                }
            }
        }
        for (uint x = 0; x < mid; ++x) {
            for (uint y = 0; y < mid; ++y) {     // copy 1st quadrant to 2,3,4
                val = pupil[mid + y][mid + x];
                if (val > 0) {
                    if (x) {
                        pupil[mid + y][mid - x] = val;
                        area += val;
                    }
                    if (y) {
                        pupil[mid - y][mid + x] = val;
                        area += val;
                    }
                    if (x && y) {
                        pupil[mid - y][mid - x] = val;
                        area += val;
                    }
                }
            }
        }

        return area;

    }

    
}


PupilInfo::PupilInfo( string filename, uint16_t pixels )
    : nPixels(pixels), pupilRadius(0), filename(filename) {

}


PupilInfo::PupilInfo( uint16_t pixels, double pupilRadius )
    : nPixels(pixels), pupilRadius(pupilRadius), filename("") {

}

      
uint64_t PupilInfo::size( void ) const {
    static uint64_t sz = sizeof(uint16_t) + sizeof(double) + 1;
    sz += filename.length();
    return sz;
}


uint64_t PupilInfo::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = pack(ptr,nPixels);
    count += pack(ptr+count,pupilRadius);
    count += pack(ptr+count,filename);
    return count;
}


uint64_t PupilInfo::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = unpack(ptr,nPixels,swap_endian);
    count += unpack(ptr+count,pupilRadius,swap_endian);
    count += unpack(ptr+count,filename,swap_endian);
    return count;
}


bool PupilInfo::operator<(const PupilInfo& rhs) const {
    if(filename == rhs.filename) {
        if(nPixels == rhs.nPixels) {
            return (pupilRadius < rhs.pupilRadius);
        } else return (nPixels < rhs.nPixels);
    } else return (filename < rhs.filename);
}


void Pupil::calculatePupilSize (double &frequencyCutoff, double &pupilRadiusInPixels, uint16_t &nPupilPixels, double wavelength, uint32_t nPixels, double telescopeDiameter, double arcSecsPerPixel) {
    static double radians_per_arcsec = M_PI / (180.0 * 3600.0);         // (2.0*PI)/(360.0*3600.0)
    double radians_per_pixel = arcSecsPerPixel * radians_per_arcsec;
    double q_number = wavelength / (radians_per_pixel * telescopeDiameter);
    frequencyCutoff = (double) nPixels / q_number;
    nPupilPixels = nPixels >> 2;
    pupilRadiusInPixels = frequencyCutoff / 2.0;                   // telescope radius in pupil pixels...
    if (nPupilPixels < pupilRadiusInPixels) {            // this should only be needed for oversampled images
        uint16_t goodsizes[] = { 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45, 48, 50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128, 135, 144 };
        for (int i = 0; (nPupilPixels = max (goodsizes[i], nPupilPixels)) < pupilRadiusInPixels; ++i);     // find right size
    }
    nPupilPixels <<= 1;
}


Pupil::Pupil( uint16_t pixels, double pupilRadius )
    : nPixels(pixels), radius(pupilRadius), area(0){

}

Pupil::Pupil(Pupil&& rhs) : redux::util::Array<double>(std::move(reinterpret_cast<redux::util::Array<double>&>(rhs))),
    nPixels(std::move(rhs.nPixels)), radius(std::move(rhs.radius)),
    area(std::move(rhs.area)), pupilSupport(std::move(rhs.pupilSupport)),
    otfSupport(std::move(rhs.otfSupport)), pupilInOTF(std::move(rhs.pupilInOTF)) {

}


Pupil::Pupil(const Pupil& rhs) : redux::util::Array<double>(reinterpret_cast<const redux::util::Array<double>&>(rhs)),
    nPixels(rhs.nPixels), radius(rhs.radius), area(rhs.area),
    pupilSupport(rhs.pupilSupport), otfSupport(rhs.otfSupport), pupilInOTF(rhs.pupilInOTF)  {
    
}

            uint16_t nPixels;
            double radius;
            double area;

uint64_t Pupil::size( void ) const {
    uint64_t sz = Array<double>::size();
    sz += sizeof(nPixels) + sizeof(radius) + sizeof(area);
    sz += pupilSupport.size()*sizeof(size_t) + sizeof(uint64_t);
    sz += otfSupport.size()*sizeof(size_t) + sizeof(uint64_t);
    sz += pupilInOTF.size()*2*sizeof(size_t) + sizeof(uint64_t);
    return sz;
}


uint64_t Pupil::pack( char* data ) const {
    using redux::util::pack;
    uint64_t count = Array<double>::pack(data);
    count += pack(data+count,nPixels);
    count += pack(data+count,radius);
    count += pack(data+count,area);
    count += pack(data+count,pupilSupport);
    count += pack(data+count,otfSupport);
    count += pack(data+count,(uint64_t)pupilInOTF.size());
    for( const auto& it: pupilInOTF ) {
        count += pack(data+count,it.first);
        count += pack(data+count,it.second);
    }
    return count;
}


uint64_t Pupil::unpack( const char* data, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Array<double>::unpack(data,swap_endian);
    count += unpack(data+count,radius,swap_endian);
    count += unpack(data+count,nPixels,swap_endian);
    count += unpack(data+count,area,swap_endian);
    count += unpack(data+count,pupilSupport,swap_endian);
    count += unpack(data+count,otfSupport,swap_endian);
    uint64_t tmp;
    count += unpack(data+count,tmp,swap_endian);
    pupilInOTF.resize(tmp);
    for( auto& it: pupilInOTF ) {
        count += unpack(data+count,it.first,swap_endian);
        count += unpack(data+count,it.second,swap_endian);
    }
    return count;
}


bool Pupil::load( const string& filename, uint16_t pixels ) {
    
    if ( bfs::exists( bfs::path( filename ) ) ) {
        redux::file::readFile( filename, *this );
        if( nDimensions() != 2 || dimSize(0) != pixels || dimSize(1) != pixels ) {    // mismatch
            LOG_WARN << "File: " << filename << " does not match this " << pixels << "x" << pixels << " Pupil."
            << printArray(dimensions(),"  filedims");
            clear();
        } else {
            *this *= 1.1;
            *this -= 0.05;
            // TODO rescale file to right size
            nPixels = pixels;
            radius = 0;
            normalize();
            generateSupport(1E-9);                         // TODO: tweak or make into a config parameter?
            return true;
        }
    }
    return false;
    
}


void Pupil::generate( uint16_t pixels, double pupilRadius ) {
    
    nPixels = pixels;
    radius = pupilRadius;
    resize( nPixels, nPixels );
    auto ptr = reshape( nPixels, nPixels );        // returns a 2D shared_ptr
    area = makePupil( ptr.get(), nPixels, radius );
    
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
    for (auto & it: subOTF) {                                           // find the indices where the pupil is > threshold
        if( it > threshold ) {
            pupilSupport.push_back(cnt);
            size_t otfOffset = (cnt / nPixels) * otfPixels + (cnt % nPixels) + (nPixels / 2) + (nPixels / 2) * otfPixels;
            pupilInOTF.push_back(make_pair (cnt, otfOffset));
            area += it;
        }
        cnt++;
    }
    
    FourierTransform::autocorrelate(OTF);                               // auto-correlate the pupil to generate the support of the OTF.

    LOG_DEBUG << "Generated pupilSupport with " << pupilSupport.size() << " elements. (full size = " << nElements() << ", area = " << area << " )";
    double* tmpPtr = OTF.get();
    for (size_t index = 0; index < OTF.nElements(); ++index) {           // map indices where the OTF-mask (auto-correlated pupil-mask) is non-zero.
        if (fabs(tmpPtr[index]) > threshold) {
            otfSupport.push_back(index);
        }
    }
    LOG_DEBUG << "Generated otfSupport with " << otfSupport.size() << " elements. (full size = " << OTF.nElements() << ", area ~ " << (4.0*area) << " )";
}


void Pupil::normalize( void ) {
    
    ArrayStats stats;
    stats.getMinMaxMean(*this);

    if( stats.min != 0.0 || stats.max != 1.0 ) {
        LOG_WARN << "The pupil will be naively re-scaled to the interval [0,1], if this is not what you want, ask someone to implement a better normalization.";
    }

    // FIXME:  decide better normalization scheme to allow max values < 1 and min values > 0 (i.e. realistic transmission of an aperture)
    *this -= stats.min;
    if( stats.min == stats.max ) {
        LOG_WARN << "The pupil is a constant value.";
    } else {
        *this *= 1.0/(stats.max-stats.min);
    }
    
    

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

