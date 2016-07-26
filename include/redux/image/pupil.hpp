#ifndef REDUX_IMAGE_PUPIL_HPP
#define REDUX_IMAGE_PUPIL_HPP

#include "redux/util/array.hpp"

#include <mutex>
#include <vector>

namespace redux {

    namespace image {
        
        struct PupilInfo {
            
            uint16_t nPixels;
            double pupilRadius;
            std::string filename;
            
            PupilInfo( std::string, uint16_t pixels=0 );
            PupilInfo( uint16_t pixels, double pupilRadius );
            
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            bool operator<( const PupilInfo& rhs ) const;
            
        };


        struct Pupil : public redux::util::Array<double> {
            
            static void calculatePupilSize(double &, double &, uint16_t&, double, uint32_t, double, double );
            
            Pupil(void) {};
            Pupil( uint16_t pixels, double pupilRadius );
            Pupil(Pupil&& rhs);
            Pupil(const Pupil& rhs);
            
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            bool load( const std::string& filename, uint16_t pixels );
            
            void generate( uint16_t pixels, double pupilRadius );
            void generateSupport(double threshold=0);                           //!< Gets the indices of elements in the pupil/otf which are >threshold
            void normalize( void );                                             //!< Scale pupil to the interval [0,1]
            
            Pupil& operator=( const Pupil& rhs );
            bool operator<(const Pupil& rhs) const;                             //!< So that the Pupil-struct can be stored in comparative containers (set/map)

            uint16_t nPixels;
            double radius;
            double area;
            
            std::vector<size_t> pupilSupport, otfSupport;                       //!< The indices to the support of the pupil/otf (i.e. elements greater than some threshold)
            std::vector<std::pair<size_t,size_t>> pupilInOTF;                   //!< Maps the pupil-support into OTF-space (which is (2*nPixels,2*nPixels))
            std::mutex mtx;

        };


    }   // image

}   // redux


#endif  // REDUX_IMAGE_PUPIL_HPP
