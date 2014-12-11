#ifndef REDUX_IMAGE_STATISTICS_HPP
#define REDUX_IMAGE_STATISTICS_HPP

//#include "redux/image/image.hpp"
#include "redux/util/array.hpp"
#include "redux/image/fouriertransform.hpp"

namespace redux {

    namespace image {

        enum StatType { ST_VALUES=1, ST_RMS, ST_NOISE=4, ST_ALL=7 };

        struct Statistics {
            typedef std::shared_ptr<Statistics> Ptr;

            Statistics() : clip(-1), cutoff(-1) {}
            template <typename T> void getMinMaxMean(const redux::util::Array<T>& data);
            template <typename T> void getRmsStddev(const redux::util::Array<T>& data, double mean);
            template <typename T> void getNoise(const redux::util::Array<T>& data);
            void getNoise(const redux::image::FourierTransform& ft);
            
            template <typename T> void getStats(const redux::util::Array<T>& data, int flags=ST_ALL);
            template <typename T> void getStats(uint32_t borderClip, const redux::util::Array<T>& data, int flags=ST_ALL);
            
            static size_t size( void );
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );

            int clip;
            double cutoff;
            double min, max, median;
            double mean, rms, stddev;
            double noise, noiseRMS;
            uint8_t noiseType;              // flag indicating noise statistics (not used atm.)

        };

    }   // image

}   // redux


#endif  // REDUX_IMAGE_STATISTICS_HPP
