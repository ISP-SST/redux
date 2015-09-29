#ifndef REDUX_UTIL_ARRAYSTATS_HPP
#define REDUX_UTIL_ARRAYSTATS_HPP

#include "redux/util/array.hpp"
#include "redux/image/fouriertransform.hpp"

namespace redux {

    namespace util {

        enum StatType { ST_VALUES=1, ST_RMS, ST_NOISE=4, ST_ALL=7 };

        struct ArrayStats {
            typedef std::shared_ptr<ArrayStats> Ptr;

            ArrayStats() : clip(-1), cutoff(-1), min(0), max(0), median(0), sum(0), norm(0),
                           mean(0), rms(0), stddev(0), noise(0), noiseRMS(0) {}
            template <typename T> void getMinMaxMean(const T* data, size_t count);
            template <typename T> void getMinMaxMean(const redux::util::Array<T>& data) {
                if (data.dense()) getMinMaxMean(data.ptr(),data.nElements());
                else {
                    redux::util::Array<T> tmp;
                    data.copy(tmp);
                    getMinMaxMean(tmp.get(),data.nElements());
                }
            }

            template <typename T> void getRmsStddev(const T* data, size_t count);
            template <typename T> void getRmsStddev(const redux::util::Array<T>& data) {
                if (data.dense()) getRmsStddev(data.ptr(),data.nElements());
                else {
                    redux::util::Array<T> tmp;
                    data.copy(tmp);
                    getRmsStddev(tmp.get(),data.nElements());
                }
            }
            template <typename T> void getNoise(const redux::util::Array<T>& data, int smooth=0);
            void getNoise(const redux::image::FourierTransform& ft);
            
            template <typename T> void getStats(const T* data, size_t count, int flags=ST_ALL);
            template <typename T> void getStats(const redux::util::Array<T>& data, int flags=ST_ALL);
            template <typename T> void getStats(uint32_t borderClip, const redux::util::Array<T>& data, int flags=ST_ALL);
            
            static size_t size( void );
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );

            int clip;
            double cutoff;
            double min, max, median;
            double sum, norm, mean, rms, stddev;
            double noise, noiseRMS;

        };

    }   // image

}   // redux


#endif  // REDUX_UTIL_ARRAYSTATS_HPP
