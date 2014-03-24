#ifndef REDUX_IMAGE_STATISTICS_HPP
#define REDUX_IMAGE_STATISTICS_HPP

//#include "redux/image/image.hpp"
#include "redux/util/array.hpp"

namespace redux {

    namespace image {

        enum StatType { ST_VALUES=1, ST_RMS, ST_NOISE=4, ST_ALL=7 };

        template <typename T>
        struct Statistics  {
            typedef typename std::shared_ptr<Statistics> Ptr;

            void getStats(const redux::util::Array<T>& data, int flags=ST_ALL);
            void getStats(uint32_t borderClip, const redux::util::Array<T>& data, int flags=ST_ALL);
            
            T min, max, median;
            double mean, rms, stddev;
            double noisePower, noiseRMS;
            uint8_t noiseType;              // flag indicating noise statistics (not used atm.)

        };

    }   // image

}   // redux


#endif  // REDUX_IMAGE_STATISTICS_HPP
