#ifndef QIITEST_IMAGE_UTILS_HPP__
#define QIITEST_IMAGE_UTILS_HPP__

#include "redux/image/image.hpp"
#include "redux/imglib/imgorg.h"

namespace testsuite {

    namespace image {

        redux::image::Image allocImage(redux::imglib::imgType type, int w, int h);

        void fillWithRandomData(redux::image::Image & image, int min = 0, int max = 0);

        bool equivalent(const redux::image::Image & image1, const redux::image::Image & image2);

    }

}

#endif //QIITEST_IMAGE_UTILS_HPP__
