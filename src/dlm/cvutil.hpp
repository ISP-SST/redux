#ifndef REDUX_DLM_CVUTIL_HPP
#define REDUX_DLM_CVUTIL_HPP

#ifdef REDUX_WITH_OPENCV

#include <opencv2/core/core.hpp>

#include <stdio.h>          // has to be included before idl_export.h, otherwise FILE is undefined.
#include <idl_export.h>

namespace redux {

    // access array as a cv::Mat. NB: data is not copied.
    cv::Mat arrayToMat( const IDL_VPTR& array, int verbose=0 );

    // returns a float Mat with values scaled to [0, 1.0]
    cv::Mat getImgAsGrayFloat( const IDL_VPTR& array, int verbose=0 );

}

#endif  // REDUX_WITH_OPENCV

#endif  // REDUX_DLM_CVUTIL_HPP
