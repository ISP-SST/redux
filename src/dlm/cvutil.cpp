#ifdef REDUX_WITH_OPENCV

#include "cvutil.hpp"

//#include "rdx.hpp"
/*
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/imgproc/imgproc.hpp>*/

#include <iostream>

using namespace cv;

/*
using namespace redux::file;
using namespace redux::image;
using namespace redux::util;

using namespace std;
namespace bfs = boost::filesystem;*/

Mat redux::arrayToMat (const IDL_VPTR& in, int verbose) {

    try {

        IDL_ENSURE_ARRAY (in);

        int nDims = in->value.arr->n_dim;
        std::vector<int> dims;
        std::copy (in->value.arr->dim, in->value.arr->dim + nDims, back_inserter (dims));
        std::reverse (dims.begin(), dims.end());

        int type = 0;

        switch (in->type) {
            case IDL_TYP_BYTE:   type |= CV_8U; break;
            case IDL_TYP_INT:    type |= CV_16S; break;
            case IDL_TYP_LONG:   type |= CV_32S; break;
            case IDL_TYP_FLOAT:  type |= CV_32F; break;
            case IDL_TYP_DOUBLE: type |= CV_64F; break;
            default:
                if (verbose) std::cerr << "Unsupported data type. " << (int) in->type << std::endl;
                return Mat();
        }

        Mat img (nDims, dims.data(), CV_MAKETYPE (type, 1), in->value.arr->data);

        return std::move (img);


    } catch (cv::Exception& e) {
        if (verbose) std::cerr << "OpenCV error: " << e.msg << std::endl;
    }


    return Mat();

}


// returns a float Mat with values scaled to [0, 1.0]
Mat redux::getImgAsGrayFloat (const IDL_VPTR& in, int verbose) {

    try {
        Mat img = arrayToMat (in, verbose);
        Mat img2;
        img.convertTo (img2, CV_32FC1);
        double minValue, maxValue;
        cv::minMaxLoc (img2, &minValue, &maxValue);
        if (verbose > 1) std::cout << "getImgAsGrayFloat:  minValue = " << minValue <<  "  maxValue = " << maxValue << std::endl;
        img2 = (img2 - minValue) / (maxValue - minValue);
        return std::move (img2);
    } catch (cv::Exception& e) {
        if (verbose) std::cerr << "OpenCV error: " << e.msg << std::endl;
    }

    return Mat();

}


#endif


