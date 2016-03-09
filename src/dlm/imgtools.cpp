#include "imgtools.hpp"

#include <redux/util/stringutil.hpp>
#include <redux/util/opencv.hpp>


#include <stdio.h>
#include <iostream>
#include <set>


#ifdef REDUX_WITH_OPENCV
#    include <opencv2/core/core.hpp>
#    include <opencv2/features2d/features2d.hpp>
#    include <opencv2/calib3d/calib3d.hpp>
#    include <opencv2/imgproc/imgproc.hpp>
using namespace cv;
#endif



using namespace std;

using namespace redux::util;

namespace {

    static IDL_MEMINT dims3x3[] = { 3, 3 };

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT help;
        IDL_INT by_distance;
        float eps;
        IDL_INT margin;
        IDL_INT max_dist;
        float max_scale;
        IDL_INT max_points;
        IDL_INT niter;
        IDL_INT nrefpoints;
        IDL_INT show;
        float threshold;
        IDL_INT verbose;
        IDL_VPTR h_init;
        IDL_VPTR points;
    } KW_RESULT;

    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "BY_DIST",    IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (by_distance) },
        { (char*) "EPS",        IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (eps) },
        { (char*) "H_INIT",     IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (h_init) },
        { (char*) "HELP",       IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (help) },
        { (char*) "MARGIN",     IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (margin) },
        { (char*) "MAX_DIST",   IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (max_dist) },
        { (char*) "MAX_POINTS", IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (max_points) },
        { (char*) "MAX_SCALE",  IDL_TYP_FLOAT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (max_scale) },
        { (char*) "NITER",      IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (niter) },
        { (char*) "NREF",       IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (nrefpoints) },
        { (char*) "POINTS",     IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (points) },
        { (char*) "SHOW",       IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (show) },
        { (char*) "THRESHOLD",  IDL_TYP_FLOAT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (threshold) },
        { (char*) "VERBOSE",    IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (verbose) },
        { NULL }
    };


#ifdef REDUX_WITH_OPENCV
    Mat arrayToMat (const IDL_VPTR& in, int verbose = 0) {

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
                    if (verbose) cout << "Unsupported data type. " << (int) in->type << endl;
                    return Mat();
            }

            Mat img (nDims, dims.data(), CV_MAKETYPE (type, 1), in->value.arr->data);

            return std::move (img);


        } catch (cv::Exception& e) {
            if (verbose) cout << "OpenCV error: " << e.msg << endl;
        }


        return Mat();

    }


    Mat getImgAsGrayFloat (const IDL_VPTR& in, int verbose = 0) {

        try {
            Mat img = arrayToMat (in, verbose);
            Mat img2;
            img.convertTo (img2, CV_32FC1);
            double minValue, maxValue;
            cv::minMaxLoc (img2, &minValue, &maxValue);
            img2 = (img2 - minValue) / (maxValue - minValue);
            return std::move (img2);
        } catch (cv::Exception& e) {
            if (verbose) cout << "OpenCV error: " << e.msg << endl;
        }

        return Mat();

    }


    bool transformCheck (const Mat& trans, double maxValue = 10000, double maxScaleDiff = 1.5) {
        if (!checkRange (trans, true, 0, -maxValue, maxValue)) return false;

        double det = fabs (determinant (trans (cv::Rect_<int> (0, 0, 2, 2))));
        double minScale = std::max( 0.1, (1.0/maxScaleDiff) );
        if (det < minScale || det > maxScaleDiff ) {
            return false;
        }

        return true;
    }


    // find approximate (affine) transforms by matching 3 central points in the object, to a reference + its nNeighbours nearest points
    vector<Mat> getInitializations (const vector<KeyPoint>& obj, const vector<KeyPoint>& target,
                                    double maxValue = 10000, double maxScaleDiff=1.5) {

        size_t nPairs = 3;      // 3 points needed for getAffineTransform()

        vector<Mat> initializations;

        if (obj.size() < nPairs || target.size() < nPairs) {    // not enough input keypoints
            return std::move (initializations);
        }

        vector<size_t> objInd (obj.size());
        vector<size_t> targetInd (target.size());
        std::iota (objInd.begin(), objInd.end(), 0);
        std::iota (targetInd.begin(), targetInd.end(), 0);

        std::set< map<size_t, size_t> > combinations;

        do {
            std::sort (targetInd.begin(), targetInd.end());
            do {
                map<size_t, size_t> tmp;
                for (size_t i = 0; i < nPairs; ++i) {
                    tmp[objInd[i]] = targetInd[i];
                    if( i ) {       // distances should be in the scale range
                        double objDist = norm( obj[objInd[i]].pt - obj[objInd[0]].pt );
                        double ratio = norm( target[targetInd[i]].pt - target[targetInd[0]].pt ) / objDist;
                        if( ratio > maxScaleDiff || ratio < 1.0/maxScaleDiff ) {
                            continue;
                        }
                    }
                }
                if( tmp.size() != nPairs ) continue;
                combinations.insert (tmp);
            } while (next_permutation (targetInd.begin(), targetInd.end()));
        } while (next_permutation (objInd.begin(), objInd.end()));

        vector<Point2f> objPoints (nPairs);
        vector<Point2f> targetPoints (nPairs);

        for (auto & c : combinations) {
            size_t i = 0;
            for (auto & p : c) {
                objPoints[i] = obj[ p.first ].pt;
                targetPoints[i] = target[ p.second ].pt;
                i++;
            }
            Mat trans = getAffineTransform (objPoints, targetPoints);
            if (transformCheck (trans, maxValue, maxScaleDiff)) {
                Mat H_init = Mat::eye(3, 3, CV_32F);            // unit
                trans.copyTo(H_init(cv::Rect_<int> (0, 0, 3, 2))); // copy (affine) initialization to the top 2 rows
                initializations.push_back(H_init);
            }
        }

        return std::move (initializations);

    }

    double metric (const Mat& H, const vector<Point2f>& obj, const vector<Point2f>& scene, double* maxVal = nullptr) {
        double ret (0), mx (0);
        size_t nP = std::min(obj.size(),scene.size());
        vector<Point2f> mappedPoints (obj.size());
        perspectiveTransform (obj, mappedPoints, H);
        for (size_t i=0; i<nP; ++i) {
            double dist = norm( mappedPoints[i] - scene[i] );
            ret += dist;
            mx = std::max(mx, dist);
        }
        if (maxVal) *maxVal = mx;
        return ret / nP;
    }

    vector<DMatch> matchNearest (const vector<KeyPoint>& kp1, const vector<KeyPoint>& kp2, const Mat& H, double maxDistance = 30) {

        vector<DMatch> matches;
        vector<Point2f> points1, points2, mappedPoints;

        KeyPoint::convert (kp1, points1);
        KeyPoint::convert (kp2, points2);

        size_t nPoints = points1.size();
        mappedPoints.resize (nPoints);
        std::pair<int, double> nearest;
        perspectiveTransform (points1, mappedPoints, H);

        for (size_t i = 0; i < nPoints; ++i) {
            nearest = make_pair (-1, 1E12);

            for (size_t j = 0; j < points2.size(); ++j) {
                double dist = norm( mappedPoints[i] - points2[j]);
                if (fabs (dist) < nearest.second) {
                    nearest.second = dist;
                    nearest.first = j;
                }
            }

            if (nearest.first >= 0 && nearest.second < maxDistance) {
                matches.push_back (DMatch (i, nearest.first, nearest.second));
            }
        }

        nPoints = points2.size();
        mappedPoints.resize (nPoints);
        perspectiveTransform (points2, mappedPoints, H.inv());

        for (size_t i = 0; i < nPoints; ++i) {
            nearest = make_pair (-1, 1E12);

            for (size_t j = 0; j < points1.size(); ++j) {
                double dist = norm(mappedPoints[i] - points1[j]);
                if (fabs (dist) < nearest.second) {
                    nearest.second = dist;
                    nearest.first = j;
                }
            }

            if (nearest.first >= 0 && nearest.second < maxDistance) {
                matches.push_back (DMatch (nearest.first, i, nearest.second));
            }
        }

        std::sort (matches.begin(), matches.end(),
        [] (const DMatch & a, const DMatch & b) {
            if (a.queryIdx == b.queryIdx) return a.trainIdx < b.trainIdx;

            return a.queryIdx < b.queryIdx;
        });

        auto last = std::unique (matches.begin(), matches.end(),
        [] (const DMatch & a, const DMatch & b) { return (a.queryIdx == b.queryIdx && a.trainIdx == b.trainIdx); });

        matches.erase (last, matches.end());

        std::sort (matches.begin(), matches.end());

        return std::move (matches);

    }
#endif

}



IDL_VPTR redux::img_align (int argc, IDL_VPTR* argv, char* argk) {

#ifndef REDUX_WITH_OPENCV
    cerr << "img_align: redux has to be re-compiled with OpenCV enabled to be able to use this function." << endl;
    return IDL_GettmpInt(0);
#else
    IDL_VPTR ret;
    IDL_MakeTempArray (IDL_TYP_FLOAT, 2, dims3x3, IDL_ARR_INI_NOP, &ret);
    Mat retMat = arrayToMat (ret);

    if (argc < 2) {
        cerr << "img_align takes 2 images as input." << endl;
        return ret;
    }

    KW_RESULT kw;
    kw.by_distance = 0;
    kw.eps = 1E-3;
    kw.help = 0;
    kw.margin = 30;
    kw.max_dist = 80;
    kw.max_scale = 1.1;
    kw.max_points = -1;
    kw.niter = 100;
    kw.nrefpoints = 4;
    kw.show = 0;
    kw.verbose = 0;
    kw.threshold = 0;
    (void) IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    IDL_VPTR img1_in = argv[0];
    IDL_VPTR img2_in = argv[1];
    IDL_ENSURE_SIMPLE (img1_in);
    IDL_ENSURE_SIMPLE (img2_in);

    Mat raw1 = arrayToMat(img1_in, kw.verbose);
    Mat raw2 = arrayToMat(img2_in, kw.verbose);
    Mat imgFloat1 = getImgAsGrayFloat (img1_in, kw.verbose);
    Mat imgFloat2 = getImgAsGrayFloat (img2_in, kw.verbose);

    cv::SimpleBlobDetector::Params params;
    params.minThreshold = 50;
    params.thresholdStep = 10;
    params.maxThreshold = 240;
    params.filterByColor = false;
    params.filterByInertia = false;
    params.filterByConvexity = false;
    params.filterByArea = false;

    try {

        Mat result (imgFloat1.rows, imgFloat1.cols, CV_8UC3);
        Mat imgByte1 (imgFloat1.rows, imgFloat1.cols, CV_8UC1);
        Mat imgByte2 (imgFloat2.rows, imgFloat2.cols, CV_8UC1);
        
        Mat H, H_init;

        if (kw.threshold > 0.0 && kw.threshold < 1.0) {
            if (kw.verbose > 1) {
                cout << "Applying threshold " << kw.threshold << " to input images." << endl;
            }
            params.minThreshold = kw.threshold*255;
            // apply hard/noise threshold at threshold/4 to make the cross-correlation more distinct.
            threshold (imgFloat1, imgFloat1, kw.threshold/4, 0, THRESH_TOZERO);
            threshold (imgFloat2, imgFloat2, kw.threshold/4, 0, THRESH_TOZERO);
        }

        imgFloat1.convertTo (imgByte1, CV_8UC1, 255);
        imgFloat2.convertTo (imgByte2, CV_8UC1, 255);

        double maxTransformationValue = std::max(imgFloat1.cols, imgFloat1.rows) + 100;

        SimpleBlobDetector detector (params);
        vector<KeyPoint> keypoints1, keypoints2;
        
        detector.detect (imgByte1, keypoints1);
        detector.detect (imgByte2, keypoints2);
        Point2f last(imgFloat1.cols-kw.margin, imgFloat1.rows-kw.margin);
        keypoints1.erase( std::remove_if(keypoints1.begin(), keypoints1.end(),
                                         [&kw,&last](const KeyPoint& kp) {
                                             if ( std::isfinite( norm(kp.pt) ) ) {
                                                 if( (kp.pt.x > kw.margin) && (kp.pt.y > kw.margin) &&
                                                     (kp.pt.x < last.x) && (kp.pt.y < last.y)
                                                 ) return false;
                                             }
                                             return true;
                                        }),
                          keypoints1.end());
        
        keypoints2.erase( std::remove_if(keypoints2.begin(), keypoints2.end(),
                                         [&kw,&last](const KeyPoint& kp) {
                                             if ( std::isfinite( norm(kp.pt) ) ) {
                                                 if( (kp.pt.x > kw.margin) && (kp.pt.y > kw.margin) &&
                                                     (kp.pt.x < last.x) && (kp.pt.y < last.y)
                                                 ) return false;
                                             }
                                             return true;
                                        }),
                          keypoints2.end());

        if ( kw.verbose > 1 || keypoints1.size() < 3 || keypoints2.size() < 3 ) {
            cout << "Detected " << keypoints1.size() << " keypoints in image 1 and " << keypoints2.size() << " keypoints in image 2" << endl;
        }
        
        if( keypoints1.size() < 3 || keypoints2.size() < 3 ) {
            cout << "Not enough points for fitting." << endl;
            return ret;
        }
        
        std::sort (keypoints1.begin(), keypoints1.end(),
        [] (const KeyPoint & a, const KeyPoint & b) {
            return (a.size > b.size);
        });
        std::sort (keypoints2.begin(), keypoints2.end(),
        [] (const KeyPoint & a, const KeyPoint & b) {
            return (a.size > b.size);
        });
    
        size_t nLargest = std::min<size_t>(keypoints1.size(),3*kw.nrefpoints);
        nLargest = std::min(keypoints2.size(),nLargest);
        vector<KeyPoint> largestKP1 (keypoints1.begin(), keypoints1.begin() + nLargest);
        vector<KeyPoint> largestKP2 (keypoints2.begin(), keypoints2.begin() + nLargest);
      
        Point2f mid (imgFloat1.cols / 2, imgFloat1.rows / 2);
        std::sort (largestKP1.begin(), largestKP1.end(),
        [mid] (const KeyPoint & a, const KeyPoint & b) {
            return (norm(a.pt - mid) < norm(b.pt - mid));
        });
        std::sort (largestKP2.begin(), largestKP2.end(),
        [mid] (const KeyPoint & a, const KeyPoint & b) {
            return (norm(a.pt - mid) < norm(b.pt - mid));
        });

//         if( kw.show > 1 ) {
//             drawKeypoints(imgByte1, largestKP1, result, Scalar(0,0,255) );
//             imshow( "Good Matches & Object detection", result );
//             waitKey(0); 
//             drawKeypoints(imgByte2, largestKP2, result, Scalar(0,0,255) );
//             imshow( "Good Matches & Object detection", result );
//             waitKey(0); 
//         }
      
        largestKP1.resize(kw.nrefpoints);
        largestKP2.resize(kw.nrefpoints);
        
//         if( kw.show ) {
//             drawKeypoints(imgByte1, largestKP1, result, Scalar(0,0,255) );
//             imshow( "Good Matches & Object detection", result );
//             waitKey(0); 
//             drawKeypoints(imgByte2, largestKP2, result, Scalar(0,0,255) );
//             imshow( "Good Matches & Object detection", result );
//             waitKey(0); 
//         }
      
        std::vector<Mat> initializations;
        if( kw.h_init ) {
            if( kw.h_init->type == IDL_TYP_UNDEF ) {
                IDL_VPTR tmp;
                IDL_MakeTempArray (IDL_TYP_FLOAT, 2, dims3x3, IDL_ARR_INI_ZERO, &tmp);
                IDL_VarCopy( tmp, kw.h_init );
            }
            H_init = arrayToMat(kw.h_init);
            if( transformCheck( H_init, maxTransformationValue, kw.max_scale ) ) {
                initializations.push_back( H_init );
            }
        }

        if( initializations.empty() ) {
            initializations = getInitializations (largestKP1, largestKP2, maxTransformationValue, kw.max_scale);
        }
        if ( kw.verbose > 1 && initializations.size() > 1 ) {
            cout << "Matching the " << kw.nrefpoints << " largest keypoints in the images." << endl;
            cout << "Restricting the transformation to have shifts smaller than " << maxTransformationValue
                 << " and a maximal scaling of " << kw.max_scale << " gave " << initializations.size()
                 << " valid permutations." << endl;
        }

        std::map<double, Mat> results;

        for (auto & i : initializations) {
            try {
                double correlation = findTransformECC (imgByte1, imgByte2, i, MOTION_HOMOGRAPHY,
                                                       TermCriteria (TermCriteria::COUNT + TermCriteria::EPS,
                                                                     kw.niter, kw.eps));

                if (correlation > 0) {
                    results.insert (make_pair (correlation, i));
                    if (fabs (1.0 - correlation) < kw.eps) break;
                }
            } catch (cv::Exception& e) { }
        }

        H = results.rbegin()->second;

        vector< DMatch > matches = matchNearest (keypoints1, keypoints2, H, kw.max_dist);

        if (kw.verbose) {
            if (kw.verbose > 1) {
                vector<double> cc;
                for (auto r : results) cc.push_back (r.first);
                std::sort(cc.rbegin(),cc.rend());
                cout << printArray (cc, "correlations") << endl;
                cout << matches.size() << " pairs matched using a max_distance of " << kw.max_dist << " for pairing." << endl;
            }
        }

        if (kw.max_points < 0) kw.max_points = keypoints1.size() * 0.85;   // by default, do refinement by dropping the weakest 15%

        if (kw.max_points > 3) {     // at least 4 points needed for findHomography
            //vector< DMatch > refineMatches;
            if ( size_t(kw.max_points) < keypoints1.size() ) {
                if (kw.by_distance) {         // sort according to distance from image centre
                    Point2f mid (imgFloat1.cols / 2, imgFloat1.rows / 2);
                    std::sort (keypoints1.begin(), keypoints1.end(),
                    [mid] (const KeyPoint & a, const KeyPoint & b) {
                        return (norm(a.pt - mid) > norm(b.pt - mid));
                    });
                } 
                keypoints1.resize (kw.max_points);
                matches = matchNearest (keypoints1, keypoints2, H, kw.max_dist);
            } //else refineMatches = matches;

            vector<Point2f> obj, scene;

//             if( kw.show > 1 ) {
//                 drawKeypoints(imgByte1, keypoints1, result, Scalar(0,0,255) );
//                 imshow( "Good Matches & Object detection", result );
//                 waitKey(0); 
//             }
        
            if (kw.verbose > 1) {
                cout << "Using " << matches.size() << " pairs to refine the fit." << endl;
            }
            for (auto & m : matches) {
                obj.push_back (keypoints1[ m.queryIdx ].pt);
                scene.push_back (keypoints2[ m.trainIdx ].pt);
            }

            double max1, max2;
            double metric1 = metric (H, obj, scene, &max1);
            H = findHomography (obj, scene, CV_LMEDS, 0.5);       //  CV_RANSAC  CV_LMEDS  0
            double metric2 = metric (H, obj, scene, &max2);

            if (kw.verbose > 1) {
                cout << "   -> error (avg,max): (" << metric1 << "," << max1 << ") -> (" << metric2 << "," << max2 << ")" << endl;
            }
        }

        if (kw.verbose) {
            cout << "Transformation matrix:\n" << H << endl;
        }
        
        H.assignTo( retMat, retMat.type() );
        
        if( kw.h_init ) {
            H.assignTo( H_init, H_init.type() );
        }
        
        if( kw.points && matches.size()) {
            IDL_VPTR points;
            IDL_MEMINT dims[] = { 2, 2, static_cast<IDL_MEMINT>(matches.size()) };    // x/y, im#, nMatches 
            IDL_MakeTempArray( IDL_TYP_FLOAT, 3, dims, IDL_ARR_INI_ZERO, &points );
            IDL_VarCopy( points, kw.points );
            Mat pointsMat = arrayToMat( kw.points );
            int mCount(0);
            for ( auto & m : matches ) {
                pointsMat.at<float>(mCount, 0, 0 ) = keypoints1[ m.queryIdx ].pt.x;
                pointsMat.at<float>(mCount, 0, 1 ) = keypoints1[ m.queryIdx ].pt.y;
                pointsMat.at<float>(mCount, 1, 0 ) = keypoints2[ m.trainIdx ].pt.x;
                pointsMat.at<float>(mCount++, 1, 1 ) = keypoints2[ m.trainIdx ].pt.y;
            }
        }

        return ret;

    } catch (cv::Exception& e) {
        cout << "OpenCV error: " << e.msg << endl;
    }

    return ret;
    
#endif
    
}


IDL_VPTR redux::img_project (int argc, IDL_VPTR* argv, char* argk) {

#ifndef REDUX_WITH_OPENCV
    cerr << "img_project: redux has to be re-compiled with OpenCV enabled to be able to use this function." << endl;
    return IDL_GettmpInt(0);
#else
    if (argc != 2) {
        cerr << "img_project needs 2 arguments. A 3x3 or 2x3 transformation matrix, and an image." << endl;
        return IDL_GettmpInt (0);
    }

    KW_RESULT kw;
    kw.help = 0;
    kw.nrefpoints = 4;
    kw.verbose = 0;
    kw.threshold = 0.05;
    (void) IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    IDL_VPTR H_in = argv[0];
    IDL_VPTR img_in = argv[1];
    IDL_ENSURE_SIMPLE (img_in);
    IDL_ENSURE_ARRAY (img_in);
    IDL_ENSURE_SIMPLE (H_in);
    IDL_ENSURE_ARRAY (H_in);

    const Mat img = arrayToMat (img_in);
    const Mat H = arrayToMat (H_in);


    IDL_VPTR out;
    IDL_MakeTempArray (img_in->type, img_in->value.arr->n_dim,  img_in->value.arr->dim, IDL_ARR_INI_NOP, &out);

    Mat outImg = arrayToMat (out);
    img.copyTo (outImg);
    Size dsize = outImg.size();
    int flags = INTER_LINEAR;
    int borderMode = BORDER_CONSTANT;
    const Scalar borderValue = Scalar();

    if (H.cols == 3 && H.rows == 3) {
        if (kw.verbose) cout << "Applying perspective transform."  << endl;
        warpPerspective (img, outImg, H, dsize, flags, borderMode, borderValue);
    } else if (H.cols == 3 && H.rows == 2) {
        if (kw.verbose) cout << "Applying affine transform." << endl;
        warpAffine (img, outImg, H, dsize, flags, borderMode, borderValue);
    } else {
        cerr << "img_project: Transformation matrix has to have 2 or 3 rows." << endl;
        return IDL_GettmpInt (0);
    }

    return out;
    
#endif

}
