#include "imgtools.hpp"

#include "idlutil.hpp"

#include <redux/file/fileio.hpp>
#include <redux/file/fileana.hpp>
#include <redux/image/image.hpp>
#include <redux/image/descatter.hpp>
#include <redux/image/fouriertransform.hpp>
#include <redux/image/utils.hpp>
#include <redux/util/array.hpp>
#include <redux/util/stringutil.hpp>
#include <redux/util/arraystats.hpp>

#include <atomic>
#include <future>
#include <iostream>
#include <numeric>
#include <set>
#include <thread>

#include <boost/asio.hpp>
#include <boost/thread/thread.hpp>
#include <boost/filesystem.hpp>

#ifdef REDUX_WITH_OPENCV
#    include "cvutil.hpp"
#    include <redux/util/opencv.hpp>
#    include <opencv2/core/core.hpp>
#    include <opencv2/features2d/features2d.hpp>
#    include <opencv2/calib3d/calib3d.hpp>
#    include <opencv2/imgproc/imgproc.hpp>
using namespace cv;
#endif

using namespace redux::file;
using namespace redux::image;
using namespace redux::util;
using namespace redux;

using namespace std;
namespace bfs = boost::filesystem;


namespace {

    static char nl[] = "\n";

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT help;
        IDL_INT by_distance;
        float eps;
        IDL_INT margin;
        IDL_INT max_dist;
        IDL_INT max_shift;
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

    
#ifdef REDUX_WITH_OPENCV
    
    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_MEMINT dims3x3[] = { 3, 3 };
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "BY_DIST",    IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (by_distance) },
        { (char*) "EPS",        IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (eps) },
        { (char*) "HELP",       IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (help) },
        { (char*) "H_INIT",     IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (h_init) },
        { (char*) "MARGIN",     IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (margin) },
        { (char*) "MAX_DIST",   IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (max_dist) },
        { (char*) "MAX_POINTS", IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (max_points) },
        { (char*) "MAX_SCALE",  IDL_TYP_FLOAT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (max_scale) },
        { (char*) "MAX_SHIFT",  IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (max_shift) },
        { (char*) "NITER",      IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (niter) },
        { (char*) "NREF",       IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (nrefpoints) },
        { (char*) "POINTS",     IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (points) },
        { (char*) "SHOW",       IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (show) },
        { (char*) "THRESHOLD",  IDL_TYP_FLOAT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (threshold) },
        { (char*) "VERBOSE",    IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (verbose) },
        { NULL }
    };
    
    bool transformCheck (const Mat& trans, const Size& imgSize, double maxValue = 10000, double maxScaleDiff = 1.5) {

        if( (trans.at<double>(0,0) > 0 && trans.at<double>(0,2) > maxValue) ||
            (trans.at<double>(0,0) < 0 && abs(trans.at<double>(0,2)-imgSize.width) > maxValue)
        ) return false;
        if( (trans.at<double>(1,1) > 0 && trans.at<double>(1,2) > maxValue) ||
            (trans.at<double>(1,1) < 0 && abs(trans.at<double>(1,2)-imgSize.height) > maxValue)
        ) return false;

        double det = fabs (determinant (trans (cv::Rect_<int> (0, 0, 2, 2))));
        double minScale = std::max( 0.1, (1.0/maxScaleDiff) );
        if (det < minScale || det > maxScaleDiff ) {
            return false;
        }

        return true;
    }


    // find approximate (affine) transforms by matching 3 central points in the object, to a reference + its nNeighbours nearest points
    vector<Mat> getInitializations (const vector<KeyPoint>& obj, const vector<KeyPoint>& target, const Size& imgSize,
                                    double maxShift = 10000, double maxScaleDiff=1.5) {

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
                    Point2f oPoint = obj[objInd[i]].pt;
                    Point2f tPoint = target[targetInd[i]].pt;
                    double refDist = norm( oPoint - tPoint );
                    refDist = min( refDist, norm( oPoint - Point2f(tPoint.x, imgSize.height-tPoint.y) ));
                    refDist = min( refDist, norm( oPoint - Point2f(imgSize.width-tPoint.x, tPoint.y) ));
                    refDist = min( refDist, norm( oPoint - Point2f(imgSize.width-tPoint.x, imgSize.height-tPoint.y) ));
                    if( refDist > maxShift ) continue;      // image-shift is too large
                    if( i ) {       // relative distances should be in the scale range
                        double objDist = norm( obj[objInd[i]].pt - obj[objInd[0]].pt );
                        double ratio = norm( target[targetInd[i]].pt - target[targetInd[0]].pt ) / objDist;
                        if( ratio > maxScaleDiff || ratio < 1.0/maxScaleDiff ) {
                            continue;
                        }
                    }
                    if( tmp.size() != nPairs ) continue;
                    combinations.insert (tmp);
                }
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
            if (transformCheck (trans, imgSize, maxShift, maxScaleDiff)) {
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
    KW_RESULT kw;
    kw.by_distance = 0;
    kw.eps = 1E-3;
    kw.help = 0;
    kw.margin = 30;
    kw.max_dist = 80;
    kw.max_scale = 1.04;
    kw.max_points = -1;
    kw.max_shift = 200;
    kw.niter = 100;
    kw.nrefpoints = 4;
    kw.show = 0;
    kw.verbose = 0;
    kw.threshold = 0;
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    IDL_VPTR ret;
    IDL_MakeTempArray (IDL_TYP_FLOAT, 2, dims3x3, IDL_ARR_INI_NOP, &ret);
    Mat retMat = arrayToMat (ret);

    if (nPlainArgs < 2) {
        cerr << "img_align takes 2 images as input." << endl;
        return ret;
    }

    IDL_VPTR img1_in = argv[0];
    IDL_VPTR img2_in = argv[1];
    IDL_ENSURE_SIMPLE (img1_in);
    IDL_ENSURE_SIMPLE (img2_in);

    Mat raw1 = arrayToMat(img1_in, kw.verbose);
    Mat raw2 = arrayToMat(img2_in, kw.verbose);
    Mat imgFloat1 = getImgAsGrayFloat (img1_in, kw.verbose);
    Mat imgFloat2 = getImgAsGrayFloat (img2_in, kw.verbose);
    
    cv::Size imgSize = imgFloat1.size();
    if( imgSize != imgFloat2.size() ) {
        cerr << "img_align: input images must be of the same size." << endl;
        return ret;
    }

    cv::SimpleBlobDetector::Params params;
    params.minThreshold = 50;
    params.thresholdStep = 10;
    params.maxThreshold = 240;
    params.filterByColor = false;
    params.filterByInertia = false;
    params.filterByConvexity = false;
    params.filterByArea = false;

    try {

        Mat result (imgSize, CV_8UC3);
        Mat imgByte1 (imgSize, CV_8UC1);
        Mat imgByte2 (imgSize, CV_8UC1);
        
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
            Mat H_initd;
            H_init.assignTo(H_initd,CV_64F);
            if( transformCheck( H_initd, imgSize, kw.max_shift, kw.max_scale ) ) {
                initializations.push_back( H_init );
            }
        }

        if( initializations.empty() ) {
            initializations = getInitializations (largestKP1, largestKP2, imgSize, kw.max_shift, kw.max_scale);
        }
        if ( kw.verbose > 1 && initializations.size() > 1 ) {
            cout << "Matching the " << kw.nrefpoints << " largest keypoints in the images." << endl;
            cout << "Restricting the transformation to have shifts smaller than " << kw.max_shift
                 << " and a maximal scaling of " << kw.max_scale << " gave " << initializations.size()
                 << " valid permutations." << endl;
        }

        std::map<double, Mat> results;
        std::vector<std::future<double>> correlations;
        std::vector<double> correlations2(initializations.size());
        TermCriteria term_crit = TermCriteria(TermCriteria::COUNT + TermCriteria::EPS, kw.niter, kw.eps);
        auto fit_func = std::bind( findTransformECC, std::ref(imgByte1), std::ref(imgByte2), std::placeholders::_1,
            MOTION_HOMOGRAPHY, term_crit, noArray()
        );
        
        boost::asio::io_service service;
        boost::thread_group pool;
        std::mutex mtx;
        
        for ( size_t i=0; i<initializations.size(); ++i) {
            service.post([&,i](){
                try {
                    const Mat& init = initializations[i];
                    double correlation = fit_func( init );
                    if ( correlation > 0 ) {    // negative correlation means findTransformECC failed.
                        std::unique_lock<std::mutex> lock(mtx);
                        results.insert( make_pair(correlation, init) );
                    }
                } catch (cv::Exception& e) { }
            });
        }
        
        for(uint16_t t = 0; t < thread::hardware_concurrency(); ++t) {
            pool.create_thread(boost::bind(&boost::asio::io_service::run, &service));
        }
        pool.join_all();
        
        if( results.empty() ) {
            cout << "Failed to match detected pinholes, try to increase nref (current: "
                 << kw.nrefpoints << ")  or adjust threshold (current: "
                 << kw.threshold << ")." << endl;
            return ret;
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
            cout << "Scaling: " << fabs(determinant (H (cv::Rect_<int> (0, 0, 2, 2)))) << endl;
            cv::Point2f shift( H.at<double>(0,2), H.at<double>(1,2) );
            if( H.at<double>(0,0) < 0 ) {
                shift.x -= imgSize.width;
            }
            if( H.at<double>(1,1) < 0 ) {
                shift.y -= imgSize.height;
            }
            cout << "Shift: " << shift << endl;
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

    } catch( const cv::Exception& e ) {
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

    KW_RESULT kw;
    kw.help = 0;
    kw.nrefpoints = 4;
    kw.verbose = 0;
    kw.threshold = 0.05;
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if (nPlainArgs != 2) {
        cerr << "img_project needs 2 arguments. A 3x3 or 2x3 transformation matrix, and an image." << endl;
        return IDL_GettmpInt (0);
    }

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
    int flags = INTER_CUBIC;                  // INTER_LINEAR, INTER_CUBIC, INTER_AREA, INTER_LANCZOS4
    int borderMode = BORDER_CONSTANT;
    const Scalar borderValue = Scalar();

    try {
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
    } catch( const cv::Exception& e ) {
        cout << "OpenCV error: " << e.msg << endl;
    }

    return out;
    
#endif

}


IDL_VPTR redux::img_remap (int argc, IDL_VPTR* argv, char* argk) {

#ifndef REDUX_WITH_OPENCV
    cerr << "img_remap: redux has to be re-compiled with OpenCV enabled to be able to use this function." << endl;
    return IDL_GettmpInt(0);
#else
    KW_RESULT kw;
    kw.help = 0;
    kw.nrefpoints = 4;
    kw.verbose = 0;
    kw.threshold = 0.05;
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if (nPlainArgs != 3) {
        cerr << "img_remap needs 3 arguments. An image and 2 coordinate maps." << endl;
        return IDL_GettmpInt (0);
    }

    IDL_VPTR img_in = argv[0];
    IDL_ENSURE_SIMPLE( img_in );
    IDL_ENSURE_ARRAY( img_in );
    
    IDL_VPTR map_in = argv[1];
    IDL_ENSURE_SIMPLE( map_in );
    IDL_ENSURE_ARRAY( map_in );
    
    IDL_VPTR map_in2 = argv[2];
    IDL_ENSURE_SIMPLE( map_in2 );
    IDL_ENSURE_ARRAY( map_in2 );

    if( img_in-> value.arr->n_dim != 2 ) {
        cerr << "img_remap: only supports a single 2D image at the moment." << endl;
        return IDL_GettmpInt (0);
    }
    
    const Mat img = arrayToMat( img_in );
    const Mat cmap = arrayToMat( map_in );
    const Mat cmap2 = arrayToMat( map_in2 );

cerr << "img_remap: " << printArray(img_in->value.arr->dim,img_in->value.arr->n_dim,"img") << endl;
cerr << "img_remap: " << printArray(map_in->value.arr->dim,map_in->value.arr->n_dim,"cmap") << endl;
cerr << "img_remap: " << printArray(map_in->value.arr->dim,map_in->value.arr->n_dim,"cmap2") << endl;

    IDL_VPTR out;
    IDL_MakeTempArray( img_in->type, 2, map_in->value.arr->dim, IDL_ARR_INI_NOP, &out );
    Mat outImg = arrayToMat( out );
    
    //Size dsize = outImg.size();
    int interpolationMethod = INTER_LINEAR;
    int borderMode = BORDER_CONSTANT;
    const Scalar borderValue = Scalar();

    try {
        if( true ) {
            /*if (kw.verbose)*/ cout << "Remapping 1."  << endl;
            remap( img, outImg, cmap, cmap2, interpolationMethod, borderMode, borderValue);
        } else {
            /*if (kw.verbose)*/ cout << "Remapping 2." << endl;
            //warpAffine (img, outImg, cmap, dsize, interpolationMethod, borderMode, borderValue);
        }
    } catch( const cv::Exception& e ) {
        cout << "OpenCV error: " << e.msg << endl;
    }

    return out;
    
#endif

}


IDL_VPTR rdx_find_shift(int argc, IDL_VPTR* argv, char* argk) {

#ifndef REDUX_WITH_OPENCV
    cerr << "rdx_find_shift: redux has to be re-compiled with OpenCV enabled to be able to use this function." << endl;
    return IDL_GettmpInt(0);
#else
    
    IDL_VPTR ret;
    IDL_MEMINT dims[] = {3,2};
    IDL_MakeTempArray( IDL_TYP_FLOAT, 2, dims, IDL_ARR_INI_NOP, &ret );
    Mat retMat = arrayToMat( ret );

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

    IDL_VPTR ref = argv[0];
    IDL_ENSURE_SIMPLE(ref);
    IDL_ENSURE_ARRAY(ref);
    IDL_VPTR shifted = argv[1];
    IDL_ENSURE_SIMPLE (shifted);
    IDL_ENSURE_ARRAY(shifted);
    
    Mat grayRef = getImgAsGrayFloat( ref, kw.verbose );
    Mat grayShifted = getImgAsGrayFloat( shifted, kw.verbose );

    try {
        
        const int warp_mode = MOTION_TRANSLATION; // MOTION_EUCLIDEAN;
 
        Mat warp_matrix = Mat::eye(2, 3, CV_32F);

        int number_of_iterations = 5000;
        double termination_eps = 1e-10;
        TermCriteria criteria(TermCriteria::COUNT+TermCriteria::EPS, number_of_iterations, termination_eps);

        findTransformECC( grayRef, grayShifted, warp_matrix, warp_mode, criteria );
 
        warp_matrix.assignTo( retMat, retMat.type() );

    } catch( const cv::Exception& e ) {
        cout << "OpenCV error: " << e.msg << endl;
    }

    return ret;
    
#endif
    
}

namespace {
    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD;
        IDL_INT help;
        UCHAR nthreads;
        IDL_INT verbose;
        float thres;
        IDL_VPTR mask;
    } KW_FILLPIX;


    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR kw_fillpix_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "HELP",           IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_FILLPIX,help) },
        { (char*) "MASK",           IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_FILLPIX,mask) },
        { (char*) "NTHREADS",       IDL_TYP_BYTE,  1, 0,                      0, (char*) IDL_KW_OFFSETOF2(KW_FILLPIX,nthreads) },
        { (char*) "THRESHOLD",      IDL_TYP_FLOAT, 1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_FILLPIX,thres) },
        { (char*) "VERBOSE",        IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_FILLPIX,verbose) },
        { NULL }
    };
}

string rdx_fillpix_info( int lvl ) {
    string ret = "RDX_FILLPIX";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"        ");          // newline if lvl>1
        ret += "   Syntax:   out = rdx_fillpix(image(s), /KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      MASK                Mask with non-zero values where filling is supposed to be performed.\n"
                    "      NTHREADS            Number of threads.\n"
                    "      THRESHOLD           Value below which a pixel qualifies for filling. (default: 0)\n"
                    "      VERBOSE             Verbosity, default is 0 (only error output).\n";
        }
    } else ret += "\n";
    return ret;
}



IDL_VPTR rdx_fillpix( int argc, IDL_VPTR* argv, char* argk ) {

    IDL_VPTR images = argv[0];
    IDL_ENSURE_SIMPLE(images);
    IDL_ENSURE_ARRAY(images);
    
    shared_ptr<uint8_t> maskData;
    shared_ptr<uint8_t*> mask2D;
    
    KW_FILLPIX kw;
    kw.nthreads = std::thread::hardware_concurrency();
    kw.mask = nullptr;
    (void) IDL_KWProcessByOffset( argc, argv, argk, kw_fillpix_pars, (IDL_VPTR*)0, 255, &kw );

    if( kw.help ) {
        cout << rdx_fillpix_info(2) << endl;
        return IDL_GettmpInt(0);
    }
        
    kw.nthreads = max<UCHAR>(1, min<UCHAR>(kw.nthreads, thread::hardware_concurrency()));
   
    //int nDims = images->value.arr->n_dim;
    /*size_t nImages = 1;
    if( nDims == 3 ) {
        nImages = images->value.arr->dim[2];
    }*/
    IDL_MEMINT xSize = images->value.arr->dim[1];
    IDL_MEMINT ySize = images->value.arr->dim[0];
    size_t nPixels = xSize*ySize;
    //size_t frameSize = nPixels*images->value.arr->elt_len;

    if( kw.mask ) {
        IDL_ENSURE_SIMPLE(kw.mask);
        IDL_ENSURE_ARRAY(kw.mask);
        if( kw.mask->value.arr->n_dim != 2 ) {
            cout << "rdx_fillpix: mask has to be a 2D array." << endl;
            return IDL_GettmpInt(0);
        }
        if( (kw.mask->value.arr->dim[0] != ySize) || (kw.mask->value.arr->dim[1] != xSize) ) {
            cout << "rdx_fillpix: mask size does not match image(s)." << endl;
            return IDL_GettmpInt(0);
        }
        maskData = castOrCopy<uint8_t>( kw.mask );
        mask2D = reshapeArray( maskData.get(), ySize, xSize );
    } else {
        mask2D = sharedArray<uint8_t>( ySize, xSize );
        uint8_t* mPtr = *(mask2D.get());
        memset( mPtr, 0, nPixels );
        auto imData = castOrCopy<double>( images );
        double* dPtr = imData.get();
        for( size_t n=0; n<nPixels; ++n ) {
            if( dPtr[n] <= kw.thres ) mPtr[n] = 1;
        }
    }
    
    UCHAR dataType = images->type;
    IDL_VPTR ret;
    char* retData = IDL_MakeTempArray( dataType, images->value.arr->n_dim, images->value.arr->dim, IDL_ARR_INI_NOP, &ret );
    memcpy( retData, images->value.arr->data, images->value.arr->n_elts*images->value.arr->elt_len );

    
    try {
            
        switch( dataType ) {
            case( IDL_TYP_BYTE ): {
                auto data = reinterpret_cast<UCHAR*>( retData );
                auto images2D = reshapeArray( data, ySize, xSize );
                fillPixels( images2D.get(), ySize, xSize, mask2D.get() );
                break;
            }
            case( IDL_TYP_INT ): {
                auto data = reinterpret_cast<IDL_INT*>( retData );
                auto images2D = reshapeArray( data, ySize, xSize );
                fillPixels( images2D.get(), ySize, xSize, mask2D.get() );
                break;
            }
            case( IDL_TYP_LONG ): {
                auto data = reinterpret_cast<IDL_LONG*>( retData );
                auto images2D = reshapeArray( data, ySize, xSize );
                fillPixels( images2D.get(), ySize, xSize, mask2D.get() );
                break;
            }
            case( IDL_TYP_FLOAT ): {
                auto data = reinterpret_cast<float*>( retData );
                auto images2D = reshapeArray( data, ySize, xSize );
                fillPixels( images2D.get(), ySize, xSize, mask2D.get() );
                break;
            }
            case( IDL_TYP_DOUBLE ): {
                auto data = reinterpret_cast<double*>( retData );
                auto images2D = reshapeArray( data, ySize, xSize );
                fillPixels( images2D.get(), ySize, xSize, mask2D.get() );
                break;
            }
            default: ;
        }

    } catch( const exception& e ) {
        cout << "rdx_fillpix: unhandled exception: " << e.what() << endl;
    }
    
    return ret;
    
}



namespace {
    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT help;
        UCHAR nthreads;
        IDL_VPTR times;
        IDL_INT verbose;
    } LOAD_FILES_KW;


    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR load_files_kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "HELP",             IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(LOAD_FILES_KW,help) },
        { (char*) "NTHREADS",         IDL_TYP_BYTE,  1, 0,                      0, (char*) IDL_KW_OFFSETOF2(LOAD_FILES_KW,nthreads) },
        { (char*) "TIMES",            IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(LOAD_FILES_KW,times) },
        { (char*) "VERBOSE",          IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(LOAD_FILES_KW,verbose) },
        { NULL }
    };
}

string load_files_info( int lvl ) {
    string ret = "RDX_LOADFILES";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"      ");          // newline if lvl>1
        ret += "   Syntax:   out = rdx_loadfiles(file_list, /KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      NTHREADS            Number of threads.\n"
                    "      TIMES               (output) Array with timestamps from headers.\n"
                    "      VERBOSE             Verbosity, default is 0 (only error output).\n";
        }
    } else ret += "\n";
    return ret;
}


IDL_VPTR load_files( int argc, IDL_VPTR* argv, char* argk ) {
    
    LOAD_FILES_KW kw;
    kw.nthreads = std::thread::hardware_concurrency();
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, load_files_kw_pars, (IDL_VPTR*)0, 255, &kw);
        
    if( nPlainArgs < 1 ) {
        return IDL_GettmpInt(0);
    }
    
    IDL_VPTR filelist = argv[0];
    IDL_VPTR ret;
    IDL_ENSURE_STRING( filelist )
    IDL_ENSURE_SIMPLE( filelist );
    
    try {
        
        vector<string> filenames;
        if ( !(filelist->flags & IDL_V_ARR) ) {
            bfs::path fn( string(filelist->value.str.s) );
            if( bfs::is_regular_file(fn) ) {
                filenames.push_back( fn.string() );
            } else return IDL_GettmpInt(0);
        } else {
            IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>(filelist->value.arr->data);
            for( int i=0; i<filelist->value.arr->n_elts; ++i) {
                bfs::path fn( string(strptr[i].s) );
                if( bfs::is_regular_file(fn) ) {
                    filenames.push_back( fn.string() );
                }
            }
        }
        
        size_t nImages = filenames.size();
        if( nImages == 0 ) {
            return IDL_GettmpInt(0);
        }

        kw.nthreads = max<UCHAR>(1, min<UCHAR>(kw.nthreads, thread::hardware_concurrency()));

        if( kw.help ) {
            cout << load_files_info(2) << endl;
            return IDL_GettmpInt(0);
        }
        
        // Get size etc from first image.
        std::shared_ptr<redux::file::FileMeta> meta = redux::file::getMeta( filenames[0] );
        if( !meta || (meta->nDims() != 2) ) {        // Only allow images for now
            cout << "rdx_loadfiles: Failed to get meta, or not 2D input." << endl;
            return IDL_GettmpInt(0);
        }
        
        string statusString;
        if( kw.verbose ) {
            statusString = "Loading " + to_string(nImages)
            + " files using " +to_string((int)kw.nthreads) + string(" thread") + ((kw.nthreads>1)?"s.":".");
            cout << statusString << ((kw.verbose == 1)?"\n":"") << flush;
        }
        
        unique_ptr<double[]> times;
        double* timesPtr = nullptr;
        if( kw.times && nImages > 1 ) {
            times.reset( new double[nImages] );
            timesPtr = times.get();
        }
    
        IDL_MEMINT dims[] = { static_cast<IDL_MEMINT>(meta->dimSize(1)),
                              static_cast<IDL_MEMINT>(meta->dimSize(0)),
                              static_cast<IDL_MEMINT>(nImages) };

        char* rawPtr;
        if( nImages == 1 ) {
            rawPtr = IDL_MakeTempArray( meta->getIDLType(), 2, dims, IDL_ARR_INI_NOP, &ret );
            readFile( filenames[0], rawPtr, meta );
            if( kw.times ) {
                IDL_ALLTYPES tmp;
                tmp.d = meta->getAverageTime().time_of_day().total_nanoseconds()*1E-9;
                IDL_StoreScalar( kw.times, IDL_TYP_DOUBLE, &tmp );
            }
            return ret;
        } else {
            if( kw.times ) {
                IDL_VPTR tmp;
                (void)IDL_MakeTempArray( IDL_TYP_DOUBLE, 1, dims+2, IDL_ARR_INI_ZERO, &tmp );
                IDL_VarCopy( tmp, kw.times );
                timesPtr = (double*)kw.times->value.arr->data;
            }
            rawPtr = IDL_MakeTempArray( meta->getIDLType(), 3, dims, IDL_ARR_INI_NOP, &ret );
        }
        
        size_t frameSize = meta->dataSize();
        
        atomic<size_t> nLoaded(0);
        postLoadCallback postLoad = [=,&nLoaded]( char* data, size_t i, std::shared_ptr<redux::file::FileMeta>& meta ) {
            if( timesPtr ) timesPtr[ i ] = meta->getAverageTime().time_of_day().total_nanoseconds()*1E-9;
            size_t nl = nLoaded++;
            if( kw.verbose > 1 ) printProgress( statusString, (nl*100.0/(nImages-1)));
        };
        
        loadFiles( filenames, rawPtr, frameSize, kw.nthreads, postLoad );

        if( kw.verbose > 1 ) {
            cout << endl;
        }

    } catch (const exception& e ) {
        cout << "rdx_loadfiles: unhandled exception: " << e.what() << endl;
        return IDL_GettmpInt(0);
    }
    
    return ret;    

}


typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    IDL_INT preserve_type;
    IDL_INT smooth;
    IDL_INT larger;
    IDL_INT invert;
    float threshold;
    IDL_INT verbose;
    //IDL_VPTR time;
   // IDL_STRING split_chars;
} MAKE_MASK_KW;


// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR make_mask_kw_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "HELP",           IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MAKE_MASK_KW,help) },
    { (char*) "INVERT",         IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MAKE_MASK_KW,invert) },
    { (char*) "LARGER",         IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MAKE_MASK_KW,larger) },
    { (char*) "PRESERVE_TYPE",  IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MAKE_MASK_KW,preserve_type) },
    { (char*) "SMOOTH",         IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MAKE_MASK_KW,smooth) },
    { (char*) "THRESHOLD",      IDL_TYP_FLOAT, 1,           0, 0, (char*) IDL_KW_OFFSETOF2(MAKE_MASK_KW,threshold) },
    { (char*) "VERBOSE",        IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MAKE_MASK_KW,verbose) },
    { NULL }
};


string make_mask_info( int lvl ) {
    string ret = "RDX_MAKE_MASK";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"      ");          // newline if lvl>1
        ret += "   Syntax:   mask = rdx_make_mask(input,/KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      INVERT              Invert resulting mask.\n"
                    "      LARGER              Filter away features larger than SMOOTH (default is to filter away smaller features).\n"
                    "      PRESERVE_TYPE       Make output the same type as input instead of BYTE.\n"
                    "      SMOOTH              Smoothing range (default = 0, used as parameter for morph_open/morph_close).\n"
                    "      THRESHOLD           Level of initial thresholding of the input. (default=0.01  [0.0-1.0])\n"
                    "      VERBOSE             Verbosity, default is 0 (only error output).\n";
        }
    } else ret += "\n";
    return ret;
}


IDL_VPTR rdx_make_mask( int argc, IDL_VPTR* argv, char* argk ) {
    
    MAKE_MASK_KW kw;
    kw.threshold = 0;
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, make_mask_kw_pars, (IDL_VPTR*)0, 255, &kw );

    if( nPlainArgs < 1 ) {
        return IDL_GettmpInt(0);
    }
    
    if( kw.help ) {
        cout << make_mask_info(2) << endl;
        return IDL_GettmpInt(0);
    }
    
#ifndef REDUX_WITH_OPENCV
    cerr << "rdx_make_mask: redux has to be re-compiled with OpenCV enabled to be able to use this function." << endl;
    return IDL_GettmpInt(0);
#else

    try {
        
        IDL_VPTR input = argv[0];
        IDL_ENSURE_ARRAY( input );
        IDL_ENSURE_SIMPLE( input );
        
        IDL_VPTR mask;
        if ( kw.preserve_type ) {
            IDL_MakeTempArray( input->type, input->value.arr->n_dim,  input->value.arr->dim, IDL_ARR_INI_NOP, &mask );
        } else {
            IDL_MakeTempArray( IDL_TYP_BYTE, input->value.arr->n_dim,  input->value.arr->dim, IDL_ARR_INI_NOP, &mask );
        }
        
        Mat inputMat = redux::arrayToMat( input );
        Mat maskMat = redux::arrayToMat( mask );
        
        cv::make_mask( inputMat, maskMat, kw.threshold, kw.smooth, kw.larger, kw.invert );
            
        return mask;
        
    } catch (const exception& e ) {
        cout << "rdx_make_mask: unhandled exception: " << e.what() << endl;
        return IDL_GettmpInt(0);
    }
    
    
#endif

}


typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    IDL_INT verbose;
} MAKE_WIN_KW;


// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR make_win_kw_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "HELP",           IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MAKE_WIN_KW,help) },
    { (char*) "VERBOSE",        IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MAKE_WIN_KW,verbose) },
    { NULL }
};


string make_win_info( int lvl ) {
    string ret = "RDX_MAKE_WINDOW";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"      ");          // newline if lvl>1
        ret += "   Syntax:   win = rdx_make_window(img_size, blend_region, /KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      VERBOSE             Verbosity, default is 0 (only error output).\n";
        }
    } else ret += "\n";
    return ret;
}


IDL_VPTR rdx_make_win( int argc, IDL_VPTR* argv, char* argk ) {
    
    MAKE_WIN_KW kw;
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, make_win_kw_pars, (IDL_VPTR*)0, 255, &kw );

    if( kw.help ) {
        cout << make_win_info(2) << endl;
        return IDL_GettmpInt(0);
    }
    
    if( nPlainArgs < 2 ) {
        cout << "rdx_make_window: needs 2 arguments: nPixels & pupil-radius (in pixels). " << endl;
        return IDL_GettmpInt (0);
    }
    
    IDL_LONG nPixels = IDL_LongScalar(argv[0]);
    IDL_LONG nBlend = IDL_LongScalar(argv[1]);
    IDL_VPTR tmp;
 
    Array<float> win(nPixels, nPixels);
    win = 1.0;
    redux::image::apodizeInPlace( win, nBlend );
        
    IDL_MEMINT dims[] = { nPixels, nPixels };
    float* tmpData = (float*)IDL_MakeTempArray( IDL_TYP_FLOAT, 2, dims, IDL_ARR_INI_NOP, &tmp );
    win.copyTo<float>(tmpData);
    
    return tmp;

}

namespace {
    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT check;
        IDL_INT help;
        UCHAR nthreads;
        IDL_INT padding;
        IDL_INT pinh_align;
        IDL_INT verbose;
        IDL_INT lun;
        float fp_thres;
        float limit;
        IDL_VPTR summed;
        IDL_VPTR nsummed;
        IDL_VPTR dark;
        IDL_VPTR gain;
        IDL_VPTR bs_gain;
        IDL_VPTR bs_psf;
        IDL_VPTR time;
        IDL_VPTR xyc;
       // IDL_STRING split_chars;
    } SI_KW;
    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR si_kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "BACKSCATTER_GAIN", IDL_TYP_UNDEF, 1, IDL_KW_VIN|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,bs_gain) },
        { (char*) "BACKSCATTER_PSF",  IDL_TYP_UNDEF, 1, IDL_KW_VIN|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,bs_psf) },
        { (char*) "CHECK",            IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(SI_KW,check) },
        { (char*) "DARK",             IDL_TYP_UNDEF, 1, IDL_KW_VIN|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,dark) },
        { (char*) "FILLPIX_THRESHOLD",IDL_TYP_FLOAT, 1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(SI_KW,fp_thres) },
        { (char*) "GAIN",             IDL_TYP_UNDEF, 1, IDL_KW_VIN|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,gain) },
        { (char*) "HELP",             IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(SI_KW,help) },
        { (char*) "LIMIT",            IDL_TYP_FLOAT, 1, 0,                      0, (char*) IDL_KW_OFFSETOF2(SI_KW,limit) },
        { (char*) "LUN",              IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(SI_KW,lun) },
        { (char*) "NSUMMED",          IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,nsummed) },
        { (char*) "NTHREADS",         IDL_TYP_BYTE,  1, 0,                      0, (char*) IDL_KW_OFFSETOF2(SI_KW,nthreads) },
        { (char*) "PADDING",          IDL_TYP_INT,   1, 0,                      0, (char*) IDL_KW_OFFSETOF2(SI_KW,padding) },
        { (char*) "PINHOLE_ALIGN",    IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(SI_KW,pinh_align) },
        { (char*) "SUMMED",           IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,summed) },
        { (char*) "TIME",             IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,time) },
        { (char*) "VERBOSE",          IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(SI_KW,verbose) },
        { (char*) "XYC",              IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,xyc) },
       // { (char*) "SPLIT_CHARS", IDL_TYP_STRING, 1, 0, 0, (char*)IDL_KW_OFFSETOF2(SI_KW,split_chars) },
        { NULL }
    };
}

string sum_images_info( int lvl ) {
    string ret = "RDX_SUMIMAGES";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"      ");          // newline if lvl>1
        ret += "   Syntax:   out = rdx_sumimages(img_array, /KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      BACKSCATTER_GAIN    Image containing the imprinted background pattern. (for descattering).\n"
                    "      BACKSCATTER_PSF     Image containing the scattering PSF. (for descattering).\n"
                    "      CHECK               Discard images which are statistically deviant.\n"
                    "      DARK                Dark field.\n"
                    "      GAIN                Gain table (inverted flat-field).\n"
                    "      LIMIT               Allowed deviation from the median. (0.0175)\n"
                    "      NSUMMED             (output) Number of images actually summed.\n"
                    "      NTHREADS            Number of threads.\n"
                    "      PADDING             Padding size for the descattering procedure. (256)\n"
                    "      PINHOLE_ALIGN       Do sub-pixel alignment before summing.\n"
                    "      SUMMED              (output) Raw sum.\n"
                    "      VERBOSE             Verbosity, default is 0 (only error output).\n"
                    "      XYC                 (output) Coordinates of align-feature and image-shifts.\n";
        }
    } else ret += "\n";
    return ret;
}

IDL_VPTR sum_files( int argc, IDL_VPTR* argv, char* argk );

IDL_VPTR sum_images( int argc, IDL_VPTR* argv, char* argk ) {
    
    IDL_VPTR images = argv[0];
    IDL_ENSURE_SIMPLE( images );
    IDL_ENSURE_ARRAY( images );
    
    if( images->type == IDL_TYP_STRING ) return sum_files( argc, argv, argk );  // looks like a file-list.

    SI_KW kw;
    kw.dark = nullptr;
    kw.gain = nullptr;
    kw.nsummed = nullptr;
    kw.summed = nullptr;
    kw.bs_gain = nullptr;
    kw.bs_psf = nullptr;
    kw.limit = 0.0175;
    kw.padding = 256;
    kw.nthreads = std::thread::hardware_concurrency();
    (void)IDL_KWProcessByOffset( argc, argv, argk, si_kw_pars, (IDL_VPTR*)0, 255, &kw );
    
    if( kw.help ) {
        cout << sum_images_info(2) << endl;
        return IDL_GettmpInt(0);
    }
    
    kw.nthreads = max<UCHAR>(1, min<UCHAR>(kw.nthreads, thread::hardware_concurrency()));
    kw.padding = max<IDL_INT>(0, min<IDL_INT>(kw.padding, 4096));   // prevent insane padding
    
    try {
        
        IDL_VPTR ret;
        shared_ptr<double> darkData, gainData, bsGainData;
        shared_ptr<uint8_t> maskData;
        shared_ptr<uint8_t*> mask2D;
        shared_ptr<fftw_complex> bsOtfData;
        double *darkPtr = nullptr;
        double *gainPtr = nullptr;
        double *bsGainPtr = nullptr;
        fftw_complex *bsOtfPtr = nullptr;
        IDL_MEMINT xCalibSize(0);
        IDL_MEMINT yCalibSize(0);
        
        if( kw.dark ) {
            IDL_ENSURE_SIMPLE( kw.dark );
            IDL_ENSURE_ARRAY( kw.dark );
            if( kw.dark->value.arr->n_dim == 2 ) {
                darkData = castOrCopy<double>( kw.dark );
                darkPtr = darkData.get();
                xCalibSize = kw.dark->value.arr->dim[1];
                yCalibSize = kw.dark->value.arr->dim[0];
            } else cout << "dark must be a 2D image." << endl;
        }
        if( kw.gain ) {
            IDL_ENSURE_SIMPLE( kw.gain );
            IDL_ENSURE_ARRAY( kw.gain );
            if( kw.gain->value.arr->n_dim == 2 ) {
                if( (xCalibSize && (xCalibSize != kw.gain->value.arr->dim[1])) ||
                    (yCalibSize && (yCalibSize != kw.gain->value.arr->dim[0])) ) {
                    cout << "Dark & Gain files have different sizes." << endl;
                    return IDL_GettmpInt(0);
                }
                gainData = castOrCopy<double>( kw.gain  );
                gainPtr = gainData.get();
                xCalibSize = kw.gain->value.arr->dim[1];
                yCalibSize = kw.gain->value.arr->dim[0];
                maskData.reset( new uint8_t[xCalibSize*yCalibSize] );
                mask2D = reshapeArray( maskData.get(), yCalibSize, xCalibSize );
                redux::image::make_mask( gainPtr, maskData.get(), yCalibSize, xCalibSize, 0, 5, true, true ); // filter away larger features than ~5 pixels and invert
            } else cout << "gain must be a 2D image." << endl;
        }
        if( kw.bs_gain && kw.bs_psf ) {
            IDL_ENSURE_SIMPLE( kw.bs_gain );
            IDL_ENSURE_ARRAY( kw.bs_gain );
            IDL_ENSURE_SIMPLE( kw.bs_psf );
            IDL_ENSURE_ARRAY( kw.bs_psf );
            if( darkData && gainData ) {
                if( kw.bs_gain->value.arr->n_dim==2 && kw.bs_psf->value.arr->n_dim==2 ) {
                    IDL_MEMINT bsXsize = kw.bs_gain->value.arr->dim[1];
                    IDL_MEMINT bsYsize = kw.bs_gain->value.arr->dim[0];
                    if( bsXsize != kw.bs_psf->value.arr->dim[1] || bsYsize != kw.bs_psf->value.arr->dim[0]) {
                        cout << "Backscatter Gain & PSF have different sizes." << endl;
                        return IDL_GettmpInt(0);
                    }
                    if( (xCalibSize && (xCalibSize != bsXsize)) ||
                        (yCalibSize && (yCalibSize != bsYsize)) ) {
                        cout << "Backscatter Gain & PSF have different size than the dark/gain." << endl;
                        return IDL_GettmpInt(0);
                    }
                    size_t xSizePadded = bsXsize + 2*kw.padding;
                    size_t ySizePadded = bsYsize + 2*kw.padding;
                    size_t dataSize = xSizePadded*ySizePadded;

                    //bsGainData.reset( fftw_alloc_real(dataSize), fftw_free );
                    bsGainData.reset( (double*)fftw_malloc(dataSize*sizeof(double)), fftw_free );
                    
                    bsGainPtr = bsGainData.get();
                    memset( bsGainPtr, 0, dataSize*sizeof(double) );
                    copyInto( kw.bs_gain, bsGainPtr, ySizePadded, xSizePadded, kw.padding, kw.padding );
                    
                    shared_ptr<double> bsPsfData( (double*)fftw_malloc(dataSize*sizeof(double)), fftw_free );
                    double* psfPtr = bsPsfData.get();
                    memset( psfPtr, 0, dataSize*sizeof(double) );
                    copyInto( kw.bs_psf, psfPtr, ySizePadded, xSizePadded, kw.padding, kw.padding );

                    // normalize psf
                    double sum = 0;
                    for( size_t i=0; i<dataSize; ++i ) sum += psfPtr[i];
                    sum = 1.0/(sum*dataSize);
                    for( size_t i=0; i<dataSize; ++i ) psfPtr[i] *= sum;
                    
                    FourierTransform::reorder( psfPtr, ySizePadded, xSizePadded );
                    
                    FourierTransform::Plan::Ptr plan = FourierTransform::Plan::get( ySizePadded, xSizePadded, FourierTransform::Plan::R2C, 1 );
                    //bsOtfData.reset( fftw_alloc_complex( ySizePadded*(xSizePadded/2+1) ), fftw_free );
                    bsOtfData.reset( (fftw_complex*)fftw_malloc(ySizePadded*(xSizePadded/2+1)*sizeof(fftw_complex)), fftw_free );
                    bsOtfPtr = bsOtfData.get();
                    plan->forward( bsPsfData.get(), bsOtfPtr );
                    
                } else cout << "backscatter_gain/psf must be a 2D images." << endl;
            } else cout << "Backscatter correction requires both dark & gain to be present." << endl;
        } else if( kw.bs_gain || kw.bs_psf ) cout << "Both backscatter_gain and backscatter_psf must be specified." << endl;
        
        size_t nImages(0);
        IDL_MEMINT xSize(0);
        IDL_MEMINT ySize(0);
        size_t nPixels(0);
        size_t frameSize(0);
        
        UCHAR dataType = images->type;
        UCHAR* dataPtr = images->value.arr->data;
        int nDims = images->value.arr->n_dim;
        if( nDims == 3 ) {
            nImages = images->value.arr->dim[2];
            xSize = images->value.arr->dim[1];
            ySize = images->value.arr->dim[0];
            nPixels = xSize*ySize;
            frameSize = nPixels*images->value.arr->elt_len;
        } else {
            cerr << "rdx_sumimages: Only 3D image-cubes supported." << endl;
            return IDL_GettmpInt(0);
        }
        
        if( kw.check && nImages < 3 ) {
            cerr << "rdx_sumimages: Not enough statistics, skipping check." << endl;
            kw.check = 0;
        }

        size_t xSizePadded = xSize + 2*kw.padding;
        size_t ySizePadded = ySize + 2*kw.padding;
        
        string statusString;
        if( kw.verbose ) {
            statusString = (kw.check?"Checking and summing ":"Summing ") + to_string(nImages)
            + " images using " +to_string((int)kw.nthreads) + string(" thread") + ((kw.nthreads>1)?"s.":".");
            cout << statusString << ((kw.verbose == 1)?"\n":"") << flush;
        }
        
        IDL_MEMINT dims[] = { ySize, xSize }; 
        double* summedData = (double*)IDL_MakeTempArray( IDL_TYP_DOUBLE, 2, dims, IDL_ARR_INI_NOP, &ret ); //IDL_ARR_INI_ZERO
        
        unique_ptr<double[]> checked;
        unique_ptr<double[]> sums( new double [ nPixels*kw.nthreads ] );
        unique_ptr<double[]> tmp( new double [ nPixels*kw.nthreads ] );
        
        shared_ptr<float*> shiftsData;
        if( kw.pinh_align ) {
            shiftsData = sharedArray<float>( nImages, 2 );
        }
        float** shifts = shiftsData.get();
        
        double* sumPtr = sums.get();
        double* tmpPtr = tmp.get();
        memset( sumPtr, 0, nPixels*kw.nthreads*sizeof(double) );
        memset( summedData, 0, nPixels*sizeof(double) );
        double* checkedPtr = nullptr;
        
        
        atomic<size_t> imgIndex(0);
        atomic<size_t> threadIndex(0);
        std::vector<std::thread> threads;
        if( kw.check ) {
            checked.reset( new double [ 2*nImages ] );
            checkedPtr = checked.get();
            for( UCHAR t=0; t<kw.nthreads; ++t ) {      // calculate averages
                threads.push_back( std::thread(
                    [&](){
                        size_t myImgIndex;
                        while( (myImgIndex=imgIndex.fetch_add(1)) < nImages ) {
                            bool hasInf;
                            checkedPtr[myImgIndex] = getMinMaxMean( dataPtr+myImgIndex*frameSize, nPixels, dataType, 0, 0, &hasInf );
                            if( hasInf ) checkedPtr[myImgIndex] = std::numeric_limits<double>::infinity();
                        }
                    }));
            }
            for( auto& th : threads ) th.join();
            memcpy( checkedPtr+nImages, checkedPtr, nImages*sizeof(double) );
            nth_element( checkedPtr, checkedPtr+nImages/2, checkedPtr+nImages );        // median
            double tmean = kw.limit*checkedPtr[ nImages/2 ];
            memcpy( checkedPtr, checkedPtr+nImages, nImages*sizeof(double) );
            for( size_t i=0; i<nImages; ++i ) {
                size_t offset = std::min( std::max(nImages+i-1, nImages), 2*nImages-3 );       // restrict to array
                checkedPtr[i] = (abs(checkedPtr[i]-medianOf3(checkedPtr+offset)) <= tmean);
            }
            imgIndex = 0;
            threadIndex = 0;
        }


        
        atomic<size_t> nSummed(0);
        mutex mtx;
        
#ifdef REDUX_WITH_OPENCV
        Mat refImg;
        promise<cv::Rect> subImgROI;
        shared_future<cv::Rect> fut = subImgROI.get_future();

        const int warp_mode = MOTION_TRANSLATION; // MOTION_EUCLIDEAN;
 
        int number_of_iterations = 200;
        double termination_eps = 1e-10;
        TermCriteria criteria(TermCriteria::COUNT+TermCriteria::EPS, number_of_iterations, termination_eps);

        int flags = INTER_CUBIC;                  // INTER_LINEAR, INTER_CUBIC, INTER_AREA, INTER_LANCZOS4
        int borderMode = BORDER_CONSTANT;
        const Scalar borderValue = Scalar();
#endif
        
        auto sumFunc = [&]( size_t imgIndex, double* mySumPtr, double* myTmpPtr, UCHAR* myDataPtr) {
            try {
                if( kw.pinh_align ) {
                    
                    if( darkPtr || gainPtr ) {
                        applyDarkAndGain( myDataPtr, myTmpPtr, darkPtr, gainPtr, nPixels, dataType );
                    } else {
                        copyToRaw( myDataPtr, myTmpPtr, nPixels, dataType );
                    }
                    if( bsOtfPtr ) {
                        redux::image::descatter( myTmpPtr, ySize, xSize, bsGainPtr, bsOtfPtr, ySizePadded, xSizePadded, 1, 50, 1E-8 );
                    }
                    
                    if ( mask2D ) {
                        shared_ptr<double*> tmp2D = reshapeArray( myTmpPtr, ySize, xSize );
                        fillPixels( tmp2D.get(), (size_t)ySize, (size_t)xSize, mask2D.get() );
                    }
                    
#ifdef REDUX_WITH_OPENCV
                    Mat cvImg( ySize, xSize, CV_64FC1, myTmpPtr );
                    cv::Rect roi;
                    if( !imgIndex ) {
                        
                        int margin = 100;
                        int refSize = 99;

                        roi = cv::Rect( cv::Point(margin,margin), cv::Point(ySize-margin, xSize-margin) );
                        Mat clippedImg( cvImg, roi );
                        
                        double minVal, maxVal;
                        cv::Point maxLoc;
                        minMaxLoc( clippedImg, &minVal, &maxVal, nullptr, &maxLoc );
                        maxLoc += cv::Point(margin,margin);
                        
                        shifts[0][0] = maxLoc.x;    // save location as shift for the reference image.
                        shifts[0][1] = maxLoc.y;

                        roi = cv::Rect( maxLoc-cv::Point(refSize/2,refSize/2), cv::Size(refSize,refSize) );
                        
                        Mat subImg( cvImg, roi );
                        refImg.create( refSize, refSize, CV_32FC1 );
                        subImg.convertTo( refImg, CV_32F );

                        threshold( refImg, refImg, minVal+0.2*(maxVal-minVal), 0, THRESH_TOZERO );
                         
                        int nonZ = cv::countNonZero(refImg);
                        while( cv::countNonZero(refImg) == nonZ ) {
                             refImg.adjustROI( -1, -1, -1, -1 );
                             refSize -= 2;
                        }
                        refImg.adjustROI( 2, 2, 2, 2 );
                        refSize += 4;
                        refImg = refImg.clone();

                        roi = cv::Rect( maxLoc-cv::Point(refSize/2,refSize/2), cv::Size(refSize,refSize) );
                        Mat( cvImg, roi ).convertTo( refImg, CV_32F );

                        subImgROI.set_value( roi );
                      
                         
                    } else {

                        roi = fut.get();
                        Mat subImgFloat( roi.size(), CV_32FC1 );
                        Mat( cvImg, roi ).convertTo( subImgFloat, CV_32F );

                        Mat warp_matrix = Mat::eye( 2, 3, CV_32F );
                        findTransformECC( subImgFloat, refImg, warp_matrix, warp_mode, criteria );

                        shifts[imgIndex][0] = warp_matrix.at<float>(0,2);
                        shifts[imgIndex][1] = warp_matrix.at<float>(1,2);

                        double* slaskPtr = myTmpPtr+nPixels*kw.nthreads;
                        Mat slask( ySize, xSize, CV_64FC1, slaskPtr );
                        warpAffine( cvImg, slask, warp_matrix, cvImg.size(), flags, borderMode, borderValue );

                        memcpy( myTmpPtr, slaskPtr, nPixels*sizeof(double) );
                    }
#endif                        
                    
                    
                    for( size_t n=0; n<nPixels; ++n ) mySumPtr[n] += myTmpPtr[n];
                    
                } else {
                    addToRaw( myDataPtr, mySumPtr, nPixels, dataType );
                }
                
                size_t ns = nSummed++;
                if( kw.verbose > 1 ) printProgress( statusString, (ns*100.0/(nImages-1)));
            } catch( const exception& e ) {
                cout << "rdx_sumimages:sumFunc: unhandled exception: " << e.what() << endl;
            }
        };

        threads.clear();
        for( UCHAR t=0; t<kw.nthreads; ++t ) {
            threads.push_back( std::thread(
                [&](){
                    size_t myImgIndex;
                    size_t myThreadIndex = threadIndex.fetch_add(1);
                    double* mySumPtr = sumPtr+myThreadIndex*nPixels;
                    double* myTmpPtr = tmpPtr+myThreadIndex*nPixels;
                    UCHAR* myDataPtr = dataPtr+myThreadIndex*frameSize;
                    while( (myImgIndex=imgIndex.fetch_add(1)) < nImages ) {
                        if( !checkedPtr || checkedPtr[myImgIndex] > 0 ) {
                            sumFunc( myImgIndex, mySumPtr, myTmpPtr, reinterpret_cast<UCHAR*>(myDataPtr) );
                        } else {
                            if ( kw.lun ) {
                                string tmp = "Image #" + to_string(myImgIndex);
                                IDL_PoutRaw( kw.lun, (char*)tmp.c_str(), tmp.length() );
                                IDL_PoutRaw( kw.lun, nl, 1 );
                            }
                        }
                    }
                    std::unique_lock<mutex> lock(mtx);
                    for( size_t i=0; i<nPixels; ++i ) summedData[i] += mySumPtr[i];
                }));
        }
        for (auto& th : threads) th.join();
        
#ifdef REDUX_WITH_OPENCV
        if( kw.pinh_align > 1 ) {
            Mat cvImg( ySize, xSize, CV_64FC1, summedData );
            Mat warp_matrix = Mat::eye( 2, 3, CV_32F );
            for( size_t n=1; n<nSummed; ++n ) {
                warp_matrix.at<float>(0,2) -= shifts[n][0];
                warp_matrix.at<float>(1,2) -= shifts[n][1];
            }
            warp_matrix.at<float>(0,2) /= nSummed;
            warp_matrix.at<float>(1,2) /= nSummed;
            Mat slask( ySize, xSize, CV_64FC1, tmpPtr );
            warpAffine( cvImg, slask, warp_matrix, cvImg.size(), flags, borderMode, borderValue );
            memcpy( summedData, tmpPtr, nPixels*sizeof(double) );
        }
#endif

        if( kw.verbose > 1 ) {
            printProgress( statusString, 100.0 );
            if( nSummed == nImages ) {
                cout << "  All ok." << endl;
            } else {
                cout << "  Check failed for " << (nImages-nSummed) << " files." << endl;
            }
        }
        
        if( kw.summed ) {
            IDL_VPTR tmpSummed;
            double* tmpData = (double*)IDL_MakeTempArray( IDL_TYP_DOUBLE, 2, dims, IDL_ARR_INI_NOP, &tmpSummed );
            memcpy( tmpData, summedData, nPixels*sizeof(double));
            IDL_VarCopy( tmpSummed, kw.summed );
        }
        
        if( kw.nsummed ) {
            IDL_VPTR tmpNS = IDL_GettmpLong( nSummed );
            IDL_VarCopy( tmpNS, kw.nsummed );
        }
        
        if( kw.pinh_align && kw.xyc ) {
            IDL_VPTR tmp;
            IDL_MEMINT dims[] = { 2, static_cast<IDL_MEMINT>(nImages) }; 
            float* tmpData = (float*)IDL_MakeTempArray( IDL_TYP_FLOAT, 2, dims, IDL_ARR_INI_ZERO, &tmp ); //IDL_ARR_INI_ZERO
            memcpy( tmpData, *shifts, 2*nImages*sizeof(float));
            IDL_VarCopy( tmp, kw.xyc );
        }
        
        
        if( nSummed ) {
            double nImg_inv = 1.0/nSummed;
            for( size_t i=0; i<nPixels; ++i ) summedData[i] *= nImg_inv;
        }
        
        if( !kw.pinh_align ) {

            applyDarkAndGain( summedData, summedData, darkData.get(), gainData.get(), nPixels );
            
            if( bsOtfPtr ) {
                redux::image::descatter( summedData, ySize, xSize, bsGainPtr, bsOtfPtr, ySizePadded, xSizePadded, kw.nthreads, 50, 1E-8 );
            }
            
            if ( mask2D ) {
                shared_ptr<double*> summed2D = reshapeArray( summedData, ySize, xSize );
                fillPixels( summed2D.get(), (size_t)ySize, (size_t)xSize, mask2D.get() );
            }
        }
        
        return ret;
        
    } catch (const exception& e ) {
        cout << "rdx_sumimages: unhandled exception: " << e.what() << endl;
        return IDL_GettmpInt(0);
    }

    
}


string sum_files_info( int lvl ) {
    string ret = "RDX_SUMFILES";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"       ");          // newline if lvl>1
        ret += "   Syntax:   out = rdx_sumfiles(file_list, /KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      BACKSCATTER_GAIN    Image containing the imprinted background pattern. (for descattering).\n"
                    "      BACKSCATTER_PSF     Image containing the scattering PSF. (for descattering).\n"
                    "      CHECK               Discard images which are statistically deviant.\n"
                    "      DARK                Dark field.\n"
                    "      GAIN                Gain table (inverted flat-field).\n"
                    "      LIMIT               Allowed deviation from the median. (0.0175)\n"
                    "      LUN                 IDL file unit (id) where discarded files will be logged.\n"
                    "      NSUMMED             (output) Number of images actually summed.\n"
                    "      NTHREADS            Number of threads.\n"
                    "      PADDING             Padding size for the descattering procedure. (256)\n"
                    "      PINHOLE_ALIGN       Do sub-pixel alignment before summing.\n"
                    "      SUMMED              (output) Raw sum.\n"
                    "      TIME                (output) Average timestamp from file-headers.\n"
                    "      VERBOSE             Verbosity, default is 0 (only error output).\n"
                    "      XYC                 (output) Coordinates of align-feature and image-shifts.\n";
        }
    } else ret += "\n";
    return ret;
}


IDL_VPTR sum_files( int argc, IDL_VPTR* argv, char* argk ) {
    
    IDL_VPTR filenames = argv[0];
    IDL_ENSURE_SIMPLE( filenames );
    
    if( filenames->type != IDL_TYP_STRING ) return sum_images( argc, argv, argk );

#ifndef REDUX_WITH_OPENCV
    cerr << "rdx_sumfiles: redux has to be re-compiled with OpenCV enabled to be able to use this function." << endl;
    return IDL_GettmpInt(0);
#else
    
    SI_KW kw;
    kw.dark = nullptr;
    kw.gain = nullptr;
    kw.nsummed = nullptr;
    kw.summed = nullptr;
    kw.bs_gain = nullptr;
    kw.bs_psf = nullptr;
    kw.limit = 0.0175;
    kw.padding = 256;
    kw.nthreads = std::thread::hardware_concurrency();
    (void)IDL_KWProcessByOffset( argc, argv, argk, si_kw_pars, (IDL_VPTR*)0, 255, &kw );
    
    if( kw.help ) {
        cout << sum_files_info(2) << endl;
        return IDL_GettmpInt(0);
    }
    
    kw.nthreads = max<UCHAR>(1, min<UCHAR>(kw.nthreads, thread::hardware_concurrency()));
    kw.padding = max<IDL_INT>(0, min<IDL_INT>(kw.padding, 4096));   // prevent insane padding
    
    try {
        
        IDL_VPTR ret;
        shared_ptr<double> darkData, gainData, bsGainData;
        shared_ptr<uint8_t> maskData;
        shared_ptr<uint8_t*> mask2D;
        shared_ptr<fftw_complex> bsOtfData;
        double *darkPtr = nullptr;
        double *gainPtr = nullptr;
        double *bsGainPtr = nullptr;
        fftw_complex *bsOtfPtr = nullptr;
        IDL_MEMINT xCalibSize(0);
        IDL_MEMINT yCalibSize(0);
        
        if( kw.dark ) {
            IDL_ENSURE_SIMPLE( kw.dark );
            IDL_ENSURE_ARRAY( kw.dark );
            if( kw.dark->value.arr->n_dim == 2 ) {
                darkData = castOrCopy<double>( kw.dark );
                darkPtr = darkData.get();
                xCalibSize = kw.dark->value.arr->dim[1];
                yCalibSize = kw.dark->value.arr->dim[0];
            } else cout << "dark must be a 2D image." << endl;
        }
        if( kw.gain ) {
            IDL_ENSURE_SIMPLE( kw.gain );
            IDL_ENSURE_ARRAY( kw.gain );
            if( kw.gain->value.arr->n_dim == 2 ) {
                if( (xCalibSize && (xCalibSize != kw.gain->value.arr->dim[1])) ||
                    (yCalibSize && (yCalibSize != kw.gain->value.arr->dim[0])) ) {
                    cout << "Dark & Gain files have different sizes." << endl;
                    return IDL_GettmpInt(0);
                }
                gainData = castOrCopy<double>( kw.gain  );
                gainPtr = gainData.get();
                xCalibSize = kw.gain->value.arr->dim[1];
                yCalibSize = kw.gain->value.arr->dim[0];
                maskData.reset( new uint8_t[xCalibSize*yCalibSize] );
                mask2D = reshapeArray( maskData.get(), yCalibSize, xCalibSize );
                redux::image::make_mask( gainPtr, maskData.get(), yCalibSize, xCalibSize, 0, 5, true, true ); // filter away larger features than ~5 pixels and invert
            } else cout << "gain must be a 2D image." << endl;
        }
        if( kw.bs_gain && kw.bs_psf ) {
            IDL_ENSURE_SIMPLE( kw.bs_gain );
            IDL_ENSURE_ARRAY( kw.bs_gain );
            IDL_ENSURE_SIMPLE( kw.bs_psf );
            IDL_ENSURE_ARRAY( kw.bs_psf );
            if( darkData && gainData ) {
                if( kw.bs_gain->value.arr->n_dim==2 && kw.bs_psf->value.arr->n_dim==2 ) {
                    IDL_MEMINT bsXsize = kw.bs_gain->value.arr->dim[1];
                    IDL_MEMINT bsYsize = kw.bs_gain->value.arr->dim[0];
                    if( bsXsize != kw.bs_psf->value.arr->dim[1] || bsYsize != kw.bs_psf->value.arr->dim[0]) {
                        cout << "Backscatter Gain & PSF have different sizes." << endl;
                        return IDL_GettmpInt(0);
                    }
                    if( (xCalibSize && (xCalibSize != bsXsize)) ||
                        (yCalibSize && (yCalibSize != bsYsize)) ) {
                        cout << "Backscatter Gain & PSF have different size than the dark/gain." << endl;
                        return IDL_GettmpInt(0);
                    }
                    size_t xSizePadded = bsXsize + 2*kw.padding;
                    size_t ySizePadded = bsYsize + 2*kw.padding;
                    size_t dataSize = xSizePadded*ySizePadded;

                    bsGainData.reset( fftw_alloc_real(dataSize), fftw_free );
                    bsGainPtr = bsGainData.get();
                    memset( bsGainPtr, 0, dataSize*sizeof(double) );
                    copyInto( kw.bs_gain, bsGainPtr, ySizePadded, xSizePadded, kw.padding, kw.padding );
                    
                    shared_ptr<double> bsPsfData( fftw_alloc_real(dataSize), fftw_free );
                    double* psfPtr = bsPsfData.get();
                    memset( psfPtr, 0, dataSize*sizeof(double) );
                    copyInto( kw.bs_psf, psfPtr, ySizePadded, xSizePadded, kw.padding, kw.padding );

                    // normalize psf
                    double sum = 0;
                    for( size_t i=0; i<dataSize; ++i ) sum += psfPtr[i];
                    sum = 1.0/(sum*dataSize);
                    for( size_t i=0; i<dataSize; ++i ) psfPtr[i] *= sum;
                    
                    FourierTransform::reorder( psfPtr, ySizePadded, xSizePadded );
                    
                    FourierTransform::Plan::Ptr plan = FourierTransform::Plan::get( ySizePadded, xSizePadded, FourierTransform::Plan::R2C, 1 );
                    bsOtfData.reset( fftw_alloc_complex( ySizePadded*(xSizePadded/2+1) ), fftw_free );
                    bsOtfPtr = bsOtfData.get();
                    plan->forward( bsPsfData.get(), bsOtfPtr );
                    
                } else cout << "backscatter_gain/psf must be a 2D images." << endl;
            } else cout << "Backscatter correction requires both dark & gain to be present." << endl;
        } else if( kw.bs_gain || kw.bs_psf ) cout << "Both backscatter_gain and backscatter_psf must be specified." << endl;
        
        vector<string> existingFiles;
        if ( !(filenames->flags & IDL_V_ARR) ) {
            bfs::path fn( string(filenames->value.str.s) );
            if( bfs::is_regular_file(fn) ) {
                existingFiles.push_back( fn.string() );
            } else return IDL_GettmpInt(0);
        } else {
            IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>(filenames->value.arr->data);
            for( int i=0; i<filenames->value.arr->n_elts; ++i ) {
                bfs::path fn( string(strptr[i].s) );
                if( bfs::is_regular_file(fn) ) {
                    existingFiles.push_back( fn.string() );
                }
            }
        }
        size_t nImages = existingFiles.size();
        kw.nthreads = static_cast<UCHAR>( min<size_t>(kw.nthreads, nImages) );

        if( !nImages ) { 
            cout << "rdx_sumfiles: No input files." << endl;
            return IDL_GettmpInt(0);
        }
        // Get size etc from first image.
        std::shared_ptr<redux::file::FileMeta> meta = redux::file::getMeta( existingFiles[0] );
        if( !meta || (meta->nDims() != 2) ) {        // Only allow 2D images
            cout << "rdx_sumfiles: Failed to get meta, or not 2D input." << endl;
            return IDL_GettmpInt(0);
        }
        
        IDL_MEMINT xSize = meta->dimSize(1);
        IDL_MEMINT ySize = meta->dimSize(0);
        size_t nPixels = xSize*ySize;
        size_t frameSize = meta->dataSize();
        UCHAR dataType = meta->getIDLType();
        
        size_t xSizePadded = xSize + 2*kw.padding;
        size_t ySizePadded = ySize + 2*kw.padding;
        
        vector<boost::posix_time::time_duration> times;
        if( kw.time ) {
            times.resize( nImages );
        }
        
        string statusString;
        if( kw.verbose ) {
            statusString = (kw.check?"Checking and summing ":"Summing ") + to_string(nImages)
            + " files using " +to_string((int)kw.nthreads) + string(" thread") + ((kw.nthreads>1)?"s.":".");
            cout << statusString << ((kw.verbose == 1)?"\n":"") << flush;
        }
        
        IDL_MEMINT dims[] = { ySize, xSize }; 
        double* summedData = (double*)IDL_MakeTempArray( IDL_TYP_DOUBLE, 2, dims, IDL_ARR_INI_NOP, &ret ); //IDL_ARR_INI_ZERO
        
        unique_ptr<char[]> loadBuffer( new char[ kw.nthreads*frameSize ] );
        char* loadPtr = loadBuffer.get();
        
        unique_ptr<double[]> checked;
        unique_ptr<double[]> sums( new double [ nPixels*kw.nthreads ] );
        unique_ptr<double[]> tmp( new double [ 2*nPixels*kw.nthreads ] );
        
        shared_ptr<float*> shiftsData;
        float** shifts = nullptr;
        if( kw.pinh_align ) {
            shiftsData = sharedArray<float>( nImages, 2 );
            shifts = shiftsData.get();
        }
        
        double* sumPtr = sums.get();
        double* tmpPtr = tmp.get();
        memset( sumPtr, 0, nPixels*kw.nthreads*sizeof(double) );
        memset( summedData, 0, nPixels*sizeof(double) );
        double* checkedPtr = nullptr;
        
        if( kw.check ) {
            if( nImages < 3 ) {
                cerr << "rdx_sumimages: Not enough statistics, skipping check." << endl;
                kw.check = 0;
            } else {
                checked.reset( new double [ 2*nImages ] );
                checkedPtr = checked.get();
            }
        }
        
        atomic<size_t> imgIndex(0);
        atomic<size_t> threadIndex(0);
        atomic<size_t> nSummed(0);
        
        mutex mtx;

        Mat refImg;
        promise<cv::Rect> subImgROI;
        shared_future<cv::Rect> fut = subImgROI.get_future();

        const int warp_mode = MOTION_TRANSLATION; // MOTION_EUCLIDEAN;
 
        int number_of_iterations = 200;
        double termination_eps = 1e-10;
        TermCriteria criteria(TermCriteria::COUNT+TermCriteria::EPS, number_of_iterations, termination_eps);

        int flags = INTER_CUBIC;                  // INTER_LINEAR, INTER_CUBIC, INTER_AREA, INTER_LANCZOS4
        int borderMode = BORDER_CONSTANT;
        const Scalar borderValue = Scalar();

        auto sumFunc = [&]( size_t imgIndex, double* mySumPtr, double* myTmpPtr, UCHAR* myLoadPtr ) {

                if( checkedPtr  ) {
                    bool hasInf;
                    checkedPtr[imgIndex] = getMinMaxMean( myLoadPtr, nPixels, dataType, 0, 0, &hasInf );
                    if( hasInf ) {
                        checkedPtr[imgIndex] = std::numeric_limits<double>::infinity();
                        return;
                    }
                }
                
                if( kw.pinh_align ) {
                    
                    if( darkPtr || gainPtr ) {
                        applyDarkAndGain( myLoadPtr, myTmpPtr, darkPtr, gainPtr, nPixels, dataType );
                    } else {
                        copyToRaw( myLoadPtr, myTmpPtr, nPixels, dataType );
                    }
                    if( bsOtfPtr ) {
                        redux::image::descatter( myTmpPtr, ySize, xSize, bsGainPtr, bsOtfPtr, ySizePadded, xSizePadded, 1, 50, 1E-8 );
                    }
                    
                    if ( mask2D ) {
                        shared_ptr<double*> tmp2D = reshapeArray( myTmpPtr, ySize, xSize );
                        fillPixels( tmp2D.get(), (size_t)ySize, (size_t)xSize, mask2D.get() );
                    }
                    
                    Mat cvImg( ySize, xSize, CV_64FC1, myTmpPtr );
                    cv::Rect roi;
                    if( !imgIndex ) {
                        
                        int margin = 100;
                        int refSize = 99;

                        roi = cv::Rect( cv::Point(margin,margin), cv::Point(ySize-margin, xSize-margin) );
                        Mat clippedImg( cvImg, roi );
                        
                        double minVal, maxVal;
                        cv::Point maxLoc;
                        minMaxLoc( clippedImg, &minVal, &maxVal, nullptr, &maxLoc );
                        maxLoc += cv::Point(margin,margin);
                        
                        shifts[0][0] = maxLoc.x;    // save location as shift for the reference image.
                        shifts[0][1] = maxLoc.y;

                        roi = cv::Rect( maxLoc-cv::Point(refSize/2,refSize/2), cv::Size(refSize,refSize) );
                        
                        Mat subImg( cvImg, roi );
                        refImg.create( refSize, refSize, CV_32FC1 );
                        subImg.convertTo( refImg, CV_32F );

                        threshold( refImg, refImg, minVal+0.2*(maxVal-minVal), 0, THRESH_TOZERO );
                         
                        int nonZ = cv::countNonZero(refImg);
                        while( cv::countNonZero(refImg) == nonZ ) {
                             refImg.adjustROI( -1, -1, -1, -1 );
                             refSize -= 2;
                        }
                        refImg.adjustROI( 2, 2, 2, 2 );
                        refSize += 4;
                        refImg = refImg.clone();

                        roi = cv::Rect( maxLoc-cv::Point(refSize/2,refSize/2), cv::Size(refSize,refSize) );
                        Mat( cvImg, roi ).convertTo( refImg, CV_32F );

                        subImgROI.set_value( roi );
                         
                    } else {

                        roi = fut.get();
                        Mat subImgFloat( roi.size(), CV_32FC1 );
                        Mat( cvImg, roi ).convertTo( subImgFloat, CV_32F );

                        Mat warp_matrix = Mat::eye( 2, 3, CV_32F );
                        findTransformECC( subImgFloat, refImg, warp_matrix, warp_mode, criteria );

                        shifts[imgIndex][0] = warp_matrix.at<float>(0,2);
                        shifts[imgIndex][1] = warp_matrix.at<float>(1,2);

                        double* slaskPtr = myTmpPtr+nPixels*kw.nthreads;
                        Mat slask( ySize, xSize, CV_64FC1, slaskPtr );
                        warpAffine( cvImg, slask, warp_matrix, cvImg.size(), flags, borderMode, borderValue );

                        memcpy( myTmpPtr, slaskPtr, nPixels*sizeof(double) );
                        
                    }
                    
                    for( size_t n=0; n<nPixels; ++n ) mySumPtr[n] += myTmpPtr[n];

                } else {
                    addToRaw( myLoadPtr, mySumPtr, nPixels, dataType );
                }
                
        };

        std::vector<std::thread> threads;
        for( UCHAR t=0; t<kw.nthreads; ++t ) {
            threads.push_back( std::thread(
                [&](){
                    size_t myImgIndex;
                    size_t myThreadIndex = threadIndex.fetch_add(1);
                    double* mySumPtr = sumPtr+myThreadIndex*nPixels;
                    double* myTmpPtr = tmpPtr+myThreadIndex*nPixels;
                    char* myLoadPtr = loadPtr+myThreadIndex*frameSize;
                    shared_ptr<redux::file::FileMeta> myMeta;
                    while( (myImgIndex=imgIndex.fetch_add(1)) < nImages ) {
                        try {
                            readFile( existingFiles[myImgIndex], myLoadPtr, myMeta );
                            sumFunc( myImgIndex, mySumPtr, myTmpPtr, reinterpret_cast<UCHAR*>(myLoadPtr) );
                            if( kw.time ) {
                                 times[myImgIndex] = myMeta->getAverageTime().time_of_day();
                            }
                           size_t ns = nSummed++;
                            if( kw.verbose > 1 ) printProgress( statusString, (ns*100.0/(nImages-1)));
                        } catch( const exception& e ) {
                            cout << "rdx_sumfiles: Failed to load file: " << existingFiles[myImgIndex] << "  Reason: " << e.what() << endl;
                        }
                    }
                    std::unique_lock<mutex> lock(mtx);
                    for( size_t i=0; i<nPixels; ++i ) summedData[i] += mySumPtr[i];
                }));
        }
        for (auto& th : threads) th.join();
        
        size_t nDiscarded(0);
        if( checkedPtr ) {
            set<size_t> discarded;
            memcpy( checkedPtr+nImages, checkedPtr, nImages*sizeof(double) );
            nth_element( checkedPtr, checkedPtr+nImages/2, checkedPtr+nImages );        // median
            double tmean = kw.limit*checkedPtr[ nImages/2 ];
            memcpy( checkedPtr, checkedPtr+nImages, nImages*sizeof(double) );
            for( size_t i=0; i<nImages; ++i ) {
                if( isfinite(checkedPtr[i]) ) {
                    size_t offset = std::min( std::max(nImages+i-1, nImages), 2*nImages-3 );       // restrict to array
                    checkedPtr[i] = (abs(checkedPtr[i]-medianOf3(checkedPtr+offset)) <= tmean);
                    if ( checkedPtr[i] == 0 ) {
                        if ( kw.lun ) {
                            IDL_PoutRaw( kw.lun, (char*)existingFiles[i].c_str(), existingFiles[i].length() );
                            IDL_PoutRaw( kw.lun, nl, 1 );
                        }
                        ++nDiscarded;
                    }
                } else {
                    IDL_PoutRaw( kw.lun, (char*)existingFiles[i].c_str(), existingFiles[i].length() );
                    IDL_PoutRaw( kw.lun, nl, 1 );
                    ++nDiscarded;
                }
            }
            if( nDiscarded ) {
                imgIndex = 0;
                threadIndex = 0;
                threads.clear();
                for( UCHAR t=0; t<kw.nthreads; ++t ) {
                    threads.push_back( std::thread(
                        [&](){
                            size_t myImgIndex;
                            size_t myThreadIndex = threadIndex.fetch_add(1);
                            double* mySumPtr = sumPtr+myThreadIndex*nPixels;
                            memset( mySumPtr, 0, nPixels*sizeof(double) );
                            double* myTmpPtr = tmpPtr+myThreadIndex*nPixels;
                            char* myLoadPtr = loadPtr+myThreadIndex*frameSize;
                            shared_ptr<redux::file::FileMeta> myMeta;
                            while( (myImgIndex=imgIndex.fetch_add(1)) < nImages ) {
                                if( checkedPtr[myImgIndex] == 0 ) {     // subtract discarded images
                                    try {
                                        readFile( existingFiles[myImgIndex], loadPtr, myMeta );
                                        if( kw.time ) times[myImgIndex] = boost::posix_time::time_duration(0,0,0);
                                        sumFunc( myImgIndex, mySumPtr, myTmpPtr, reinterpret_cast<UCHAR*>(myLoadPtr) );
                                        --nSummed;
                                        if( shifts ) {
                                            shifts[myImgIndex][0] = 0;
                                            shifts[myImgIndex][1] = 0;
                                        }
                                    } catch( const exception& e ) {
                                        cout << "rdx_sumfiles: Failed to load file: " << existingFiles[myImgIndex] << "  Reason: " << e.what() << endl;
                                    }
                                } 
                            }
                            std::unique_lock<mutex> lock(mtx);
                            for( size_t i=0; i<nPixels; ++i ) summedData[i] -= mySumPtr[i];
                        }));
                }
                for (auto& th : threads) th.join();
            }
        }
        
        if( kw.pinh_align > 1 ) {
            Mat cvImg( ySize, xSize, CV_64FC1, summedData );
            Mat warp_matrix = Mat::eye( 2, 3, CV_32F );
            for( size_t n=1; n<nSummed; ++n ) {
                warp_matrix.at<float>(0,2) -= shifts[n][0];
                warp_matrix.at<float>(1,2) -= shifts[n][1];
            }
            warp_matrix.at<float>(0,2) /= nSummed;
            warp_matrix.at<float>(1,2) /= nSummed;
            Mat slask( ySize, xSize, CV_64FC1, tmpPtr );
            warpAffine( cvImg, slask, warp_matrix, cvImg.size(), flags, borderMode, borderValue );
            memcpy( summedData, tmpPtr, nPixels*sizeof(double) );
        }
        
        if( kw.verbose > 1 ) {
            printProgress( statusString, 100.0 );
            if( nDiscarded ) {
                cout << "  Check failed for " << nDiscarded << " file" << ((nDiscarded>1)?"s.":".") << endl;
            } else {
                cout << "  All ok." << endl;
            }
        }
        
        if( kw.summed ) {
            IDL_VPTR tmpSummed;
            double* tmpData = (double*)IDL_MakeTempArray( IDL_TYP_DOUBLE, 2, dims, IDL_ARR_INI_NOP, &tmpSummed );
            memcpy( tmpData, summedData, nPixels*sizeof(double));
            IDL_VarCopy( tmpSummed, kw.summed );
        }
        
        if( kw.nsummed ) {
            IDL_VPTR tmpNS = IDL_GettmpLong( nSummed );
            IDL_VarCopy( tmpNS, kw.nsummed );
        }
        
        if( kw.pinh_align && kw.xyc ) {
            IDL_VPTR tmp;
            IDL_MEMINT dims[] = { 2, static_cast<IDL_MEMINT>(nImages) }; 
            float* tmpData = (float*)IDL_MakeTempArray( IDL_TYP_FLOAT, 2, dims, IDL_ARR_INI_ZERO, &tmp ); //IDL_ARR_INI_ZERO
            memcpy( tmpData, *shifts, 2*nImages*sizeof(float));
            IDL_VarCopy( tmp, kw.xyc );
        }
        
        if( kw.time ) {
            boost::posix_time::time_duration avgTime(0,0,0);
            for( auto& t: times ) {
                avgTime += t;
            }
            if( nSummed ) {
                avgTime /= nSummed;
            }
            string tStr = boost::posix_time::to_simple_string(avgTime);
            IDL_VPTR tmpTimeString = IDL_StrToSTRING( (char*)tStr.c_str() );
            IDL_VarCopy( tmpTimeString, kw.time );
        }

        if( nSummed ) {
            double nImg_inv = 1.0/nSummed;
            for( size_t i=0; i<nPixels; ++i ) summedData[i] *= nImg_inv;
        }

        if( !kw.pinh_align ) {

            applyDarkAndGain( summedData, summedData, darkData.get(), gainData.get(), nPixels );
            
            if( bsOtfPtr ) {
                redux::image::descatter( summedData, ySize, xSize, bsGainPtr, bsOtfPtr, ySizePadded, xSizePadded, kw.nthreads, 50, 1E-8 );
            }
            
            if ( mask2D ) {
                shared_ptr<double*> summed2D = reshapeArray( summedData, ySize, xSize );
                fillPixels( summed2D.get(), (size_t)ySize, (size_t)xSize, mask2D.get() );
            }
        }
        
        return ret;
        
    } catch (const exception& e ) {
        cout << "rdx_sumfiles: unhandled exception: " << e.what() << endl;
        return IDL_GettmpInt(0);
    }

#endif

}


/*typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; // Must be first entry in structure
    IDL_INT help;
    float range;
    IDL_INT telea;
    IDL_INT verbose;
    //IDL_VPTR time;
   // IDL_STRING split_chars;
} INP_KW;


// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR inp_kw_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "HELP",           IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(INP_KW,help) },
    { (char*) "RANGE",          IDL_TYP_FLOAT, 1,           0, 0, (char*) IDL_KW_OFFSETOF2(INP_KW,range) },
    { (char*) "TELEA",          IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(INP_KW,telea) },
//    { (char*) "TIME",             IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(INP_KW,time) },
    { (char*) "VERBOSE",        IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(INP_KW,verbose) },
   // { (char*) "SPLIT_CHARS", IDL_TYP_STRING, 1, 0, 0, (char*)IDL_KW_OFFSETOF2(SI_KW,split_chars) },
    { NULL }
};
*/

//#include <opencv2/photo/photo.hpp>

IDL_VPTR redux::inpaint( int argc, IDL_VPTR* argv, char* argk ) {
    
#ifndef REDUX_WITH_OPENCV
    cerr << "inpaint: redux has to be re-compiled with OpenCV enabled to be able to use this function." << endl;
    return IDL_GettmpInt(0);
#else
    cerr << "inpaint: redux has to be re-compiled with OpenCV enabled to be able to use this function." << endl;
    return IDL_GettmpInt(0);

#ifdef sdfhshsfgjdjdgj
    
    INP_KW kw;
    kw.range = 1;
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, inp_kw_pars, (IDL_VPTR*)0, 255, &kw );

    if( nPlainArgs < 2 ) {
        cout << "inpaint: Needs 2 inputs, image and mask." << endl;
        return IDL_GettmpInt(0);
    }
    
    int type = cv::INPAINT_NS;
    if( kw.help ) {
        cout << "inpaint: No help available." << endl;
        return IDL_GettmpInt(0);
    }
    
    if( kw.range < 0 ) {
        if ( kw.verbose ) cout << "inpaint: Nothing to do." << endl;
        return argv[0];
    }
    
    if( kw.telea ) {
        type = CV_INPAINT_TELEA;
    }
    
    IDL_VPTR img = argv[0];
    IDL_ENSURE_ARRAY( img );
    IDL_ENSURE_SIMPLE( img );
    
    IDL_VPTR mask = argv[1];
    IDL_ENSURE_ARRAY( mask );
    IDL_ENSURE_SIMPLE( mask );
    
    Mat imgMat = arrayToMat( img );
    Mat imgFloat = getImgAsGrayFloat( img );
    Mat maskMat = arrayToMat( mask );
    
    Mat byteMask( maskMat.rows, maskMat.cols, CV_8UC1 );
    maskMat.convertTo( byteMask, CV_8UC1, 1 );

    IDL_VPTR ret;
    IDL_MakeTempArray( img->type, img->value.arr->n_dim, img->value.arr->dim, IDL_ARR_INI_NOP, &ret );
    Mat retMat = arrayToMat( ret );
    
//     IDL_MEMINT dims[] = { 3, 1024, 1024 };
//     int dims2[] = { 1024, 1024 };
//     char* data = IDL_MakeTempArray( IDL_TYP_BYTE, 3, dims, IDL_ARR_INI_NOP, &ret );
//    Mat tmpMat( imgMat.size(), CV_MAKETYPE(CV_8U, 3) );
    
    double minValue, maxValue;
    cv::minMaxLoc( imgMat, &minValue, &maxValue );
    
/*    typedef Vec<uint8_t, 3> VT;
    std::transform( tmpMat.begin<VT>(), tmpMat.end<VT>(), imgFloat.begin<float>(), tmpMat.begin<VT>(),
                    [=](const VT& a, const float& b){
                        float scaled = (b-minValue)/(maxValue-minValue) * std::numeric_limits<uint32_t>::max();
                        uint32_t scaled_int = static_cast<uint32_t>(scaled);
                        uint8_t* ptr = reinterpret_cast<uint8_t*>(&scaled_int);
                        if( REDUX_BYTE_ORDER == REDUX_BIG_ENDIAN ) ptr++;
                        return VT(*ptr++,*ptr++,*ptr++);
                    }  );

*/
    cout << "Blaha: minValue=" << minValue << "  maxValue=" << maxValue << endl;
    //retMat = maxValue;
//    cvtColor( imgMat, retMat, CV_GRAY2RGB );

    double alpha = (255.0-0.0)/(maxValue-minValue);
    double beta = 0.0 - minValue*alpha;
    
    Mat byteImg, result;
    imgMat.convertTo( byteImg, CV_8UC1, alpha, beta );
    double minValue2, maxValue2;
    cv::minMaxLoc( byteImg, &minValue2, &maxValue2 );
    
    cout << "Blaha: minValue2=" << minValue2 << "  maxValue2=" << maxValue2 << endl;

    cv::inpaint( byteImg, byteMask, result, kw.range, type );
    
    cv::minMaxLoc( result, &minValue2, &maxValue2 );
    
    cout << "Blaha: minValue3=" << minValue2 << "  maxValue3=" << maxValue2 << endl;
    
    result.convertTo( retMat, retMat.type(), 1.0/alpha, minValue );
    
    cv::minMaxLoc( retMat, &minValue2, &maxValue2 );
    
    cout << "Blaha: minValue4=" << minValue2 << "  maxValue4=" << maxValue2 << endl;
    
//    return ret;
    
//    byteMask = 1 - byteMask;
    
    //Mat blaha = byteMask( cv::Rect(300,300,500,500) );
    //blaha = 1;
    byteMask = 1 - byteMask;
    imgMat.copyTo( retMat, byteMask );


/*    cv::inpaint( tmpMat.clone(), byteMask, tmpMat, kw.range, type );
    
    std::transform( tmpMat.begin<VT>(), tmpMat.end<VT>(), imgFloat.begin<float>(), imgFloat.begin<float>(),
                    [=](const VT& a, const float& b){
                        uint32_t scaled_int;
                        uint8_t* ptr = reinterpret_cast<uint8_t*>(&scaled_int);
                        if( REDUX_BYTE_ORDER == REDUX_BIG_ENDIAN ) ptr++;
                        *ptr++ = a[0];
                        *ptr++ = a[1];
                        *ptr++ = a[2];
                        float scaled = (b-minValue)/(maxValue-minValue) * std::numeric_limits<uint32_t>::max();
                        return scaled_int*(maxValue-minValue)/std::numeric_limits<uint32_t>::max() + minValue;
                    }  );
    
    imgFloat.assignTo( retMat, retMat.type() );
*/

    return ret;
    
#endif
#endif

}

namespace {
    static int dummy RDX_UNUSED =
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_find_shift, (char*)"RDX_FIND_SHIFT", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1 ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_fillpix, (char*)"RDX_FILLPIX", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, rdx_fillpix_info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)load_files, (char*)"RDX_LOADFILES", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, load_files_info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_make_mask,  (char*)"RDX_MAKE_MASK",  0, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, make_mask_info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_make_win,  (char*)"RDX_MAKE_WINDOW",  2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, make_win_info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)sum_images, (char*)"RDX_SUMIMAGES", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, sum_images_info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)sum_files,  (char*)"RDX_SUMFILES",  1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, sum_files_info );
}
