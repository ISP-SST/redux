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
#include <redux/util/progresswatch.hpp>

#include <atomic>
#include <future>
#include <iostream>
#include <numeric>
#include <set>
#include <thread>

#include <boost/asio.hpp>
#include <boost/thread/thread.hpp>
#include <boost/filesystem.hpp>

#define DEBUG_DLM_

#ifdef RDX_WITH_OPENCV
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
        IDL_INT by_size;
        float eps;
        IDL_INT margin;
        IDL_INT max_dist;
        IDL_INT max_shift;
        float max_scale;
        IDL_INT max_points;
        IDL_INT niter;
        IDL_INT nrefpoints;
        IDL_INT orientation;
        IDL_INT show;
        IDL_INT smooth;
        IDL_INT verbose;
        IDL_VPTR h_init;
        IDL_VPTR status;
        IDL_VPTR threshold;
        IDL_VPTR points;
    } KW_RESULT;

    
#ifdef RDX_WITH_OPENCV
    
    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_MEMINT dims3x3[] = { 3, 3 };
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "BY_SIZE",    IDL_TYP_INT,   1, 0,           0, (char*) IDL_KW_OFFSETOF (by_size) },
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
        { (char*) "ORIENTATION",IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (orientation) },
        { (char*) "POINTS",     IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (points) },
        { (char*) "SHOW",       IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (show) },
        { (char*) "SMOOTH",     IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (smooth) },
        { (char*) "STATUS",     IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (status) },
        { (char*) "THRESHOLD",  IDL_TYP_UNDEF, 1, IDL_KW_VIN|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (threshold) },
        { (char*) "VERBOSE",    IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (verbose) },
        { NULL }
    };
    
    bool transformCheck (const Mat& trans, const Size& imgSize1, const Size& imgSize2, double approximateScale,
                         double maxShift=10000, double maxScaleDiff=1.5, int orientation=0 ) {

        // check that orthogonality is approximately preserved
        double ortho = trans.at<float>(0,0)*trans.at<float>(0,1) + trans.at<float>(1,0)*trans.at<float>(1,1);
        if( fabs(ortho) > 0.05 ) {
            return false;
        }
        
        double det = determinant( trans (cv::Rect_<int> (0, 0, 2, 2)) );
        // check that the determinant has the right sign (i.e. if the matrix has the right orientaion)
        if( orientation && (orientation*det<0) ) return false;
        
        // check that the scaling is close enough to the one measured from pinhole distances.
        double scl = sqrt( fabs(det) );
        if (scl < approximateScale/maxScaleDiff || scl > approximateScale*maxScaleDiff ) {
            return false;
        }
        
        // check that the midpoint of img1 is mapped close enough to the midpoint of img2
        vector<Point2f> refMid(1, Point2f(imgSize1.width/2,imgSize1.height/2));
        vector<Point2f> mappedRefMid(1);
        perspectiveTransform( refMid, mappedRefMid, trans );
        double shift = norm( mappedRefMid[0] - Point2f(imgSize2.width/2,imgSize2.height/2) );
        if( shift > maxShift ) {
            //cout << "   Discarding: shift = " << shift << " maxShift = " << maxShift << endl;
            return false;
        }

        return true;
    }

    vector<int> getNearestIndex( const vector<Point2f>& pts, const vector<KeyPoint>& kps, double maxDistance=30 ) {

        int nPoints = pts.size();
        int nKP = kps.size();
        vector<int> ret(nPoints,-1);
        
        std::pair<int, double> nearest;

        for( int i=0; i<nPoints; ++i ) {
            nearest = make_pair (-1, 1E12);
            for (int j=0; j<nKP; ++j) {
                double dist = norm( pts[i] - kps[j].pt);
                if( fabs (dist) < nearest.second ) {
                    nearest.second = dist;
                    nearest.first = j;
                }
            }
            if (nearest.first >= 0 && nearest.second < maxDistance) {
                ret[i] = nearest.first;
            }
        }

        return std::move (ret);

    }

    // find approximate (affine) transforms by matching 3 central points in the object, to a reference + its nNeighbours nearest points
    vector<Mat> getInitializations( const vector<KeyPoint>& ref, const vector<KeyPoint>& target, const Size& imgSize1, const Size& imgSize2,
                                    double approximateScale, double maxShift = 10000, double maxScaleDiff=1.5, int orientation=0) {

        vector<Mat> initializations;
        size_t nRP = ref.size();
        size_t nTP = target.size();
        
        if( nRP < 3 || nTP < 3 ) {    // not enough input keypoints
            return std::move(initializations);
        }
        
        double minScale = approximateScale/maxScaleDiff;
        double maxScale = approximateScale*maxScaleDiff;

        vector<Point2f> refPoints(3);
        vector<Point2f> targetPoints(3);
        vector<Point2f> allRefPoints,mappedRefPoints;

        KeyPoint::convert(ref, allRefPoints);
        mappedRefPoints.resize( nRP );

        std::map<vector<int>, Mat> maps;

        for( size_t r0=0; r0<nRP-2; ++r0 ) {
        for( size_t r1=r0+1; r1<nRP-1; ++r1 ) {
        for( size_t r2=r1+1; r2<nRP; ++r2 ) {
            refPoints[0] = ref[r0].pt;
            refPoints[1] = ref[r1].pt;
            refPoints[2] = ref[r2].pt;
            double norm01 = norm( refPoints[0] - refPoints[1] );
            double norm02 = norm( refPoints[0] - refPoints[2] );
            double norm12 = norm( refPoints[1] - refPoints[2] );
            for( size_t i=0; i<nTP; ++i ) {
                for( size_t j=0; j<nTP; ++j ) {
                    if( j == i ) continue;
                    double normij = norm( target[i].pt - target[j].pt );
                    if( (normij < norm01*minScale) || (normij > norm01*maxScale) ) continue;
                    for( size_t k=0; k<nTP; ++k ) {
                        if( k==i || k==j ) continue;
                        double normik = norm( target[i].pt - target[k].pt );
                        if( (normik < norm02*minScale) || (normik > norm02*maxScale) ) continue;
                        double normjk = norm( target[j].pt - target[k].pt );
                        if( (normjk < norm12*minScale) || (normjk > norm12*maxScale) ) continue;
                        targetPoints[0] = target[i].pt;
                        targetPoints[1] = target[j].pt;
                        targetPoints[2] = target[k].pt;
                        Mat trans = getAffineTransform( refPoints, targetPoints );
                        Mat H_init = Mat::eye(3, 3, CV_32F);            // unit
                        trans.copyTo(H_init(cv::Rect_<int> (0, 0, 3, 2))); // copy (affine) initialization to the top 2 rows
                        if (transformCheck (H_init, imgSize1, imgSize2, approximateScale, maxShift, maxScaleDiff, orientation)) {
                            perspectiveTransform( allRefPoints, mappedRefPoints, H_init );
                            maps.emplace( getNearestIndex(mappedRefPoints,target), H_init);
                        }
                    }
                }
            }
        }}}
        
        for( auto &it: maps ) initializations.push_back( it.second );

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

    vector<DMatch> matchNearest (const vector<KeyPoint>& kp1, const vector<KeyPoint>& kp2, const Mat& H, double maxDistance=30, double scale=1.0) {

        vector<DMatch> matches;
        vector<Point2f> points1, points2, mappedPoints;

        KeyPoint::convert (kp1, points1);
        KeyPoint::convert (kp2, points2);

        size_t nPoints = points1.size();
        mappedPoints.resize(nPoints);
        perspectiveTransform( points1, mappedPoints, H );
        for( size_t i=0; i<nPoints; ++i ) {
            double nearestDistance = 1E12;
            int nearestIndex = -1;
            for( size_t j=0; j<points2.size(); ++j ) {
                double dist = norm( mappedPoints[i] - points2[j] );
                if( dist < nearestDistance) {
                    nearestDistance = dist;
                    nearestIndex = j;
                }
            }
            if( nearestIndex >= 0 && (nearestDistance < maxDistance*scale) ) {
                matches.push_back( DMatch(i, nearestIndex, nearestDistance) );
            }
        }

        nPoints = points2.size();
        mappedPoints.resize(nPoints);
        perspectiveTransform( points2, mappedPoints, H.inv() );
        for( size_t i=0; i<nPoints; ++i ) {
            double nearestDistance = 1E12;
            int nearestIndex = -1;
            for( size_t j=0; j<points1.size(); ++j ) {
                double dist = norm( mappedPoints[i] - points1[j] );
                if( dist < nearestDistance ) {
                    nearestDistance = dist;
                    nearestIndex = j;
                }
            }
            if (nearestIndex >= 0 && nearestDistance < maxDistance) {
                matches.push_back( DMatch(nearestIndex, i, nearestDistance) );
            }
        }

        std::sort( matches.begin(), matches.end());
        set<int> allQueryInd;
        set<int> allTrainInd;
        matches.erase( std::remove_if( matches.begin(), matches.end(),
                       [&]( const DMatch& m ) {
                            if( allQueryInd.count(m.queryIdx) || allTrainInd.count(m.trainIdx) ) return true;
                            allQueryInd.insert( m.queryIdx );
                            allTrainInd.insert( m.trainIdx );
                            return false;
                        }), matches.end() );
        
        return std::move (matches);

    }
#endif

}


IDL_VPTR redux::img_align (int argc, IDL_VPTR* argv, char* argk) {
	
	static_assert( RDX_WITH_OPENCV == 1, "redux has to be compiled with OpenCV support");

    KW_RESULT kw;
    kw.by_size = 0;
    kw.eps = 1E-3;
    kw.help = 0;
    kw.margin = 30;
    kw.max_dist = 80;
    kw.max_scale = 1.04;
    kw.max_points = -1;
    kw.max_shift = 200;
    kw.niter = 30;
    kw.nrefpoints = 4;
    kw.threshold = 0;
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if( kw.status ) {
        IDL_ALLTYPES tmp;
        tmp.i = 0;
        IDL_StoreScalar( kw.status, IDL_TYP_INT, &tmp );
    }
    
    IDL_VPTR ret;
    IDL_MakeTempArray (IDL_TYP_FLOAT, 2, dims3x3, IDL_ARR_INI_ZERO, &ret);
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

    Mat imgIn1 = getImgAsGrayFloat( img1_in, kw.verbose );
    Mat imgIn2 = getImgAsGrayFloat( img2_in, kw.verbose );
    
    cv::Size imgSize1 = imgIn1.size();
    cv::Size imgSize2 = imgIn2.size();
    int maxSize = max( imgIn1.size().height, imgIn1.size().width );
    maxSize = max( maxSize, imgIn2.size().height );
    maxSize = max( maxSize, imgIn2.size().width );
    
    cv::Size imgSize(maxSize, maxSize);   // The smallest square that fits both images

    std::vector<float> thresholds;
    
    cv::SimpleBlobDetector::Params params;
    params.thresholdStep = 5;
    params.maxThreshold = 255;
    params.minDistBetweenBlobs = 20;
    params.filterByColor = false;
    params.filterByInertia = false;
    params.filterByConvexity = false;
    params.filterByArea = true;
    params.minArea = 10;
    params.maxArea = 150;
    

    try {

        Mat imgByte( imgSize, CV_8UC1 );
        Mat imgByte1( imgSize1, CV_8UC1 );
        Mat imgByte2( imgSize2, CV_8UC1 );
        Mat tmp1( imgSize1, CV_8UC1 );
        Mat tmp2( imgSize2, CV_8UC1 );
        Mat H, H_init;

#ifdef DEBUG_DLM_
        Mat result( imgSize, CV_8UC3 );
#endif

        imgIn1.convertTo( imgByte1, CV_8UC1, 255 );
        imgIn2.convertTo( imgByte2, CV_8UC1, 255 );
        
        kw.margin = std::max<IDL_INT>( kw.margin, 0 );
        int m2 = 2*kw.margin;
        Rect ROI1( kw.margin, kw.margin, imgSize1.width-m2, imgSize1.height-m2 );
        Rect ROI2( kw.margin, kw.margin, imgSize2.width-m2, imgSize2.height-m2 );
        if( kw.margin > 0 ) {
            if( m2 >= imgSize1.width || m2 >= imgSize1.height ||
                m2 >= imgSize2.width || m2 >= imgSize2.height ) {
                cerr << "img_align: margin is too large." << endl;
                return ret;
            }
            if ( kw.verbose > 1 ) cout << "img_align: applying margin: " << kw.margin << endl;
            Mat mask1( imgByte1.size(), CV_8UC1, Scalar::all(1) );
            mask1(ROI1).setTo( Scalar::all(0) );
            imgByte1.setTo( Scalar::all(0), mask1 );
            Mat mask2( imgByte2.size(), CV_8UC1, Scalar::all(1) );
            mask2(ROI2).setTo( Scalar::all(0) );
            imgByte2.setTo( Scalar::all(0), mask2 );
        }
        
        if( kw.threshold ) {
            thresholds = getAsVector<float>(kw.threshold);
            if( thresholds.size() ) {
                thresholds.resize( 2, thresholds[0] );  // if only 1 value, it will be copied
                thresholds[0] = std::min( std::max(thresholds[0], 0.0f), 1.0f );
                thresholds[1] = std::min( std::max(thresholds[1], 0.0f), 1.0f );
                if ( kw.verbose > 1 ) cout << "img_align: applying " << printArray(thresholds,"threshold") << endl;
                threshold( imgByte1, imgByte1, thresholds[0]*255, 0, THRESH_TOZERO );
                threshold( imgByte1, imgByte1, thresholds[1]*255, 0, THRESH_TOZERO );
            }
        }
        
        if( kw.smooth ) {
            if( kw.smooth%2 == 0 ) kw.smooth--;     // size has to be odd in the Gaussian filtering
            if ( kw.verbose > 1 ) cout << "img_align: applying Gaussian filter of size " << kw.smooth << endl;
            GaussianBlur( imgByte1, imgByte1, Size(kw.smooth,kw.smooth), 0, 0 );
            GaussianBlur( imgByte2, imgByte2, Size(kw.smooth,kw.smooth), 0, 0 );
        }
        
        double otsu1 = threshold( imgByte1, tmp1, 0, 0, THRESH_TOZERO|THRESH_OTSU );
        double otsu2 = threshold( imgByte2, tmp2, 0, 0, THRESH_TOZERO|THRESH_OTSU );
        
        if ( kw.verbose > 1 ) cout << "img_align: Detection cutoffs: [" << otsu1 << "," << otsu2 << "]" << endl;
        
        vector<KeyPoint> keypoints1, keypoints2;
        
        params.minThreshold = otsu1/4;

        Ptr<SimpleBlobDetector> detector;
        
#if CV_MAJOR_VERSION >= 3
        detector = SimpleBlobDetector::create(params);
#else
        detector = Ptr<SimpleBlobDetector>( new SimpleBlobDetector(params) );
#endif
        detector->detect( imgByte1, keypoints1 );
        
        params.minThreshold = otsu2/4;
#if CV_MAJOR_VERSION >= 3
        detector = SimpleBlobDetector::create(params);
#else
        detector = Ptr<SimpleBlobDetector>( new SimpleBlobDetector(params) );
#endif
        detector->detect( imgByte2, keypoints2 );
        
        //imgIn1.convertTo( imgByte1, CV_8UC1, 255 );
        //imgIn2.convertTo( imgByte2, CV_8UC1, 255 );
        
        copyMakeBorder( imgByte1, imgByte, 0, maxSize-imgSize1.height, 0, maxSize-imgSize1.width, BORDER_CONSTANT, Scalar::all(0) );
        imgByte.copyTo(imgByte1);
        copyMakeBorder( imgByte2, imgByte, 0, maxSize-imgSize2.height, 0, maxSize-imgSize2.width, BORDER_CONSTANT, Scalar::all(0) );
        imgByte.copyTo(imgByte2);
    
        Point2f mid1(imgIn1.cols/2, imgIn1.rows/2);
        Point2f mid2(imgIn2.cols/2, imgIn2.rows/2);

        if ( kw.verbose > 1 ) {
            cout << "Detected " << keypoints1.size() << " keypoints in image 1 and " << keypoints2.size() << " keypoints in image 2" << endl;
        }
        
        if( keypoints1.size() < 3 || keypoints2.size() < 3 ) {
            if( kw.status ) {
                IDL_ALLTYPES tmp;
                tmp.i = -1;
                IDL_StoreScalar( kw.status, IDL_TYP_INT, &tmp );
            } else {
                cout << "Not enough points for fitting, try to adjust threshold (" <<
                    printArray(thresholds,"current") << ")." << endl;
            }
            return ret;
        }
        
        vector<float> nearestNeighbours1;
        for( KeyPoint& kp: keypoints1 ) {   // calculate something roughly proportional to the power of each peak
            double nearestNeighbour = numeric_limits<float>::max();
            for( KeyPoint& tmpkp: keypoints1 ) {
                if( &tmpkp != &kp ) {
                    nearestNeighbour = min(nearestNeighbour,cv::norm(kp.pt - tmpkp.pt));
                }
            }
            nearestNeighbour = min(nearestNeighbour,cv::norm(kp.pt-mid1));
            nearestNeighbours1.push_back(nearestNeighbour);
            float val = imgByte1.at<uchar>( kp.pt );
            kp.response = kp.size * val;
        }
        std::nth_element( nearestNeighbours1.begin(), nearestNeighbours1.begin() + nearestNeighbours1.size()/2, nearestNeighbours1.end() );
        float medianNN1 = *( nearestNeighbours1.begin() + nearestNeighbours1.size() / 2 );

        vector<float> nearestNeighbours2;
        for( KeyPoint& kp: keypoints2 ) {
            double nearestNeighbour = numeric_limits<float>::max();
            for( KeyPoint& tmpkp: keypoints2 ) {
                if( &tmpkp != &kp ) {
                    nearestNeighbour = min(nearestNeighbour,cv::norm(kp.pt - tmpkp.pt));
                }
            }
            nearestNeighbours2.push_back(nearestNeighbour);
            float val = imgByte2.at<uchar>( kp.pt );
            kp.response = kp.size * val;
        }
        std::nth_element( nearestNeighbours2.begin(), nearestNeighbours2.begin() + nearestNeighbours2.size()/2, nearestNeighbours2.end() );
        float medianNN2 = *( nearestNeighbours2.begin() + nearestNeighbours2.size() / 2 );

        float approximateScale = medianNN2/medianNN1;
        
        // sort w.r.t. something proportional to the power of each peak
        std::sort (keypoints1.begin(), keypoints1.end(),
            [] (const KeyPoint& a, const KeyPoint& b) {
                return (a.response > b.response);
        });
        std::sort (keypoints2.begin(), keypoints2.end(),
            [] (const KeyPoint& a, const KeyPoint& b) {
                return (a.response > b.response);
        });
        
        vector<KeyPoint> selectedKP1( keypoints1.begin(), keypoints1.end() );
        vector<KeyPoint> selectedKP2( keypoints2.begin(), keypoints2.end() );
        if( selectedKP1.size() > size_t(3*kw.nrefpoints) ) selectedKP1.resize(3*kw.nrefpoints);
        if( selectedKP2.size() > size_t(3*kw.nrefpoints) ) selectedKP2.resize(3*kw.nrefpoints);
        
        // sort w.r.t. distance from center
        std::sort (selectedKP1.begin(), selectedKP1.end(),
            [&] (const KeyPoint& a, const KeyPoint& b) {
                return (cv::norm(a.pt - mid1) < cv::norm(b.pt - mid1));
        });
        std::sort (selectedKP2.begin(), selectedKP2.end(),
            [&] (const KeyPoint& a, const KeyPoint& b) {
                return (cv::norm(a.pt - mid2) < cv::norm(b.pt - mid2));
        });
        if( selectedKP1.size() > size_t(kw.nrefpoints) ) selectedKP1.resize(kw.nrefpoints);
        if( selectedKP2.size() > size_t(kw.nrefpoints+2) ) selectedKP2.resize(kw.nrefpoints+2);
        
        std::vector<Mat> initializations;
        if( kw.h_init ) {
            if( kw.h_init->type == IDL_TYP_UNDEF ) {
                IDL_VPTR tmp;
                IDL_MakeTempArray (IDL_TYP_FLOAT, 2, dims3x3, IDL_ARR_INI_ZERO, &tmp);
                IDL_VarCopy( tmp, kw.h_init );
            }
            H_init = arrayToMat( kw.h_init );
            Mat H_initf;
            H_init.assignTo(H_initf, CV_32F);
            if( transformCheck( H_initf, imgSize1, imgSize2, approximateScale, kw.max_shift, kw.max_scale, kw.orientation ) ) {
                initializations.push_back( H_initf );
            }
        }

        if( initializations.empty() ) {
            initializations = getInitializations( selectedKP1, selectedKP2, imgSize1, imgSize2, approximateScale, kw.max_shift, kw.max_scale, kw.orientation);
        }

        vector< DMatch > matches;
        if ( initializations.size() ) {
            
            if ( kw.verbose > 1 ) {
                cout << "Matching the " << kw.nrefpoints << " largest keypoints in the images." << endl;
                cout << "Restricting the transformation to have shifts smaller than " << kw.max_shift
                    << " and a maximal scaling of " << (approximateScale*kw.max_scale) << " gave " << initializations.size()
                    << " valid permutations." << endl;
            }

            size_t nMatches(0);
            for( auto& h: initializations ) {       // find the mapping which can pair the most keypoints.
                vector< DMatch > tmpmatches = matchNearest( keypoints1, keypoints2, h, medianNN1/3, approximateScale );
                vector<Point2f> obj, scene;
                for( auto & m: tmpmatches ) {
                    obj.push_back( keypoints1[ m.queryIdx ].pt );
                    scene.push_back( keypoints2[ m.trainIdx ].pt );
                }
                Mat HH = findHomography( obj, scene, CV_LMEDS );
                HH.assignTo( h, CV_32F );
                size_t nm = tmpmatches.size();
                h.at<float>(2,2) = static_cast<float>(nm);  // temporarily store nMatches in the map
                if( nm > nMatches ) {
                    matches = std::move(tmpmatches);
                    nMatches = nm;
                }
            }
            std::sort (initializations.begin(), initializations.end(),  // sort by number of succesfully paired pinholes.
                [] (const Mat& a, const Mat& b) {
                    return (a.at<float>(2,2) > b.at<float>(2,2));
            });
            for( auto& h: initializations ) {
                h.at<float>(2,2) = 1.0;  // reset the temporary value 
            }
            if( initializations.size() > 8 ) initializations.resize( 8 );   // keep the 8 maps which could pair the most pinholes
        }
        
        // Create a mask from the best mapping found above.
        imgByte.setTo( Scalar::all(0) );
        for( auto& m: matches ) {
            KeyPoint& kp = keypoints2[ m.trainIdx ];
            circle( imgByte, kp.pt, 3*kp.size, Scalar(255), -1 );
        }
        
        std::map<double, Mat> results;
        TermCriteria term_crit = TermCriteria(TermCriteria::COUNT + TermCriteria::EPS, kw.niter, kw.eps);
        auto fit_func = std::bind( findTransformECC, std::ref(imgByte1), std::ref(imgByte2), std::placeholders::_1,
            MOTION_HOMOGRAPHY, term_crit, imgByte );
        
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
            if( kw.status ) {
                IDL_ALLTYPES tmp;
                tmp.i = -2;
                IDL_StoreScalar( kw.status, IDL_TYP_INT, &tmp );
            } else {
                cout << "Failed to match detected pinholes, try to increase nref (current: "
                    << kw.nrefpoints << ")." << endl;
            }
            return ret;
        }
        
        H = results.rbegin()->second;
        matches = matchNearest( keypoints1, keypoints2, H, medianNN1/3, approximateScale );
        
        if( kw.verbose > 1 ) {
            vector<double> cc;
            for (auto r : results) cc.push_back (r.first);
            std::sort(cc.rbegin(),cc.rend());
            cout << printArray (cc, "correlations") << endl;
            cout << matches.size() << " pairs matched using a max_distance of " << kw.max_dist << " for pairing." << endl;
        }

        if( kw.max_points < 0 ) kw.max_points = 25;   // by default, do refinement by using the strongest 25 pinholes

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

        if( matches.size() > 3 ) {     // at least 4 points needed for findHomography
            
            if( kw.max_points > 3 ) {

                if( kw.by_size ) {                              // sort by size*intensity
                    std::sort( matches.begin(), matches.end(),
                        [&] (const DMatch& a, const DMatch& b) {
                            return (keypoints1[ a.queryIdx ].response*keypoints2[ a.trainIdx ].response >
                                    keypoints1[ b.queryIdx ].response*keypoints2[ b.trainIdx ].response);
                    });
                } else {                                        // sort according to distance from image centre
                    std::sort( matches.begin(), matches.end(),
                        [&] (const DMatch& a, const DMatch& b) {
                            return (norm(keypoints1[ a.queryIdx ].pt - mid1) > norm(keypoints1[ b.queryIdx ].pt - mid1));
                    });
                }

                if( size_t(kw.max_points) < matches.size() ) matches.resize( kw.max_points );
            }
        
            if( kw.verbose > 1 ) {
                cout << "Using " << matches.size() << " pairs to refine the fit." << endl;
            }

            vector<Point2f> obj, scene;
            for( auto & m: matches ) {
                obj.push_back (keypoints1[ m.queryIdx ].pt);
                scene.push_back (keypoints2[ m.trainIdx ].pt);
            }

            double max1, max2;
            double metric1 = metric (H, obj, scene, &max1);
            H = findHomography( obj, scene, CV_LMEDS );
            double metric2 = metric (H, obj, scene, &max2);
         
            if (kw.verbose > 1) {
                cout << "   -> error (avg,max): (" << metric1 << "," << max1 << ") -> (" << metric2 << "," << max2 << ")" << endl;
            }
        }

        if (kw.verbose) {
            cout << "Transformation matrix:\n" << H << endl;
            cout << "Scale: " << sqrt(fabs(determinant (H (cv::Rect_<int> (0, 0, 2, 2))))) << endl;
            vector<Point2f> refMid(1, Point2f(imgSize1.width/2,imgSize1.height/2));
            vector<Point2f> mappedRefMid(1);
            perspectiveTransform( refMid, mappedRefMid, H );
            Point2f diff =  mappedRefMid[0] - Point2f(imgSize2.width/2,imgSize2.height/2);
            double shift = norm( diff );
            cout << "Shift: " << shift << " " << diff << endl;
         }
        
        H.assignTo( retMat, retMat.type() );
        
        if( kw.h_init ) {
            H.assignTo( H_init, H_init.type() );
        }
        
        return ret;

    } catch( const cv::Exception& e ) {
        cout << "OpenCV error: " << e.msg << endl;
    }

    return ret;

    
}


IDL_VPTR redux::img_project (int argc, IDL_VPTR* argv, char* argk) {

#ifndef RDX_WITH_OPENCV
    cerr << "img_project: redux has to be re-compiled with OpenCV enabled to be able to use this function." << endl;
    return IDL_GettmpInt(0);
#else

    KW_RESULT kw;
    kw.help = 0;
    kw.nrefpoints = 4;
    kw.verbose = 0;
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

    const Mat img = arrayToMat(img_in);
    const Mat H = arrayToMat(H_in);
    
    Size2f inSize = img.size();

    vector<Point2f> corners({{0,0}, {inSize.width,0}, {0,inSize.height}, {inSize.width,inSize.height}});
    vector<Point2f>  mappedCorners;
    perspectiveTransform( corners, mappedCorners, H );
    int maxX = mappedCorners[0].x;
    int minX = mappedCorners[0].x;
    int maxY = mappedCorners[0].y;
    int minY = mappedCorners[0].y;
    for( auto& c: mappedCorners ) {
        maxX = max<int>(maxX, c.x+0.5);
        minX = min<int>(minX, c.x+0.5);
        maxY = max<int>(maxY, c.y+0.5);
        minY = min<int>(minY, c.y+0.5);
    }

    IDL_VPTR out;
    IDL_ARRAY_DIM outDims = { maxX-minX+1, maxY-minY+1, 0, 0, 0, 0, 0, 0 };
    IDL_MakeTempArray (img_in->type, img_in->value.arr->n_dim, outDims, IDL_ARR_INI_NOP, &out);

    Mat outImg = arrayToMat(out);
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

#ifndef RDX_WITH_OPENCV
    cerr << "img_remap: redux has to be re-compiled with OpenCV enabled to be able to use this function." << endl;
    return IDL_GettmpInt(0);
#else
    KW_RESULT kw;
    kw.help = 0;
    kw.nrefpoints = 4;
    kw.verbose = 0;
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

#ifndef RDX_WITH_OPENCV
    cerr << "rdx_find_shift: redux has to be re-compiled with OpenCV enabled to be able to use this function." << endl;
    return IDL_GettmpInt(0);
#else
    
    IDL_VPTR ret;
    IDL_MEMINT dims[] = {3,2};
    IDL_MakeTempArray( IDL_TYP_FLOAT, 2, dims, IDL_ARR_INI_NOP, &ret );
    Mat retMat = arrayToMat( ret );

    KW_RESULT kw;
    kw.by_size = 0;
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
        cout << "OpenCv error: " << e.msg << endl;
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
    
#ifndef RDX_WITH_OPENCV
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



namespace apz {
                
    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT blend;
        IDL_INT help;
        IDL_INT margin;
        IDL_INT verbose;
    } KW_RESULT;


    // NOTE:  The keywords MUST be listed in alphabetical order !!
    IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "BLEND",          IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(blend) },
        { (char*) "HELP",           IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(help) },
        { (char*) "MARGIN",         IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(margin) },
        { (char*) "VERBOSE",        IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF(verbose) },
        { NULL }
    };


    string make_win_info( int lvl ) {
        string ret = "RDX_MAKE_WINDOW";
        if( lvl > 0 ) {
            ret += ((lvl > 1)?"\n":"    ");          // newline if lvl>1
            ret += "   Syntax:   win = rdx_make_window( size, /KEYWORDS )\n";
            ret += ((lvl > 1)?"":"                   ");
            ret += "             win = rdx_make_window( columnns, rows, /KEYWORDS )\n";
            if( lvl > 1 ) {
                ret +=  "   Accepted Keywords:\n"
                        "      HELP                Display this info.\n"
                        "      BLEND               Size of smoothing area. (default=12.5% of image width/height).\n"
                        "      MARGIN              Border size. (default = 0, i.e. smoothing all the way to the edge).\n"
                        "      VERBOSE             Verbosity, default is 0 (only error output).\n";
            }
        } else ret += "\n";
        return ret;
    }


}


IDL_VPTR rdx_make_win( int argc, IDL_VPTR* argv, char* argk ) {
    
    apz::KW_RESULT kw;
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, apz::kw_pars, (IDL_VPTR*)0, 255, &kw );

    if( kw.help || nPlainArgs < 1 ) {
        cout << apz::make_win_info(2) << endl;
        return IDL_GettmpInt(0);
    }
    
    IDL_LONG nCols = IDL_LongScalar(argv[0]);
    IDL_LONG nRows = nCols;
    if( nPlainArgs > 1 ) {
        nRows = IDL_LongScalar(argv[1]);
    }

    if( kw.blend == 0 ) {
        kw.blend = (nRows + nCols) / 16;
    }
    
    IDL_VPTR tmp;
    IDL_MEMINT dims[] = { nCols, nRows };
    float* tmpData = (float*)IDL_MakeTempArray( IDL_TYP_FLOAT, 2, dims, IDL_ARR_INI_NOP, &tmp );
    
    std::fill( tmpData, tmpData+nCols*nRows, 1.0 );
    
    shared_ptr<float*> tmpArr = reshapeArray( tmpData, nRows, nCols );
    apodizeInPlace( tmpArr.get(), nRows, nCols, kw.blend, kw.blend, kw.margin, kw.margin );
    
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
        IDL_INT filter;
        float fp_thres;
        float limit;
        IDL_VPTR framenumbers;
        IDL_VPTR discarded;
        IDL_VPTR summed;
        IDL_VPTR nsummed;
        IDL_VPTR dark;
        IDL_VPTR gain;
        IDL_VPTR bs_gain;
        IDL_VPTR bs_psf;
        IDL_VPTR time_beg;
        IDL_VPTR time_end;
        IDL_VPTR time_avg;
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
        { (char*) "DISCARDED",        IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,discarded) },
        { (char*) "FILLPIX_THRESHOLD",IDL_TYP_FLOAT, 1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(SI_KW,fp_thres) },
        { (char*) "FILTER",           IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(SI_KW,filter) },
        { (char*) "FRAMENUMBERS",     IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,framenumbers) },
        { (char*) "GAIN",             IDL_TYP_UNDEF, 1, IDL_KW_VIN|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,gain) },
        { (char*) "HELP",             IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(SI_KW,help) },
        { (char*) "LIMIT",            IDL_TYP_FLOAT, 1, 0,                      0, (char*) IDL_KW_OFFSETOF2(SI_KW,limit) },
        { (char*) "LUN",              IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(SI_KW,lun) },
        { (char*) "NSUMMED",          IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,nsummed) },
        { (char*) "NTHREADS",         IDL_TYP_BYTE,  1, 0,                      0, (char*) IDL_KW_OFFSETOF2(SI_KW,nthreads) },
        { (char*) "PADDING",          IDL_TYP_INT,   1, 0,                      0, (char*) IDL_KW_OFFSETOF2(SI_KW,padding) },
        { (char*) "PINHOLE_ALIGN",    IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(SI_KW,pinh_align) },
        { (char*) "SUMMED",           IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,summed) },
        { (char*) "TIME_AVG",         IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,time_avg) },
        { (char*) "TIME_BEG",         IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,time_beg) },
        { (char*) "TIME_END",         IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(SI_KW,time_end) },
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
                    "      FILTER              Size of median-filter to be applied to means before checking. (3)\n"
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
    
    if( kw.filter<2 && kw.check ) {
        kw.filter = 3;
    }
    
    kw.nthreads = max<UCHAR>(1, min<UCHAR>(kw.nthreads, thread::hardware_concurrency()));
    kw.padding = max<IDL_INT>(0, min<IDL_INT>(kw.padding, 4096));   // prevent insane padding
    
    try {
        
        IDL_VPTR ret;
        shared_ptr<double> darkData, gainData, bsGainData;
        shared_ptr<uint8_t> maskData;
        shared_ptr<uint8_t*> mask2D;
        shared_ptr<double> bsOtfData;
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
                make_mask( gainPtr, maskData.get(), yCalibSize, xCalibSize, 0, 5, true, true ); // filter away larger features than ~5 pixels and invert
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
                    bsOtfData.reset( (double*)fftw_malloc(ySizePadded*(xSizePadded/2+1)*sizeof(fftw_complex)), fftw_free );
                    bsOtfPtr = reinterpret_cast<fftw_complex*>(bsOtfData.get());
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
        
        if( kw.filter > 0 ) {
            kw.check = 1;
        }
        
        if( kw.check && nImages < 3 ) {
            cerr << "rdx_sumimages: Not enough statistics, skipping check." << endl;
            kw.check = kw.filter = 0;
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
        
#ifdef RDX_WITH_OPENCV
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
                        descatter( myTmpPtr, ySize, xSize, bsGainPtr, bsOtfPtr, ySizePadded, xSizePadded, 1, 50, 1E-8 );
                    }
                    
                    if ( mask2D ) {
                        shared_ptr<double*> tmp2D = reshapeArray( myTmpPtr, ySize, xSize );
                        fillPixels( tmp2D.get(), (size_t)ySize, (size_t)xSize, mask2D.get() );
                    }
                    
#ifdef RDX_WITH_OPENCV
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
        
#ifdef RDX_WITH_OPENCV
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
                descatter( summedData, ySize, xSize, bsGainPtr, bsOtfPtr, ySizePadded, xSizePadded, kw.nthreads, 50, 1E-8 );
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
                    "      DISCARDED           (output) Framenumbers of the discarded frames.\n"
                    "      FILTER              Size of median-filter to be applied to means before checking. (3)\n"
                    "      FRAMENUMBERS        (output) Framenumbers of the summed frames.\n"
                    "      GAIN                Gain table (inverted flat-field).\n"
                    "      LIMIT               Allowed deviation from the median. (0.0175)\n"
                    "      LUN                 IDL file unit (id) where discarded files will be logged.\n"
                    "      NSUMMED             (output) Number of images actually summed.\n"
                    "      NTHREADS            Number of threads.\n"
                    "      PADDING             Padding size for the descattering procedure. (256)\n"
                    "      PINHOLE_ALIGN       Do sub-pixel alignment before summing.\n"
                    "      SUMMED              (output) Raw sum.\n"
                    "      TIME_BEG            (output) Begin-time (for the first frame).\n"
                    "      TIME_END            (output) End-time (for the last frame).\n"
                    "      TIME_AVG            (output) Average timestamp from file-headers.\n"
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

#ifndef RDX_WITH_OPENCV
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
    
    if( kw.filter<2 && kw.check ) {
        kw.filter = 3;
    }
    
    kw.nthreads = max<UCHAR>(1, min<UCHAR>(kw.nthreads, thread::hardware_concurrency()));
    kw.padding = max<IDL_INT>(0, min<IDL_INT>(kw.padding, 4096));   // prevent insane padding

    try {
        
        IDL_VPTR ret;
        shared_ptr<double> darkData, gainData, bsGainData;
        shared_ptr<uint8_t> maskData;
        shared_ptr<uint8_t*> mask2D;
        shared_ptr<double> bsOtfData;
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
                make_mask( gainPtr, maskData.get(), yCalibSize, xCalibSize, 0, 5, true, true ); // filter away larger features than ~5 pixels and invert
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
                    bsOtfData.reset( (double*)fftw_malloc( ySizePadded*(xSizePadded/2+1)*sizeof(fftw_complex) ), fftw_free );
                    bsOtfPtr = reinterpret_cast<fftw_complex*>(bsOtfData.get());
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

        size_t nFiles = existingFiles.size();
        //kw.nthreads = static_cast<UCHAR>( min<size_t>(kw.nthreads, nFiles) );

        if( !nFiles ) { 
            cout << "rdx_sumfiles: No input files." << endl;
            return IDL_GettmpInt(0);
        }
        // Get size etc from first image.
        std::shared_ptr<redux::file::FileMeta> meta = redux::file::getMeta( existingFiles[0] );
        if( !meta ) {        // Only allow 2D images
            cout << "rdx_sumfiles: Failed to get metadata." << endl;
            return IDL_GettmpInt(0);
        }
        size_t nDims = meta->nDims();
        IDL_MEMINT xSize(0);
        IDL_MEMINT ySize(0);
        size_t nTotalFrames( nFiles );
        vector<size_t> nFrames( nFiles, 1 );
        if( nDims == 2 ) {
            xSize = meta->dimSize(1);
            ySize = meta->dimSize(0);
        } else if( nDims == 3 ) {
            xSize = meta->dimSize(2);
            ySize = meta->dimSize(1);
        } else {                                // Only allow 2D & 3D files
            cout << "rdx_sumfiles: Only 2D and 3D files supported." << endl;
            return IDL_GettmpInt(0);
        }
        
        size_t nPixels = xSize*ySize;
        size_t frameSize = nPixels * meta->elementSize();
        UCHAR dataType = meta->getIDLType();
        size_t xSizePadded = xSize + 2*kw.padding;
        size_t ySizePadded = ySize + 2*kw.padding;
        
        bool done(false);
        mutex mtx;
        boost::asio::io_service ioService;
        std::shared_ptr<boost::asio::io_service::work> workLoop( new boost::asio::io_service::work(ioService) );
        boost::thread_group pool;
        for( uint16_t t=0; t < kw.nthreads; ++t ) {
            pool.create_thread( [&](){
                while( !done ) {
                    try {
                        ioService.run();
                    } catch( exception& e ) {
                        cerr << "Exception in summing thread: " << e.what() << endl;
                    } catch( ... ) {
                        cerr << "Unhandled exception in thread." << endl;
                    }
                }
            });
        }
        
        ProgressWatch progWatch;
        if( nDims == 3 ) {
            nFrames.clear();
            nTotalFrames = 0;
            progWatch.set( nFiles );
            for( auto &fn: existingFiles ) {
                ioService.post([&](){
                    auto tmpMeta = redux::file::getMeta( fn );
                    if( meta ) {
                        size_t frames = tmpMeta->dimSize(0);
                        unique_lock<mutex> lock(mtx);
                        nFrames.push_back( frames );
                        nTotalFrames += frames;
                    }
                    ++progWatch;
                });
            }
            progWatch.wait();
        }
        
        size_t maxFileSize = *std::max_element( nFrames.begin(), nFrames.end() ) * frameSize;
        atomic<size_t> nSummed(0);
        atomic<size_t> frameIndex(0);
        
        vector<int32_t> frameNumbers( nTotalFrames );
        vector<bpx::ptime> time_beg;
        vector<bpx::ptime> time_end;
        if( kw.time_beg || kw.time_end || kw.time_avg ) {
            time_beg.resize( nTotalFrames );
            time_end.resize( nTotalFrames );
        }
        
        string statusString;
        if( kw.verbose ) {
            statusString = (kw.check?"Checking and summing ":"Summing ") + to_string(nTotalFrames)
            + " frames using " +to_string((int)kw.nthreads) + string(" thread") + ((kw.nthreads>1)?"s.":":");
            cout << statusString << ((kw.verbose == 1)?"\n":"") << flush;
        }
        
        IDL_MEMINT dims[] = { xSize, ySize }; 
        double* summedData = (double*)IDL_MakeTempArray( IDL_TYP_DOUBLE, 2, dims, IDL_ARR_INI_ZERO, &ret ); //IDL_ARR_INI_ZERO
        
        unique_ptr<double[]> checked;
        unique_ptr<double[]> sums( new double [ nPixels*kw.nthreads ] );
        unique_ptr<double[]> tmp( new double [ 2*nPixels*kw.nthreads ] );
        
        shared_ptr<float*> shiftsData;
        float** shifts = nullptr;
        if( kw.pinh_align ) {
            shiftsData = sharedArray<float>( nTotalFrames, 2 );
            shifts = shiftsData.get();
        }
        
        if( kw.filter > 1 ) {
            kw.check = 1;
        }
        
        double* sumPtr = sums.get();
        double* tmpPtr = tmp.get();
        memset( sumPtr, 0, nPixels*kw.nthreads*sizeof(double) );
        double* checkedPtr = nullptr;
        if( kw.check ) {
            if( nTotalFrames < 3 ) {
                cerr << "rdx_sumimages: Not enough statistics, skipping check." << endl;
                kw.check = kw.filter = 0;
            } else if( kw.pinh_align ) {
                cerr << "rdx_sumimages: Checking together with aligning is not yet supported, skipping check." << endl;
                kw.check = kw.filter = 0;
            } else {
                checked.reset( new double [ 2*nTotalFrames ] );
                checkedPtr = checked.get();
            }
        }
        atomic<size_t> threadIndex(0);

        Mat refImg;
        promise<cv::Rect> subImgROI;
        shared_future<cv::Rect> fut = subImgROI.get_future();

        const int warp_mode = MOTION_TRANSLATION; // MOTION_EUCLIDEAN;
 
        int number_of_iterations = 50;
        double termination_eps = 1e-4;
        TermCriteria criteria(TermCriteria::COUNT+TermCriteria::EPS, number_of_iterations, termination_eps);

        int flags = INTER_CUBIC;                  // INTER_LINEAR, INTER_CUBIC, INTER_AREA, INTER_LANCZOS4
        int borderMode = BORDER_CONSTANT;
        const Scalar borderValue = Scalar();

        thread_local int myThreadIndex(-1);
        atomic<bool> refLoaded(false);
        
        function<void(size_t,shared_ptr<char>,size_t)> 
            sumFunc = [&]( size_t frameIndex, shared_ptr<char> threadBuffer, size_t frameOffset ) {
            
                // for pinh_align we need the reference (first) image before the rest, so re-post other images until we get it.
                if( kw.pinh_align && frameIndex && !refLoaded.load()  ) {
                    ioService.post( std::bind(sumFunc,frameIndex,threadBuffer, frameOffset) );
                    return;
                }
                
                if( myThreadIndex == -1 ) myThreadIndex = threadIndex.fetch_add(1);
                
                double* threadSum = sumPtr+myThreadIndex*nPixels;
                double* threadTmp = tmpPtr+myThreadIndex*nPixels;
                UCHAR* framePtr = reinterpret_cast<UCHAR*>( threadBuffer.get() + frameOffset );
                
                if( checkedPtr  ) {
                    bool hasInf;
                    checkedPtr[frameIndex] = getMinMaxMean( framePtr, nPixels, dataType, 0, 0, &hasInf );
                    if( hasInf ) {
                        checkedPtr[frameIndex] = std::numeric_limits<double>::infinity();
                        ++progWatch;
                        return;
                    }
                }
                
                if( kw.pinh_align ) {
                    
                    if( !frameIndex ) refLoaded.store(true);
                    
                    try {
                            
                        if( darkPtr || gainPtr ) {
                            applyDarkAndGain( framePtr, threadTmp, darkPtr, gainPtr, nPixels, dataType );
                        } else {
                            copyToRaw( framePtr, threadTmp, nPixels, dataType );
                        }
                        if( bsOtfPtr ) {
                            descatter( threadTmp, ySize, xSize, bsGainPtr, bsOtfPtr, ySizePadded, xSizePadded, 1, 50, 1E-8 );
                        }
                        
                        if ( mask2D ) {
                            shared_ptr<double*> tmp2D = reshapeArray( threadTmp, ySize, xSize );
                            fillPixels( tmp2D.get(), (size_t)ySize, (size_t)xSize, mask2D.get() );
                        }
                        
                        Mat cvImg( ySize, xSize, CV_64FC1, threadTmp );
                        cv::Rect roi;
                        if( !frameIndex ) {
                            
                            int margin = 100;
                            int refSize = 29;

                            roi = cv::Rect( cv::Point(margin, margin), cv::Point(xSize-margin, ySize-margin) );
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

//                             int nonZ = cv::countNonZero(refImg);
//                             while( cv::countNonZero(refImg) == nonZ ) {
//                                  refImg.adjustROI( -1, -1, -1, -1 );
//                                  refSize -= 2;
//                             }
//                             refImg.adjustROI( 2, 2, 2, 2 );
//                             refSize += 4;
//                             refImg = refImg.clone();
// 
//                             roi = cv::Rect( maxLoc-cv::Point(refSize/2,refSize/2), cv::Size(refSize,refSize) );
                            Mat( cvImg, roi ).convertTo( refImg, CV_32F );

                            subImgROI.set_value( roi );
                             
                        } else {

                            roi = fut.get();

                            Mat subImgFloat( roi.size(), CV_32FC1 );
                            double* slaskPtr = threadTmp+nPixels*kw.nthreads;
                            Mat slask( ySize, xSize, CV_64FC1, slaskPtr );
                            
                            Mat( cvImg, roi ).convertTo( subImgFloat, CV_32F );

                            Mat warp_matrix = Mat::eye( 2, 3, CV_32F );
                            findTransformECC( subImgFloat, refImg, warp_matrix, warp_mode, criteria );

                            shifts[frameIndex][0] = warp_matrix.at<float>(0,2);
                            shifts[frameIndex][1] = warp_matrix.at<float>(1,2);

                            warpAffine( cvImg, slask, warp_matrix, cvImg.size(), flags, borderMode, borderValue );

                            memcpy( threadTmp, slaskPtr, nPixels*sizeof(double) );
                            
                        }
                        
                        for( size_t n=0; n<nPixels; ++n ) threadSum[n] += threadTmp[n];
                        
                    } catch( const cv::Exception& e ) {
                        nSummed--;
                        if( checkedPtr  ) {
                            checkedPtr[frameIndex] = std::numeric_limits<double>::infinity();
                        }
                    }

                } else {
                    addToRaw( framePtr, threadSum, nPixels, dataType );
                }
                
                size_t ns = nSummed++;
                if( kw.verbose > 1 ) printProgress( statusString, (ns*100.0/(nTotalFrames-1)));
                ++progWatch;
        };

        progWatch.set( nFiles+nTotalFrames );
        size_t frameCount(0);
        for( size_t i=0; i<existingFiles.size(); ++i ) {
            ioService.post([&,i,frameCount](){
                shared_ptr<char> threadBuffer( new char[ maxFileSize ], []( char*& p ) { delete[] p; } );
                shared_ptr<redux::file::FileMeta> threadMeta;
                try {
                    readFile( existingFiles[i], threadBuffer.get(), threadMeta );
                    if( threadMeta ) {
                        if( !time_beg.empty() ) {
                            vector<bpx::ptime> startTimes = threadMeta->getStartTimes();
                            bpx::time_duration expTime = threadMeta->getExposureTime();
                            if( startTimes.size() == nFrames[i] ) {
                                std::copy( startTimes.begin(), startTimes.end(), time_beg.begin()+frameCount);
                                std::transform( startTimes.begin(), startTimes.end(), time_end.begin()+frameCount,
                                    [expTime](const bpx::ptime& pt){ return pt+expTime; }
                                );
                            }
                        }
                        if( kw.framenumbers ) {
                          vector<size_t> fn = threadMeta->getFrameNumbers();
                          if( fn.size() == nFrames[i] ) {
                              if( fn.front() ) {
                                  std::copy( fn.begin(), fn.end(), frameNumbers.begin()+frameCount );
                              } else {
                                  std::transform( fn.begin(), fn.end(), frameNumbers.begin()+frameCount,
                                      [frameCount](const size_t& a){ return a+frameCount; }
                                  );
                              }
                          }
                        }
                    }
                    size_t bufferOffset(0);
                    for( size_t j=0; j<nFrames[i]; ++j ) {
                        //sumFunc( frameCount+j, threadBuffer, bufferOffset);
                        ioService.post( std::bind(sumFunc,frameCount+j,threadBuffer, bufferOffset) );
                        bufferOffset += frameSize;
                    }

                } catch( const exception& e ) {
                    cout << "rdx_sumfiles: Failed to load file: " << existingFiles[i] << "  Reason: " << e.what() << endl;
                    progWatch.decreaseTarget( nFrames[i] );
                }
                ++progWatch;
            });
            frameCount += nFrames[i];
        }

        progWatch.wait();

        for( size_t t=0; t<kw.nthreads; ++t ) {
            std::transform( summedData, summedData+nPixels, sumPtr+t*nPixels, summedData, std::plus<double>() );
        }

        
        size_t nDiscarded(0);
        if( checkedPtr ) {
            set<size_t> discarded;
            memcpy( checkedPtr+nTotalFrames, checkedPtr, nTotalFrames*sizeof(double) );
            nth_element( checkedPtr, checkedPtr+nTotalFrames/2, checkedPtr+nTotalFrames );        // total median
            double tmean = kw.limit*checkedPtr[ nTotalFrames/2 ];
            memcpy( checkedPtr, checkedPtr+nTotalFrames, nTotalFrames*sizeof(double) );
            median_filter( checkedPtr, nTotalFrames, kw.filter );
            for( size_t i=0; i<nTotalFrames; ++i ) {
                if( isfinite(checkedPtr[i]) ) {
                    checkedPtr[i] = (abs(checkedPtr[i]-checkedPtr[nTotalFrames+i]) <= tmean);
                    if ( checkedPtr[i] > 0 ) continue;      // check passed
                } else checkedPtr[i] = 0;
                ++nDiscarded;
            }
            if( nDiscarded ) {
                memset( sumPtr, 0, nPixels*kw.nthreads*sizeof(double) );
                frameIndex = 0;
                threadIndex = 0;
                progWatch.set( 0 );
                size_t frameCount(0);
                for( size_t i=0; i<existingFiles.size(); ++i ) {
                    bool allOk(true);
                    for( size_t j=0; j<nFrames[i]; ++j ) {
                        if( checkedPtr[frameCount+j] == 0 ) {
                            if( kw.lun ) {
                                string msg = existingFiles[i];
                                if( nFrames[i] > 1 ) msg += ", frame # " + to_string(j);
                                IDL_PoutRaw( kw.lun, (char*)msg.c_str(), msg.length() );
                                IDL_PoutRaw( kw.lun, nl, 1 );
                            }
                            allOk = false;
                        }
                    }
                    if( !allOk ) {
                        progWatch.increaseTarget(1);
                        ioService.post([&,i,frameCount](){
                            shared_ptr<char> threadBuffer( new char[ maxFileSize ], []( char*& p ) { delete[] p; } );
                            shared_ptr<redux::file::FileMeta> threadMeta;
                            try {
                                readFile( existingFiles[0], threadBuffer.get(), threadMeta );
                                size_t bufferOffset(0);
                                for( size_t j=0; j<nFrames[i]; ++j ) {
                                    if( checkedPtr[frameCount+j] == 0 ) {     // subtract discarded images
                                        progWatch.increaseTarget(1);
                                        ioService.post( std::bind(sumFunc, frameCount+j, threadBuffer, bufferOffset) );
                                        --nSummed;
                                        frameNumbers[frameCount+j] *= -1;
                                        if( shifts ) {
                                            shifts[frameCount+j][0] = 0;
                                            shifts[frameCount+j][1] = 0;
                                        }
                                        if( !time_beg.empty() ) {
                                            time_beg[frameCount+j] = time_end[frameCount+j] = bpx::ptime();
                                        }
                                    }
                                    bufferOffset += frameSize;
                                }
                            } catch( const exception& e ) {
                                cout << "rdx_sumfiles: Failed to load file: " << existingFiles[i] << "  Reason: " << e.what() << endl;
                            }
                            ++progWatch;
                        });
                    }
                    frameCount += nFrames[i];
                }

                progWatch.wait();

                for( size_t t=0; t<kw.nthreads; ++t ) {
                    std::transform( summedData, summedData+nPixels, sumPtr+t*nPixels, summedData, std::minus<double>() );
                }

            }
        }
        
        done = true;
        workLoop.reset();
        ioService.stop();
        pool.join_all();

        vector<int32_t> discarded;
        if( kw.discarded || kw.framenumbers ) {
            for( size_t i=0; i<frameNumbers.size(); ) {
                int32_t val = frameNumbers[i];
                if( val < 0 ) {
                    discarded.push_back(-val);
                    frameNumbers.erase( frameNumbers.begin()+i );
                } else {
                    i++;
                }
            }
        }
        
        nDiscarded = discarded.size();
        nSummed = frameNumbers.size();
        if( nSummed+nDiscarded != nTotalFrames ) {
            cout << "  FrameCount mismatch: " << nSummed << "+" << nDiscarded << " != " << nTotalFrames << endl;
        }
        
        if( kw.pinh_align ) {
            Mat cvImg( ySize, xSize, CV_64FC1, summedData );
            Mat warp_matrix = Mat::eye( 2, 3, CV_32F );
            for( size_t n=1; n<nTotalFrames; ++n ) {
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
                cout << "  Check failed for " << nDiscarded << " frame" << ((nDiscarded>1)?"s.":".") << endl;
            } else {
                cout << "  All ok." << endl;
            }
        }
        
        if( kw.summed ) {
            IDL_VPTR tmp;
            double* tmpData = (double*)IDL_MakeTempArray( IDL_TYP_DOUBLE, 2, dims, IDL_ARR_INI_NOP, &tmp );
            memcpy( tmpData, summedData, nPixels*sizeof(double));
            IDL_VarCopy( tmp, kw.summed );
        }

        if( nDiscarded && kw.discarded ) {
            IDL_VPTR tmp;
            IDL_MEMINT nD = discarded.size();
            IDL_MEMINT dims[] = { nD }; 
            int32_t* tmpData = (int32_t*)IDL_MakeTempArray( IDL_TYP_LONG, 1, dims, IDL_ARR_INI_NOP, &tmp );
            memcpy( tmpData, discarded.data(), nD*sizeof(int32_t));
            IDL_VarCopy( tmp, kw.discarded );
        }
        
        if( kw.framenumbers ) {
            IDL_VPTR tmp;
            IDL_MEMINT nFN = frameNumbers.size();
            IDL_MEMINT dims[] = { nFN }; 
            int32_t* tmpData = (int32_t*)IDL_MakeTempArray( IDL_TYP_LONG, 1, dims, IDL_ARR_INI_NOP, &tmp );
            memcpy( tmpData, frameNumbers.data(), nFN*sizeof(int32_t));
            IDL_VarCopy( tmp, kw.framenumbers );
        }

        if( kw.nsummed ) {
            IDL_VPTR tmpNS = IDL_GettmpLong( nSummed );
            IDL_VarCopy( tmpNS, kw.nsummed );
        }
        
        if( kw.pinh_align && kw.xyc ) {
            IDL_VPTR tmp;
            IDL_MEMINT dims[] = { 2, static_cast<IDL_MEMINT>(nTotalFrames) }; 
            float* tmpData = (float*)IDL_MakeTempArray( IDL_TYP_FLOAT, 2, dims, IDL_ARR_INI_ZERO, &tmp ); //IDL_ARR_INI_ZERO
            memcpy( tmpData, *shifts, 2*nTotalFrames*sizeof(float));
            IDL_VarCopy( tmp, kw.xyc );
        }
        
        if( !time_beg.empty() ) {
            time_beg.erase( std::remove_if( time_beg.begin(), time_beg.end(),
                                       [](const bpx::ptime& a){ return a.is_special(); } ), time_beg.end() );
            time_end.erase( std::remove_if( time_end.begin(), time_end.end(),
                                       [](const bpx::ptime& a){ return a.is_special(); } ), time_end.end() );
            std::sort( time_beg.begin(), time_beg.end() );
            std::sort( time_end.begin(), time_end.end() );
        }
        
        if( kw.time_beg ) {
            string tStr = bpx::to_simple_string(time_beg.begin()->time_of_day());
            IDL_VPTR tmpTimeString = IDL_StrToSTRING( (char*)tStr.c_str() );
            IDL_VarCopy( tmpTimeString, kw.time_beg );
        }

        if( kw.time_end ) {
            string tStr = bpx::to_simple_string(time_end.rbegin()->time_of_day());
            IDL_VPTR tmpTimeString = IDL_StrToSTRING( (char*)tStr.c_str() );
            IDL_VarCopy( tmpTimeString, kw.time_end );
        }

        if( kw.time_avg ) {
            bpx::time_duration avgTime(0,0,0);
            time_beg.insert( time_beg.end(), time_end.begin(), time_end.end() );
            size_t tCount(0);
            for( bpx::ptime& t: time_beg ) {
                avgTime += t.time_of_day();
                tCount++;
            }
            if( tCount ) {
                avgTime /= tCount;
            }
            string tStr = bpx::to_simple_string(avgTime);
            IDL_VPTR tmpTimeString = IDL_StrToSTRING( (char*)tStr.c_str() );
            IDL_VarCopy( tmpTimeString, kw.time_avg );
        }

        if( nSummed ) {
            double nImg_inv = 1.0/nSummed;
            for( size_t i=0; i<nPixels; ++i ) summedData[i] *= nImg_inv;
        }

        if( !kw.pinh_align ) {

            applyDarkAndGain( summedData, summedData, darkData.get(), gainData.get(), nPixels );
            
            if( bsOtfPtr ) {
                descatter( summedData, ySize, xSize, bsGainPtr, bsOtfPtr, ySizePadded, xSizePadded, kw.nthreads, 50, 1E-8 );
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
    
#ifndef RDX_WITH_OPENCV
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
                        if( RDX_BYTE_ORDER == RDX_BIG_ENDIAN ) ptr++;
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
                        if( RDX_BYTE_ORDER == RDX_BIG_ENDIAN ) ptr++;
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


namespace mz {
    
    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT blend;
        IDL_INT clip;
        IDL_INT crop;
        IDL_INT help;
        IDL_INT margin;
        IDL_INT transpose;
        IDL_INT verbose;
    } KW_RESULT;

    IDL_KW_PAR kw_pars[] = {   IDL_KW_FAST_SCAN,       // NOTE:  The keywords MUST be listed in alphabetical order !!
        { ( char* ) "BLEND",         IDL_TYP_INT, 1,           0,                 0, ( char* ) IDL_KW_OFFSETOF ( blend ) },
        { ( char* ) "CLIP",          IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( clip ) },
        { ( char* ) "CROP",          IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( crop ) },
        { ( char* ) "HELP",          IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( help ) },
        { ( char* ) "MARGIN",        IDL_TYP_INT, 1,           0,                 0, ( char* ) IDL_KW_OFFSETOF ( margin ) },
        { ( char* ) "TRANSPOSE",     IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( transpose ) },
        { ( char* ) "VERBOSE",       IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( verbose ) },
        { NULL }
    };
        
    string info( int lvl ) {
        
        string ret = "RDX_MOZAIC";
        if( lvl > 0 ) {
            ret += ((lvl > 1)?"\n":"         ");          // newline if lvl>1
            ret += "   Syntax:   img = rdx_mozaic(patches, x_pos, y_pos, /KEYWORDS)\n";
            if( lvl > 1 ) {
//                 ret +=  "   Combine patches into an image.\n"
//                         "   patches has to be an array of 2D images and the positions has to have the same number of elements as patches.";
                ret +=  "   Accepted Keywords:\n"
                        "      BLEND               Size of smoothing region. (default = PatchSize/8)\n"
                        "      CLIP                Remove empty border after mozaic.\n"
                        "      CROP                Crop away the fixed boundary.\n"
                        "      HELP                Display this info.\n"
                        "      MARGIN              Ignore outermost m pixels in each patch (default = PatchSize/8)\n"
                        "      TRANSPOSE           Transpose each patch before mozaic, and transpose the final image back afterwards.\n"
                        "      VERBOSE             Verbosity, default is 0 (only error output).\n";
            }
        } else ret += "\n";

        return ret;
        
    }

}


IDL_VPTR rdx_mozaic( int argc, IDL_VPTR *argv, char *argk ) {

    mz::KW_RESULT kw;
    kw.blend = kw.margin = -1;
    IDL_KWProcessByOffset( argc, argv, argk, mz::kw_pars, (IDL_VPTR*)0, 255, &kw );

    if( kw.help ) {
        cout << mz::info(2) << endl;
        return IDL_GettmpInt(0);
    }
    
    IDL_VPTR patches = argv[0];
    IDL_ENSURE_SIMPLE( patches );
    IDL_ENSURE_ARRAY( patches );
    
    if( patches->value.arr->n_dim < 2 ) {
        cout << "Patches must be >= 2D." << endl;
        cout << mz::info(2) << endl;
        return IDL_GettmpInt(0);
    }
    
    size_t pCols  = patches->value.arr->dim[0];
    size_t pRows  = patches->value.arr->dim[1];
    size_t nPatches = 1;
    if( patches->value.arr->n_dim > 2 ) {
        nPatches = patches->value.arr->dim[2];
        for( UCHAR i=3; i<patches->value.arr->n_dim; ++i ) {
            nPatches *= patches->value.arr->dim[i];
        }
    }
    
    vector<int32_t> posX = getAsVector<int32_t>( argv[1] );
    vector<int32_t> posY = getAsVector<int32_t>( argv[2] );
    size_t nPosX = posX.size();
    size_t nPosY = posY.size();
    
    if( nPosX != nPatches || nPosY != nPatches ) {
        cout << "Number of positions must match number of patches." << endl;
        cout << mz::info(2) << endl;
        return IDL_GettmpInt(0);
    }
    
   
    int32_t maxPosX(0);
    int32_t minPosX( std::numeric_limits<int32_t >::max() );
    for( auto& i: posX ) {
        if( i > maxPosX ) maxPosX = i;
        if( i < minPosX ) minPosX = i;
    }
    maxPosX -= minPosX;
    for( auto& i: posX ) i -= minPosX;
    
    int32_t maxPosY(0);
    int32_t minPosY( std::numeric_limits<int32_t >::max() );
    for( auto& i: posY ) {
        if( i > maxPosY ) maxPosY = i;
        if( i < minPosY ) minPosY = i;
    }
    maxPosY -= minPosY;
    for( auto& i: posY ) i -= minPosY;
    
    if( kw.margin < 0 ) {
        kw.margin = std::max(pCols,pRows)/8; // number of pixels to cut from the edges of each path,
    }
    if( kw.blend < 0 ) {
        kw.blend = (std::min(pCols,pRows)-2*kw.margin)/3;
    }
    

    if ( kw.margin > int(pCols/2) || kw.margin > int(pRows/2) ) {
        cout << "Margin is too big, nothing will be left." << endl;
        return IDL_GettmpInt(0);
    }
    
    size_t imgCols = maxPosX+pCols;
    size_t imgRows = maxPosY+pRows;

    if ( kw.verbose > 0 ) {
        cout << "Mozaic:  nPatches = " << nPatches << endl;
        cout << "          imgSize = (" << imgCols << "," << imgRows << ")" << endl;
        cout << "        patchSize = (" << pCols << "," << pRows << ")" << endl;
        if ( kw.verbose > 1 ) {
            cout << "      " << printArray(posX,"patchPositionsX") << endl;
            cout << "      " << printArray(posY,"patchPositionsY") << endl;
        }
        if( kw.crop ) cout << "             crop = YES" << endl;
        else cout << "             clip = " << ( kw.clip?"YES":"NO" ) << endl;
        cout << "           margin = " << kw.margin << endl;
        cout << "            blend = " << kw.blend << endl;
        cout << "        transpose = " << ( kw.transpose?"YES":"NO" ) << endl;
    }

    UCHAR dataType = patches->type;
    IDL_VPTR ret;
    
    try {
        
        double** tmpImg = redux::util::newArray<double>( imgRows, imgCols );
    
        switch( dataType ) {
            case( IDL_TYP_BYTE ): {
                auto tmpPatches = reshapeArray( reinterpret_cast<const UCHAR*>(patches->value.arr->data), nPatches, pRows, pCols );
                mozaic( tmpImg, imgRows, imgCols, tmpPatches.get(), nPatches, pRows, pCols, posY.data(), posX.data(), kw.blend, kw.margin, kw.transpose );
                break;
            }
            case( IDL_TYP_INT ): {
                auto tmpPatches = reshapeArray( reinterpret_cast<const IDL_INT*>(patches->value.arr->data), nPatches, pRows, pCols );
                mozaic( tmpImg, imgRows, imgCols, tmpPatches.get(), nPatches, pRows, pCols, posY.data(), posX.data(), kw.blend, kw.margin, kw.transpose );
                break;
            }
            case( IDL_TYP_LONG ): {
                auto tmpPatches = reshapeArray( reinterpret_cast<const IDL_LONG*>(patches->value.arr->data), nPatches, pRows, pCols );
                mozaic( tmpImg, imgRows, imgCols, tmpPatches.get(), nPatches, pRows, pCols, posY.data(), posX.data(), kw.blend, kw.margin, kw.transpose );
                break;
            }
            case( IDL_TYP_FLOAT ): {
                auto tmpPatches = reshapeArray( reinterpret_cast<const float*>(patches->value.arr->data), nPatches, pRows, pCols );
                mozaic( tmpImg, imgRows, imgCols, tmpPatches.get(), nPatches, pRows, pCols, posY.data(), posX.data(), kw.blend, kw.margin, kw.transpose );
                break;
            }
            case( IDL_TYP_DOUBLE ): {
                auto tmpPatches = reshapeArray( reinterpret_cast<const double*>(patches->value.arr->data), nPatches, pRows, pCols );
                mozaic( tmpImg, imgRows, imgCols, tmpPatches.get(), nPatches, pRows, pCols, posY.data(), posX.data(), kw.blend, kw.margin, kw.transpose );
                break;
            }
            default: ;
        }
       
        size_t rm = std::max(pCols,pRows)/8; // hardcoded margin matching the one calculated in momfbd_read
        if( kw.crop && (imgRows>2*rm) && (imgCols>2*rm) ) {
            imgRows -= 2*rm;
            imgCols -= 2*rm;
            for( size_t i(0); i<imgRows; ++i ) {
                memcpy( *tmpImg+i*imgCols, tmpImg[rm+i]+rm, imgCols*sizeof(double) );
            }
        } else if( kw.clip ) {
            redux::image::img_trim( tmpImg, imgRows, imgCols, 1E-15 );
        }
        
        redux::util::transpose( *tmpImg, imgRows, imgCols );
        
        IDL_MEMINT dims[] = { static_cast<IDL_LONG64>(imgRows), static_cast<IDL_LONG64>(imgCols) };
        char* retData = IDL_MakeTempArray( dataType, 2, dims, IDL_ARR_INI_ZERO, &ret );
        copyToIDL( *tmpImg, reinterpret_cast<UCHAR*>(retData), imgRows*imgCols, dataType );
        
        delArray( tmpImg );
        
    } catch( const exception& e ) {
        cout << "rdx_mozaic: unhandled exception: " << e.what() << endl;
    }
    
    return ret;

}


namespace {
    static int dummy RDX_UNUSED =
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_find_shift, (char*)"RDX_FIND_SHIFT", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1 ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_fillpix, (char*)"RDX_FILLPIX", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, rdx_fillpix_info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_make_mask,  (char*)"RDX_MAKE_MASK",  0, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, make_mask_info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_make_win,  (char*)"RDX_MAKE_WINDOW",  1, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, apz::make_win_info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)sum_images, (char*)"RDX_SUMIMAGES", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, sum_images_info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_mozaic, (char*)"RDX_MOZAIC", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, mz::info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)sum_files,  (char*)"RDX_SUMFILES",  1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, sum_files_info );
}
