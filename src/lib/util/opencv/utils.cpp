#include <redux/util/opencv.hpp>


#include <iostream>
using namespace std;
using namespace cv;



void cv::make_mask( InputArray image, InputOutputArray mask, double thres, int smooth, bool filterLarger, bool invert) {
    
    Mat src = image.getMat();
    Mat dst = mask.getMat();
    
    CV_Assert( !src.empty() );
    CV_Assert( !dst.empty() );

    Mat graySrc;
    src.convertTo( graySrc, CV_32FC1 );
    
    threshold( graySrc, graySrc, thres, 1, THRESH_BINARY ); // THRESH_BINARY,THRESH_TOZERO
    
    double minValue, maxValue;
    minMaxLoc( graySrc, &minValue, &maxValue );
    if( maxValue == minValue ) {        // uniform input, no need to filter.
        if( invert ) {
            graySrc = 1 - graySrc;
        }
        graySrc.assignTo( dst, dst.type() );
        return;
    }

    if( smooth < 2 && !invert ) {         // if no filtering was requested, return 
        graySrc.assignTo( dst, dst.type() );
        return;
    }    
    
    Mat filteredMask( graySrc.rows, graySrc.cols, CV_8UC1 );
    graySrc.convertTo( filteredMask, CV_8UC1, 1 );
    
    if( smooth > 1 ) {
        Mat originalMask = filteredMask.clone();
        // shape:  enum { MORPH_RECT, MORPH_CROSS, MORPH_ELLIPSE };
        // operations: enum { MORPH_ERODE, MORPH_DILATE, MORPH_OPEN, MORPH_CLOSE, MORPH_GRADIENT, MORPH_TOPHAT, MORPH_BLACKHAT };
        Point anchor( smooth/2, smooth/2 );
        Mat element = getStructuringElement( MORPH_RECT, Size( smooth, smooth ), anchor );
        morphologyEx( filteredMask, filteredMask, MORPH_CLOSE, element, anchor, 1, BORDER_REFLECT_101 );
        if( !(smooth & 1) ) anchor -= Point( 1, 1 );    // move anchor-point for the reverse transform, so that asymmetry-shifts cancel
        morphologyEx( filteredMask, filteredMask, MORPH_OPEN, element, anchor, 1, BORDER_REFLECT_101 ); 
        bitwise_xor( filteredMask, originalMask, filteredMask );
        if( filterLarger ) {
            bitwise_or( 1-filteredMask, originalMask, filteredMask );
        } else {
            bitwise_xor( filteredMask, originalMask, filteredMask );
        }
    }
    
    if( invert ) {
        filteredMask = 1 - filteredMask;
    }
    
    filteredMask.assignTo( dst, dst.type() );

}

