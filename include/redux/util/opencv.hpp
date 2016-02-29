#ifndef REDUX_UTIL_OPENCV_HPP
#define REDUX_UTIL_OPENCV_HPP

#ifdef REDUX_WITH_OPENCV

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

namespace cv {

enum
{
    MOTION_TRANSLATION = 0,
    MOTION_EUCLIDEAN   = 1,
    MOTION_AFFINE      = 2,
    MOTION_HOMOGRAPHY  = 3
};

/** @brief Finds the geometric transform (warp) between two images in terms of the ECC criterion @cite EP08 .

@param templateImage single-channel template image; CV_8U or CV_32F array.
@param inputImage single-channel input image which should be warped with the final warpMatrix in
order to provide an image similar to templateImage, same type as temlateImage.
@param warpMatrix floating-point \f$2\times 3\f$ or \f$3\times 3\f$ mapping matrix (warp).
@param motionType parameter, specifying the type of motion:
 -   **MOTION_TRANSLATION** sets a translational motion model; warpMatrix is \f$2\times 3\f$ with
     the first \f$2\times 2\f$ part being the unity matrix and the rest two parameters being
     estimated.
 -   **MOTION_EUCLIDEAN** sets a Euclidean (rigid) transformation as motion model; three
     parameters are estimated; warpMatrix is \f$2\times 3\f$.
 -   **MOTION_AFFINE** sets an affine motion model (DEFAULT); six parameters are estimated;
     warpMatrix is \f$2\times 3\f$.
 -   **MOTION_HOMOGRAPHY** sets a homography as a motion model; eight parameters are
     estimated;\`warpMatrix\` is \f$3\times 3\f$.
@param criteria parameter, specifying the termination criteria of the ECC algorithm;
criteria.epsilon defines the threshold of the increment in the correlation coefficient between two
iterations (a negative criteria.epsilon makes criteria.maxcount the only termination criterion).
Default values are shown in the declaration above.
@param inputMask An optional mask to indicate valid values of inputImage.

The function estimates the optimum transformation (warpMatrix) with respect to ECC criterion
(@cite EP08), that is

\f[\texttt{warpMatrix} = \texttt{warpMatrix} = \arg\max_{W} \texttt{ECC}(\texttt{templateImage}(x,y),\texttt{inputImage}(x',y'))\f]

where

\f[\begin{bmatrix} x' \\ y' \end{bmatrix} = W \cdot \begin{bmatrix} x \\ y \\ 1 \end{bmatrix}\f]

(the equation holds with homogeneous coordinates for homography). It returns the final enhanced
correlation coefficient, that is the correlation coefficient between the template image and the
final warped input image. When a \f$3\times 3\f$ matrix is given with motionType =0, 1 or 2, the third
row is ignored.

Unlike findHomography and estimateRigidTransform, the function findTransformECC implements an
area-based alignment that builds on intensity similarities. In essence, the function updates the
initial transformation that roughly aligns the images. If this information is missing, the identity
warp (unity matrix) should be given as input. Note that if images undergo strong
displacements/rotations, an initial transformation that roughly aligns the images is necessary
(e.g., a simple euclidean/similarity transform that allows for the images showing the same image
content approximately). Use inverse warping in the second image to take an image close to the first
one, i.e. use the flag WARP_INVERSE_MAP with warpAffine or warpPerspective. See also the OpenCV
sample image_alignment.cpp that demonstrates the use of the function. Note that the function throws
an exception if algorithm does not converges.

@sa
estimateRigidTransform, findHomography
 */
double findTransformECC( InputArray templateImage, InputArray inputImage,
                         InputOutputArray warpMatrix, int motionType = MOTION_AFFINE,
                         TermCriteria criteria = TermCriteria(TermCriteria::COUNT+TermCriteria::EPS, 50, 0.001),
                         InputArray inputMask = noArray());


} // cv

#endif          // if WITH_OPENCV

#endif      // REDUX_UTIL_OPENCV_HPP
