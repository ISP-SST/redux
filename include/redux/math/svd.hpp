#ifndef REDUX_MATH_SVD_HPP
#define REDUX_MATH_SVD_HPP

#include <opencv/cv.h>

namespace redux {

    namespace math {

        /*! Singlular value decomposition
         *  If A is (m x n) then svd produces a vector S (min(m,n) x 1) with (non-negative) singular values in decreasing order,
         *  and unitary matrices U (m x m) and V (n x n) so that A = U*(I*S)*V^T
         */
        template<typename T>
        void svd2(T* A, int rows, int cols, T* S, T* U=nullptr, T* V=nullptr) {
            
            CvMat A_ = cvMat(rows, cols, cv::DataType<T>::type, A);
            CvMat S_ = cvMat(std::min(rows, cols), 1, cv::DataType<T>::type, S);
            if( U && V ) {
                CvMat U_ = cvMat(rows, rows, cv::DataType<T>::type, U); 
                CvMat V_ = cvMat(cols, cols, cv::DataType<T>::type, V); 
                cvSVD( &A_, &S_, &U_, &V_ );
            } else {
                cvSVD( &A_, &S_ );
            }

        }
        
        void svd(double* A_U, int rows, int cols, double* S, double* V);
        
     
        

    }

}

#endif      //  REDUX_MATH_SVD_HPP
