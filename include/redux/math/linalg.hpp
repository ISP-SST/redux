#ifndef REDUX_MATH_LINALG_HPP
#define REDUX_MATH_LINALG_HPP

#include "redux/util/array.hpp"

#include <gsl/gsl_permutation.h>

namespace redux {

    namespace math {

        /*! Singlular value decomposition
         *  If A is (m x n) then svd produces a vector S (min(m,n) x 1) with (non-negative) singular values in decreasing order,
         *  and unitary matrices U (m x m) and V (n x n) so that A = U*(I*S)*V^T
         */
 /*       template<typename T>
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
 */       
        void normalize_columns( double* A, size_t rows, size_t cols );
        void normalize_columns( redux::util::Array<double>& A );
        
        void normalize_rows( double* A, size_t rows, size_t cols );
        void normalize_rows( redux::util::Array<double>& A );
        
        void gauss_jordan( double* A, size_t rows, size_t cols );
        void gauss_jordan( redux::util::Array<double>& A );
        
        void lu_decomp( double* A, int rows, int cols );
        void lu_decomp( redux::util::Array<double>& A);
        
        void qr_decomp(const double* A, int rows, int cols, double* Q, double* R);
        void qr_decomp(const redux::util::Array<double>& A, redux::util::Array<double>& Q, redux::util::Array<double>& R);
        
        void qr_decomp_pivot(const double* A, int rows, int cols, double* Q, double* R, gsl_permutation* p);
        void qr_decomp_pivot(const double* A, int rows, int cols, double* Q, double* R);
        void qr_decomp_pivot(const redux::util::Array<double>& A, redux::util::Array<double>& Q, redux::util::Array<double>& R, gsl_permutation* p);
        
        void cod_decomp(const double* A, int rows, int cols, double* Q, double* R, double* Z);
        void cod_decomp(const redux::util::Array<double>& A, redux::util::Array<double>& Q, redux::util::Array<double>& R, redux::util::Array<double>& Z);
        
        void svd(double* A_U, int rows, int cols, double* S, double* V);
        
     
        
 /*       template<typename T>
        void svd(T* data, int rows, int cols, T* sigma, T* u, T* v) {
            
            cv::Mat data_(rows, cols, cv::DataType<T>::type, data); // does not copy or free
            cv::Mat sigma_(rows, cols, cv::DataType<T>::type, sigma);
            cv::Mat u_(rows, rows, cv::DataType<T>::type, u);
            cv::Mat v_(cols, cols, cv::DataType<T>::type, v);
            
//             cv::Mat data_(rows, cols, CV_64F, data); // does not copy or free
//             cv::Mat sigma_(std::min(rows, cols), 1, CV_64F, sigma, CV_AUTOSTEP);
//             cv::Mat u_(rows, rows, CV_64F, u, CV_AUTOSTEP);
//             cv::Mat v_(cols, cols, CV_64F, v, CV_AUTOSTEP);
//             
            //cvSVD( &data_, &sigma_, &u_, &v_ );
            cv::SVD::compute( data_, sigma_, u_, v_ );
            //cvSVD( &data_, &sigma_, &u_, &v_ );

        }
*/

    }

}

#endif      //  REDUX_MATH_LINALG_HPP
